/* Copyright (C) 2012-2017 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */
#include <NTL/BasicThreadPool.h>

#include "recryption.h"
#include "EncryptedArray.h"
#include "EvalMap.h"
#include "powerful.h"
#include "CtPtrs.h"
#include "intraSlot.h"
#include "norms.h"
#include "sample.h"
#include "debugging.h"

NTL_CLIENT

long thinRecrypt_initial_level=0;

#define PRINT_LEVELS

/*constexpr*/ double RecryptData::magicConst = 1.0;

/************* Some local functions *************/
static void x2iInSlots(ZZX& poly, long i,
		       vector<ZZX>& xVec, const EncryptedArray& ea);


static inline double pow(long a, long b) {return pow(double(a), double(b));}

RecryptData::~RecryptData()
{
  if (alMod!=NULL)     delete alMod;
  if (ea!=NULL)        delete ea;
  if (firstMap!=NULL)  delete firstMap;
  if (secondMap!=NULL) delete secondMap;
  if (p2dConv!=NULL)   delete p2dConv;
}



/** We want to get the smallest value of e-e', subject to a
 * few constraints: For RecryptData::magicConst, an exponent e,
 * a norm-t secret key, and plaintext space mod p^r, we need to
 * find integers a and b so that:
 *
 *    (1) (4a+ 8p^r)(t+1)*magicConst <= q  = p^e +1
 *    (2) (4b+ 8)(t+1) * magicConst <= q-4= p^e -3
 *
 * Then e' is the largest exponent such that p^{e'} <= 2(a+b).
 *
 * Note that if we let e,e' tend to infinity and set a=b=p^{e'}/4,
 * then the two constraints above degenerate to
 *
 *    2(t+1)*magicConst < p^{e-e'}
 *
 * so the smallest value of e-e' that we can hope for is
 *
 *    e-e' = ceiling( log_p( 2(t+1)*magicConst ) )
 *
 * The setAE procedure tries to find a setting of e,e',a (and
 * b = p^{e'}/2 -a) that satisfies the constraints (1,2) and
 * yeilds the smallest e-e', so long as e is "not too big".
 * Specifically, it minimizes e-e' while ensuring p^e < 2^{30},
 * if possible (and failing that it will just satisfy the
 * constraints (1,2)).
 *
 * Once e, e', a are set, it copmutes and returns the largest 
 * Hamming-weight for the key for which constraints (1,2) still hold.
 * NOTE: it returns the Hamming weight, *not* the size of the key.
 *       The size can be computed by calling the function
 *       sampleHWtBoundedEffectiveBound(context, weight)
 */
long RecryptData::setAE(long& a, long& e, long& ePrime,
                    const FHEcontext& context, long t)
{
  if (t<=0) t = RecryptData::defSkHwt;
  double bound = (t+1) * RecryptData::magicConst * 4;
  long p = context.zMStar.getP();
  //  long r = context().alMod.getR();
  long p2r = context.alMod.getPPowR();

  // We must have p^e > 8 p^r (t+1)*magicCons
  double eMin = ceil(log(1+ 2*bound*p2r)/log(p));

  // make sure that p^e for this smallest e is single-precision
  assert(eMin*log(p) < log(NTL_SP_BOUND));

  // The loop below tries to find e such that p^e < 2^30
  a = 0;
  e = eMin;
  ePrime = 0;
  long eMinusEprime = e; // want to minimize e-e'
  long eTry = e;
  long b = 0;
  do {
    double p2e = pow(p,eTry);

    // Solve for the largest a,b satisfying constraints (1,2)
    long aTry = floor( (p2e+1)/bound ) - 2*p2r;
    long bTry = floor( (p2e-3)/bound ) - 2;
    assert(aTry > 0 && bTry > 0);     // sanity check
    aTry -= (aTry % p2r); // reduce to nearest multiply of p^r

    long ePrimeTry = floor( log(2*(aTry+bTry))/log(p) );
    if (eTry - ePrimeTry < eMinusEprime) {
      e = eTry;
      ePrime = ePrimeTry;
      eMinusEprime = e - ePrime;

      // Set the value of a, reduced by a bit if possible
      long pToEprimeOver2 = floor(pow(p,ePrime)/2.0);
      b = bTry;
      a = pToEprimeOver2 - bTry;
      a -= (a % p2r);
      if (a<0) a=0;
    }
  } while ((++eTry)*log(p) <= log(1L<<30));  

  // Try to increase t, while maintaining constraints (1,2)
  double boundA = (4*a + 8*p2r)*RecryptData::magicConst;
  double boundB = (4*b + 8)*RecryptData::magicConst;
  long q = ceil(pow(p,e))+1;
  while (t++) {
    if (t>512 || boundA*t > q || boundB*t > q-4)
      break;
  }
  t--;

  long hwt;  // Find the largest Hamming weight that yeilds size <= t
  for (hwt=256; hwt>0; hwt--) {
    if (sampleHWtBoundedEffectiveBound(context, hwt)<=t)
      break;
  }  
#ifdef DEBUG_PRINTOUT
  cerr << "RecryptData::setAE(): e="<<e<<", e'="<<ePrime
       << ", a="<<a<<", sk-hwt="<<hwt<<" (size="<<t<<")\n";
#endif
  return hwt;
}


bool RecryptData::operator==(const RecryptData& other) const
{
  if (mvec != other.mvec) return false;
  if (hwt != other.hwt) return false;

  return true;
}



long fhe_disable_fat_boot = 0;

// The main method
void RecryptData::init(const FHEcontext& context, const Vec<long>& mvec_,
		       long t, bool build_cache_, bool minimal)
{
  if (alMod != NULL) { // were we called for a second time?
    cerr << "@Warning: multiple calls to RecryptData::init\n";
    return;
  }
  assert(computeProd(mvec_) == (long)context.zMStar.getM()); // sanity check

  // Record the arguments to this function
  mvec = mvec_;
  build_cache = build_cache_;

  hwt = setAE(a, e, ePrime, context, t);
  long p = context.zMStar.getP();
  long r = context.alMod.getR();

  // First part of Bootstrapping works wrt plaintext space p^{r'}
  alMod = new PAlgebraMod(context.zMStar, e-ePrime+r);
  ea = new EncryptedArray(context, *alMod);
         // Polynomial defaults to F0, PAlgebraMod explicitly given

  p2dConv = new PowerfulDCRT(context, mvec);

  if (fhe_disable_fat_boot) return;

  // Initialize the linear polynomial for unpacking the slots
  zz_pBak bak; bak.save(); ea->getAlMod().restoreContext();
  long nslots = ea->size();
  long d = ea->getDegree();

  const Mat<zz_p>& CBi=ea->getDerived(PA_zz_p()).getNormalBasisMatrixInverse();

  vector<ZZX> LM;
  LM.resize(d);
  for (long i = 0; i < d; i++) // prepare the linear polynomial
    LM[i] = rep(CBi[i][0]);

  vector<ZZX> C; 
  ea->buildLinPolyCoeffs(C, LM); // "build" the linear polynomial

  unpackSlotEncoding.resize(d);  // encode the coefficients

  for (long j = 0; j < d; j++) {
    vector<ZZX> v(nslots);
    for (long k = 0; k < nslots; k++) v[k] = C[j];
    ea->encode(unpackSlotEncoding[j], v);
  }
  firstMap = new EvalMap(*ea, minimal, mvec, true, build_cache);
  secondMap = new EvalMap(*context.ea, minimal, mvec, false, build_cache);
}

/********************************************************************/
/********************************************************************/

#ifdef DEBUG_PRINTOUT
#include "debugging.h"
long printFlag = FLAG_PRINT_VEC;
#endif

// Extract digits from fully packed slots
void extractDigitsPacked(Ctxt& ctxt, long botHigh, long r, long ePrime,
			 const vector<ZZX>& unpackSlotEncoding);

// Extract digits from unpacked slots
void extractDigitsThin(Ctxt& ctxt, long botHigh, long r, long ePrime);

static
long makeDivisible(vec_ZZ& vec, long p2e, long p2r, long q, long a, 
                   double& U_norm, const PAlgebra& palg);

 
// bootstrap a ciphertext to reduce noise
void FHEPubKey::reCrypt(Ctxt &ctxt)
{
  FHE_TIMER_START;

  // Some sanity checks for dummy ciphertext
  long ptxtSpace = ctxt.getPtxtSpace();
  if (ctxt.isEmpty()) return;
  if (ctxt.parts.size()==1 && ctxt.parts[0].skHandle.isOne()) {
    // Dummy encryption, just ensure that it is reduced mod p
    ZZX poly = to_ZZX(ctxt.parts[0]);
    for (long i=0; i<poly.rep.length(); i++)
      poly[i] = to_ZZ( rem(poly[i],ptxtSpace) );
    poly.normalize();
    ctxt.DummyEncrypt(poly);
    return;
  }

  assert(recryptKeyID>=0); // check that we have bootstrapping data

  long p = getContext().zMStar.getP();
  long r = getContext().alMod.getR();
  long p2r = getContext().alMod.getPPowR();

  long intFactor = ctxt.intFactor;

  // the bootstrapping key is encrypted relative to plaintext space p^{e-e'+r}.
  const RecryptData& rcData = getContext().rcData;
  long e = rcData.e;
  long ePrime = rcData.ePrime;
  long p2ePrime = power_long(p,ePrime);
  long q = power_long(p,e)+1;
  assert(e>=r);

#ifdef DEBUG_PRINTOUT
  cerr << "reCrypt: p="<<p<<", r="<<r<<", e="<<e<<" ePrime="<<ePrime
       << ", q="<<q<<endl;
#endif

  // can only bootstrap ciphertext with plaintext-space dividing p^r
  assert(p2r % ptxtSpace == 0);


#ifdef PRINT_LEVELS
  CheckCtxt(ctxt, "init");
#endif

  ctxt.dropSmallAndSpecialPrimes();

#ifdef PRINT_LEVELS
  CheckCtxt(ctxt, "after mod down");
#endif



  FHE_NTIMER_START(AAA_preProcess);

  // Make sure that this ciphertxt is in canonical form
  if (!ctxt.inCanonicalForm()) ctxt.reLinearize();

  // Mod-switch down if needed
  IndexSet s = ctxt.getPrimeSet() / context.specialPrimes;
  assert(s <= context.ctxtPrimes);
  if (s.card()>2) { // leave only bottom two primes
    long frst = s.first();
    long scnd = s.next(frst);
    IndexSet s2(frst,scnd);
    s.retain(s2); // retain only first two primes
  }
  ctxt.modDownToSet(s);

  // key-switch to the bootstrapping key
  ctxt.reLinearize(recryptKeyID);

  // "raw mod-switch" to the bootstrapping mosulus q=p^e+1.
  vector<ZZX> zzParts; // the mod-switched parts, in ZZX format
  double noise = ctxt.rawModSwitch(zzParts, q);

  // Add multiples of p2r and q to make the zzParts divisible by p^{e'}
  long maxU=0;
  double maxU_norm = 0;
  for (long i=0; i<(long)zzParts.size(); i++) {
    // make divisible by p^{e'}
    double U_norm;
    long newMax = makeDivisible(zzParts[i].rep, p2ePrime, p2r, q,
				rcData.a, U_norm, context.zMStar);
    zzParts[i].normalize();   // normalize after working directly on the rep
    if (maxU < newMax)  maxU = newMax;
    if (maxU_norm < U_norm)  maxU_norm = U_norm;
  }

  // Check that the estimated noise is still low
  if (noise + maxU_norm*p2r*(skBounds[recryptKeyID]+1) > q/2)
    cerr << " ******** warning: noise/q after makeDivisible = "
	 << ((noise + maxU_norm*p2r*(skBounds[recryptKeyID]+1))/q) << endl;

#ifdef PRINT_LEVELS
   if (dbgKey) {
     const RecryptData& rcData = ctxt.getContext().rcData;
     ZZX ptxt;
     rawDecrypt(ptxt, zzParts, dbgKey->sKeys[recryptKeyID], q);

     Vec<ZZ> powerful;
     ZZX ptxt_alt;
     rcData.p2dConv->ZZXtoPowerful(powerful, ptxt);
     vecRed(powerful, powerful, q, false);
     rcData.p2dConv->powerfulToZZX(ptxt_alt, powerful);

#if 1
     xdouble max_pwrfl = conv<xdouble>(largestCoeff(powerful));
     xdouble max_canon = embeddingLargestCoeff(ptxt_alt,rcData.ea->getPAlgebra());
     double noise_bnd = noise + maxU_norm*p2r*(skBounds[recryptKeyID]+1);

     cerr << "  after makeDivisible";

     cerr << ", maxU=" << maxU;
     cerr << ", maxU_norm=" << maxU_norm;

     cerr << ", max_pwrfl/q=" << (max_pwrfl/q);

     double ratio0 = log(max_canon/noise_bnd)/log(2.0);
     cerr << ", log2(max_canon/bound)=" << ratio0;
     if (ratio0 > 0) cerr << " BAD-BOUND";

     double ratio1 = log(max_pwrfl/max_canon)/log(2.0);
     cerr << ", log2(max_pwrfl/max_canon)=" << ratio1;
     if (ratio1 > 0) cerr << " BAD-BOUND";

     cerr << "\n";
#else
     cerr << "  after makeDivisible, noiseEst="
	  << (noise + maxU*p2r*(skBounds[recryptKeyID]+1))
	  << ", maxCanon="
	  << embeddingLargestCoeff(ptxt,rcData.ea->getPAlgebra());
     cerr << ",  maxCoeff="<<largestCoeff(ptxt)
	  << ", maxPowfl="<<largestCoeff(powerful)
	  << endl;
#endif
  }
#endif

  for (long i=0; i<(long)zzParts.size(); i++)
    zzParts[i] /= p2ePrime;   // divide by p^{e'}

  // Multiply the post-processed cipehrtext by the encrypted sKey
#ifdef DEBUG_PRINTOUT
  cerr << "+ Before recryption ";
  decryptAndPrint(cerr, recryptEkey, *dbgKey, *dbgEa, printFlag);
#endif

  double p0size = to_double(embeddingLargestCoeff(zzParts[0], context.zMStar));
  double p1size = to_double(embeddingLargestCoeff(zzParts[1], context.zMStar));
  // FIXME: This might be slow without Armadillo

  // NOTE: here we lose the intFactor associated with ctxt.
  // We will restore it below.
  ctxt = recryptEkey;

  ctxt.multByConstant(zzParts[1]);
  ctxt.addConstant(zzParts[0]);

#ifdef DEBUG_PRINTOUT
  cerr << "+ Before linearTrans1 ";
  decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, printFlag);
#endif
  FHE_NTIMER_STOP(AAA_preProcess);

#ifdef PRINT_LEVELS
  CheckCtxt(ctxt, "after preProcess");
#endif


  // Move the powerful-basis coefficients to the plaintext slots
  FHE_NTIMER_START(AAA_LinearTransform1);
  ctxt.getContext().rcData.firstMap->apply(ctxt);
  FHE_NTIMER_STOP(AAA_LinearTransform1);


#ifdef PRINT_LEVELS
  CheckCtxt(ctxt, "after LinearTransform1");
#endif

#ifdef DEBUG_PRINTOUT
  cerr << "+ After linearTrans1 ";
  decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, printFlag);
#endif

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)
  FHE_NTIMER_START(AAA_extractDigitsPacked);
  extractDigitsPacked(ctxt, e-ePrime, r, ePrime,
		      context.rcData.unpackSlotEncoding);
  FHE_NTIMER_STOP(AAA_extractDigitsPacked);


#ifdef PRINT_LEVELS
  CheckCtxt(ctxt, "after extractDigitsPacked");
#endif

#ifdef DEBUG_PRINTOUT
  cerr << "+ Before linearTrans2 ";
  decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, printFlag);
#endif

  // Move the slots back to powerful-basis coefficients
  FHE_NTIMER_START(AAA_LinearTransform2);
  ctxt.getContext().rcData.secondMap->apply(ctxt);
  FHE_NTIMER_STOP(AAA_LinearTransform2);


#ifdef PRINT_LEVELS
  CheckCtxt(ctxt, "after linearTransform2");
#endif

  // restore intFactor
  if (intFactor != 1)
    ctxt.intFactor = MulMod(ctxt.intFactor, intFactor, ptxtSpace);
}

/*********************************************************************/
/*********************************************************************/

// Return in poly a polynomial with X^i encoded in all the slots
static void x2iInSlots(ZZX& poly, long i,
		       vector<ZZX>& xVec, const EncryptedArray& ea)
{
  xVec.resize(ea.size());
  ZZX x2i = ZZX(i,1);
  for (long j=0; j<(long)xVec.size(); j++) xVec[j] = x2i;
  ea.encode(poly, xVec);
}

// Make every entry of vec divisible by p2e by adding/subtracting
// multiples of p2r and q, while keeping the added multiples small.
// Specifically, for q = 1 mod p2e and any a < p2e/(2*p2r), any
// integer z can be made divisible by p2e via z' = z + u*p2r + v*q,
// with |u|*p2r <= a and |v| <= p2e/2 -a.
// Returns the largest absolute values of the u's and the new entries.
static long makeDivisible(vec_ZZ& vec, long p2e, long p2r, long q, long a, 
                          double& U_norm, const PAlgebra& palg)
{
  assert(q>0 && p2e>0 && p2r>0 && a>=0
         && q % p2e == 1 && a % p2r == 0 && a*2 < p2e);
  long aa = a / p2r;

  ZZX vec_orig;  conv(vec_orig, vec);
  PolyRed(vec_orig, vec_orig, q);

#ifdef DEBUG_PRINTOUT
  ZZX vec_orig; conv(vec_orig, vec); // backup vector
  zzX uVec(INIT_SIZE, vec.length());
  zzX vVec(INIT_SIZE, vec.length());
#endif

  long maxU = 0;
  for (long i=0; i<vec.length(); i++) {
    ZZ& z = vec[i];
    long u, v;

    // What to add to z to make it divisible by p2e?
    long zMod = rem(z, p2e); // zMod is in [0,p2e-1]
    if (zMod > p2e/2) { // need to add a positive number
      zMod = p2e - zMod;
      u = zMod/p2r;
      if (u > a) u = a;
    }
    else {              // need to add a negative number
      u = -(zMod/p2r);
      if (u < -a) u = -a;
      zMod = -zMod;
    }
    v = zMod - u*p2r;
    z += u*p2r + to_ZZ(q)*v; // make z divisible by p2e

    if (rem(z,p2e) != 0) { // sanity check
      cerr << "**error: original z["<<i<<"]=" << (z-(u*p2r+to_ZZ(q)*v))
	   << std::dec << ", p^r="<<p2r << ", p^e="<<p2e << endl;
      cerr << "z' = z + "<<u<<"*p^r +"<<v<<"*q = "<<z<<endl;
      exit(1);
    }
    if (abs(u) > maxU) maxU = abs(u);
#ifdef DEBUG_PRINTOUT
    uVec[i] = u;
    vVec[i] = v;
#endif
  }

#ifdef DEBUG_PRINTOUT
  if (dbgEa) {
    const PAlgebra& palg = dbgEa->getPAlgebra();
    double U_norm = conv<double>(embeddingLargestCoeff(uVec, palg));
    double V_norm = conv<double>(embeddingLargestCoeff(vVec, palg));
    cerr << "  makeDivisible: maxU=" << (maxU*p2r)
         << ", U_norm=" << (U_norm*p2r)
         << ", V_norm=" << (V_norm*q)   << endl;
  }
#endif

  ZZX vec_new;
  conv(vec_new, vec);
  vec_new = vec_new - vec_orig;
  PolyRed(vec_new, vec_new, q);

  assert(divide(vec_new, vec_new, p2r));

  U_norm = conv<double>(embeddingLargestCoeff(vec_new, palg));

  return maxU;
}

#ifdef FHE_BOOT_THREADS

// Extract digits from fully packed slots, multithreaded version
void extractDigitsPacked(Ctxt& ctxt, long botHigh, long r, long ePrime,
			 const vector<ZZX>& unpackSlotEncoding)
{
  FHE_TIMER_START;

  // Step 1: unpack the slots of ctxt
  FHE_NTIMER_START(unpack);
  ctxt.cleanUp();

  // Apply the d automorphisms and store them in scratch area
  long d = ctxt.getContext().zMStar.getOrdP();

  vector<Ctxt> unpacked(d, Ctxt(ZeroCtxtLike, ctxt));
  { // explicit scope to force all temporaries to be released
    vector< shared_ptr<DoubleCRT> > coeff_vector;
    coeff_vector.resize(d);

    FHE_NTIMER_START(unpack1);
    for (long i = 0; i < d; i++)
      coeff_vector[i] = shared_ptr<DoubleCRT>(new 
        DoubleCRT(unpackSlotEncoding[i], ctxt.getContext(), ctxt.getPrimeSet()) );
    FHE_NTIMER_STOP(unpack1);

    FHE_NTIMER_START(unpack2);
    vector<Ctxt> frob(d, Ctxt(ZeroCtxtLike, ctxt));

    NTL_EXEC_RANGE(d, first, last)
    // FIXME: implement using hoisting!
        for (long j = first; j < last; j++) { // process jth Frobenius 
          frob[j] = ctxt;
          frob[j].frobeniusAutomorph(j);
          frob[j].cleanUp();
          // FIXME: not clear if we should call cleanUp here
        }
    NTL_EXEC_RANGE_END

    FHE_NTIMER_STOP(unpack2);

    FHE_NTIMER_START(unpack3);
    Ctxt tmp1(ZeroCtxtLike, ctxt);
    for (long i = 0; i < d; i++) {
      for (long j = 0; j < d; j++) {
        tmp1 = frob[j];
        tmp1.multByConstant(*coeff_vector[mcMod(i+j, d)]);
        unpacked[i] += tmp1;
      }
    }
    FHE_NTIMER_STOP(unpack3);
  }
  FHE_NTIMER_STOP(unpack);

#ifdef DEBUG_PRINTOUT
  cerr << "+ After unpack ";
  decryptAndPrint(cerr, unpacked[0], *dbgKey, *dbgEa, printFlag);
#endif

  NTL_EXEC_RANGE(d, first, last)
      for (long i = first; i < last; i++) {
        extractDigitsThin(unpacked[i], botHigh, r, ePrime);
      }
  NTL_EXEC_RANGE_END

#ifdef DEBUG_PRINTOUT
  cerr << "+ Before repack ";
  decryptAndPrint(cerr, unpacked[0], *dbgKey, *dbgEa, printFlag);
#endif

  // Step 3: re-pack the slots
  FHE_NTIMER_START(repack);
  const EncryptedArray& ea2 = *ctxt.getContext().ea;
  ZZX xInSlots;
  vector<ZZX> xVec(ea2.size());
  ctxt = unpacked[0];
  for (long i=1; i<d; i++) {
    x2iInSlots(xInSlots, i, xVec, ea2);
    unpacked[i].multByConstant(xInSlots);
    ctxt += unpacked[i];
  }
  FHE_NTIMER_STOP(repack);
#ifdef DEBUG_PRINTOUT
  cerr << "+ After repack ";
  decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, printFlag);
#endif
}


#else

// Extract digits from fully packed slots
void extractDigitsPacked(Ctxt& ctxt, long botHigh, long r, long ePrime,
			 const vector<ZZX>& unpackSlotEncoding)
{
  FHE_TIMER_START;

  // Step 1: unpack the slots of ctxt
  FHE_NTIMER_START(unpack);
  ctxt.cleanUp();

  // Apply the d automorphisms and store them in scratch area
  long d = ctxt.getContext().zMStar.getOrdP();

  vector<Ctxt> unpacked(d, Ctxt(ZeroCtxtLike, ctxt));
  { // explicit scope to force all temporaries to be released
    vector< shared_ptr<DoubleCRT> > coeff_vector;
    coeff_vector.resize(d);
    for (long i = 0; i < d; i++)
      coeff_vector[i] = shared_ptr<DoubleCRT>(new 
        DoubleCRT(unpackSlotEncoding[i], ctxt.getContext(), ctxt.getPrimeSet()) );
    Ctxt tmp1(ZeroCtxtLike, ctxt);
    Ctxt tmp2(ZeroCtxtLike, ctxt);

    // FIXME: implement using hoisting!
    for (long j = 0; j < d; j++) { // process jth Frobenius 
      tmp1 = ctxt;
      tmp1.frobeniusAutomorph(j);
      tmp1.cleanUp();
      // FIXME: not clear if we should call cleanUp here

      for (long i = 0; i < d; i++) {
        tmp2 = tmp1;
        tmp2.multByConstant(*coeff_vector[mcMod(i+j, d)]);
        unpacked[i] += tmp2;
      }
    }
  }
  FHE_NTIMER_STOP(unpack);

#ifdef DEBUG_PRINTOUT
  cerr << "+ After unpack ";
  decryptAndPrint(cerr, unpacked[0], *dbgKey, *dbgEa, printFlag);
  cerr << "    extracting "<<(topHigh+1)<<" digits\n";
#endif

  for (long i=0; i<(long)unpacked.size(); i++) {
    extractDigitsThin(unpacked[i], botHigh, r, ePrime); 
  }

#ifdef DEBUG_PRINTOUT
  cerr << "+ Before repack ";
  decryptAndPrint(cerr, unpacked[0], *dbgKey, *dbgEa, printFlag);
#endif

  // Step 3: re-pack the slots
  FHE_NTIMER_START(repack);
  const EncryptedArray& ea2 = *ctxt.getContext().ea;
  ZZX xInSlots;
  vector<ZZX> xVec(ea2.size());
  ctxt = unpacked[0];
  for (long i=1; i<d; i++) {
    x2iInSlots(xInSlots, i, xVec, ea2);
    unpacked[i].multByConstant(xInSlots);
    ctxt += unpacked[i];
  }
  FHE_NTIMER_STOP(repack);
}

#endif


// Use packed bootstrapping, so we can bootstrap all in just one go.
void packedRecrypt(const CtPtrs& cPtrs,
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea)
{
  FHEPubKey& pKey = (FHEPubKey&)cPtrs[0]->getPubKey();

  // Allocate temporary ciphertexts for the recryption
  int nPacked = divc(cPtrs.size(), ea.getDegree()); // ceil(totoalNum/d)
  std::vector<Ctxt> cts(nPacked, Ctxt(pKey));

  repack(CtPtrs_vectorCt(cts), cPtrs, ea);  // pack ciphertexts
  //  cout << "@"<< lsize(cts)<<std::flush;
  for (Ctxt& c: cts) {     // then recrypt them
    c.reducePtxtSpace(2);  // we only have recryption data for binary ctxt
#ifdef DEBUG_PRINTOUT
    ZZX ptxt;
    decryptAndPrint((cout<<"  before recryption "), c, *dbgKey, *dbgEa);
    dbgKey->Decrypt(ptxt, c);
    c.DummyEncrypt(ptxt);
    decryptAndPrint((cout<<"  after recryption "), c, *dbgKey, *dbgEa);
#else
    pKey.reCrypt(c);
#endif
  }
  unpack(cPtrs, CtPtrs_vectorCt(cts), ea, unpackConsts);
}

// recrypt all ctxt at level < belowLvl
void packedRecrypt(const CtPtrs& array,
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea, long belowLvl)
{
  std::vector<Ctxt*> v;
  for (long i=0; i<array.size(); i++)
    if ( array.isSet(i) && !array[i]->isEmpty()
         && array[i]->bitCapacity()<belowLvl*(array[i]->getContext().BPL()) )
      v.push_back(array[i]);
  packedRecrypt(CtPtrs_vectorPt(v), unpackConsts, ea);
}
void packedRecrypt(const CtPtrMat& m,
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea, long belowLvl)
{
  std::vector<Ctxt*> v;
  for (long i=0; i<m.size(); i++)
    for (long j=0; j<m[i].size(); j++)
      if ( m[i].isSet(j) && !m[i][j]->isEmpty()
           && m[i][j]->bitCapacity()<belowLvl*(m[i][j]->getContext().BPL()) )
        v.push_back(m[i][j]);
  packedRecrypt(CtPtrs_vectorPt(v), unpackConsts, ea);
}



//===================== Thin Bootstrapping stuff ==================

ThinRecryptData::~ThinRecryptData()
{
  if (alMod!=NULL)     delete alMod;
  if (ea!=NULL)        delete ea;
  if (coeffToSlot!=NULL)  delete coeffToSlot;
  if (slotToCoeff!=NULL) delete slotToCoeff;
}

bool ThinRecryptData::operator==(const ThinRecryptData& other) const
{
  if (mvec != other.mvec) return false;
  if (hwt != other.hwt) return false;

  return true;
}


// This code was copied from RecryptData::init, and is mostly
// the same, except for the linear-map-related stuff.
// FIXME: There is really too much code (and data!) duplication here.
void ThinRecryptData::init(const FHEcontext& context, const Vec<long>& mvec_,
		       long t, bool build_cache_, bool minimal)
{
  if (alMod != NULL) { // were we called for a second time?
    cerr << "@Warning: multiple calls to ThinRecryptData::init\n";
    return;
  }
  assert(computeProd(mvec_) == (long)context.zMStar.getM()); // sanity check

  // Record the arguments to this function
  mvec = mvec_;
  build_cache = build_cache_;

  hwt = RecryptData::setAE(a,e,ePrime,context,t);
  long p = context.zMStar.getP();
  long r = context.alMod.getR();

  // First part of Bootstrapping works wrt plaintext space p^{r'}
  alMod = new PAlgebraMod(context.zMStar, e-ePrime+r);
  ea = new EncryptedArray(context, *alMod);
         // Polynomial defaults to F0, PAlgebraMod explicitly given

  coeffToSlot = new ThinEvalMap(*ea, minimal, mvec, true, build_cache);
  slotToCoeff = new ThinEvalMap(*context.ea, minimal, mvec, false, build_cache);
}


// Extract digits from thinly packed slots


long fhe_force_chen_han = 0;

void extractDigitsThin(Ctxt& ctxt, long botHigh, long r, long ePrime)
{
  FHE_TIMER_START;

  Ctxt unpacked(ctxt);
  unpacked.cleanUp();

  vector<Ctxt> scratch;

  long p = ctxt.getContext().zMStar.getP();
  long p2r = power_long(p,r);
  long topHigh = botHigh + r-1;


  // degree Chen/Han technique is p^{bot-1}(p-1)r
  // degree of basic technique is p^{bot-1}p^r, 
  //     or p^{bot-1}p^{r-1} if p==2, r > 1, and bot+r > 2

  bool use_chen_han = false;
  if (r > 1) {
    double chen_han_cost = log(p-1) + log(r);
    double basic_cost;
    if (p == 2 && botHigh + r > 2)
       basic_cost = (r-1)*log(p);
    else
       basic_cost = r*log(p);

    //cerr << "*** basic: " << basic_cost << "\n";
    //cerr << "*** chen/han: " << chen_han_cost << "\n";


    double thresh = 1.5;
    if (p == 2) thresh = 1.75;
    // increasing thresh makes chen_han less likely to be chosen.
    // For p == 2, the basic algorithm is just squaring, 
    // and so is a bit cheaper, so we raise thresh a bit.
    // This is all a bit heuristic.

    if (basic_cost > thresh*chen_han_cost)
      use_chen_han = true;
  }

  if (fhe_force_chen_han > 0)
    use_chen_han = true;
  else if (fhe_force_chen_han < 0)
    use_chen_han = false;


  if (use_chen_han) {
    // use Chen and Han technique

    extendExtractDigits(scratch, unpacked, botHigh, r);

#if 0
    for (long i: range(scratch.size())) {
      CheckCtxt(scratch[i], "**");
    }
#endif

    for (long j = 0; j < botHigh; j++) {
      unpacked -= scratch[j];
      unpacked.divideByP();
    }

    if (p==2 && botHigh>0)   // For p==2, subtract also the previous bit
      unpacked += scratch[botHigh-1];
    unpacked.negate();

    if (r>ePrime) {          // Add in digits from the bottom part, if any
      long topLow = r-1 - ePrime;
      Ctxt tmp = scratch[topLow];
      for (long j=topLow-1; j>=0; --j) {
	tmp.multByP();
	tmp += scratch[j];
      }
      if (ePrime>0)
	tmp.multByP(ePrime); // multiply by p^e'
      unpacked += tmp;
    }
    unpacked.reducePtxtSpace(p2r); // Our plaintext space is now mod p^r

    ctxt = unpacked;
  }
  else {

    if (p==2 && r>1 && topHigh+1 > 2)
      topHigh--; // For p==2 we sometime get a bit for free

    extractDigits(scratch, unpacked, topHigh+1);

    // set upacked = -\sum_{j=botHigh}^{topHigh} scratch[j] * p^{j-botHigh}
    if (topHigh >= LONG(scratch.size())) {
      topHigh = scratch.size() -1;
      cerr << " @ suspect: not enough digits in extractDigitsPacked\n";
    }

    unpacked = scratch[topHigh];
    for (long j=topHigh-1; j>=botHigh; --j) {
      unpacked.multByP();
      unpacked += scratch[j];
    }
    if (p==2 && botHigh>0)   // For p==2, subtract also the previous bit
      unpacked += scratch[botHigh-1];
    unpacked.negate();

    if (r>ePrime) {          // Add in digits from the bottom part, if any
      long topLow = r-1 - ePrime;
      Ctxt tmp = scratch[topLow];
      for (long j=topLow-1; j>=0; --j) {
	tmp.multByP();
	tmp += scratch[j];
      }
      if (ePrime>0)
	tmp.multByP(ePrime); // multiply by p^e'
      unpacked += tmp;
    }
    unpacked.reducePtxtSpace(p2r); // Our plaintext space is now mod p^r
    ctxt = unpacked;
  }

}


// Hack to get at private fields of public key
struct FHEPubKeyHack { // The public key
  const FHEcontext& context; // The context

  //! @var Ctxt pubEncrKey
  //! The public encryption key is an encryption of 0,
  //! relative to the first secret key
  Ctxt pubEncrKey;

  std::vector<long> skHwts; // The Hamming weight of the secret keys
  std::vector<KeySwitch> keySwitching; // The key-switching matrices

  // The keySwitchMap structure contains pointers to key-switching matrices
  // for re-linearizing automorphisms. The entry keySwitchMap[i][n] contains
  // the index j such that keySwitching[j] is the first matrix one needs to
  // use when re-linearizing s_i(X^n). 
  std::vector< std::vector<long> > keySwitchMap;

  NTL::Vec<int> KS_strategy; // NTL Vec's support I/O, which is
                             // more convenient

  // bootstrapping data

  long recryptKeyID; // index of the bootstrapping key
  Ctxt recryptEkey;  // the key itself, encrypted under key #0

};

//#define PRINT_LEVELS

// bootstrap a ciphertext to reduce noise
void FHEPubKey::thinReCrypt(Ctxt &ctxt)
{
  FHE_TIMER_START;

  // Some sanity checks for dummy ciphertext
  long ptxtSpace = ctxt.getPtxtSpace();
  if (ctxt.isEmpty()) return;

  if (ctxt.parts.size()==1 && ctxt.parts[0].skHandle.isOne()) {
    // Dummy encryption, just ensure that it is reduced mod p
    ZZX poly = to_ZZX(ctxt.parts[0]);
    for (long i=0; i<poly.rep.length(); i++)
      poly[i] = to_ZZ( rem(poly[i],ptxtSpace) );
    poly.normalize();
    ctxt.DummyEncrypt(poly);
    return;
  }

  assert(recryptKeyID>=0); // check that we have bootstrapping data

  long p = ctxt.getContext().zMStar.getP();
  long r = ctxt.getContext().alMod.getR();
  long p2r = ctxt.getContext().alMod.getPPowR();

  long intFactor = ctxt.intFactor;

  const ThinRecryptData& trcData = ctxt.getContext().trcData;

  // the bootstrapping key is encrypted relative to plaintext space p^{e-e'+r}.
  long e = trcData.e;
  long ePrime = trcData.ePrime;
  long p2ePrime = power_long(p,ePrime);
  long q = power_long(p,e)+1;
  assert(e>=r);

#ifdef DEBUG_PRINTOUT
  cerr << "reCrypt: p="<<p<<", r="<<r<<", e="<<e<<" ePrime="<<ePrime
       << ", q="<<q<<endl;
#endif

  // can only bootstrap ciphertext with plaintext-space dividing p^r
  assert(p2r % ptxtSpace == 0);

#ifdef PRINT_LEVELS
  CheckCtxt(ctxt, "init");
#endif

  ctxt.dropSmallAndSpecialPrimes();

  if (thinRecrypt_initial_level) {
    // experimental code...we should drop down
    // to a reasonably small level before doing the 
    // first linear map.
    long first = context.ctxtPrimes.first();
    long last = min(context.ctxtPrimes.last(), first + thinRecrypt_initial_level - 1);
    ctxt.bringToSet(IndexSet(first, last));
  }

#ifdef PRINT_LEVELS
  CheckCtxt(ctxt, "after mod down");
#endif

  // Move the slots to powerful-basis coefficients
  FHE_NTIMER_START(AAA_slotToCoeff);
  trcData.slotToCoeff->apply(ctxt);
  FHE_NTIMER_STOP(AAA_slotToCoeff);

#ifdef PRINT_LEVELS
  CheckCtxt(ctxt, "after slotToCoeff");
#endif

  FHE_NTIMER_START(AAA_bootKeySwitch);

  // Make sure that this ciphertxt is in canonical form
  if (!ctxt.inCanonicalForm()) ctxt.reLinearize();

  // Mod-switch down if needed
  IndexSet s = ctxt.getPrimeSet() / context.specialPrimes;
  assert(s <= context.ctxtPrimes);
  if (s.card()>2) { // leave only bottom two primes
    long frst = s.first();
    long scnd = s.next(frst);
    IndexSet s2(frst,scnd);
    s.retain(s2); // retain only first two primes
  }
  ctxt.modDownToSet(s);

  // key-switch to the bootstrapping key
  ctxt.reLinearize(recryptKeyID);

  // "raw mod-switch" to the bootstrapping mosulus q=p^e+1.
  vector<ZZX> zzParts; // the mod-switched parts, in ZZX format
  double noise = ctxt.rawModSwitch(zzParts, q);

#ifdef PRINT_LEVELS

   if (dbgKey) {
     const RecryptData& rcData = ctxt.getContext().rcData;
     ZZX ptxt;
     rawDecrypt(ptxt, zzParts, dbgKey->sKeys[recryptKeyID], q);

     Vec<ZZ> powerful;
     ZZX ptxt_alt;
     rcData.p2dConv->ZZXtoPowerful(powerful, ptxt);
     vecRed(powerful, powerful, q, false);
     rcData.p2dConv->powerfulToZZX(ptxt_alt, powerful);

     xdouble max_pwrfl = conv<xdouble>(largestCoeff(powerful));
     xdouble max_canon = embeddingLargestCoeff(ptxt_alt,rcData.ea->getPAlgebra());
     double noise_bnd = noise;

     ZZX skey_poly;
     dbgKey->sKeys[recryptKeyID].toPoly(skey_poly);
     double skey_bound = skBounds[recryptKeyID];
     double skey_canon = conv<double>(embeddingLargestCoeff(skey_poly,rcData.ea->getPAlgebra()));

     cerr << "  before makeDivisible";

     cerr << ", max_canon=" << max_canon;
     cerr << ", skey_canon=" << skey_canon;

     double ratio0 = log(max_canon/noise_bnd)/log(2.0);
     cerr << ", log2(max_canon/bound)=" << ratio0;
     if (ratio0 > 0) cerr << " BAD-BOUND";

     double ratio1 = log(max_pwrfl/max_canon)/log(2.0);
     cerr << ", log2(max_pwrfl/max_canon)=" << ratio1;
     if (ratio1 > 0) cerr << " BAD-BOUND";

     double ratio2 = log(skey_canon/skey_bound)/log(2.0);
     cerr << ", log2(skey_canon/skey_bound)=" << ratio2;
     if (ratio2 > 0) cerr << " BAD-BOUND";

     cerr << "\n";
  }
#endif

  // Add multiples of p2r and q to make the zzParts divisible by p^{e'}
  assert(zzParts.size() == 2);


  ZZX U_vec[2];
  long maxU=0;
  double maxU_norm = 0;
  for (long i=0; i<(long)zzParts.size(); i++) {
    // make divisible by p^{e'}
    double U_norm;
    long newMax = makeDivisible(zzParts[i].rep, p2ePrime, p2r, q,
				trcData.a, U_norm, context.zMStar);
    zzParts[i].normalize();   // normalize after working directly on the rep
    if (maxU < newMax)  maxU = newMax;
    if (maxU_norm < U_norm)  maxU_norm = U_norm;
  }


  // Check that the estimated noise is still low
  if (noise + maxU_norm*p2r*(skBounds[recryptKeyID]+1) > q/2)
    cerr << " ******** warning: noise/q after makeDivisible = "
	 << ((noise + maxU_norm*p2r*(skBounds[recryptKeyID]+1))/q) << endl;

#ifdef PRINT_LEVELS
   if (dbgKey) {
     const RecryptData& rcData = ctxt.getContext().rcData;
     ZZX ptxt;
     rawDecrypt(ptxt, zzParts, dbgKey->sKeys[recryptKeyID], q);

     Vec<ZZ> powerful;
     ZZX ptxt_alt;
     rcData.p2dConv->ZZXtoPowerful(powerful, ptxt);
     vecRed(powerful, powerful, q, false);
     rcData.p2dConv->powerfulToZZX(ptxt_alt, powerful);


#if 1
     xdouble max_pwrfl = conv<xdouble>(largestCoeff(powerful));
     xdouble max_canon = embeddingLargestCoeff(ptxt_alt,rcData.ea->getPAlgebra());
     double noise_bnd = noise + maxU_norm*p2r*(skBounds[recryptKeyID]+1);

     cerr << "  after makeDivisible";

     cerr << ", noise=" << noise;
     cerr << ", maxU_norm=" << maxU_norm;
     cerr << ", p2r=" << p2r;
     cerr << ", sk_bnd=" << skBounds[recryptKeyID];
     cerr << ", max_canon=" << max_canon;

     cerr << ", max_pwrfl/q=" << (max_pwrfl/q);

     double ratio0 = log(max_canon/noise_bnd)/log(2.0);
     cerr << ", log2(max_canon/bound)=" << ratio0;
     if (ratio0 > 0) cerr << " BAD-BOUND";

     double ratio1 = log(max_pwrfl/max_canon)/log(2.0);
     cerr << ", log2(max_pwrfl/max_canon)=" << ratio1;
     if (ratio1 > 0) cerr << " BAD-BOUND";

     cerr << "\n";
#else
     cerr << "  after makeDivisible, noiseEst="
	  << (noise + maxU*p2r*(skBounds[recryptKeyID]+1))
	  << ", maxCanon="
	  << embeddingLargestCoeff(ptxt,rcData.ea->getPAlgebra());
     cerr << ",  maxCoeff="<<largestCoeff(ptxt)
	  << ", maxPowfl="<<largestCoeff(powerful)
	  << endl;
#endif
  }
#endif

  for (long i=0; i<(long)zzParts.size(); i++)
    zzParts[i] /= p2ePrime;   // divide by p^{e'}

  // Multiply the post-processed cipehrtext by the encrypted sKey
#ifdef DEBUG_PRINTOUT
  cerr << "+ Before recryption ";
  decryptAndPrint(cerr, recryptEkey, *dbgKey, *dbgEa, printFlag);
#endif

  double p0size = to_double(embeddingLargestCoeff(zzParts[0], context.zMStar));
  double p1size = to_double(embeddingLargestCoeff(zzParts[1], context.zMStar));
  // FIXME: This might be slow without Armadillo

  // NOTE: here we lose the intFactor associated with ctxt.
  // We will restore it below.
  ctxt = recryptEkey;


  ctxt.multByConstant(zzParts[1], p1size);
  ctxt.addConstant(zzParts[0], p0size);

#ifdef DEBUG_PRINTOUT
  cerr << "+ Before linearTrans1 ";
  decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, printFlag);
#endif
  FHE_NTIMER_STOP(AAA_bootKeySwitch);

#ifdef PRINT_LEVELS
   CheckCtxt(ctxt, "after bootKeySwitch");
#endif

  // Move the powerful-basis coefficients to the plaintext slots
  FHE_NTIMER_START(AAA_coeffToSlot);
  trcData.coeffToSlot->apply(ctxt);
  FHE_NTIMER_STOP(AAA_coeffToSlot);


#ifdef PRINT_LEVELS
   CheckCtxt(ctxt, "after coeffToSlot");
#endif

#ifdef DEBUG_PRINTOUT
  cerr << "+ After linearTrans1 ";
  decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, printFlag);
#endif

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)
  FHE_NTIMER_START(AAA_extractDigitsThin);
  extractDigitsThin(ctxt, e-ePrime, r, ePrime);
  FHE_NTIMER_STOP(AAA_extractDigitsThin);


#ifdef PRINT_LEVELS
   CheckCtxt(ctxt, "after extractDigitsThin");
#endif

#ifdef DEBUG_PRINTOUT
  cerr << "+ Before linearTrans2 ";
  decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, printFlag);
#endif

  // restore intFactor
  if (intFactor != 1)
    ctxt.intFactor = MulMod(ctxt.intFactor, intFactor, ptxtSpace);
}



