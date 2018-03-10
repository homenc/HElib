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


/************* Some local functions *************/
static void x2iInSlots(ZZX& poly, long i,
		       vector<ZZX>& xVec, const EncryptedArray& ea);

// Make every entry of vec divisible by p^e by adding/subtracting multiples
// of p^r and q, while keeping the added multiples small. 
template<class VecInt>
long makeDivisible(VecInt& vec, long p2e, long p2r, long q, double alpha);
static inline double pow(long a, long b) {return pow(double(a), double(b));}

RecryptData::~RecryptData()
{
  if (alMod!=NULL)     delete alMod;
  if (ea!=NULL)        delete ea;
  if (firstMap!=NULL)  delete firstMap;
  if (secondMap!=NULL) delete secondMap;
  if (p2dConv!=NULL)   delete p2dConv;
}


/** Computing the recryption parameters
 *
 * To get the smallest possible value of e-e', the params need to satisfy:
 *  (p^e +1)/4 =>
 *       max { (t+1)( 1+ (alpha/2)*(p^e/p^{ceil(log_p(t+2))}) ) + noise      }
 *           { (t+1)( 1+ ((1-alpha)/2)*(p^e/p^{ceil(log_p(t+2))}) +p^r/2) +1 },
 *
 * where noise is taken to be twice the mod-switching additive term, namely
 * noise = p^r *sqrt((t+1)*phi(m)/3). Denoting rho=(t+1)/p^{ceil(log_p(t+2))}
 * (and ignoring fome +1 terms), this is equivalent to:
 *
 *   p^e > max { 4(t+noise)/(1-2*alpha*rho), 2(t+1)p^r/(1-2(1-alpha)rho) }.
 *
 * We first compute the optimal value for alpha (which must be in [0,1]),
 * that makes the two terms in the max{...} as close as possible, and
 * then compute the smallest value of e satisfying this constraint.
 *
 * If this value is too big then we try again with e-e' one larger,
 * which means that rho is a factor of p smaller.
 */

// Some convenience functions
static double lowerBound1(long p, long r, long ePrime, long t,
			  double alpha, double noise)
{
  return (t+1)*(1+ alpha*pow(p,r+ePrime-1)/2)+noise;
}
static double lowerBound2(long p, long r, long ePrime, long t, double alpha)
{
  return (t+1)*(1+ (1-alpha)*pow(p,r+ePrime-1)/2 + pow(p,r)/2)+1;
}

static void setAlphaE(double& alpha, long& e, double rho, double gamma,
		      double noise, double logp, long p2r, long t)
{
  alpha = (1 +gamma*(2*rho-1))/(2*rho*(1+gamma));
  if (alpha<0) alpha=0;
  else if (alpha>1) alpha=1;

  if (alpha<1) {
    double ratio = 4*(t+noise)/(1-2*alpha*rho);
    e = floor(1+ log(ratio)/logp);
  }
  else
    e = floor(1+ log(2*(t+1)*p2r)/logp);
}

bool RecryptData::operator==(const RecryptData& other) const
{
  if (mvec != other.mvec) return false;
  if (hwt != other.hwt) return false;
  if (conservative != other.conservative) return false;

  return true;
}



// The main method
void RecryptData::init(const FHEcontext& context, const Vec<long>& mvec_,
		       long t, bool consFlag, bool build_cache_, bool minimal)
{
  if (alMod != NULL) { // were we called for a second time?
    cerr << "@Warning: multiple calls to RecryptData::init\n";
    return;
  }
  assert(computeProd(mvec_) == (long)context.zMStar.getM()); // sanity check

  // Record the arguments to this function
  mvec = mvec_;
  conservative = consFlag;
  build_cache = build_cache_;

  if (t <= 0) t = defSkHwt+1; // recryption key Hwt
  hwt = t;
  long p = context.zMStar.getP();
  long phim = context.zMStar.getPhiM();
  long r = context.alMod.getR();
  long p2r = context.alMod.getPPowR();
  double logp = log((double)p);

  double noise = p2r * sqrt((t+1)*phim/3.0);
  double gamma = 2*(t+noise)/((t+1)*p2r); // ratio between numerators

  long logT = ceil(log((double)(t+2))/logp); // ceil(log_p(t+2))
  double rho = (t+1)/pow(p,logT);

  if (!conservative) {   // try alpha, e with this "aggresive" setting
    setAlphaE(alpha, e, rho, gamma, noise, logp, p2r, t);
    ePrime = e -r +1 -logT;

    // If e is too large, try again with rho/p instead of rho
    long bound = (1L << (context.bitsPerLevel-1)); // halfSizePrime/2
    if (pow(p,e) > bound) { // try the conservative setting instead
      cerr << "* p^e="<<pow(p,e)<<" is too big (bound="<<bound<<")\n";
      conservative = true;
    }
  }
  if (conservative) { // set alpha, e with a "conservative" rho/p
    setAlphaE(alpha, e, rho/p, gamma, noise, logp, p2r, t);
    ePrime = e -r -logT;
  }

  // Compute highest key-Hamming-weight that still works (not more than 256)
  double qOver4 = (pow(p,e)+1)/4;
  for (t-=10; qOver4>=lowerBound2(p,r,ePrime,t,alpha)
	 &&  qOver4>=lowerBound1(p,r,ePrime,t,alpha,noise) && t<257; t++);
  skHwt = t-1;

  // First part of Bootstrapping works wrt plaintext space p^{r'}
  alMod = new PAlgebraMod(context.zMStar, e-ePrime+r);
  ea = new EncryptedArray(context, *alMod);
         // Polynomial defaults to F0, PAlgebraMod explicitly given


  p2dConv = new PowerfulDCRT(context, mvec);

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

  // the bootstrapping key is encrypted relative to plaintext space p^{e-e'+r}.
  long e = getContext().rcData.e;
  long ePrime = getContext().rcData.ePrime;
  long p2ePrime = power_long(p,ePrime);
  long q = power_long(p,e)+1;
  assert(e>=r);

#ifdef DEBUG_PRINTOUT
  cerr << "reCrypt: p="<<p<<", r="<<r<<", e="<<e<<" ePrime="<<ePrime
       << ", q="<<q<<endl;
#endif

  // can only bootstrap ciphertext with plaintext-space dividing p^r
  assert(p2r % ptxtSpace == 0);

  FHE_NTIMER_START(preProcess);

  // Make sure that this ciphertxt is in canonical form
  if (!ctxt.inCanonicalForm()) ctxt.reLinearize();

  // Mod-switch down if needed
  IndexSet s = ctxt.getPrimeSet() / getContext().specialPrimes; // set minus
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
  noise = sqrt(noise);

  // Add multiples of p2r and q to make the zzParts divisible by p^{e'}
  long maxU=0;
  for (long i=0; i<(long)zzParts.size(); i++) {
    // make divisible by p^{e'}
    long newMax = makeDivisible(zzParts[i].rep, p2ePrime, p2r, q,
				getContext().rcData.alpha);
    zzParts[i].normalize();   // normalize after working directly on the rep
    if (maxU < newMax)  maxU = newMax;
  }

  // Check that the estimated noise is still low
  if (noise + maxU*p2r*(skHwts[recryptKeyID]+1) > q/2) 
    cerr << " * noise/q after makeDivisible = "
	 << ((noise + maxU*p2r*(skHwts[recryptKeyID]+1))/q) << endl;

  for (long i=0; i<(long)zzParts.size(); i++)
    zzParts[i] /= p2ePrime;   // divide by p^{e'}

  // Multiply the post-processed cipehrtext by the encrypted sKey
#ifdef DEBUG_PRINTOUT
  cerr << "+ Before recryption ";
  decryptAndPrint(cerr, recryptEkey, *dbgKey, *dbgEa, printFlag);
#endif

  double p0size = to_double(coeffsL2Norm(zzParts[0]));
  double p1size = to_double(coeffsL2Norm(zzParts[1]));
  ctxt = recryptEkey;
  ctxt.multByConstant(zzParts[1], p1size*p1size);
  ctxt.addConstant(zzParts[0], p0size*p0size);

#ifdef DEBUG_PRINTOUT
  cerr << "+ Before linearTrans1 ";
  decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, printFlag);
#endif
  FHE_NTIMER_STOP(preProcess);

  // Move the powerful-basis coefficients to the plaintext slots
  FHE_NTIMER_START(LinearTransform1);
  ctxt.getContext().rcData.firstMap->apply(ctxt);
  FHE_NTIMER_STOP(LinearTransform1);

#ifdef DEBUG_PRINTOUT
  cerr << "+ After linearTrans1 ";
  decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, printFlag);
#endif

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)
  extractDigitsPacked(ctxt, e-ePrime, r, ePrime,
		      context.rcData.unpackSlotEncoding);

#ifdef DEBUG_PRINTOUT
  cerr << "+ Before linearTrans2 ";
  decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, printFlag);
#endif

  // Move the slots back to powerful-basis coefficients
  FHE_NTIMER_START(LinearTransform2);
  ctxt.getContext().rcData.secondMap->apply(ctxt);
  FHE_NTIMER_STOP(LinearTransform2);
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

// Make every entry of vec divisible by p^e by adding/subtracting multiples
// of p^r and q, while keeping the added multiples small. Specifically, for
// e>=2, an integer z can be made divisible by p^e via
//   z' = z + u * p^r + v*p^r * q,
// with
//   |u|<=ceil(alpha p^{r(e-1)}/2) and |v|<=0.5+floor(beta p^{r(e-1)}/2),
// for any alpha+beta=1. We assume that r<e and that q-1 is divisible by p^e.
// Returns the largest absolute values of the u's and the new entries.
//
// This code is more general than we need, for bootstrapping we will always
// have e>r.
template<class VecInt>
long makeDivisible(VecInt& vec, long p2e, long p2r, long q, double alpha)
{
  assert(((p2e % p2r == 0) && (q % p2e == 1)) ||
	 ((p2r % p2e == 0) && (q % p2r == 1)));

  long maxU =0;
  ZZ maxZ;
  for (long i=0; i<vec.length(); i++) {
    ZZ z, z2; conv(z, vec[i]);
    long u=0, v=0;

    long zMod1=0, zMod2=0;
    if (p2r < p2e && alpha>0) {
      zMod1 = rem(z,p2r);
      if (zMod1 > p2r/2) zMod1 -= p2r; // map to the symmetric interval

      // make z divisible by p^r by adding a multiple of q
      z2 = z - to_ZZ(zMod1)*q;
      zMod2 = rem(z2,p2e); // z mod p^e, still divisible by p^r
      if (zMod2 > p2e/2) zMod2 -= p2e; // map to the symmetric interval
      zMod2 /= -p2r; // now z+ p^r*zMod2=0 (mod p^e) and |zMod2|<=p^{r(e-1)}/2

      u = ceil(alpha * zMod2);
      v = zMod2 - u; // = floor((1-alpha) * zMod2)
      z = z2 + u*p2r + to_ZZ(q)*v*p2r;
    }
    else { // r >= e or alpha==0, use only mulitples of q
      zMod1 = rem(z,p2e);
      if (zMod1 > p2e/2) zMod1 -= p2e; // map to the symmetric interval
      z -= to_ZZ(zMod1) * q;
    }
    if (abs(u) > maxU) maxU = abs(u);
    if (abs(z) > maxZ) maxZ = abs(z);

    if (rem(z,p2e) != 0) { // sanity check
      cerr << "**error: original z["<<i<<"]=" << vec[i]
	   << std::dec << ", p^r="<<p2r << ", p^e="<<p2e << endl;
      cerr << "z' = z - "<<zMod1<<"*q = "<< z2<<endl;
      cerr << "z''=z' +" <<u<<"*p^r +"<<v<<"*p^r*q = "<<z<<endl;
      exit(1);
    }
    conv(vec[i], z); // convert back to native format
  }
  return maxU;
}
// explicit instantiation for vec_ZZ and vec_long
template long makeDivisible<vec_ZZ>(vec_ZZ& v, long p2e,
				    long p2r, long q, double alpha);
// template long makeDivisible<vec_long>(vec_long& v, long p2e,
// 				      long p2r, long q, double alpha);

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

  // Step 2: extract the digits top-1,...,0 from the slots of unpacked[i]
  long p = ctxt.getContext().zMStar.getP();
  long p2r = power_long(p,r);
  long topHigh = botHigh + r-1;

#ifdef DEBUG_PRINTOUT
  cerr << "+ After unpack ";
  decryptAndPrint(cerr, unpacked[0], *dbgKey, *dbgEa, printFlag);
  cerr << "    extracting "<<(topHigh+1)<<" digits\n";
#endif

  if (p==2 && r>2)
    topHigh--; // For p==2 we sometime get a bit for free

  FHE_NTIMER_START(extractDigits);

  NTL_EXEC_RANGE(d, first, last)
      for (long i = first; i < last; i++) {
        vector<Ctxt> scratch;
    
        if (topHigh<=0) { // extracting LSB = no-op
          scratch.assign(1,unpacked[i]);
        } else {          // extract digits topHigh...0, store them in scratch
          extractDigits(scratch, unpacked[i], topHigh+1);
        }

        // set upacked[i] = -\sum_{j=botHigh}^{topHigh} scratch[j] * p^{j-botHigh}
        if (topHigh >= (long)scratch.size()) {
          topHigh = scratch.size() -1;
          cerr << " @ suspect: not enough digits in extractDigitsPacked\n";
        }
    
        unpacked[i] = scratch[topHigh];
        for (long j=topHigh-1; j>=botHigh; --j) {
          unpacked[i].multByP();
          unpacked[i] += scratch[j];
        }
        if (p==2 && botHigh>0) {   // For p==2, subtract also the previous bit
          //cerr << scratch.size() << " " <<  botHigh-1 << "\n";
          unpacked.at(i) += scratch.at(botHigh-1);
        }
        unpacked[i].negate();
    
        if (r>ePrime) {          // Add in digits from the bottom part, if any
          long topLow = r-1 - ePrime;
          Ctxt tmp = scratch[topLow];
          for (long j=topLow-1; j>=0; --j) {
    	tmp.multByP();
    	tmp += scratch[j];
          }
          if (ePrime>0)
    	tmp.multByP(ePrime); // multiply by p^e'
          unpacked[i] += tmp;
        }
        unpacked[i].reducePtxtSpace(p2r); // Our plaintext space is now mod p^r
      }
  NTL_EXEC_RANGE_END

  FHE_NTIMER_STOP(extractDigits);

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

  vector<Ctxt> scratch; // used below 
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

  // Step 2: extract the digits top-1,...,0 from the slots of unpacked[i]
  long p = ctxt.getContext().zMStar.getP();
  long p2r = power_long(p,r);
  long topHigh = botHigh + r-1;

#ifdef DEBUG_PRINTOUT
  cerr << "+ After unpack ";
  decryptAndPrint(cerr, unpacked[0], *dbgKey, *dbgEa, printFlag);
  cerr << "    extracting "<<(topHigh+1)<<" digits\n";
#endif

  if (p==2 && r>2)
    topHigh--; // For p==2 we sometime get a bit for free

  FHE_NTIMER_START(extractDigits);
  for (long i=0; i<(long)unpacked.size(); i++) {
    if (topHigh<=0) { // extracting LSB = no-op
      scratch.assign(1,unpacked[i]);
    } else {          // extract digits topHigh...0, store them in scratch
      extractDigits(scratch, unpacked[i], topHigh+1);
    }

    // set upacked[i] = -\sum_{j=botHigh}^{topHigh} scratch[j] * p^{j-botHigh}
    if (topHigh >= (long)scratch.size()) {
      topHigh = scratch.size() -1;
      cerr << " @ suspect: not enough digits in extractDigitsPacked\n";
    }

    unpacked[i] = scratch[topHigh];
    for (long j=topHigh-1; j>=botHigh; --j) {
      unpacked[i].multByP();
      unpacked[i] += scratch[j];
    }
    if (p==2 && botHigh>0)   // For p==2, subtract also the previous bit
      unpacked[i] += scratch[botHigh-1];
    unpacked[i].negate();

    if (r>ePrime) {          // Add in digits from the bottom part, if any
      long topLow = r-1 - ePrime;
      Ctxt tmp = scratch[topLow];
      for (long j=topLow-1; j>=0; --j) {
	tmp.multByP();
	tmp += scratch[j];
      }
      if (ePrime>0)
	tmp.multByP(ePrime); // multiply by p^e'
      unpacked[i] += tmp;
    }
    unpacked[i].reducePtxtSpace(p2r); // Our plaintext space is now mod p^r
  }
  FHE_NTIMER_STOP(extractDigits);

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
         && array[i]->findBaseLevel()<belowLvl )
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
           && m[i][j]->findBaseLevel()<belowLvl )
        v.push_back(m[i][j]);
  packedRecrypt(CtPtrs_vectorPt(v), unpackConsts, ea);
}
