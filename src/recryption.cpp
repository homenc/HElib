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


#ifdef DEBUG_PRINTOUT
#include "debugging.h"
long printFlag = FLAG_PRINT_VEC;
#endif

/************************ Some local functions ***********************/
/*********************************************************************/
static void
printSizesPowerful(const vector<ZZX>& zzParts, const DoubleCRT& sKey,
                   const RecryptData& rcData, long q, double noise);

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

#if 1
static long makeDivisible(vec_ZZ& vec, long p2e, long p2r, long q, long a, 
                          double& U_norm, const PAlgebra& palg, PowerfulDCRT* p2dConv = 0)
{
  //OLD: assert(q>0 && p2e>0 && p2r>0 && a>=0 && q % p2e == 1 && a % p2r == 0 && a*2 < p2e);
  helib::assertTrue<helib::InvalidArgument>(q > 0l, "q must be positive");
  helib::assertTrue<helib::InvalidArgument>(p2e > 0l, "p2e must be positive");
  helib::assertTrue<helib::InvalidArgument>(p2r > 0l, "p2r must be positive");
  helib::assertTrue<helib::InvalidArgument>(a >= 0l, "a must be non-negative");

  helib::assertEq<helib::InvalidArgument>(q % p2e, 1l, "q must equal 1 modulo p2e");
  helib::assertEq<helib::InvalidArgument>(a % p2r, 0l, "p2r must divide a");

  helib::assertTrue<helib::InvalidArgument>(a * 2 < p2e, "a must be less than half of p2e");

  long aa = a / p2r;


#ifdef DEBUG_PRINTOUT
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
      if (u > aa) u = aa;
    }
    else {              // need to add a negative number
      u = -(zMod/p2r);
      if (u < -aa) u = -aa;
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
    double V_norm;
    if (p2dConv) {
       ZZX poly;
       p2dConv->powerfulToZZX(poly, conv<Vec<ZZ>>(uVec));
       U_norm = conv<double>(embeddingLargestCoeff(poly, palg)) *p2r;
       p2dConv->powerfulToZZX(poly, conv<Vec<ZZ>>(vVec));
       V_norm = conv<double>(embeddingLargestCoeff(poly, palg));
    }
    else {
      U_norm = conv<double>(embeddingLargestCoeff(uVec, palg)) *p2r;
      V_norm = conv<double>(embeddingLargestCoeff(vVec, palg));
    }
    cerr << "  makeDivisible: maxU=" << (maxU*p2r)
         << ", U_norm=" << U_norm
         << ", V_norm=" << V_norm << endl;
  }
#endif

  return maxU;
}
#else
// experimental, randomized version
static long makeDivisible(vec_ZZ& vec, long p2e, long p2r, long q, long a, 
                          double& U_norm, const PAlgebra& palg)
{
  //OLD: assert(q>0 && p2e>0 && p2r>0 && a>=0 && q % p2e == 1 && a % p2r == 0 && a*2 < p2e);
  helib::assertTrue<helib::InvalidArgument>(q > 0l, "q must be positive");
  helib::assertTrue<helib::InvalidArgument>(p2e > 0l, "p2e must be positive");
  helib::assertTrue<helib::InvalidArgument>(p2r > 0l, "p2r must be positive");
  helib::assertTrue<helib::InvalidArgument>(a >= 0l, "a must be non-negative");

  helib::assertEq<helib::InvalidArgument>(q % p2e, 1l, "q must equal 1 modulo p2e");
  helib::assertEq<helib::InvalidArgument>(a % p2r, 0, "p2r must divide a");

  helib::assertTrue<helib::InvalidArgument>(a * 2 < p2e, "a must be less than half of p2e");
  long aa = a / p2r;


  vec_ZZ orig_vec = vec;
  double sum_norm;
  long maxU;

  zzX uVec(INIT_SIZE, vec.length());
  zzX vVec(INIT_SIZE, vec.length());

  for (long trials = 0; trials < 10; trials++) {
    vec_ZZ try_vec = orig_vec;

    long maxU1 = 0;
    for (long i=0; i<vec.length(); i++) {
      ZZ& z = try_vec[i];
      long u, v;

      // What to add to z to make it divisible by p2e?
      long zMod = rem(z, p2e); // zMod is in [0,p2e-1]
      if (zMod > p2e/2) { // need to add a positive number
	zMod = p2e - zMod;
	u = zMod/p2r;
	if (u > aa) u = aa;
        if (u > 0) {
          long ran = RandomBnd(u+1); //0..u
          u -= ran;
        }
      }
      else {              // need to add a negative number
	u = -(zMod/p2r);
	if (u < -aa) u = -aa;
        if (-u > 0) {
          long ran = RandomBnd(-u+1); //0..|u|
          u += ran;
        }
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
      if (abs(u) > maxU1) maxU1 = abs(u);

      uVec[i] = u;
      vVec[i] = v;
    }

    double U_norm1 = conv<double>(embeddingLargestCoeff(uVec, palg)) *p2r;
    double V_norm1 = conv<double>(embeddingLargestCoeff(vVec, palg));
    double sum_norm1 = U_norm1 + V_norm1;
#ifdef DEBUG_PRINTOUT
      cerr << "  makeDivisible: maxU1=" << (maxU1*p2r)
	   << ", U_norm1=" << U_norm1
	   << ", V_norm1=" << V_norm1 
           << ", sum=" << (U_norm1+V_norm1) 
           << "\n";
#endif

    if (trials == 0 || sum_norm1 < sum_norm) {
      sum_norm = sum_norm1;
      maxU = maxU1;
      U_norm = U_norm1;
      vec = try_vec;
    }
  }

  return maxU;
}

#endif

static inline double pow(long a, long b) {return pow(double(a), double(b));}
/*********************************************************************/
/*********************************************************************/

RecryptData::~RecryptData()
{
  if (alMod!=NULL)     delete alMod;
  if (ea!=NULL)        delete ea;
  if (firstMap!=NULL)  delete firstMap;
  if (secondMap!=NULL) delete secondMap;
  if (p2dConv!=NULL)   delete p2dConv;
}


/**
 * Summary of Appendix A from https://ia.cr/2014/873 (version from 2019):
 * Assume that we already chosen a, e, e' and t (which induces the
 * secret-key size tau).
 * 
 * Going into the recryption procedure after "raw mod-switching", we
 * have a ciphertext (c0,c1) where the ci's are "random modulo q" in
 * their powerful basis. Denoting x = c0+c1*s (without mod-q reduction),
 * then |x|< |c0|+|c1*s|< q + B*||c1*s||, where |X| is powerful-basis
 * norm, ||X|| is canonical embedding norm, and B is some bound on the
 * ratio between the two.
 * 
 * Let tau be our bound on the secret key canonical-mebedding norm,
 * and we think of c1 as having random coefficients in [+-q/2] (say in
 * the powerful basis). Then ||c1|| < A*q whp (for some other bound A),
 * and therefore ||c1*s||< A*q*tau. Hence we get |x| < q*(1+B*A*tau).
 * The quantity A*B for this ring is recorded as cM in the PAlgebra,
 * so we have |x|/q < 1 + cM*tau < (1+tau)*cM.
 *
 * We also assume that the "noise term" after mod-q reduction is bounded
 * by |[x]_q| < 2*p^r*(1+tau)*cM (this expression is twice the added
 * noise term from mod-switching).
 * 
 * After makeDivisible relative to e' and a (with a divisible by p^r),
 * and b = p^e'/2 -a, we have a ciphertext (c0',c1') s.t.
 *     x' = c0'+c1'*s = x+p^r(u0+u1*s)+(v0+v1*s),
 * where |u0|,|u1|<a and |v0|,|v1|<b. It follows from the above that
 * 
 *   |x'|/q  < (2+b)(1+tau)*cM, and
 *   |[x']_q|< p^r(2+a)(1+tau)*cM
 * 
 * To be able to use the Lemma 5.1 from https://ia.cr/2014/873, we
 * need to have |x'|/q + |[x']_q| <= (q-1)/2 = p^e/2. Using the bounds
 * from above, a sufficient condition for this is
 * 
 *    (1) (p^{e'}/2 + 2(p^r+1))(tau+1)*cM <= (q-1)/2  = p^e/2
 *
 * (This is Equation (9) in Appendix A of https://ia.cr/2014/873,
 * but note that the a here is a*p^r there.)
 * 
 * Note that as we let e,e' tend to infinity the constraint above
 * degenerates to (tau+1)*cM < p^{e-e'}, so the smallest value
 * of e-e' that we can hope for is
 *
 *    (2) e-e' = 1 + floor( log_p( (tau+1)*cM) )
 *
 * The setAE procedure tries to minimize e-e' subject to (1), and
 * in addition subject to the constraint that e is "not too big".
 * Specifically, it tries to ensure p^e<2^{30}, and failing that it
 * uses the smallest e for which 2(p^r+1)(tau+1)*cM*2 <= p^e, and the
 * largest e' for that value of e.
 *
 * Once e,e' are set, it splits p^{e'}/2=a+b with a,b about equal and
 * a divisible by p^r. Then it computes and returns the largest Hamming
 * weight for the key (that implies the norm tau) for which constraint
 * (1) still holds.
 *
 * NOTE: setAE returns the Hamming weight, *not* the norm tau. The norm
 * can be computed from the weight using sampleHWtBoundedEffectiveBound.
 **/
long RecryptData::setAE(long& a, long& e, long& ePrime,
                    const FHEcontext& context, long targetWeight)
{
  bool default_target=false;
  if (targetWeight<=0) {
    targetWeight = RecryptData::defSkHwt;
    default_target=true;
  }

  double skSize = sampleHWtBoundedEffectiveBound(context, targetWeight);

  double magicConst = context.zMStar.get_cM();
  double factor = (skSize+1) * magicConst * 2;
  long p = context.zMStar.getP();
  long p2r = context.alMod.getPPowR();
  long frstTerm = 2*(p2r+1);
  double logp = log(p);    // log p

  // Start with the smallest e s.t. p^e > (2p^r+1)(skSize+1)*magicCons
  ePrime = 0;
  e = ceil( log(frstTerm*factor) /logp );
  //OLD: assert(e*logp < log(NTL_SP_BOUND));
  helib::assertTrue(e*logp < log(NTL_SP_BOUND), "p^e for this smallest e must be single precision");
  // p^e for this smallest e must be single precision

  // Loop to try and find better solutions, that still satisfy
  //          (p^{e'}/2 + 2*p^r+1) * (s+1)*magicConst*2 <= p^e
  //          <=> p^{e'}/2 + 2*p^r +1 <= p^e/factor
  //          <=> p^{e'} <=  2*(p^e/factor - 2*p^r -1)
  // stop when you hit p^e > 2^30
  long eTry = e;
  long eMinusEprime = e; // want to minimize e-e'
  do {
    double p2e = pow(p,eTry);

    // Solve for the largest e' satisfying constraints (1)
    // p^{e'} <=  2*(p^e/factor - 2*p^r -1)
    long ePrimeTry = floor( log(2*(p2e/factor -frstTerm)) / logp);
    if (eTry - ePrimeTry < eMinusEprime) {
      e = eTry;
      ePrime = ePrimeTry;
      eMinusEprime = e - ePrime;
    }
  } while ((++eTry)*log(p) <= log(1L<<30));  

  // Split p^{e'}/2 into a+b with a divisible by p^r
  long pToEprimeOver2 = floor(pow(p,ePrime)/2.0);
  a = 2*pToEprimeOver2 /5; // empirically we want a little below a half
  a -= (a % p2r);
 
  if (default_target) {
    // Try to increase t, while maintaining constraint (1)
    double bound = (frstTerm + pToEprimeOver2)*magicConst * 2;
    long p2e = pow(p,e);
    while (targetWeight++) {
      skSize = sampleHWtBoundedEffectiveBound(context, targetWeight);
      if (targetWeight>256 || bound*(skSize+1) >= p2e)
	break;
    }
    targetWeight--;
  }

#ifdef DEBUG_PRINTOUT
  long b = pToEprimeOver2 - a;
  cerr << "RecryptData::setAE(): e="<<e<<", e'="<<ePrime
       << ", a="<<a<<", b="<<b<<", sk-hwt="<<targetWeight
       << " (size="<<skSize<<"), magicConst="<<magicConst<<endl;
#endif
  return targetWeight;
}


bool RecryptData::operator==(const RecryptData& other) const
{
  if (mvec != other.mvec) return false;
  if (skHwt != other.skHwt) return false;

  return true;
}



// The main method
void RecryptData::init(const FHEcontext& context, const Vec<long>& mvec_,
                  bool enableThick, long t, bool build_cache_, bool minimal)
{
  if (alMod != NULL) { // were we called for a second time?
    cerr << "@Warning: multiple calls to RecryptData::init\n";
    return;
  }
  helib::assertEq(computeProd(mvec_), (long)context.zMStar.getM(), "Cyclotomic polynomial mismatch"); // sanity check

  // Record the arguments to this function
  mvec = mvec_;
  build_cache = build_cache_;

  skHwt = setAE(a, e, ePrime, context, t);
  long p = context.zMStar.getP();
  long r = context.alMod.getR();

  // First part of Bootstrapping works wrt plaintext space p^{r'}
  alMod = new PAlgebraMod(context.zMStar, e-ePrime+r);
  ea = new EncryptedArray(context, *alMod);
         // Polynomial defaults to F0, PAlgebraMod explicitly given

  p2dConv = new PowerfulDCRT(context, mvec);

  if (!enableThick) return;

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

// Extract digits from fully packed slots
void extractDigitsPacked(Ctxt& ctxt, long botHigh, long r, long ePrime,
			 const vector<ZZX>& unpackSlotEncoding);

// Extract digits from unpacked slots
void extractDigitsThin(Ctxt& ctxt, long botHigh, long r, long ePrime);

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

  //OLD: assert(recryptKeyID>=0); // check that we have bootstrapping data
  helib::assertTrue(recryptKeyID>=0l, "No bootstrapping data");

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
  //OLD: assert(e>=r);
  helib::assertTrue(e>=r, "rcData.e must be at least alMod.r");

#ifdef DEBUG_PRINTOUT
  cerr << "reCrypt: p="<<p<<", r="<<r<<", e="<<e<<" ePrime="<<ePrime
       << ", q="<<q<<endl;
  CheckCtxt(ctxt, "init");
#endif

  // can only bootstrap ciphertext with plaintext-space dividing p^r
  //OLD: assert(p2r % ptxtSpace == 0);
  helib::assertEq(p2r % ptxtSpace, 0l, "ptxtSpace must divide p^r when bootstrapping");

  ctxt.dropSmallAndSpecialPrimes();

#ifdef DEBUG_PRINTOUT
  CheckCtxt(ctxt, "after mod down");
#endif


  FHE_NTIMER_START(AAA_preProcess);

  // Make sure that this ciphertxt is in canonical form
  if (!ctxt.inCanonicalForm()) ctxt.reLinearize();

  // Mod-switch down if needed
  IndexSet s = ctxt.getPrimeSet() / context.specialPrimes;
  //OLD: assert(s <= context.ctxtPrimes);
  helib::assertTrue(s <= context.ctxtPrimes, "Not enough room to mod down when bootstrapping");
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

#ifdef DEBUG_PRINTOUT
  if (dbgKey) {
    cerr << "  before makeDivisible (recryption modulus q="<<q
         << "), noise_bnd=" << noise<<endl;
    printSizesPowerful(zzParts, dbgKey->sKeys[recryptKeyID],
                       ctxt.getContext().rcData, q, noise);
  }
#endif

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
#ifdef DEBUG_PRINTOUT
  double newNoise = noise + maxU_norm*p2r*(skBounds[recryptKeyID]+1);
  cerr << "  after makeDivisible, maxU=" << maxU
       << ", maxU_norm="<<maxU_norm<<", p2r="<<p2r
       << ", noise_bnd="<<newNoise<<", sk_bnd="<< skBounds[recryptKeyID]
       << endl;
   if (dbgKey)
     printSizesPowerful(zzParts, dbgKey->sKeys[recryptKeyID],
                        ctxt.getContext().rcData, q, newNoise);
#endif

  for (long i=0; i<(long)zzParts.size(); i++)
    zzParts[i] /= p2ePrime;   // divide by p^{e'}

  // Multiply the post-processed cipehrtext by the encrypted sKey

  double p0size = to_double(embeddingLargestCoeff(zzParts[0], context.zMStar));
  double p1size = to_double(embeddingLargestCoeff(zzParts[1], context.zMStar));
  // FIXME: This might be slow without Armadillo

  // NOTE: here we lose the intFactor associated with ctxt.
  // We will restore it below.
  ctxt = recryptEkey;

  ctxt.multByConstant(zzParts[1]);
  ctxt.addConstant(zzParts[0]);

#ifdef DEBUG_PRINTOUT
  CheckCtxt(ctxt, "after preProcess");
#endif
  FHE_NTIMER_STOP(AAA_preProcess);

  // Move the powerful-basis coefficients to the plaintext slots
  FHE_NTIMER_START(AAA_LinearTransform1);
  ctxt.getContext().rcData.firstMap->apply(ctxt);
  FHE_NTIMER_STOP(AAA_LinearTransform1);

#ifdef DEBUG_PRINTOUT
  CheckCtxt(ctxt, "after LinearTransform1");
#endif

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)
  FHE_NTIMER_START(AAA_extractDigitsPacked);
  extractDigitsPacked(ctxt, e-ePrime, r, ePrime,
		      context.rcData.unpackSlotEncoding);
  FHE_NTIMER_STOP(AAA_extractDigitsPacked);


#ifdef DEBUG_PRINTOUT
  CheckCtxt(ctxt, "after extractDigitsPacked");
#endif

  // Move the slots back to powerful-basis coefficients
  FHE_NTIMER_START(AAA_LinearTransform2);
  ctxt.getContext().rcData.secondMap->apply(ctxt);
  FHE_NTIMER_STOP(AAA_LinearTransform2);


#ifdef DEBUG_PRINTOUT
  CheckCtxt(ctxt, "after linearTransform2");
#endif

  // restore intFactor
  if (intFactor != 1)
    ctxt.intFactor = MulMod(ctxt.intFactor, intFactor, ptxtSpace);
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

  //#ifdef DEBUG_PRINTOUT
  //  CheckCtxt(unpacked[0], "after unpack");
  //#endif

  NTL_EXEC_RANGE(d, first, last)
  for (long i = first; i < last; i++) {
    extractDigitsThin(unpacked[i], botHigh, r, ePrime);
  }
  NTL_EXEC_RANGE_END

  //#ifdef DEBUG_PRINTOUT
  //CheckCtxt(unpacked[0], "before repack");
  //#endif

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
  //#ifdef DEBUG_PRINTOUT
  //CheckCtxt(ctxt, "after repack");
  //#endif
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

  //#ifdef DEBUG_PRINTOUT
  //  CheckCtxt(unpacked[0], "after unpack");
  //#endif

  for (long i=0; i<(long)unpacked.size(); i++) {
    extractDigitsThin(unpacked[i], botHigh, r, ePrime); 
  }

  //#ifdef DEBUG_PRINTOUT
  //  CheckCtxt(unpacked[0], "before repack");
  //#endif

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
    pKey.reCrypt(c);
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
  if (coeffToSlot!=NULL)  delete coeffToSlot;
  if (slotToCoeff!=NULL) delete slotToCoeff;
}


// This code was copied from RecryptData::init, and is mostly
// the same, except for the linear-map-related stuff.
// FIXME: There is really too much code (and data!) duplication here.
void ThinRecryptData::init(const FHEcontext& context, const Vec<long>& mvec_,
                      bool alsoThick, long t, bool build_cache_, bool minimal)
{
  RecryptData::init(context, mvec_, alsoThick, t, build_cache_, minimal);
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

  //OLD: assert(recryptKeyID>=0); // check that we have bootstrapping data
  helib::assertTrue(recryptKeyID>=0l, "Bootstrapping data not present");

  long p = ctxt.getContext().zMStar.getP();
  long r = ctxt.getContext().alMod.getR();
  long p2r = ctxt.getContext().alMod.getPPowR();

  long intFactor = ctxt.intFactor;

  const ThinRecryptData& trcData = ctxt.getContext().rcData;

  // the bootstrapping key is encrypted relative to plaintext space p^{e-e'+r}.
  long e = trcData.e;
  long ePrime = trcData.ePrime;
  long p2ePrime = power_long(p,ePrime);
  long q = power_long(p,e)+1;
  //OLD: assert(e>=r);
  helib::assertTrue(e>=r, "trcData.e must be at least alMod.r");

  // can only bootstrap ciphertext with plaintext-space dividing p^r
  //OLD: assert(p2r % ptxtSpace == 0);
  helib::assertEq(p2r % ptxtSpace, 0l, "ptxtSpace must divide p^r when thin bootstrapping");

#ifdef DEBUG_PRINTOUT
  CheckCtxt(ctxt, "init");
#endif

  ctxt.dropSmallAndSpecialPrimes();

#ifdef DROP_BEFORE_THIN_RECRYPT
  // experimental code...we should drop down to a reasonably low level
  // before doing the first linear map.
  long first = context.ctxtPrimes.first();
  long last = min(context.ctxtPrimes.last(),
                  first + thinRecrypt_initial_level - 1);
  ctxt.bringToSet(IndexSet(first, last));
#endif

#ifdef DEBUG_PRINTOUT
  CheckCtxt(ctxt, "after mod down");
#endif

  // Move the slots to powerful-basis coefficients
  FHE_NTIMER_START(AAA_slotToCoeff);
  trcData.slotToCoeff->apply(ctxt);
  FHE_NTIMER_STOP(AAA_slotToCoeff);

#ifdef DEBUG_PRINTOUT
  CheckCtxt(ctxt, "after slotToCoeff");
#endif

  FHE_NTIMER_START(AAA_bootKeySwitch);

  // Make sure that this ciphertxt is in canonical form
  if (!ctxt.inCanonicalForm()) ctxt.reLinearize();

  // Mod-switch down if needed
  IndexSet s = ctxt.getPrimeSet() / context.specialPrimes;
  //OLD: assert(s <= context.ctxtPrimes);
  helib::assertTrue(s <= context.ctxtPrimes,  "Not enough room to mod down when thin bootstrapping");
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
  //OLD: assert(zzParts.size() == 2);
  helib::assertEq(zzParts.size(), (std::size_t)2, "Exactly 2 parts required for mod-switching in thin bootstrapping");

#ifdef DEBUG_PRINTOUT
  if (dbgKey) {
    cerr << "  before makeDivisible (recryption modulus q="<<q
         << "), noise_bnd=" << noise<<endl;
    printSizesPowerful(zzParts, dbgKey->sKeys[recryptKeyID],
                       ctxt.getContext().rcData, q, noise);
  }
#endif

#if 1

  // Add multiples of p2r and q to make the zzParts divisible by p^{e'}
  double maxU_norm = 0;
  long maxU=0;
  for (long i=0; i<(long)zzParts.size(); i++) {
    // make divisible by p^{e'}
    double U_norm; 

    long newMax = makeDivisible(zzParts[i].rep, p2ePrime, p2r, q,
				trcData.a, U_norm, context.zMStar);

    zzParts[i].normalize(); // normalize after working directly on the rep
    if (maxU < newMax)  maxU = newMax;
    if (maxU_norm < U_norm)  maxU_norm = U_norm;
  }
#ifdef DEBUG_PRINTOUT
  double newNoise = noise + maxU_norm*(skBounds[recryptKeyID]+1);
  cerr << "  after makeDivisible, maxU=" << maxU
       << ", maxU_norm="<<maxU_norm<<", p2r="<<p2r
       << ", noise_bnd="<<newNoise<<", sk_bnd="<< skBounds[recryptKeyID]
       << endl;
   if (dbgKey)
     printSizesPowerful(zzParts, dbgKey->sKeys[recryptKeyID],
                        ctxt.getContext().rcData, q, newNoise);
#endif

  for (long i=0; i<(long)zzParts.size(); i++)
    zzParts[i] /= p2ePrime;   // divide by p^{e'}

#else
  // Experimental version
  // Add multiples of p2r and q to make the zzParts divisible by p^{e'}
  double maxU_norm = 0;
  long maxU=0;
  for (long i=0; i<(long)zzParts.size(); i++) {
    // make divisible by p^{e'}
    double U_norm; 

    Vec<ZZ> pwrfl;
    trcData.p2dConv->ZZXtoPowerful(pwrfl, zzParts[i]);

    long newMax = makeDivisible(pwrfl, p2ePrime, p2r, q,
				trcData.a, U_norm, context.zMStar);

    trcData.p2dConv->powerfulToZZX(zzParts[i], pwrfl);

    zzParts[i].normalize(); // normalize after working directly on the rep
    if (maxU < newMax)  maxU = newMax;
    if (maxU_norm < U_norm)  maxU_norm = U_norm;
  }
#ifdef DEBUG_PRINTOUT
  double newNoise = noise + maxU_norm*(skBounds[recryptKeyID]+1);
  cerr << "  after makeDivisible, maxU=" << maxU
       << ", maxU_norm="<<maxU_norm<<", p2r="<<p2r
       << ", noise_bnd="<<newNoise<<", sk_bnd="<< skBounds[recryptKeyID]
       << endl;
   if (dbgKey)
     printSizesPowerful(zzParts, dbgKey->sKeys[recryptKeyID],
                        ctxt.getContext().rcData, q, newNoise);
#endif

  for (long i=0; i<(long)zzParts.size(); i++)
    zzParts[i] /= p2ePrime;   // divide by p^{e'}

#endif

  // Multiply the post-processed cipehrtext by the encrypted sKey

  double p0size = to_double(embeddingLargestCoeff(zzParts[0], context.zMStar));
  double p1size = to_double(embeddingLargestCoeff(zzParts[1], context.zMStar));
  // FIXME: This might be slow without Armadillo

  // NOTE: here we lose the intFactor associated with ctxt.
  // We will restore it below.
  ctxt = recryptEkey;


  ctxt.multByConstant(zzParts[1], p1size);
  ctxt.addConstant(zzParts[0], p0size);

#ifdef DEBUG_PRINTOUT
   CheckCtxt(ctxt, "after bootKeySwitch");
#endif

  FHE_NTIMER_STOP(AAA_bootKeySwitch);

  // Move the powerful-basis coefficients to the plaintext slots
  FHE_NTIMER_START(AAA_coeffToSlot);
  trcData.coeffToSlot->apply(ctxt);
  FHE_NTIMER_STOP(AAA_coeffToSlot);

#ifdef DEBUG_PRINTOUT
   CheckCtxt(ctxt, "after coeffToSlot");
#endif

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)
  FHE_NTIMER_START(AAA_extractDigitsThin);
  extractDigitsThin(ctxt, e-ePrime, r, ePrime);
  FHE_NTIMER_STOP(AAA_extractDigitsThin);


#ifdef DEBUG_PRINTOUT
   CheckCtxt(ctxt, "after extractDigitsThin");
#endif

  // restore intFactor
  if (intFactor != 1)
    ctxt.intFactor = MulMod(ctxt.intFactor, intFactor, ptxtSpace);
}


static void
printSizesPowerful(const vector<ZZX>& zzParts, const DoubleCRT& sKey,
                   const RecryptData& rcData, long q, double noise)
{
  ZZX ptxt;
  //  const RecryptData& rcData = ctxt.getContext().rcData;
  const PAlgebra& palg = rcData.ea->getPAlgebra();
  rawDecrypt(ptxt, zzParts, sKey); // no mod q

  Vec<ZZ> powerful;
  rcData.p2dConv->ZZXtoPowerful(powerful, ptxt);
  xdouble max_pwrfl = conv<xdouble>(largestCoeff(powerful));
  xdouble max_canon = embeddingLargestCoeff(ptxt, palg);
  double ratio = log(max_pwrfl/max_canon)/log(2.0);
  xdouble critical_value = (max_pwrfl/q)/q;
  cerr << "                     max_pwrfl/q^2=" << ((max_pwrfl/q)/q)
       << ", log2(max_pwrfl/max_canon)=" << ratio;
  if (ratio > 0) cerr << " BAD-BOUND";
  cerr << endl;


  vecRed(powerful, powerful, q, false);
  max_pwrfl = conv<xdouble>(largestCoeff(powerful));
  rcData.p2dConv->powerfulToZZX(ptxt, powerful);
  max_canon = embeddingLargestCoeff(ptxt, palg);
  ratio = log(max_pwrfl/max_canon)/log(2.0);
  critical_value += max_pwrfl/q;
  cerr << "        after mod q, max_pwrfl/q=" << (max_pwrfl/q);
  if (critical_value > 0.5) cerr << " BAD-BOUND";
  cerr << ", log2(max_pwrfl/max_canon)=" << ratio;
  if (ratio > 0) cerr << " BAD-BOUND";


  ratio = log(max_canon/noise)/log(2.0);
  cerr << "\n        log2(max_canon/noiseEst)=" << ratio;
  if (ratio > 0) cerr << " BAD-BOUND";

  ratio = log(noise/q)/log(2.0);
  cerr << ", log2(noiseEst/q)=" << ratio;

  cerr << endl;
}
