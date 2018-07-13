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
#include "binio.h"
#include "timing.h"
#include "FHEContext.h"
#include "Ctxt.h"
#include "FHE.h"

#ifdef DEBUG_PRINTOUT
#include "debugging.h"
#endif

NTL_CLIENT;

extern int fhe_watcher;

NTL_CLIENT

void SKHandle::read(istream& str)
{
  powerOfS = read_raw_int(str);
  powerOfX = read_raw_int(str);
  secretKeyID = read_raw_int(str);
}
 
void SKHandle::write(ostream& str) const
{
  write_raw_int(str, powerOfS);
  write_raw_int(str, powerOfX);
  write_raw_int(str, secretKeyID);
}

// A hack for recording required automorphisms (see NumbTh.h)
std::set<long>* FHEglobals::automorphVals = NULL;
std::set<long>* FHEglobals::automorphVals2 = NULL;

// Dummy encryption, this procedure just encodes the plaintext in a Ctxt object
// NOTE: for now, it leaves the intFactor field of *this alone.
// This assumption is relied upon in the reCrypt() and thinReCrypt()
// routines in recryption.cpp.

void Ctxt::DummyEncrypt(const ZZX& ptxt, double size)
{
  if (isCKKS()) {
    ptxtSpace=1;
    long factor = getContext().alMod.getCx().encodeScalingFactor();
    if (size < 0.0)
      size = getContext().zMStar.getPhiM();
    noiseVar = size * fsquare(factor);
  } else { // BGV
    if (size < 0.0)
      size = getContext().zMStar.getPhiM() * fsquare(ptxtSpace/2.0);
    noiseVar = size;
  }
  primeSet = context.ctxtPrimes;
  highWaterMark = findBaseLevel();

  // A single part, with the plaintext as data and handle pointing to 1

  long f = isCKKS()?
    1 : rem(context.productOfPrimes(context.ctxtPrimes),ptxtSpace);
  if (f == 1) {
    DoubleCRT dcrt(ptxt, context, primeSet);  
    parts.assign(1, CtxtPart(dcrt));
  } else {
    ZZX tmp;
    MulMod(tmp, ptxt, f, ptxtSpace, /*positive=*/false);
    DoubleCRT dcrt(tmp, context, primeSet);  
    parts.assign(1, CtxtPart(dcrt));
  }
}


// Sanity-check: Check that prime-set is "valid", i.e. that it 
// contains either all the special primes or none of them
bool Ctxt::verifyPrimeSet() const
{
  IndexSet s = primeSet & context.specialPrimes; // special primes in primeSet
  if (!empty(s) && s!=context.specialPrimes) return false;

  s = primeSet / s;                              // ctxt primes in primeSet
  return (s.isInterval() && s.first()<=1 && !empty(s));
}


// Compute the number of digits that we need and the esitmated
// added noise from switching this ciphertext part.
static std::pair<long, NTL::xdouble>
keySwitchNoise(const CtxtPart& p, const FHEPubKey& pubKey, long pSpace)
{
  const FHEcontext& context = p.getContext();
  const PAlgebra& palg = context.zMStar;
  double ksSize2 = to_double(context.stdev) * pSpace;
  if (palg.getPow2() > 0) // power of two
    ksSize2 *= ksSize2 * palg.getPhiM() / 2.0;
  else                    // not power of two
    ksSize2 *= ksSize2 * fsquare(palg.getM()) / 2.0;
  // FIXME: Can we instead use the KeySwitch::noiseVar value?

  long nDigits = 0;
  xdouble addedNoise = to_xdouble(0.0);
  double sizeLeft = context.logOfProduct(p.getIndexSet());
  for (size_t i=0; i<context.digits.size() && sizeLeft>0.0; i++) {    
    nDigits++;

    double digitSize = context.logOfProduct(context.digits[i]);
    if (sizeLeft<digitSize) digitSize=sizeLeft;// need only part of this digit

    // Added noise due to this digit is keySwMatrixNoise^2 * |Di|^2, 
    // where |Di| is the magnitude of the digit
    addedNoise += ksSize2 * xexp(2*digitSize);

    sizeLeft -= digitSize;
  }

  // Sanity-check: make sure that the added noise is not more than the
  // special primes can handle: After dividing the added noise by the
  // product of all the special primes, it should be smaller than the
  // added noise term due to modulus switching, i.e.,
  // keySize * phi(m) * pSpace^2 / 4

  double phim = palg.getPhiM();
  double keySize = pubKey.getSKeySize(p.skHandle.getSecretKeyID());
  double logModSwitchNoise = log(keySize) 
    +2*log((double)pSpace) +log(phim) -log(4.0);
  double logKeySwitchNoise = log(addedNoise) 
    -2*context.logOfProduct(context.specialPrimes);
  assert(logKeySwitchNoise < logModSwitchNoise);

  return std::pair<long, NTL::xdouble>(nDigits,addedNoise);
}
std::pair<long, NTL::xdouble> Ctxt::computeKSNoise(long partIdx, long pSpace)
{
  if (pSpace<=1) pSpace = ptxtSpace;
  return keySwitchNoise(parts.at(partIdx), pubKey, pSpace);
}

// Multiply vector of digits by key-switching matrix and add to *this.
// It is assumed that W has at least as many b[i]'s as there are digits.
// The vector of digits is modified in place.
void Ctxt::keySwitchDigits(const KeySwitch& W, vector<DoubleCRT>& digits)
{  // An object to hold the pseudorandom ai's, note that it must be defined
  // with the maximum number of levels, else the PRG will go out of synch.
  // FIXME: This is a bug waiting to happen.
  DoubleCRT ai(context);

  // Subsequent ai's use the evolving RNG state
  RandomState state; // backup the NTL PRG seed
  NTL::SetSeed(W.prgSeed);

  // Add the columns in, one by one
  DoubleCRT tmpDCRT(context, IndexSet::emptySet());  
  for (size_t i=0; i<digits.size(); i++) {
    FHE_NTIMER_START(KS_loop);
    ai.randomize();
    tmpDCRT = digits[i];
  
    // The operations below all use the IndexSet of tmpDCRT
  
    // add digit*a[i] with a handle pointing to base of W.toKeyID
    {FHE_NTIMER_START(KS_loop_1);
     tmpDCRT.Mul(ai,  /*matchIndexSet=*/false);
    }
    {FHE_NTIMER_START(KS_loop_2);
     this->addPart(tmpDCRT, SKHandle(1,1,W.toKeyID), /*matchPrimeSet=*/true);
    }
    // add digit*b[i] with a handle pointing to one
    {FHE_NTIMER_START(KS_loop_3);
     digits[i].Mul(W.b[i], /*matchIndexSet=*/false);
    }
    {FHE_NTIMER_START(KS_loop_4);
     this->addPart(digits[i], SKHandle(), /*matchPrimeSet=*/true);
    }
  }
} // restore random state upon destruction of the RandomState, see NumbTh.h


//! @brief How many levels in the "base-set" for that ciphertext
long Ctxt::findBaseLevel() const 
{
  IndexSet s;
  findBaseSet(s);
  if (context.containsSmallPrime()) {
    if (s.contains(context.ctxtPrimes.first()))
      return 2*card(s) -1; // 1st prime is half size
    else
      return 2*card(s);
  }
  else return card(s);     // one prime per level
}

bool CtxtPart::operator==(const CtxtPart& other) const
{
  if (((DoubleCRT&)*this)!=((DoubleCRT&)other)) return false;
  
  return (skHandle==other.skHandle);
}

// Checking equality between ciphertexts. This routine performs a
// "shallow" check, comparing only pointers to ciphertext parts.
bool Ctxt::equalsTo(const Ctxt& other, bool comparePkeys) const
{
  if (&context != &other.context) return false;
  if (comparePkeys && &pubKey != &other.pubKey) return false;

  if (parts.size() != other.parts.size()) return false;
  for (size_t i=0; i<parts.size(); i++)
    if (parts[i] != other.parts[i]) return false;

  if (primeSet != other.primeSet) return false;
  if (ptxtSpace != other.ptxtSpace) return false;

  // compare noiseVar, ignoring small deviations
  if (noiseVar == 0.0) return (other.noiseVar == 0.0);
  xdouble ratio = other.noiseVar / noiseVar;
  return (ratio>0.9 && ratio<1.1);
}

// Constructor
Ctxt::Ctxt(const FHEPubKey& newPubKey, long newPtxtSpace):
  context(newPubKey.getContext()), pubKey(newPubKey), ptxtSpace(newPtxtSpace),
  noiseVar(to_xdouble(0.0))
{
  if (ptxtSpace<2) ptxtSpace = pubKey.getPtxtSpace();
  else assert (GCD(ptxtSpace, pubKey.getPtxtSpace()) > 1); // sanity check
  primeSet=context.ctxtPrimes;
  highWaterMark = findBaseLevel();
  intFactor = 1;
  ratFactor = 1.0;
}

// Constructor
Ctxt::Ctxt(ZeroCtxtLike_type, const Ctxt& ctxt):
  context(ctxt.getPubKey().getContext()), pubKey(ctxt.getPubKey()), 
  ptxtSpace(ctxt.getPtxtSpace()),
  noiseVar(to_xdouble(0.0))
{
  // same body as previous constructor
  if (ptxtSpace<2) ptxtSpace = pubKey.getPtxtSpace();
  else assert (GCD(ptxtSpace, pubKey.getPtxtSpace()) > 1); // sanity check
  primeSet=context.ctxtPrimes;
  highWaterMark = findBaseLevel();
  intFactor = 1;
  ratFactor = 1.0;
}


// A private assignment method that does not check equality of context or
// public key, this needed for example when we copy the pubEncrKey member
// between different public keys.
Ctxt& Ctxt::privateAssign(const Ctxt& other)
{
  FHE_TIMER_START;
  if (this == &other) return *this; // both point to the same object

  parts = other.parts;
  primeSet = other.primeSet;
  ptxtSpace = other.ptxtSpace;
  noiseVar  = other.noiseVar;
  highWaterMark = other.highWaterMark;
  intFactor = other.intFactor;
  ratFactor = other.ratFactor;
  return *this;
}

// explicitly multiply intFactor by e
void Ctxt::mulIntFactor(long e)
{
  if (e==1) return; // nothing to do
  intFactor = MulMod(intFactor, e, ptxtSpace);
  long bal_e = balRem(e, ptxtSpace);
  for (auto& part : parts) part *= bal_e;
  noiseVar *= fsquare(bal_e); // because every product was scaled by bal_e
}

// Ciphertext maintenance

// mod-switch up to add the primes in s \setminus primeSet, after this call we
// have s<=primeSet. s must contain either all special primes or none of them.
void Ctxt::modUpToSet(const IndexSet &s)
{
  IndexSet setDiff = s/primeSet; // set minus (primes in s but not in primeSet)
  if (empty(setDiff)) return;    // nothing to do, no primes are added

  // scale up all the parts to use also the primes in setDiff
  double f = 0.0;
  for (long i=0; i<lsize(parts); i++) {
    // addPrimesAndScale returns the log of the product of added primes,
    // all calls return the same value = log(prod. of primes in setDiff)
    f = parts[i].addPrimesAndScale(setDiff);
  }

  // The variance estimate grows by a factor of exp(f)^2 = exp(2f)
  noiseVar *= xexp(2*f);

  // If CKKS, the rational factor grows by a factor of exp(f)
  ratFactor *= xexp(f);

  primeSet.insert(setDiff); // add setDiff to primeSet
  assert(verifyPrimeSet()); // sanity-check: ensure primeSet is still valid
}

// mod-switch down to primeSet \intersect s, after this call we have
// primeSet<=s. s must contain either all special primes or none of them.
void Ctxt::modDownToSet(const IndexSet &s)
{
  FHE_TIMER_START;
  IndexSet intersection = primeSet & s;
  if (empty(intersection)) {
    cerr << "modDownToSet called from "<<primeSet<<" to "<<s<<endl;
    exit(1);
  }
  IndexSet setDiff = primeSet / intersection; // set-minus
  if (empty(setDiff)) return;    // nothing to do, removing no primes

  // Scale down all the parts: use either a simple "drop down" (just removing
  // primes, i.e., reducing the ctxt modulo the samaller modulus), or a "real
  // modulus switching" with rounding, basically whichever yeilds smaller
  // noise. Recall that we keep the invariant that a ciphertext mod Q is
  // decrypted to Q*m (mod p), so if we just "drop down" we still need to
  // multiply by (Q^{-1} mod p).

  // Get an estimate for the added noise term for modulus switching
  xdouble addedNoiseVar = modSwitchAddedNoiseVar();
  if (noiseVar*fsquare(ptxtSpace) < addedNoiseVar) { // just "drop down"
    for (size_t i=0; i<parts.size(); i++)
      parts[i].removePrimes(setDiff);       // remove the primes not in s
    long prodInv = 1;
    if (ptxtSpace>1)
      prodInv = InvMod(rem(context.productOfPrimes(setDiff),ptxtSpace), ptxtSpace);
    if (prodInv > 1) {
      for (size_t i=0; i<parts.size(); i++)
        parts[i] *= prodInv;
      noiseVar = noiseVar*fsquare(prodInv);
    }
    // cerr << "DEGENERATE DROP\n";
  } 
  else {                                       // do real mod switching
    for (size_t i=0; i<parts.size(); i++) 
      parts[i].scaleDownToSet(intersection, ptxtSpace);

    // update the noise estimate
    double f = context.logOfProduct(setDiff);
    noiseVar /= xexp(2*f);
    noiseVar += addedNoiseVar;
    ratFactor /= xexp(f); // The factor in CKKS encryption
  }
  primeSet.remove(setDiff); // remove the primes not in s
  assert(verifyPrimeSet()); // sanity-check: ensure primeSet is still valid
  FHE_TIMER_STOP;
}

void Ctxt::bringToLevel(long lvl)
{
  IndexSet target = (context.containsSmallPrime())?
    ((lvl&1)? IndexSet(0,(lvl-1)/2) : IndexSet(1,(lvl-1)/2))
    : IndexSet(0,lvl-1);
  assert(empty(target &context.specialPrimes));// no special primes in target

  modUpToSet(target);   // add any missing primes from target
  modDownToSet(target); // remove any primes not in target
}

// Modulus-switching down
void Ctxt::modDownToLevel(long lvl)
{
  long currentLvl;
  IndexSet targetSet;
  IndexSet currentSet = primeSet & context.ctxtPrimes;
  if (context.containsSmallPrime()) {
    currentLvl = 2*card(currentSet);
    if (currentSet.contains(0))
      currentLvl--;  // first prime is half the size

    if (lvl & 1) {   // odd level, includes the half-size prime
      targetSet = IndexSet(0,(lvl-1)/2);
    } else {
      targetSet = IndexSet(1,lvl/2);
    }
  }
  else {
    currentLvl = card(currentSet);
    targetSet = IndexSet(0,lvl-1);    // one prime per level
  }

  // If target is not below the current level, nothing to do
  if (lvl >= currentLvl && currentSet==primeSet) return;

  if (lvl >= currentLvl) { // just remove the special primes
    targetSet = currentSet;
  }

  // sanity-check: interval does not contain special primes
  assert(targetSet.disjointFrom(context.specialPrimes));

  // may need to mod-UP to include the smallest prime
  if (targetSet.contains(0) && !currentSet.contains(0))
    modUpToSet(targetSet); // adds the primes in targetSet / primeSet

  modDownToSet(targetSet); // removes the primes in primeSet / targetSet
}

void Ctxt::blindCtxt(const ZZX& poly)
{
  Ctxt tmp(pubKey);
  pubKey.Encrypt(tmp,poly,ptxtSpace,/*highNoise=*/true);
  *this += tmp;
  // FIXME: This implementation is not optimized, the levels in the
  //  modulus chain should be handled much better.
  // FIXME: Need to blind the intFactor too
}

// Reduce plaintext space to a divisor of the original plaintext space
void Ctxt::reducePtxtSpace(long newPtxtSpace)
{
  long g = GCD(ptxtSpace, newPtxtSpace);
  assert (g>1); // NOTE: Will trigger if called for CKKS ciphertext
  ptxtSpace = g;
  intFactor %= g;
}


// key-switch to (1,s_i), s_i is the base key with index keyID. If
// keyID<0 then re-linearize to any key for which a switching matrix exists
void Ctxt::reLinearize(long keyID)
{
  FHE_TIMER_START;
  // Special case: if *this is empty or already re-linearized then do nothing
  if (this->isEmpty() || this->inCanonicalForm(keyID)) return;
  // this->reduce();

  // To relinearize, the primeSet must be disjoint from the special primes
  if (!primeSet.disjointFrom(context.specialPrimes))
    modDownToSet(primeSet / context.specialPrimes);

  long g = ptxtSpace;
  Ctxt tmp(pubKey, ptxtSpace); // an empty ciphertext, same plaintext space
  tmp.intFactor = intFactor;   // same intFactor, too

  double logProd = context.logOfProduct(context.specialPrimes);
  tmp.noiseVar = noiseVar * xexp(2*logProd);  // The noise after mod-UP
  tmp.ratFactor = ratFactor * xexp(logProd);// CKKS factor after mod-up

  for (CtxtPart& part : parts) {
    // For a part relative to 1 or base,  only scale and add
    if (part.skHandle.isOne() || part.skHandle.isBase(keyID)) {
      part.addPrimesAndScale(context.specialPrimes);
      tmp.addPart(part, /*matchPrimeSet=*/true);
      continue;
    }
    // Look for a key-switching matrix to re-linearize this part
    const KeySwitch& W = (keyID>=0)? 
      pubKey.getKeySWmatrix(part.skHandle,keyID) :
      pubKey.getAnyKeySWmatrix(part.skHandle);

    assert(W.toKeyID>=0);      // verify that a switching matrix exists

    if (g>1) { // g==1 for CKKS, g>1 for BGV
      g = GCD(W.ptxtSpace, g); // verify that the plaintext spaces match
      assert (g>1);
      tmp.ptxtSpace = g;
    }    
    tmp.keySwitchPart(part, W); // switch this part & update noiseVar
  }
  *this = tmp;
}

void Ctxt::cleanUp()
{
  reLinearize();
  // reduce();
  if (!primeSet.disjointFrom(context.specialPrimes)) {
    modDownToSet(primeSet / context.specialPrimes);
  }
}

// Takes as arguments a key-switching matrix W = W[s'->s] and a
// ciphertext-part p relative to s', uses W to switch p relative to
// (1,s), and adds the and result to *this.
// It is assumed that the part p does not include any of the special
// primes, and that if *this is not an empty ciphertext then its
// primeSet is p.getIndexSet() \union context.specialPrimes
void Ctxt::keySwitchPart(const CtxtPart& p, const KeySwitch& W)
{
  FHE_TIMER_START;

  // no special primes in the input part
  assert(context.specialPrimes.disjointFrom(p.getIndexSet()));

  // For parts p that point to 1 or s, only scale and add
  if (p.skHandle.isOne() || p.skHandle.isBase(W.toKeyID)) { 
    CtxtPart pp = p;
    pp.addPrimesAndScale(context.specialPrimes);
    addPart(pp, /*matchPrimeSet=*/true);
    return; 
  }

  // some sanity checks
  assert(W.fromKey == p.skHandle);  // the handles must match

  // Compute the number of digits that we need and the esitmated
  // added noise from switching this ciphertext part.
  long nDigits;
  NTL::xdouble addedNoise;
  std::tie(nDigits,addedNoise)= keySwitchNoise(p, pubKey, W.ptxtSpace);

  // Break the ciphertext part into digits, if needed, and scale up these
  // digits using the special primes. This is the most expensive operation
  // during homormophic evaluation, so it should be thoroughly optimized.

  vector<DoubleCRT> polyDigits;
  p.breakIntoDigits(polyDigits, nDigits);

  // Finally we multiply the vector of digits by the key-switching matrix
  keySwitchDigits(W, polyDigits);
  noiseVar += addedNoise; // update the noise estimate
}


// Find the IndexSet such that modDown to that set of primes makes the
// additive term due to rounding into the dominant noise term 
void Ctxt::findBaseSet(IndexSet& s) const
{
  if (isCKKS()) {
    findBaseSetCKKS(s);
    return;
  }

  if (getNoiseVar()<=0.0) { // an empty ciphertext
    s = context.ctxtPrimes;
    return;
  }
  assert(verifyPrimeSet());
  bool halfSize = context.containsSmallPrime();
  double curNoise = log(getNoiseVar())/2;
  double firstNoise = context.logOfPrime(0);
  double noiseThreshold = log(modSwitchAddedNoiseVar())*0.55; 
  // FIXME: The above should have been 0.5. Making it a bit more means
  // that we will mod-switch a little less frequently, whether this is
  // a good thing needs to be tested.

  // remove special primes, if they are included in this->primeSet
  s = getPrimeSet();
  if (!s.disjointFrom(context.specialPrimes)) { 
    // scale down noise
    curNoise -= context.logOfProduct(context.specialPrimes);
    s.remove(context.specialPrimes);
  }

  /* We compare below to noiseThreshold+1 rather than to noiseThreshold
   * to make sure that if you mod-switch down to c.findBaseSet() and
   * then immediately call c.findBaseSet() again, it will not tell you
   * to mod-switch further down. Note that mod-switching adds close to
   * noiseThreshold to the scaled noise, so if the scaled noise was
   * equal to noiseThreshold then after mod-switchign you would have
   * roughly twice as much noise. Since we're mesuring the log, it means
   * that you may have as much as noiseThreshold+log(2), which we round
   * up to noiseThreshold+1 in the test below.
   */
  if (curNoise<=noiseThreshold+1) return; // no need to mod down

  // if the first prime in half size, begin by removing it
  if (halfSize && s.contains(0)) {
    curNoise -= firstNoise;
    s.remove(0);
  }

  // while noise is larger than threshold, scale down by the next prime
  while (curNoise>noiseThreshold && !empty(s)) {
    curNoise -= context.logOfPrime(s.last());
    s.remove(s.last());
  }

  // Add 1st prime if s is empty or if this does not increase noise too much
  if (empty(s) || (!s.contains(0) && curNoise+firstNoise<=noiseThreshold)) {
    s.insert(0);
    curNoise += firstNoise;
  }

  if (curNoise>noiseThreshold && log_of_ratio()>-0.5)
    cerr << "Ctxt::findBaseSet warning: already at lowest level\n";
}

// Find the IndexSet such that modDown to that set of primes makes the
// ratFactor only a bit bigger than the additive noise term due to rounding
void Ctxt::findBaseSetCKKS(IndexSet& s) const
{
  if (getNoiseVar()<=0.0) { // an empty ciphertext
    s = context.ctxtPrimes;
    return;
  }
  // check that either all specialPrimes are in, or they are all out
  assert(verifyPrimeSet());

  bool halfSize = context.containsSmallPrime();
  double first = halfSize? context.logOfPrime(0): 0.0;

  double curFactor = log(getRatFactor());
  double threshold = log(modSwitchAddedNoiseVar())
    + log(context.alMod.getPPowR()) + log(context.zMStar.getM())/2;

  // remove special primes, if they are included in this->primeSet
  s = getPrimeSet();
  if (!s.disjointFrom(context.specialPrimes)) { 
    // scale down noise
    curFactor -= context.logOfProduct(context.specialPrimes);
    s.remove(context.specialPrimes);
  }

  // if the first prime in half size, begin by removing it
  if (halfSize && s.contains(0)) {
    curFactor -= first;
    if (curFactor<threshold) return; // cannot even remove the half prime
    s.remove(0);
  }

  // while noise is larger than threshold, scale down by the next prime
  while (!empty(s)) {
    curFactor -= context.logOfPrime(s.last());
    if (curFactor + first < threshold) break; // canot remove this prime
    s.remove(s.last());
  }

  // If curNoise < threshold, add back 1st prime
  if (empty(s) || curFactor < threshold) {
    long idx = (context.ctxtPrimes / s).first(); // 1st prime not in s
    s.insert(idx);
  }

  if (log_of_ratio()>-0.5)
    cerr << "Ctxt::findBaseSetCKKS warning: already at lowest level\n";
}

/********************************************************************/
// Ciphertext arithmetic

// Add/subtract a ciphertext part to a ciphertext.
// With negative=true we subtract, otherwise we add.
void Ctxt::addPart(const DoubleCRT& part, const SKHandle& handle, 
		   bool matchPrimeSet, bool negative)
{
  FHE_TIMER_START;

  assert (&part.getContext() == &context);

  if (parts.size()==0) { // inserting 1st part 
    primeSet = part.getIndexSet();
    parts.push_back(CtxtPart(part,handle));
    if (negative) parts.back().Negate(); // not thread-safe??
  }
  else {       // adding to a ciphertext with existing parts
    if (!(part.getIndexSet() <= primeSet)) {
      // add to the the prime-set of *this, if needed (this is expensive)
      if (matchPrimeSet) {
        IndexSet setDiff = part.getIndexSet() / primeSet; // set minus
        for (size_t i=0; i<parts.size(); i++) parts[i].addPrimes(setDiff);
        primeSet.insert(setDiff);
      }
      else // this should never happen
        throw std::logic_error("part has too many primes and matchPrimeSet==false");
    }

    DoubleCRT tmp(context, IndexSet::emptySet());
    const DoubleCRT* ptr = &part;

    // mod-UP the part if needed
    IndexSet s = primeSet / part.getIndexSet();
    if (!empty(s)) { // if need to mod-UP, do it on a temporary copy
      tmp = part;
      tmp.addPrimesAndScale(s);
      ptr = &tmp;
    }
    long j = getPartIndexByHandle(handle);
    if (j>=0) { // found a matching part, add them up
      if (negative) parts[j] -= *ptr;
      else          parts[j] += *ptr;
    } else {    // no mathing part found, just append this part
      parts.push_back(CtxtPart(*ptr,handle));
      if (negative) parts.back().Negate(); // not thread-safe??
    }
  }
}

// Add a constant polynomial
void Ctxt::addConstant(const DoubleCRT& dcrt, double size)
{
  if (getContext().alMod.getTag()==PA_cx_tag) {
    addConstantCKKS(dcrt, to_xdouble(size));
    return;
  }

  // FIXME: the other addConstant variants should do the scaling
  // in the plaintext space, so as to not add noise
  

  // If the size is not given, use size = phi(m)*(ptxtSpace/2)^2
  if (size < 0.0)
      size = double(context.zMStar.getPhiM()) * fsquare(ptxtSpace/2.0);

  // Scale the constant, then add it to the part that points to one
  long f = 1;
  if (ptxtSpace > 2) {
    f = rem(context.productOfPrimes(primeSet),ptxtSpace);
    f = MulMod(intFactor, f, ptxtSpace);
    f = balRem(f, ptxtSpace);
  }

  noiseVar += size*fsquare(f);

  IndexSet delta = dcrt.getIndexSet() / primeSet; // set minus
  if (f==1 && empty(delta)) { // just add it
    addPart(dcrt, SKHandle(0,1,0));
    return;
  }

  // work with a local copy
  DoubleCRT tmp = dcrt;
  if (!empty(delta)) tmp.removePrimes(delta);
  if (f!=1)          tmp *= f;
  addPart(tmp, SKHandle(0,1,0));
}

// Add a constant polynomial
void Ctxt::addConstant(const ZZ& c)
{
  if (isCKKS()) {
    addConstantCKKS(c);
    return;
  }
  DoubleCRT dcrt(getContext(), getPrimeSet());
  long cc = rem(c, ptxtSpace); // reduce modulo plaintext space
  if (cc > ptxtSpace/2) cc -= ptxtSpace;
  dcrt = cc;

  double size = to_double(cc);

  addConstant(dcrt, fsquare(size));
}


// Add a constant polynomial for CKKS encryption. We assume that
// the constant is scaled by PAlgebraModCx::encodeScalingFactor()
void addSomePrimes(Ctxt& c);
void Ctxt::addConstantCKKS(const DoubleCRT& dcrt, xdouble size, xdouble factor)
{
  if (factor<1.0)
    conv(factor, getContext().alMod.getCx().encodeScalingFactor());

  // If the size is not given, use size = phi(m)*factor^2
  if (size < 0.0) {
    size = context.zMStar.getPhiM() * factor * factor;
  }

  xdouble ratio = floor((ratFactor/factor) +0.5); // round to integer
  double inaccuracy = abs(conv<double>(ratio*factor/ratFactor) - 1.0);

  // Check if you need to scale up to get target accuracy of 2^{-r}
  if ((inaccuracy*getContext().alMod.getPPowR()) > 1.0) {
    addSomePrimes(*this);                   // This increases ratFactor
    ratio = floor((ratFactor/factor) +0.5); // re-compute the ratio
  }
  noiseVar += size*ratio*ratio;

  ZZ intRatio = conv<ZZ>(ratio);
  IndexSet delta = dcrt.getIndexSet() / getPrimeSet(); // set minus

  if (NTL::IsOne(intRatio) && empty(delta)) { // just add it
    addPart(dcrt, SKHandle(0,1,0));
    return;
  }

  // work with a local copy
  DoubleCRT tmp = dcrt;
  if (!empty(delta)) tmp.removePrimes(delta);

  delta = getPrimeSet() / tmp.getIndexSet(); // set minus
  if (!empty(delta)) tmp.addPrimes(delta);   // that's expensive

  if (!NTL::IsOne(intRatio)) tmp *= intRatio;
  addPart(tmp, SKHandle(0,1,0));
}

void Ctxt::addConstantCKKS(const ZZX& poly, xdouble size, xdouble factor)
{
  if (factor<1.0)
    conv(factor, getContext().alMod.getCx().encodeScalingFactor());

  xdouble ratio = floor((ratFactor/factor) +0.5); // round to integer
  double inaccuracy = abs(conv<double>(ratio*factor/ratFactor) - 1.0);

  // Check if you need to scale up to get target accuracy of 2^{-r}
  if ((inaccuracy*getContext().alMod.getPPowR()) > 1.0)
    addSomePrimes(*this);                   // This increases ratFactor

  // Now call the DoubleCRT version
  addConstantCKKS(DoubleCRT(poly,context,primeSet),size,factor);
}

void Ctxt::addConstantCKKS(const ZZ& c)
{
  xdouble xc = to_xdouble(c);
  xdouble scaled = floor(ratFactor*xc +0.5); // scaled up and rounded

  DoubleCRT dcrt(getContext(), getPrimeSet());
  dcrt = to_ZZ(scaled);

  addConstantCKKS(dcrt, /*size=*/xc*xc, /*factor=*/scaled/xc);
}

// Add at least one prime to the primeSet of c
void addSomePrimes(Ctxt& c)
{
  const FHEcontext& context = c.getContext();
  IndexSet s = c.getPrimeSet();

  // Sanity check: there should be something left to add
  assert(!s.contains(context.ctxtPrimes)
         || !s.contains(context.specialPrimes));

  // Add a ctxt prime if possible
  if (!s.contains(context.ctxtPrimes)) {
    IndexSet delta = context.ctxtPrimes / s;  // set minus
    long idx = delta.first(); // We know that |delta| >= 1

    // If this is the small prime, try to add a full-size prime instead
    if (context.containsSmallPrime()
        && idx == context.ctxtPrimes.first() && delta.card()>1) 
      idx = delta.next(idx);

    s.insert(idx);
  }
  else // otherwise , insert all the special primes
    s.insert(context.specialPrimes);

  c.modUpToSet(s);
}

void Ctxt::negate()
{
  for (size_t i=0; i<parts.size(); i++) parts[i].Negate();
}

// scale up c1, c2 so they have the same factor
void Ctxt::equalizeRationalFactors(Ctxt& c1, Ctxt &c2,
                                   pair<long,long> factors)
{
  long targetPrecision = c1.getContext().alMod.getPPowR()*2;

  // if factors are given, use them
  if (factors.first>0 && factors.second>0) {
    c1.multByConstant(to_ZZ(factors.first)); // small times a
    c1.ratFactor *= factors.first;
    c2.multByConstant(to_ZZ(factors.second));  // big times b
    c2.ratFactor *= factors.second;
    assert(closeToOne(c1.ratFactor/c2.ratFactor, targetPrecision));
    cerr << "equalizeFactors using provided scaling factors ["
         << factors.first<<','<<factors.second<<endl;
    cerr << "    resulting ratFactors are ["
         << c1.ratFactor<<','<< c1.ratFactor<<"]\n";
    return;
  }
  // If factors are not given, compute them
  Ctxt& big  = (c1.ratFactor>c2.ratFactor)? c1 : c2;
  Ctxt& small= (c1.ratFactor>c2.ratFactor)? c2 : c1;

  xdouble ratio = big.ratFactor / small.ratFactor;
  if (ratio > targetPrecision) { // just scale up small
    small.multByConstant(to_ZZ(floor(ratio+0.5)));
    cerr << "equalizeFactors scaling small factor from "<<small.ratFactor
         << "to "<<small.ratFactor<<'*'<<ratio<<'=';
    small.ratFactor *= ratio;
    cerr << small.ratFactor<<"\n    large factor is "<<big.ratFactor<<endl;
    return;
  }

  // Otherwise, need to scale both big and small

  // approximate ratio as a fraction a/b
  factors = rationalApprox(to_double(ratio), targetPrecision);
  small.multByConstant(to_ZZ(factors.first)); // small times a
  small.ratFactor *= factors.first;
  big.multByConstant(to_ZZ(factors.second));  // big times b
  big.ratFactor *= factors.second;
  cerr << "equalizeFactors scalign both factor by "
       << factors.first<<','<<factors.second<<endl;
  cerr << "    resulting ratFactors are ["
         << small.ratFactor<<','<< big.ratFactor<<"]\n";
}

// Add/subtract another ciphertxt (depending on the negative flag)
void Ctxt::addCtxt(const Ctxt& other, bool negative)
{
  FHE_TIMER_START;

  // Sanity check: same context and public key
  assert (&context==&other.context && &pubKey==&other.pubKey);

  // Special case: if *this is empty then just copy other
  if (this->isEmpty()) {
    *this = other;
    if (negative) negate();
    return;
  }

  // Verify that the plaintext spaces are compatible
  if (isCKKS())
    assert(getPtxtSpace()==1 && other.getPtxtSpace()==1);
  else // BGV
    this->reducePtxtSpace(other.getPtxtSpace());

  Ctxt tmp(pubKey, other.ptxtSpace); // a temporaty empty ciphertext
  const Ctxt* other_pt = &other;

  // make other ptxtSpace match
  if (ptxtSpace != other_pt->ptxtSpace) {
    tmp = other;
    tmp.reducePtxtSpace(ptxtSpace);
    other_pt = &tmp;
  }

  // Match the prime-sets, mod-UP the arguments if needed
  IndexSet s = other_pt->primeSet / primeSet; // set-minus
  if (!empty(s)) modUpToSet(s);

  s = primeSet / other_pt->primeSet; // set-minus
  if (!empty(s)) { // need to mod-UP the other, use a temporary copy
    if (other_pt != &tmp) { tmp = other; other_pt = &tmp; }
    tmp.modUpToSet(s);
  }

  // If approximate numbers, make sure the scaling factors are the same
  if (isCKKS() && !closeToOne(ratFactor/other_pt->ratFactor,
                              getContext().alMod.getPPowR()*2)) {
    if (other_pt != &tmp) { tmp = other; other_pt = &tmp; }
    equalizeRationalFactors(*this, tmp);
  }

#if 0
  if (intFactor != other_pt->intFactor) { // harmonize factors
    long f1 = this->intFactor;
    long f2 = other_pt->intFactor;
    long ratio = MulMod(f2, InvMod(f1, ptxtSpace), ptxtSpace);
    mulIntFactor(ratio);
  } 
#else // FIXME: The part upto endif should go in a separate function?
  long e1 = 1, e2 = 1;

  if (intFactor != other_pt->intFactor) { // harmonize factors
    long f1 = this->intFactor;
    long f2 = other_pt->intFactor;
    long ratio = MulMod(f2, InvMod(f1, ptxtSpace), ptxtSpace);

    xdouble noise_min = 
      noiseVar*fsquare(balRem(ratio, ptxtSpace)) + other_pt->noiseVar;

    // set e1, e2 so that e1*f1 == e2*f2 (mod ptxtSpace),
    // minimizing the increase in noise.

    e1 = ratio; 
    // Initial choise is e1 == ratio and e2 == 1.
    // Now look for better choices

    for (long ee2 = 2; ee2 < ptxtSpace; ee2++) {
      if (GCD(ee2, ptxtSpace) == 1) {
	long ee1 = MulMod(ee2, ratio, ptxtSpace);
	xdouble noise_est = 
	  noiseVar*fsquare(balRem(ee1, ptxtSpace)) + 
          other_pt->noiseVar*fsquare(balRem(ee2, ptxtSpace));

	if (noise_est < noise_min) {
	  noise_min = noise_est;
	  e1 = ee1;
	  e2 = ee2;
	}
      }
    }
  } 
  if (e2 != 1) {
    if (other_pt != &tmp) { tmp = other; other_pt = &tmp; }
    tmp.mulIntFactor(e2);
  }
  if (e1 != 1) mulIntFactor(e1);
#endif

  // Go over the parts of other, for each one check if
  // there is a matching part in *this
  for (size_t i=0; i<other_pt->parts.size(); i++) {
    const CtxtPart& part = other_pt->parts[i];
    long j = getPartIndexByHandle(part.skHandle);
    if (j>=0) { // found a matching part, add them up
      if (negative) parts[j] -= part;
      else          parts[j] += part;
    } else {    // no mathing part found, just append this part
      parts.push_back(part);
      if (negative) parts.back().Negate(); // not thread safe??
    }
  }
  noiseVar += other_pt->noiseVar;
  highWaterMark = std::min(highWaterMark, other_pt->highWaterMark);
}

long fhe_disable_intFactor = 0;

// Create a tensor product of c1,c2. It is assumed that *this,c1,c2
// are defined relative to the same set of primes and plaintext space.
// It is also assumed that *this DOES NOT alias neither c1 nor c2.
void Ctxt::tensorProduct(const Ctxt& c1, const Ctxt& c2)
{
  clear();                // clear *this, before we start adding things to it
  primeSet = c1.primeSet; // set the correct prime-set before we begin

  long ptxtSp = c1.getPtxtSpace();

  long f = 1;

  if (ptxtSp > 2) {
    if (fhe_disable_intFactor) {
      // c1,c2 may be scaled, so multiply by the inverse scalar if needed
      f = rem(context.productOfPrimes(c1.getPrimeSet()),ptxtSp);
      if (f!=1) f = InvMod(f, ptxtSp);
    }
    else {
      long q = rem(context.productOfPrimes(c1.getPrimeSet()),ptxtSp);
      intFactor = MulMod(c1.intFactor, c2.intFactor, ptxtSp);
      intFactor = MulMod(intFactor, q, ptxtSp);
    }
  }


  // The actual tensoring
  CtxtPart tmpPart(context, IndexSet::emptySet()); // a scratch CtxtPart
  for (size_t i=0; i<c1.parts.size(); i++) {
    CtxtPart thisPart = c1.parts[i];
    if (f!=1) thisPart *= f;
    for (size_t j=0; j<c2.parts.size(); j++) {
      tmpPart = c2.parts[j];
      // What secret key will the product point to?
      if (!tmpPart.skHandle.mul(thisPart.skHandle, tmpPart.skHandle))
	Error("Ctxt::tensorProduct: cannot multiply secret-key handles");

      tmpPart *= thisPart; // The element of the tensor product

      // Check if we already have a part relative to this secret-key handle
      long k = getPartIndexByHandle(tmpPart.skHandle);
      if (k >= 0) // found a matching part
	parts[k] += tmpPart;
      else
	parts.push_back(tmpPart);
    }
  }

  /* Compute the noise estimate as c1.noiseVar * c2.noiseVar * factor
   * where the factor depends on the handles of c1,c2. Specifically,
   * if the largest powerOfS in c1,c2 are n1,n2, respectively, then we
   * have factor = ((n1+n2) choose n2).
   */
  long n1=0,  n2=0;
  for (size_t i=0; i<c1.parts.size(); i++) // get largest powerOfS in c1
    if (c1.parts[i].skHandle.getPowerOfS() > n1)
      n1 = c1.parts[i].skHandle.getPowerOfS();
  for (size_t i=0; i<c2.parts.size(); i++) // get largest powerOfS in c2
    if (c2.parts[i].skHandle.getPowerOfS() > n2)
      n2 = c2.parts[i].skHandle.getPowerOfS();

  // compute ((n1+n2) choose n2)
  long factor = 1;
  for (long i=n1+1; i<=n1+n2; i++) factor *= i;
  for (long i=n2  ; i>1     ; i--) factor /= i;

  noiseVar = c1.noiseVar * c2.noiseVar * factor * context.zMStar.get_cM();
  noiseVar *= fsquare(f);
  ratFactor = c1.ratFactor * c2.ratFactor * f;
}


Ctxt& Ctxt::operator*=(const Ctxt& other)
{
  FHE_TIMER_START;
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return  *this;

  // FIXME: what if other.isEmpty()? should we just
  // do the following?
  //   *this = other;

  // Sanity check: plaintext spaces are compatible
  if (isCKKS())
    assert(getPtxtSpace()==1 && other.getPtxtSpace()==1);
  else { // GBV
    long g = GCD(this->ptxtSpace, other.ptxtSpace);
    assert (g>1);
    this->ptxtSpace = g;
  }

  Ctxt tmpCtxt(this->pubKey, this->ptxtSpace); // a scratch ciphertext
  long lvl = findBaseLevel();
  if (lvl > highWaterMark) {
    lvl = highWaterMark;
    std::cerr << "Ctxt::operator*=: dropping level due to high-water mark\n";
  }
  if (this == &other) {  // a squaring operation
    modDownToLevel(lvl); // mod-down if needed
#ifdef DEBUG_PRINTOUT
      checkNoise(*this, *dbgKey, "modDown " + to_string(size_t(this)));
#endif
    tmpCtxt.tensorProduct(*this, other);  // compute the actual product
    tmpCtxt.noiseVar *= 2;     // a correction factor due to dependency
#ifdef DEBUG_PRINTOUT
      checkNoise(tmpCtxt, *dbgKey, "tensorProduct " + to_string(size_t(this)));
#endif
  }
  else {                // standard multiplication between two ciphertexts
    // Sanity check: same context and public key
    assert (&context==&other.context && &pubKey==&other.pubKey);

    // Match the levels, mod-DOWN the arguments if needed
    long otherLvl = other.findBaseLevel();
    if (otherLvl > other.highWaterMark) {
      otherLvl = other.highWaterMark;
      std::cerr << "Ctxt::operator*=: dropping level due to high-water mark\n";
    }
    if (isCKKS()) {
      highWaterMark = lvl-1;
      if (lvl < otherLvl)
        lvl = otherLvl; // the larger of the two
    }
    else {
      if (lvl > otherLvl) lvl = otherLvl; // the smaller of the two
      highWaterMark = lvl-1;
    }

    // mod-DOWN *this, if needed (also removes special primes, if any)
    bringToLevel(lvl);

    // mod-DOWN other, if needed
    if (primeSet!=other.primeSet){ // use temporary copy to mod-DOWN other
      Ctxt tmpCtxt1 = other;
      tmpCtxt1.bringToLevel(lvl);
      tmpCtxt.tensorProduct(*this,tmpCtxt1); // compute the actual product
    }
    else 
      tmpCtxt.tensorProduct(*this, other);   // compute the actual product
  }
  *this = tmpCtxt; // copy the result into *this

  return *this;
}

// Higher-level multiply routines that include also modulus-switching
// and re-linearization

void Ctxt::multiplyBy(const Ctxt& other)
{
  FHE_TIMER_START;
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  *this *= other;  // perform the multiplication
  reLinearize();   // re-linearize
#ifdef DEBUG_PRINTOUT
      checkNoise(*this, *dbgKey, "reLinearize " + to_string(size_t(this)));
#endif
}

void Ctxt::multiplyBy2(const Ctxt& other1, const Ctxt& other2)
{
  FHE_TIMER_START;
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  long lvl = findBaseLevel();
  long lvl1 = other1.findBaseLevel();
  long lvl2 = other2.findBaseLevel();

  if (lvl<lvl1 && lvl<lvl2){ // if both others at higher levels than this,
    Ctxt tmp = other1;       // multiply others by each other, then by this
    if (&other1 == &other2) tmp *= tmp; // squaring rather than multiplication
    else                    tmp *= other2;

    *this *= tmp;
    reLinearize(); // re-linearize after all the multiplications
    return;
  }

  const Ctxt *first, *second;
  if (lvl<lvl2 || lvl1<lvl2) { // lvl1<=lvl<lvl2 or lvl1<=lvl,lvl2
                               // multiply by other2, then by other1
    first = &other2;
    second = &other1;
  }
  else { // multiply first by other1, then by other2
    first = &other1;
    second = &other2;
  }

  if (this == second) { // handle pointer collision
    Ctxt tmp = *second;
    *this *= *first;
    *this *= tmp;
    if (this == first) // cubing operation
      noiseVar *= 3;   // a correction factor due to dependency
    else
      noiseVar *= 2;   // a correction factor due to dependency
  } else {
    *this *= *first;
    *this *= *second;
  }
  reLinearize(); // re-linearize after all the multiplications
}

// Multiply-by-constant
void Ctxt::multByConstant(const ZZ& c)
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;
  FHE_TIMER_START;

  const ZZ* cPtr = &c;
  if (isCKKS()) {
    xdouble size = to_xdouble(c);
    noiseVar *= size*size * getContext().zMStar.get_cM();
  }
  else { // BGV
    long cc = rem(c, ptxtSpace); // reduce modulo plaintext space
    if (cc > ptxtSpace/2) cc -= ptxtSpace;
    else if (cc < -ptxtSpace/2) cc += ptxtSpace;

    double size = to_double(cc);
    noiseVar *= fsquare(size) * getContext().zMStar.get_cM();

    ZZ tmp = to_ZZ(cc);
    cPtr = &tmp;
  }

  // multiply all the parts by this constant
  for (size_t i=0; i<parts.size(); i++) parts[i] *= (*cPtr);
}

// Multiply-by-constant, it is assumed that the size of this
// constant fits in a double float
void Ctxt::multByConstant(const DoubleCRT& dcrt, double size)
{
  FHE_TIMER_START;
  if (isCKKS()) {
    multByConstantCKKS(dcrt, to_xdouble(size));
    return;
  }
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  // If the size is not given, we use the default value phi(m)*ptxtSpace^2/4
  if (size < 0.0) {
    size = context.zMStar.getPhiM() * fsquare(ptxtSpace/2.0);
  }

  // multiply all the parts by this constant
  for (size_t i=0; i<parts.size(); i++) 
    parts[i].Mul(dcrt,/*matchIndexSets=*/false);

  noiseVar *= size * context.zMStar.get_cM();
}

void Ctxt::multByConstant(const ZZX& poly, double size)
{
  FHE_TIMER_START;
  if (this->isEmpty()) return;
  DoubleCRT dcrt(poly,context,primeSet);
  multByConstant(dcrt,size);
}

void Ctxt::multByConstant(const zzX& poly, double size)
{
  FHE_TIMER_START;
  if (this->isEmpty()) return;
  DoubleCRT dcrt(poly,context,primeSet);
  multByConstant(dcrt,size);
}

void Ctxt::multByConstantCKKS(const DoubleCRT& dcrt, xdouble size, ZZ factor)
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  if (IsZero(factor))
    conv(factor, getContext().alMod.getCx().encodeScalingFactor());

  // If the size is not given, use size = phi(m)*factor^2
  xdouble xfactor = to_xdouble(factor);
  if (size < 0.0)
    size = context.zMStar.getPhiM() * xfactor * xfactor;

  noiseVar *= size * context.zMStar.get_cM();
  ratFactor *= xfactor;

  // multiply all the parts by this constant
  for (size_t i=0; i<parts.size(); i++) 
    parts[i].Mul(dcrt,/*matchIndexSets=*/false);
}

// Divide a cipehrtext by 2. It is assumed that the ciphertext
// encrypts an even polynomial and has plaintext space 2^r for r>1.
// As a side-effect, the plaintext space is halved from 2^r to 2^{r-1}
// If these assumptions are not met then the result will not be a
// valid ciphertext anymore.
void Ctxt::divideBy2()
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;
  assert (ptxtSpace % 2 == 0 && ptxtSpace>2);

  // multiply all the parts by (productOfPrimes+1)/2
  ZZ twoInverse; // set to (Q+1)/2
  getContext().productOfPrimes(twoInverse, getPrimeSet());
  twoInverse += 1;
  twoInverse /= 2;
  for (size_t i=0; i<parts.size(); i++)
    parts[i] *= twoInverse;

  noiseVar /= 4;  // noise is halved by this operation
  ptxtSpace /= 2; // and so is the plaintext space
  intFactor %= ptxtSpace; // adjust intFactor
}

// Divide a cipehrtext by p, for plaintext space p^r, r>1. It is assumed
// that the ciphertext encrypts a polynomial which is zero mod p. If this
// is not the case then the result will not be a valid ciphertext anymore.
// As a side-effect, the plaintext space is reduced from p^r to p^{r-1}.
void Ctxt::divideByP()
{
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  long p = getContext().zMStar.getP();
  assert (ptxtSpace % p == 0 && ptxtSpace>p);

  // multiply all the parts by p^{-1} mod Q (Q=productOfPrimes)
  ZZ pInverse, Q;
  getContext().productOfPrimes(Q, getPrimeSet());
  InvMod(pInverse, conv<ZZ>(p), Q);
  for (size_t i=0; i<parts.size(); i++)
    parts[i] *= pInverse;

  noiseVar  /= fsquare(p);  // noise is reduced by a p factor
  ptxtSpace /= p;           // and so is the plaintext space
  intFactor %= ptxtSpace; // adjust intFactor
}

void Ctxt::automorph(long k) // Apply automorphism F(X)->F(X^k) (gcd(k,m)=1)
{
  FHE_TIMER_START;
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  // Sanity check: verify that k \in Zm*
  assert (context.zMStar.inZmStar(k));
  long m = context.zMStar.getM();

  // Apply this automorphism to all the parts
  for (size_t i=0; i<parts.size(); i++) { 
    parts[i].automorph(k);
    if (!parts[i].skHandle.isOne()) {
      parts[i].skHandle.powerOfX = MulMod(parts[i].skHandle.powerOfX,k,m);
    }
  }
  // no change in noise variance
  FHE_TIMER_STOP;
}


// Apply F(X)->F(X^k) followed by re-liearization. The automorphism is possibly
// evaluated via a sequence of steps, to ensure that we can re-linearize the
// result of every step.
void Ctxt::smartAutomorph(long k) 
{
  FHE_TIMER_START;

  // A hack: record this automorphism rather than actually performing it
  if (isSetAutomorphVals()) { // defined in NumbTh.h
    recordAutomorphVal(k);
    return;
  }
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  // Sanity check: verify that k \in Zm*
  long m = context.zMStar.getM();
  k = mcMod(k, m);
  assert (context.zMStar.inZmStar(k));

  long keyID=getKeyID();
  if (!pubKey.isReachable(k,keyID)) {// must have key-switching matrices for it
    throw std::logic_error("no key-switching matrices for k="+std::to_string(k)
                           + ", keyID="+std::to_string(keyID));
  }

  if (!inCanonicalForm(keyID)) {     // Re-linearize the input, if needed
    reLinearize(keyID);
    assert (inCanonicalForm(keyID)); // ensure that re-linearization succeeded
  }

  while (k != 1) {
    const KeySwitch& matrix = pubKey.getNextKSWmatrix(k,keyID);
    long amt = matrix.fromKey.getPowerOfX();

    // A hack: record this automorphism rather than actually performing it
    if (isSetAutomorphVals2()) { // defined in NumbTh.h
      recordAutomorphVal2(amt);
      return;
    }
    //cerr << "********* automorph " << amt << "\n";
    automorph(amt);
    reLinearize(keyID);
    k = MulMod(k, InvMod(amt,m), m);
  }
  FHE_TIMER_STOP;
}



// applies the Frobenius automorphism p^j
void Ctxt::frobeniusAutomorph(long j) 
{
  FHE_TIMER_START;
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  long m = context.zMStar.getM();
  long p = context.zMStar.getP();
  long d = context.zMStar.getOrdP();

  j = mcMod(j, d);
  long val = PowerMod(p, j, m);
  smartAutomorph(val);
  FHE_TIMER_STOP;
}


/********************************************************************/
// Utility methods

const long Ctxt::getKeyID() const
{
  for (size_t i=0; i<parts.size(); i++)
    if (!parts[i].skHandle.isOne()) return parts[i].skHandle.getSecretKeyID();

  return 0; // no part pointing to anything, return the default key
}

// Estimates the added noise variance from mod-switching down
xdouble Ctxt::modSwitchAddedNoiseVar() const
{
  xdouble addedNoise = to_xdouble(0.0);

  // incorporate the secret keys' Hamming-weight
  for (size_t i=0; i<parts.size(); i++) {
    if (parts[i].skHandle.isOne()) {
      addedNoise += 1.0;
    }
    else {
      long keyId = parts[i].skHandle.getSecretKeyID();
      long d = parts[i].skHandle.getPowerOfS();
      xdouble h, t;
      h = pubKey.getSKeySize(keyId);

      // added noise is d! h^d
      t = h;
      for (long j = 2; j <= d; j++)
        t *= h * j;

      addedNoise += t;
    }
  }
  double roundingNoise = context.zMStar.getPhiM() * context.zMStar.get_cM();
  if (getPtxtSpace()>1)
    roundingNoise *= fsquare(ptxtSpace)/2.0;

  return addedNoise * roundingNoise;
}


// void Ctxt::reduce() const
// {
//   long n = parts.size();
//   for (long i = 0; i < n; i++) parts[i].reduce();
// }

void Ctxt::write(ostream& str) const
{
  writeEyeCatcher(str, BINIO_EYE_CTXT_BEGIN);
  
  /*  Writing out in binary:
    1.  long ptxtSpace
    2.  xdouble noiseVar
    3.  IndexSet primeSet;
    4.  vector<CtxtPart> parts;
  */  
  
  write_raw_int(str, ptxtSpace);
  write_raw_int(str, highWaterMark);
  write_raw_int(str, intFactor);
  write_raw_xdouble(str, ratFactor);
  write_raw_xdouble(str, noiseVar);
  primeSet.write(str);
  write_raw_vector(str, parts);    
 
  writeEyeCatcher(str, BINIO_EYE_CTXT_END);
}

void Ctxt::read(istream& str)
{
  assert(readEyeCatcher(str, BINIO_EYE_CTXT_BEGIN)==0);
  
  ptxtSpace = read_raw_int(str);
  highWaterMark = read_raw_int(str);
  intFactor = read_raw_int(str);
  ratFactor = read_raw_xdouble(str);
  noiseVar = read_raw_xdouble(str); 
  primeSet.read(str); 
  CtxtPart blankCtxtPart(context, IndexSet::emptySet());
  read_raw_vector(str, parts, blankCtxtPart);    

  assert(readEyeCatcher(str, BINIO_EYE_CTXT_END)==0);
}

void CtxtPart::write(ostream& str)
{ 
  this->DoubleCRT::write(str); // CtxtPart is a child.
  skHandle.write(str);
}


void CtxtPart::read(istream& str)
{
  this->DoubleCRT::read(str); // CtxtPart is a child.
  skHandle.read(str);
}


istream& operator>>(istream& str, SKHandle& handle)
{
  seekPastChar(str,'['); // defined in NumbTh.cpp
  str >> handle.powerOfS;
  str >> handle.powerOfX;
  str >> handle.secretKeyID;
  seekPastChar(str,']');
  return str;
}

ostream& operator<<(ostream& str, const CtxtPart& p)
{
  return str << "[" << ((const DoubleCRT&)p) << endl 
	     << p.skHandle << "]"; 
}

istream& operator>>(istream& str, CtxtPart& p)
{
  seekPastChar(str,'['); // defined in NumbTh.cpp
  str >> (DoubleCRT&) p;
  str >> p.skHandle;
  seekPastChar(str,']');
  return str;
}

ostream& operator<<(ostream& str, const Ctxt& ctxt)
{
  str << "["<<ctxt.ptxtSpace<<" "<<ctxt.noiseVar<<" "<<ctxt.primeSet
      << ctxt.highWaterMark << " "
      << ctxt.intFactor << " " << ctxt.ratFactor << " "
      << ctxt.parts.size() << endl;
  for (size_t i=0; i<ctxt.parts.size(); i++)
    str << ctxt.parts[i] << endl;
  return str << "]";
}

istream& operator>>(istream& str, Ctxt& ctxt)
{
  seekPastChar(str,'['); // defined in NumbTh.cpp
  str >> ctxt.ptxtSpace >> ctxt.noiseVar >> ctxt.primeSet
      >> ctxt.highWaterMark >> ctxt.intFactor >> ctxt.ratFactor;
  long nParts;
  str >> nParts;
  ctxt.parts.resize(nParts, CtxtPart(ctxt.context,IndexSet::emptySet()));
  for (long i=0; i<nParts; i++) {
    str >> ctxt.parts[i];
    assert (ctxt.parts[i].getIndexSet()==ctxt.primeSet); // sanity-check
  }
  seekPastChar(str,']');
  return str;
}


void CheckCtxt(const Ctxt& c, const char* label)
{
  cerr << "  "<<label << ", level=" << c.findBaseLevel() << ", log(noise/modulus)~" << c.log_of_ratio() << ", p^r="<<c.getPtxtSpace()<<endl;
}

// The recursive incremental-product function that does the actual work
static void recursiveIncrementalProduct(Ctxt array[], long n)
{
  if (n <= 1) return; // nothing to do

  // split the array in two, first part is the highest power of two smaller
  // than n and second part is the rest

  long ell = NTL::NumBits(n-1); // 2^{l-1} <= n-1
  long n1 = 1UL << (ell-1);     // n/2 <= n1 = 2^l < n

  // Call the recursive procedure separately on the first and second parts
  recursiveIncrementalProduct(array, n1);
  recursiveIncrementalProduct(&array[n1], n-n1);

  // Multiply the last product in the 1st part into every product in the 2nd
  if (n-n1 > 1) {
    NTL_EXEC_RANGE(n-n1, first, last)
    for (long i=n1+first; i<n1+last; i++)
      array[i].multiplyBy(array[n1-1]);
    NTL_EXEC_RANGE_END
  }
  else
    for (long i=n1; i<n; i++) array[i].multiplyBy(array[n1-1]);
}

// For i=n-1...0, set v[i]=prod_{j<=i} v[j]
// This implementation uses depth log n and (nlog n)/2 products
void incrementalProduct(vector<Ctxt>& v)
{
  long n = v.size();  // how many ciphertexts do we have
  if (n > 0) recursiveIncrementalProduct(&v[0], n); // do the actual work
}


static void recursiveTotalProduct(Ctxt& out, const Ctxt array[], long n)
{
  if (n <= 3) {
    out = array[0];
    if (n == 2)      out.multiplyBy(array[1]);
    else if (n == 3) out.multiplyBy2(array[1],array[2]);
    return;
  }

  // split the array in two

  long ell = NTL::NumBits(n-1); // 2^{l-1} <= n-1
  long n1 = 1UL << (ell-1);     // n/2 <= n1 = 2^l < n

  // Call the recursive procedure separately on the first and second parts
  Ctxt out2(ZeroCtxtLike, out);
  recursiveTotalProduct(out, array, n1);
  recursiveTotalProduct(out2, &array[n1], n-n1);

  // Multiply the beginning of the two halves
  out.multiplyBy(out2);
}

// set out=prod_{i=0}^{n-1} v[j], takes depth log n and n-1 products
// out could point to v[0], but having it pointing to any other v[i]
// will make the result unpredictable.
void totalProduct(Ctxt& out, const vector<Ctxt>& v)
{
  long n = v.size();  // how many ciphertexts do we have
  if (n > 0) recursiveTotalProduct(out, &v[0], n); // do the actual work
}


// Compute the inner product of two vectors of ciphertexts, this routine uses
// the lower-level *= operator and does only one re-linearization at the end.
void innerProduct(Ctxt& result, const vector<Ctxt>& v1, const vector<Ctxt>& v2)
{
  long n = min(v1.size(), v2.size());
  if (n<=0) {
    result.clear();
    return;
  }
  result = v1[0]; result *= v2[0];
  for (long i=1; i<n; i++) {
    Ctxt tmp = v1[i];
    tmp *= v2[i];
    result += tmp;
  }
  result.reLinearize();
}

// Compute the inner product of a ciphertext vector and a constant vector
void innerProduct(Ctxt& result,
		  const vector<Ctxt>& v1, const vector<DoubleCRT>& v2)
{
  long n = min(v1.size(), v2.size());
  if (n<=0) {
    result.clear();
    return;
  }
  result = v1[0]; result.multByConstant(v2[0]);
  for (long i=1; i<n; i++) {
    Ctxt tmp = v1[i];
    tmp.multByConstant(v2[i]);
    result += tmp;
  }
}

void innerProduct(Ctxt& result,
		  const vector<Ctxt>& v1, const vector<ZZX>& v2)
{
  long n = min(v1.size(), v2.size());
  if (n<=0) {
    result.clear();
    return;
  }
  result = v1[0]; result.multByConstant(v2[0]);
  for (long i=1; i<n; i++) {
    Ctxt tmp = v1[i];
    tmp.multByConstant(v2[i]);
    result += tmp;
  }
}



// Special-purpose modulus-switching for bootstrapping.
// Mod-switch to an externally-supplied modulus. The modulus need not be in
// the moduli-chain in the context, and does not even need to be a prime.
// The ciphertext *this is not affected, instead the result is returned in
// the zzParts vector, as a vector of ZZX'es. Returns an extimate for the
// noise variance after mod-switching.
#include "powerful.h"
double Ctxt::rawModSwitch(vector<ZZX>& zzParts, long toModulus) const
{
  // Ensure that new modulus is co-prime with plaintetx space
  const long p2r = getPtxtSpace();
  assert(toModulus>1 && p2r>1 && GCD(toModulus,p2r)==1);

  // Compute the ratio between the current modulus and the new one.
  // NOTE: toModulus is a long int, so a double for the logarithms and
  //       xdouble for the ratio itself is sufficient
  xdouble ratio = xexp(log((double)toModulus)
		       - context.logOfProduct(getPrimeSet()));

  // Compute also the ratio modulo ptxtSpace
  const ZZ fromModulus = context.productOfPrimes(getPrimeSet());
  long ratioModP = MulMod(toModulus % p2r, 
			  InvMod(rem(fromModulus,p2r),p2r), p2r);

  mulmod_precon_t precon = PrepMulModPrecon(ratioModP, p2r);

  // cerr << "## converting from mod-"<<context.productOfPrimes(getPrimeSet())
  //      << " to mod-"<<toModulus<<" (ratio="<<ratio
  //      << "), ptxtSpace="<<p2r<<endl;

  // Scale and round all the integers in all the parts
  zzParts.resize(parts.size());
  const PowerfulDCRT& p2d_conv = *context.rcData.p2dConv;
  for (size_t i=0; i<parts.size(); i++) {

    Vec<ZZ> powerful;
    p2d_conv.dcrtToPowerful(powerful, parts[i]); // conver to powerful rep

    for (long j=0; j<powerful.length(); j++) {
      const ZZ& coef = powerful[j];
      long c_mod_p = MulModPrecon(rem(coef,p2r), ratioModP, p2r, precon);
      xdouble xcoef = ratio*conv<xdouble>(coef); // the scaled coefficient

      // round xcoef to an integer which is equal to c_mod_p modulo ptxtSpace
      long rounded = conv<long>(floor(xcoef));
      long r_mod_p = rounded % p2r;
      if (r_mod_p < 0) r_mod_p += p2r; // r_mod_p in [0,p-1]

      if (r_mod_p != c_mod_p) {
        long delta = SubMod(c_mod_p, r_mod_p, p2r);
	// either add delta or subtract toModulus-delta
	rounded += delta;
	if (delta > toModulus-delta) rounded -= p2r;
      }
      // SetCoeff(zzParts[i],j,rounded);
      conv(powerful[j], rounded);  // store back in the powerful vector
    }
    //  if (deg(zzParts[i])<20) cerr << "]\n   scaled poly="<<zzParts[i]<<endl;
    p2d_conv.powerfulToZZX(zzParts[i],powerful); // conver to ZZX
  }

  // Return an estimate for the noise
  return conv<double>(noiseVar*ratio*ratio + modSwitchAddedNoiseVar());
}
