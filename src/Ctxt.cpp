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
#include "CtPtrs.h"

#include "debugging.h"
#include "norms.h"

NTL_CLIENT

extern int fhe_watcher;

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
  const FHEcontext& context = getContext();
  const PAlgebra& zMStar = context.zMStar;

  if (isCKKS()) {
    ptxtSpace=1;

    if (size < 0) {
      // HEURISTIC: we assume that we can safely model the coefficients
      // of ptxt as uniformly and independently distributed over
      // [-magBound, magBound], where magBound = encodeScalingFactor
      double magBound = context.alMod.getCx().encodeScalingFactor();
      long degBound = zMStar.getPhiM();
      noiseBound = context.noiseBoundForUniform(magBound, degBound);
    }
    else 
      noiseBound = size;
  } else { // BGV
    if (size < 0) {
      // HEURISTIC: we assume that we can safely model the coefficients
      // of ptxt as uniformly and independently distributed over
      // [-magBound, magBound], where magBound = ptxtSpace/2
      double magBound = double(ptxtSpace)/2;
      long degBound = zMStar.getPhiM();
      noiseBound = context.noiseBoundForUniform(magBound, degBound);
    }
    else
      noiseBound = size;
  }

  primeSet = context.ctxtPrimes;

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

  s = primeSet & context.ctxtPrimes;   // ctxt primes in primeSet
  return s.isInterval();
}


// Compute the number of digits that we need and the esitmated
// added noise from switching this ciphertext part.
static std::pair<long, NTL::xdouble>
keySwitchNoise(const CtxtPart& p, const FHEPubKey& pubKey, const KeySwitch& ks)
{
  const FHEcontext& context = p.getContext();
  const PAlgebra& palg = context.zMStar;

  xdouble ks_bound = ks.noiseBound;

  long nDigits = 0;
  xdouble addedNoise = to_xdouble(0.0);
  double sizeLeft = context.logOfProduct(p.getIndexSet());
  for (size_t i=0; i<context.digits.size() && sizeLeft>0.0; i++) {    
    nDigits++;

    double digitSize = context.logOfProduct(context.digits[i]);
    if (sizeLeft<digitSize) digitSize=sizeLeft;// need only part of this digit

    // Added noise due to this digit is keySwMatrixNoise * |Di|, 
    // where |Di| is the magnitude of the digit
    addedNoise += ks_bound * xexp(digitSize);

    sizeLeft -= digitSize;
  }

#if 0
  // This needs to be re-thought...and/or implemented elsewhere...

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
#endif

  return std::pair<long, NTL::xdouble>(nDigits,addedNoise);
}

std::pair<long, NTL::xdouble> Ctxt::computeKSNoise(long partIdx, const KeySwitch& ks)
{
  return keySwitchNoise(parts.at(partIdx), pubKey, ks);
}

// Multiply vector of digits by key-switching matrix and add to *this.
// It is assumed that W has at least as many b[i]'s as there are digits.
// The vector of digits is modified in place.
void Ctxt::keySwitchDigits(const KeySwitch& W, vector<DoubleCRT>& digits)
{  // An object to hold the pseudorandom ai's, note that it must be defined
  // with the maximum number of levels, else the PRG will go out of synch.
  // FIXME: This is a bug waiting to happen.

  DoubleCRT ai(context, context.ctxtPrimes | context.specialPrimes);

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
  if (intFactor !=  other.intFactor)  return false;

  // compare ratFactor, ignoring small deviations
  if (ratFactor == 0.0 && other.ratFactor != 0.0) return false;
  xdouble ratio = other.ratFactor / ratFactor;
  if (ratio<0.9 && ratio>1.1) return false;

  // compare noiseBound, ignoring small deviations
  if (noiseBound == 0.0) return (other.noiseBound == 0.0);
  ratio = other.noiseBound / noiseBound;
  return (ratio>0.9 && ratio<1.1);
}

// Constructor
Ctxt::Ctxt(const FHEPubKey& newPubKey, long newPtxtSpace):
  context(newPubKey.getContext()), pubKey(newPubKey), ptxtSpace(newPtxtSpace),
  noiseBound(to_xdouble(0.0))
{
  if (ptxtSpace<2) ptxtSpace = pubKey.getPtxtSpace();
  else assert (GCD(ptxtSpace, pubKey.getPtxtSpace()) > 1); // sanity check
  primeSet=context.ctxtPrimes;
  intFactor = 1;
  ratFactor = 1.0;
}

// Constructor
Ctxt::Ctxt(ZeroCtxtLike_type, const Ctxt& ctxt):
  context(ctxt.getPubKey().getContext()), pubKey(ctxt.getPubKey()), 
  ptxtSpace(ctxt.getPtxtSpace()),
  noiseBound(to_xdouble(0.0))
{
  // same body as previous constructor
  if (ptxtSpace<2) ptxtSpace = pubKey.getPtxtSpace();
  else assert (GCD(ptxtSpace, pubKey.getPtxtSpace()) > 1); // sanity check
  primeSet=context.ctxtPrimes;
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
  noiseBound  = other.noiseBound;
  intFactor = other.intFactor;
  ratFactor = other.ratFactor;
  return *this;
}

// explicitly multiply intFactor by e, which should be
// in the interval [0, ptxtSpace)
void Ctxt::mulIntFactor(long e)
{
  if (e==1) return; // nothing to do
  intFactor = MulMod(intFactor, e, ptxtSpace);
  long bal_e = balRem(e, ptxtSpace);
  for (auto& part : parts) part *= bal_e;
  noiseBound *= abs(bal_e); // because every part was scaled by bal_e
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

  // The noise bound grows by a factor of exp(f)
  noiseBound *= xexp(f);

  // If CKKS, the rational factor grows by a factor of exp(f)
  ratFactor *= xexp(f);

  primeSet.insert(setDiff); // add setDiff to primeSet
  assert(verifyPrimeSet()); // sanity-check: ensure primeSet is still valid
}


void Ctxt::bringToSet(const IndexSet& s) 
{
  auto cap = capacity();
  if (cap<1) {
    std::cerr << "Ctxt::bringToSet called with capacity="<<cap
	 << ", likely decryption error\n";
  }
  if (empty(s)) { // If empty, use a singleton with 1st ctxt prime
    IndexSet tmp(getContext().ctxtPrimes.first());
    modUpToSet(tmp);
    modDownToSet(tmp);
    if (cap>=1)
      std::cerr << "Ctxt::bringToSet called with empty set and capacity="
	<<cap<<", this is likely a bug\n";
  }
  else {
    modUpToSet(s);
    modDownToSet(s);
  }
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
  xdouble addedNoiseBound = modSwitchAddedNoiseBound();

  // For approximate nums, make sure that scaling factor is large enough
  // HERE!!
  if (1 && isCKKS()) {
    // Factor after mod-switching is ratFactor/(prod_{i\in setDiff} qi),
    // it must be larger than addedNoiseBound by at leasr getPPowR()
    double extraFactor = log(addedNoiseBound) + log(context.alMod.getPPowR())
      - ( log(ratFactor) - getContext().logOfProduct(setDiff) );
    // If factor is too small, scale up before mod-down
    if (extraFactor > 0) {
      xdouble xf = ceil(xexp(extraFactor));
      multByConstant(conv<ZZ>(xf)); // Increases noiseBound
      //cout << "*** multByConstant=" << conv<ZZ>(xf) << "\n";
      ratFactor *= xf;              // Up the factor accordingly
    }
  }

  // Special case that never happens, but it worth checking anyways:
  // If the current noise is smaller than the added noise term from
  // mod-switching, then there is no point is actually doing it. In
  // this case just drop the extra primes and stay with the same noise
  // (but a smaller modulus).
  // The reason this never happens is it doesn't make any sense, it is
  // only making the noise/modulus ratio worse without gaining anything.
  // But since Ctxt::modDownToSet is a public method, then maybe some
  // crazy application will call it under these circumstances, so we
  // should check for this condition.

  if (noiseBound < addedNoiseBound) { // just "drop down"
    for (size_t i=0; i<parts.size(); i++)
      parts[i].removePrimes(setDiff);       // remove the primes not in s
    long prodInv = 1;
    if (ptxtSpace>1)
      prodInv = InvMod(rem(context.productOfPrimes(setDiff),ptxtSpace), ptxtSpace);
    if (prodInv > 1) {
      for (size_t i=0; i<parts.size(); i++)
        parts[i] *= prodInv;
      noiseBound *= prodInv;
    }
    cerr << "Ctxt::modDownToSet: DEGENERATE DROP\n";
  } 
  else {                                       // do real mod switching
    for (size_t i=0; i<parts.size(); i++) 
      parts[i].scaleDownToSet(intersection, ptxtSpace);

    // update the noise estimate
    double f = context.logOfProduct(setDiff);
    noiseBound /= xexp(f);
    noiseBound += addedNoiseBound;
    ratFactor /= xexp(f); // The factor in CKKS encryption
  }
  primeSet.remove(setDiff); // remove the primes not in s
  assert(verifyPrimeSet()); // sanity-check: ensure primeSet is still valid
  FHE_TIMER_STOP;
}

void Ctxt::blindCtxt(const ZZX& poly)
{
  Ctxt tmp(pubKey);
  pubKey.Encrypt(tmp,poly,ptxtSpace,/*highNoise=*/true);
  *this += tmp;
  // FIXME: Need to blind the intFactor too
  // FIXME: highNoise does not work for CKKS
  // FIXME: This implementation is not optimized, the levels in the
  //        modulus chain should be handled much better.
}

// Reduce plaintext space to a divisor of the original plaintext space
void Ctxt::reducePtxtSpace(long newPtxtSpace)
{
  long g = GCD(ptxtSpace, newPtxtSpace);
  assert (g>1); // NOTE: Will trigger if called for CKKS ciphertext
  ptxtSpace = g;
  intFactor %= g;
}

// Drop sll smallPrimes and specialPrimes, adding ctxtPrimes
// as necessary to ensure that the scaled noise is above the
// modulus-switching added noise term.
void Ctxt::dropSmallAndSpecialPrimes()
{
  if (primeSet.disjointFrom(context.smallPrimes)) {
    // nothing to do except drop the special primes, if any
    modDownToSet(context.ctxtPrimes);
  }
  else {
    // we will be dropping some smallPrimes, and we need to figure
    // out how much we have to compensate with other ctxtPrimes

    // The target set contains only the ctxtPrimes and its size
    IndexSet target = primeSet & context.ctxtPrimes;
    double log_target = context.logOfProduct(target);

    // Compute the set of dropped primes and its total size
    IndexSet dropping= primeSet / target;
    double log_dropping = context.logOfProduct(dropping);

    // Below we ensure that the scaled ctxt is not too small
    double log_modswitch_noise = log(modSwitchAddedNoiseBound());
    double log_noise = (getNoiseBound()<=0.0)? -DBL_MAX : log(getNoiseBound());
    double log_compensation = 0;

    // For CKKS, try to ensure that the scaling factor remains larger
    // than the mod-switch added noise by a factor of getPPowR()
    if (isCKKS()) {
      double log_bound = log_modswitch_noise + log(context.alMod.getPPowR());
      double log_rf = log(getRatFactor())  // log(factor) after scaling
                      + context.logOfProduct(target) - logOfPrimeSet();
      if (log_rf < log_bound) {
        IndexSet candidates = context.ctxtPrimes / target;
        for (long i: candidates) {
          target.insert(i);
          log_compensation += context.logOfPrime(i);
          if (log_rf + log_compensation >= log_bound)
            break;
        }
      }
    }

    // In either BGV or CKKS, try to ensure that the scaled noise
    // remains not much smaller than the mod-switch added noise.
    // This is done so as not to waste too much apacity.

    log_modswitch_noise += 3*log(2.0); // 3 bits of elbow room
    if (log_noise -log_dropping +log_compensation < log_modswitch_noise) {
      // For CKKS, if we arrived here it means that ratFactor > noiseBound
      //   This is probably a security bug, it can only happen when we
      //   encrypt zero (and tell the encryption routine to use ptxtSize=0)
      if (isCKKS())
        cerr << __func__
             << ": CKKS with ratFactor>noiseBound, encrypting zero?\n";

      IndexSet candidates = context.ctxtPrimes / target;
      for (long i: candidates) {
         target.insert(i);
         log_compensation += context.logOfPrime(i);
         if (log_noise -log_dropping +log_compensation >= log_modswitch_noise)
           break;
      }
    }

    // Finally mod-switch to the right target set
    bringToSet(target);
  }
}


// key-switch to (1,s_i), s_i is the base key with index keyID. If
// keyID<0 then re-linearize to any key for which a switching matrix exists
void Ctxt::reLinearize(long keyID)
{
  FHE_TIMER_START;
  // Special case: if *this is empty or already re-linearized then do nothing
  if (this->isEmpty() || this->inCanonicalForm(keyID)) return;
  // this->reduce();

  dropSmallAndSpecialPrimes();

  long g = ptxtSpace;
  Ctxt tmp(pubKey, ptxtSpace); // an empty ciphertext, same plaintext space
  tmp.intFactor = intFactor;   // same intFactor, too

  double logProd = context.logOfProduct(context.specialPrimes);
  tmp.noiseBound = noiseBound * xexp(logProd);  // The noise after mod-UP
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
    tmp.keySwitchPart(part, W); // switch this part & update noiseBound
  }
  *this = tmp;
}

void Ctxt::cleanUp()
{
  reLinearize();
  // reduce();
  if (!primeSet.disjointFrom(context.specialPrimes)
      || !primeSet.disjointFrom(context.smallPrimes)) {
    dropSmallAndSpecialPrimes();
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
  std::tie(nDigits,addedNoise)= keySwitchNoise(p, pubKey, W);

  // Break the ciphertext part into digits, if needed, and scale up these
  // digits using the special primes. This is the most expensive operation
  // during homormophic evaluation, so it should be thoroughly optimized.

  vector<DoubleCRT> polyDigits;
  p.breakIntoDigits(polyDigits, nDigits);

  // Finally we multiply the vector of digits by the key-switching matrix
  keySwitchDigits(W, polyDigits);
  noiseBound += addedNoise; // update the noise estimate
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
        for (size_t i=0; i<parts.size(); i++) {
           Warning("addPrimes called in addPart");
           parts[i].addPrimes(setDiff);
        }
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

  // If the size is not given, we use a bound based on the assumption
  // that the coefficients are uniformly and independently distributed
  // over [-ptxtSpace/2, ptxtSpace/2]
  if (size < 0.0)
      size = context.noiseBoundForUniform(double(ptxtSpace)/2.0, context.zMStar.getPhiM());

  // Scale the constant, then add it to the part that points to one
  long f = 1;
  if (ptxtSpace > 2) {
    f = rem(context.productOfPrimes(primeSet),ptxtSpace);
    f = MulMod(intFactor, f, ptxtSpace);
    f = balRem(f, ptxtSpace);
  }

  noiseBound += size*abs(f);

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

  addConstant(dcrt, size);
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
    size = context.noiseBoundForUniform(factor, getContext().zMStar.getPhiM());
  }

  xdouble ratio = floor((ratFactor/factor) +0.5); // round to integer
  double inaccuracy = abs(conv<double>(ratio*factor/ratFactor) - 1.0);

  // Check if you need to scale up to get target accuracy of 2^{-r}
  if ((inaccuracy*getContext().alMod.getPPowR()) > 1.0) {
    addSomePrimes(*this);                   // This increases ratFactor
    ratio = floor((ratFactor/factor) +0.5); // re-compute the ratio
  }

  noiseBound += size*ratio;

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
  // just call the DoubleCRT version
  addConstantCKKS(DoubleCRT(poly,context,primeSet),size,factor);
}

void Ctxt::addConstantCKKS(const ZZ& c)
{
  xdouble xc = to_xdouble(c);
  xdouble scaled = floor(ratFactor*xc +0.5); // scaled up and rounded

  DoubleCRT dcrt(getContext(), getPrimeSet());
  dcrt = to_ZZ(scaled);

  addConstantCKKS(dcrt, /*size=*/xc, /*factor=*/scaled/xc);
}

// Add the rational constant num.first / num.second
void Ctxt::addConstantCKKS(std::pair<long,long> num)
{
#if 1
  // Check if you need to scale up to get target accuracy of 2^{-r}
  xdouble xb = to_xdouble(num.second);        // denominator

  xdouble ratio = floor((ratFactor/xb) +0.5); // round to integer
  double inaccuracy = abs(conv<double>(ratio*xb/ratFactor) - 1.0);
  if ((inaccuracy*getContext().alMod.getPPowR()) > 1.0)
    addSomePrimes(*this); // This increases ratFactor

  // scaled up and round the numerator
  xdouble scaled = floor(num.first*ratFactor/xb +0.5);

  DoubleCRT dcrt(getContext(), getPrimeSet());
  dcrt = to_ZZ(scaled);

  addConstantCKKS(dcrt, /*size=*/scaled, /*factor=*/ratFactor);
#else
  // simpler alternative?
  DoubleCRT dcrt(getContext(), getPrimeSet());
  dcrt = to_ZZ(num.first);
  addConstantCKKS(dcrt, /*size=*/xdouble(num.first), /*factor=*/xdouble(num.second));
  
#endif
}

// Add at least one prime to the primeSet of c
void addSomePrimes(Ctxt& c)
{
  const FHEcontext& context = c.getContext();
  IndexSet s = c.getPrimeSet();

  // Sanity check: there should be something left to add
  assert(s != context.allPrimes());

  // Add a ctxt prime if possible
  if (!s.contains(context.ctxtPrimes)) {
    IndexSet delta = context.ctxtPrimes / s;  // set minus
    long idx = delta.first(); // We know that |delta| >= 1

    s.insert(idx);
  }
  // else, add a small prime if possible
  else if (!s.contains(context.smallPrimes)) {
    IndexSet delta = context.smallPrimes / s;  // set minus
    long idx = delta.first(); // We know that |delta| >= 1

    s.insert(idx);
  }
  else // otherwise, insert all the special primes
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

#ifdef DEBUG_PRINTOUT
    cerr << "equalizeFactors using provided scaling factors ["
         << factors.first<<','<<factors.second<<"]\n";
    cerr << "    resulting ratFactors are ["
         << c1.ratFactor<<','<< c1.ratFactor<<"]\n";
#endif
    return;
  }
  // If factors are not given, compute them
  Ctxt& big  = (c1.ratFactor>c2.ratFactor)? c1 : c2;
  Ctxt& small= (c1.ratFactor>c2.ratFactor)? c2 : c1;

  xdouble ratio = big.ratFactor / small.ratFactor;
  if (ratio > targetPrecision) { // just scale up small
    small.multByConstant(to_ZZ(floor(ratio+0.5)));

#ifdef DEBUG_PRINTOUT
    cerr << "equalizeFactors scaling small factor from "<<small.ratFactor
         << " to "<<small.ratFactor<<'*'<<ratio<<" = "
         << (small.ratFactor*ratio) << endl;
    cerr << "    large factor is "<<big.ratFactor<<endl;
#endif
    small.ratFactor *= ratio;
    return;
  }

  // Otherwise, need to scale both big and small

  // approximate ratio as a fraction a/b
  factors = rationalApprox(to_double(ratio), targetPrecision);
  small.multByConstant(to_ZZ(factors.first)); // small times a
  small.ratFactor *= factors.first;
  big.multByConstant(to_ZZ(factors.second));  // big times b
  big.ratFactor *= factors.second;
#ifdef DEBUG_PRINTOUT
  cerr << "equalizeFactors scaling both factor by ["
       << factors.first<<','<<factors.second<<"]\n";
  cerr << "    resulting ratFactors are ["
       << small.ratFactor<<','<< big.ratFactor<<"]\n";
#endif
}

static xdouble 
NoiseNorm(xdouble noise1, xdouble noise2, long e1, long e2, long p)
{
  return noise1*balRem(e1, p) + noise2*balRem(e2, p);
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

  Ctxt tmp(pubKey, other.ptxtSpace); // a temporary empty ciphertext
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

  long e1 = 1, e2 = 1;

  if (intFactor != other_pt->intFactor) { // harmonize factors
    long f1 = intFactor;
    long f2 = other_pt->intFactor;
    // set e1, e2 so that e1*f1 == e2*f2 (mod ptxtSpace),
    // minimizing the increase in noise.

    long ratio = MulMod(f2, InvMod(f1, ptxtSpace), ptxtSpace); // f2/f1
    // so equivalently, we want e1 = e2*ratio (mod ptxtSpace)

    xdouble noise1 = noiseBound;
    xdouble noise2 = other_pt->noiseBound;

    // now we run the extended Euclidean on (ptxtSpace, ratio)
    // to generate pairs (r_i, t_i) such that r_i = t_i*ratio (mod ptxtSpace).

    long r0 = ptxtSpace, t0 = 0;
    long r1 = ratio,     t1 = 1;

    long e1_best = r1,   e2_best = t1;
    xdouble noise_best = NoiseNorm(noise1, noise2, e1_best, e2_best, ptxtSpace);

    long p = context.zMStar.getP();

    while (r1 != 0) {
       long q = r0/r1;
       long r2 = r0 % r1;
       long t2 = t0 - t1*q;
       r0 = r1; r1 = r2;
       t0 = t1; t1 = t2;

       long e1_try = mcMod(r1, ptxtSpace), e2_try = mcMod(t1, ptxtSpace);
       if (e1_try % p != 0) {
	 xdouble noise_try = NoiseNorm(noise1, noise2, e1_try, e2_try, ptxtSpace);
	 if (noise_try < noise_best) {
	    e1_best = e1_try;
	    e2_best = e2_try;
	    noise_best = noise_try;
	 }
       }
    }

    e1 = e1_best;
    e2 = e2_best;

    assert(MulMod(e1, f1, ptxtSpace) == MulMod(e2, f2, ptxtSpace));
    assert(GCD(e1, ptxtSpace) == 1 && GCD(e2, ptxtSpace) == 1);
  } 


  if (e2 != 1) {
    if (other_pt != &tmp) { tmp = other; other_pt = &tmp; }
    tmp.mulIntFactor(e2);
  }
  if (e1 != 1) mulIntFactor(e1);

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
  noiseBound += other_pt->noiseBound;
}

//long fhe_disable_intFactor = 0;

// Create a tensor product of c1,c2. It is assumed that *this,c1,c2
// are defined relative to the same set of primes and plaintext space.
// It is also assumed that *this DOES NOT alias neither c1 nor c2.
void Ctxt::tensorProduct(const Ctxt& c1, const Ctxt& c2)
{
  clear();                // clear *this, before we start adding things to it
  primeSet = c1.primeSet; // set the correct prime-set before we begin

  long ptxtSp = c1.getPtxtSpace();

  if (ptxtSp > 2) { // BGV, handle the integer factor
      long q = rem(context.productOfPrimes(c1.getPrimeSet()),ptxtSp);
      intFactor = MulMod(c1.intFactor, c2.intFactor, ptxtSp);
      intFactor = MulMod(intFactor, q, ptxtSp);
  }

  // The actual tensoring
  CtxtPart tmpPart(context, IndexSet::emptySet()); // a scratch CtxtPart
  for (long i: range(c1.parts.size())) { 
    CtxtPart thisPart = c1.parts[i];
    //    if (f!=1) thisPart *= f;
    for (long j: range(c2.parts.size())) { 
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

  // Compute the noise estimate as c1.noiseBound * c2.noiseBound

  noiseBound = c1.noiseBound * c2.noiseBound;
  ratFactor = c1.ratFactor * c2.ratFactor;
}


void computeIntervalForMul(double& lo, double& hi, const Ctxt& ctxt1, const Ctxt& ctxt2)
{
  const double slack = 5*log(2.0);
  // FIXME: 5 bits of slack...could be something more dynamic

  const FHEcontext& context = ctxt1.getContext();

  double cap1 = ctxt1.capacity();
  double cap2 = ctxt2.capacity();

  double adn = log(ctxt1.modSwitchAddedNoiseBound());
  // should be the same for both ciphertexts

  double safety = 1*log(2.0); // 1 bits of safety

  hi = min(cap1, cap2) + adn - safety;
  lo = hi - slack;


  // FIXME: this is a bit hackish...

  // The idea is that for a given ctxt with modulus q and noise
  // bound n, we want to mod switch to a new modulus q' such that
  // n/(q/q') \approx AddedNoiseBound.
  // Taking logs, this is the same as saying that
  // log(q') \approx adn + (log(q) - log(n)) = adn + ctxt.capacity();

  // Right now, we just set hi to the minimum for both ciphertexts,
  // and set lo a few bits lower, so that we have some flexibility
  // in finding an efficient dropping strategy.

  // It may be worthwhile to experiment with other strategies,
  // such as setting hi = max(cap1, cap2) + adn - safety,
  // and lo = min(hi - slack, min(cap1, cap2) + adn - safety.

  // HERE!!
  if (1 && ctxt1.isCKKS()) { // ensure large enough scaling factor
    double lvl1 = ctxt1.logOfPrimeSet();
    double rf1 = log(ctxt1.getRatFactor());
    double nrf1 = rf1 - (lvl1-lo); // log of ratFactor after scaling

    double lvl2 = ctxt2.logOfPrimeSet();
    double rf2 = log(ctxt2.getRatFactor());
    double nrf2 = rf2 - (lvl2-lo); // log of ratFactor after scaling

    double nrf = min(nrf1, nrf2);
    // increase lo as necessary to ensure that nrf is larger
    // than the modswitch added noise by at least getPPowR()

    double prec = log(ctxt1.getContext().alMod.getPPowR());
    //cout << "*** computeIntervalForMul:\n";
    //cout << "\t log(modSwAddNoise)="<<adn<<endl;
    //cout << "\t log(q)="<<lvl1<<endl;
    //cout << "\t log(noiseBound)="<<(lvl1-cap1)<<endl;
    //cout << "\t log(factor)="<<rf1<<endl;
    //cout << "\t log(prec)="<<prec<<endl;
    if (nrf < adn+prec) {  
      //cout << "\t lo: "<<lo;
      lo += adn +prec -nrf;
      //cout << " -> "<<lo<<endl;
      //cout << "\t hi: "<<hi;
      hi = max(hi, lo + 1); // ensure that hi is a little bigger than lo
      //cout << " -> "<<hi<<endl;
    }
    else {
      //cout << "\t lo: "<<lo<<endl;
      //cout << "\t hi: "<<hi<<endl;
    }
  }
}

void computeIntervalForSqr(double& lo, double& hi, const Ctxt& ctxt)
{
  computeIntervalForMul(lo, hi, ctxt, ctxt);
}

double Ctxt::naturalSize() const
{
  double lo, hi;
  computeIntervalForSqr(lo, hi, *this);
  return hi;
}

IndexSet Ctxt::naturalPrimeSet() const
{
  double lo, hi;
  computeIntervalForSqr(lo, hi, *this);

  IndexSet retval = context.modSizes.getSet4Size(lo, hi, primeSet, isCKKS());
  return retval;
}


// This is essentially operator*=, but with an extra parameter
void Ctxt::multLowLvl(const Ctxt& other_orig, bool destructive)
{
  FHE_TIMER_START;

  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) 
    return;

  if (other_orig.isEmpty()) {
    *this = other_orig;
    return;
  }


  assert(isCKKS() == other_orig.isCKKS());
  assert(&context==&other_orig.context && &pubKey==&other_orig.pubKey);
  assert(!isCKKS() || (getPtxtSpace() == 1 && other_orig.getPtxtSpace() == 1));

  Ctxt* other_pt = nullptr;
  unique_ptr<Ctxt> ct; // scratch space if needed
  if (this == &other_orig) { // squaring
    IndexSet nat = naturalPrimeSet();
    bringToSet(nat); // drop to the "natural" primeSet
    if (dbgKey) { // HERE
      cerr << "*** after bringToSet,    noise/estNoise= "
	   << realToEstimatedNoise(*this, *dbgKey)
           << ", capacity= " << this->bitCapacity()
	   << "\n";
    }
    other_pt = this;
  }
  else { // real multiplication

    // If this is a non-destructive call, make a copy of other
    if (destructive)
      other_pt = (Ctxt*) &other_orig; // cast away const-ness
    else {  // work with a copy
      ct.reset(new Ctxt(other_orig)); // make a copy
      other_pt = ct.get();            // point to it
    }

    // equalize plaintext spaces
    if (!isCKKS()) {
      long g = GCD(ptxtSpace, other_pt->ptxtSpace);
      assert (g>1);
      ptxtSpace = other_pt->ptxtSpace = g;
    }

    // Compute commonPrimeSet, which defines the modulus q of the product

    // To do this, we first compute an interval [lo, hi] in which
    // log(q) should lie in order to properly manage noise growth
    double lo, hi;
    computeIntervalForMul(lo, hi, *this, *other_pt);

    // We then compute commonPrimeSet in a way that minimizes
    // the computational cost of dropping to it
    IndexSet commonPrimeSet = 
      context.modSizes.getSet4Size(lo, hi,
                                   primeSet, other_pt->primeSet, isCKKS());

    // drop the prime sets of *this and other
    bringToSet(commonPrimeSet);
    other_pt->bringToSet(commonPrimeSet);
  }

  // Perform the actual tensor product
  Ctxt tmpCtxt(pubKey, ptxtSpace);

  tmpCtxt.tensorProduct(*this, *other_pt);
  *this = tmpCtxt;

  if (dbgKey) { // HERE
    cerr << "*** after TensorProduct, noise/estNoise= "
         << realToEstimatedNoise(*this, *dbgKey)
           << ", capacity= " << this->bitCapacity()
         << "\n";
  }
}


// Higher-level multiply routines that include also modulus-switching
// and re-linearization

void Ctxt::multiplyBy(const Ctxt& other)
{
  FHE_TIMER_START;
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  if (other.isEmpty()) {
    *this = other;
    return;
  }

  if (dbgKey) { // HERE
    cerr << "*** before multiplyBy,   noise/estNoise= "
         << realToEstimatedNoise(*this, *dbgKey)
           << ", capacity= " << this->bitCapacity()
         << "\n";
  }

  *this *= other;  // perform the multiplication

  reLinearize();   // re-linearize

  if (dbgKey) { // HERE
    cerr << "*** after reLinearize,   noise/estNoise= "
         << realToEstimatedNoise(*this, *dbgKey)
         << ", capacity= " << this->bitCapacity() << endl;
  }

#ifdef DEBUG_PRINTOUT
      checkNoise(*this, *dbgKey, "reLinearize " + to_string(size_t(this)));
#endif
}

void Ctxt::multiplyBy2(const Ctxt& other1, const Ctxt& other2)
{
  FHE_TIMER_START;
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  if (other1.isEmpty()) {
    *this = other1;
    return;
  }

  if (other2.isEmpty()) {
    *this = other2;
    return;
  }

  long cap = capacity();
  long cap1 = other1.capacity();
  long cap2 = other2.capacity();

  if (cap<cap1 && cap<cap2){ // if both others at higher levels than this,
    Ctxt tmp = other1;       // multiply others by each other, then by this
    if (&other1 == &other2) tmp *= tmp; // squaring rather than multiplication
    else                    tmp *= other2;

    *this *= tmp;
    reLinearize(); // re-linearize after all the multiplications
    return;
  }

  const Ctxt *first, *second;
  if (cap<cap2 || cap1<cap2) { // cap1<=cap<cap2 or cap1<=cap,cap2
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

  ZZ c_copy = c;

  if (isCKKS()) {
    xdouble size = fabs(to_xdouble(c));
    noiseBound *= size;
  }
  else { // BGV
    long cc = balRem(rem(c, ptxtSpace), ptxtSpace); // reduce modulo plaintext space
    noiseBound *= abs(cc);
    c_copy = cc;
  }

  // multiply all the parts by this constant
  for (long i: range(parts.size())) parts[i] *= c_copy;
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

  // If the size is not given, we use the default value coreesponding to 
  // uniform dist'n on [-ptxtSpace/2, ptxtSpace/2].
  if (size < 0.0) {
    size = context.noiseBoundForUniform(double(ptxtSpace)/2.0,
                                        getContext().zMStar.getPhiM());
  }

  // multiply all the parts by this constant
  for (long i: range(parts.size()))
    parts[i].Mul(dcrt,/*matchIndexSets=*/false);

  noiseBound *= size;
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
    size = context.noiseBoundForUniform(xfactor, context.zMStar.getPhiM());

  noiseBound *= size;
  ratFactor *= xfactor;

  // multiply all the parts by this constant
  for (long i: range(parts.size()))
    parts[i].Mul(dcrt,/*matchIndexSets=*/false);
}

void Ctxt::multByConstantCKKS(std::pair<long,long> num)
{
  multByConstant(to_ZZ(num.first)); // multiply by numerator
  ratFactor *= num.second;    // increase the scaling factor  
}

void Ctxt::multByConstantCKKS(double x)
{
  xdouble target = ratFactor/x;
  if (target < getContext().alMod.getCx().encodeScalingFactor())
    multByConstantCKKS(rationalApprox(x)); // "actual multiplication"
  else
    ratFactor = target;             // just adjust the scaling factor
}


// Divide a cipehrtext by 2. It is assumed that the ciphertext
// encrypts an even polynomial and has plaintext space 2^r for r>1.
// As a side-effect, the plaintext space is halved from 2^r to 2^{r-1}
// If these assumptions are not met then the result will not be a
// valid ciphertext anymore.

// FIXME: is this still needed/used?
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
  for (long i: range(parts.size()))
    parts[i] *= twoInverse;

  noiseBound /= 2;  // noise is halved by this operation
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
  for (long i: range(parts.size()))
    parts[i] *= pInverse;

  noiseBound  /= p;  // noise is reduced by a p factor
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
  for (long i: range(parts.size())) { 
    parts[i].automorph(k);
    if (!parts[i].skHandle.isOne()) {
      parts[i].skHandle.powerOfX = MulMod(parts[i].skHandle.powerOfX,k,m);
    }
  }
  // no change in noise bound
  FHE_TIMER_STOP;
}
void Ctxt::complexConj() //  Complex conjugate, same as automorph(m-1)
{
  FHE_TIMER_START;
  // Special case: if *this is empty then do nothing
  if (this->isEmpty()) return;

  // Apply this automorphism to all the parts
  for (long i: range(parts.size())) { 
    parts[i].complexConj();
    if (!parts[i].skHandle.isOne()) {
      parts[i].skHandle.powerOfX
        = context.zMStar.getM() - parts[i].skHandle.powerOfX;
    }
  } // no change in noise bound
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
  if (this->isEmpty() || j==0) return;

  if (isCKKS()) { // For CKKS compute complex conjugate
    if (j&1) complexConj(); // If j is even do nothing
  }
  else {          // For BGV compute frobenius
    long m = context.zMStar.getM();
    long p = context.zMStar.getP();
    long d = context.zMStar.getOrdP();

    j = mcMod(j, d);
    long val = PowerMod(p, j, m);
    smartAutomorph(val);
  }
}


/********************************************************************/
// Utility methods

const long Ctxt::getKeyID() const
{
  for (long i: range(parts.size()))
    if (!parts[i].skHandle.isOne()) return parts[i].skHandle.getSecretKeyID();

  return 0; // no part pointing to anything, return the default key
}

// Estimates the added noise bound from mod-switching down
xdouble Ctxt::modSwitchAddedNoiseBound() const
{
  xdouble addedNoise = to_xdouble(0.0);

  // incorporate the secret keys' Hamming-weight
  for (long i: range(parts.size())) { 
    if (parts[i].skHandle.isOne()) {
      addedNoise += 1.0;
    }
    else {
      long keyId = parts[i].skHandle.getSecretKeyID();
      long d = parts[i].skHandle.getPowerOfS();
      xdouble h = conv<xdouble>(pubKey.getSKeyBound(keyId));

      addedNoise += NTL::power(h, d);
    }
  }

  // FIXME-NOW: not sure if this is right when isCKKS()
  double magBound;
  magBound = double(ptxtSpace)/2.0;

  double roundingNoise = context.noiseBoundForUniform(magBound,
                                                      context.zMStar.getPhiM());
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
    2.  xdouble noiseBound
    3.  IndexSet primeSet;
    4.  vector<CtxtPart> parts;
  */  
  
  write_raw_int(str, ptxtSpace);
  write_raw_int(str, intFactor);
  write_raw_xdouble(str, ratFactor);
  write_raw_xdouble(str, noiseBound);
  primeSet.write(str);
  write_raw_vector(str, parts);
 
  writeEyeCatcher(str, BINIO_EYE_CTXT_END);
}

void Ctxt::read(istream& str)
{
  assert(readEyeCatcher(str, BINIO_EYE_CTXT_BEGIN)==0);
  
  ptxtSpace = read_raw_int(str);
  intFactor = read_raw_int(str);
  ratFactor = read_raw_xdouble(str);
  noiseBound = read_raw_xdouble(str);
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
  str << "["<<ctxt.ptxtSpace<<" "<<ctxt.noiseBound<<" "<<ctxt.primeSet
      << ctxt.intFactor << " " << ctxt.ratFactor << " "
      << ctxt.parts.size() << endl;
  for (long i: range(ctxt.parts.size()))
    str << ctxt.parts[i] << endl;
  return str << "]";
}

istream& operator>>(istream& str, Ctxt& ctxt)
{
  seekPastChar(str,'['); // defined in NumbTh.cpp
  str >> ctxt.ptxtSpace >> ctxt.noiseBound >> ctxt.primeSet
      >> ctxt.intFactor >> ctxt.ratFactor;
  long nParts;
  str >> nParts;
  ctxt.parts.resize(nParts, CtxtPart(ctxt.context,IndexSet::emptySet()));
  for (long i: range(nParts)) {
    str >> ctxt.parts[i];
    assert (ctxt.parts[i].getIndexSet()==ctxt.primeSet); // sanity-check
  }
  seekPastChar(str,']');
  return str;
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
void innerProduct(Ctxt& result, const CtPtrs& v1, const CtPtrs& v2)
{
  long n = min(v1.size(), v2.size());
  if (n<=0) {
    result.clear();
    return;
  }
  result = *v1[0]; result *= *v2[0];
  for (long i=1; i<n; i++) {
    Ctxt tmp = *v1[i];
    tmp *= *v2[i];
    result += tmp;
  }
  result.reLinearize();
}
void innerProduct(Ctxt& result, const vector<Ctxt>& v1,
                  const vector<Ctxt>& v2)
{
  innerProduct(result, CtPtrs_vectorCt((vector<Ctxt>&)v1),
               CtPtrs_vectorCt((vector<Ctxt>&)v2));
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
// the zzParts vector, as a vector of ZZX'es. 
// Returns an extimate for the noise bound after mod-switching.

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
  double scaledNoise = conv<double>(noiseBound*ratio);
  double addedNoise = conv<double>(modSwitchAddedNoiseBound());
#ifdef DEBUG_PRINTOUT
  cerr << "## Ctxt::rawModSwitch: converting from mod-"
       << context.productOfPrimes(getPrimeSet())
       << " to mod-"<<toModulus<<" (ratio="<<ratio
       << "), ptxtSpace="<<p2r<<endl;
  cerr << "             scaledNoise="<< scaledNoise
       << ", addedNoise="<<addedNoise<<endl;
#endif
  return scaledNoise + addedNoise;
  // NOTE: technically, modSwitchAddedNoise bound assumes rounding is
  // done in the polynomial basis, rather than the powerful basis,
  // but the same bounds are still valid
}
