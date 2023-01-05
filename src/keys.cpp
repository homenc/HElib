/* Copyright (C) 2012-2020 IBM Corp.
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
 
/* Copyright (C) 2022 Intel Corporation
* SPDX-License-Identifier: Apache-2.0
*
* Added functionallity for separating the SK, PK, and key switching matrices.
*/

#include <queue>

#include <helib/keys.h>
#include <helib/timing.h>
#include <helib/EncryptedArray.h>
#include <helib/Ptxt.h>
#include "binio.h"
#include <helib/sample.h>
#include <helib/norms.h>
#include <helib/apiAttributes.h>
#include <helib/fhe_stats.h>
#include <helib/log.h>
#include "internal_symbols.h" // DECRYPT_ON_PWFL_BASIS

#include "io.h"

namespace helib {

/******** Utility function to generate RLWE instances *********/

// Assumes that c1 is already chosen by the caller
double RLWE1(DoubleCRT& c0, const DoubleCRT& c1, const DoubleCRT& s, long p)
// Returns a high-probability bound on the L-infty norm
// of the canonical embedding of the decryption of (c0, c1) w/r/to s
{
  // Used with p=1 for CKKS, p>=2 for BGV
  assertTrue<InvalidArgument>(
      p > 0,
      "Cannot generate RLWE instance with nonpositive p"); // Used with p=1 for
                                                           // CKKS, p>=2 for BGV
  const Context& context = s.getContext();
  const PAlgebra& palg = context.getZMStar();

  // choose a short error e
  double stdev = to_double(context.getStdev());
  if (palg.getPow2() == 0) // not power of two
    stdev *= sqrt(palg.getM());
  double bound = c0.sampleGaussianBounded(stdev);

  // Set c0 =  p*e - c1*s.
  // It is assumed that c0,c1 are defined with respect to the same set of
  // primes, but s may be defined relative to a different set. Either way
  // the primes for of c0,c1 are unchanged.
  if (p > 1) {
    c0 *= p;
    bound *= p;
  }

  DoubleCRT tmp(c1);
  tmp.Mul(s, /*matchIndexSets=*/false); // multiply but don't mod-up
  c0 -= tmp;

  return bound;
}

// Choose random c0,c1 such that c0+s*c1 = p*e for a short e
// Returns the variance of the noise canonical-embedding entries
double RLWE(DoubleCRT& c0,
            DoubleCRT& c1,
            const DoubleCRT& s,
            long p,
            NTL::ZZ* prgSeed)
{
  // choose c1 at random (using prgSeed if not nullptr)
  c1.randomize(prgSeed);
  return RLWE1(c0, c1, s, p);
}

/******************** PubKey implementation **********************/
/********************************************************************/
// Computes the keySwitchMap pointers, using breadth-first search (BFS)

PubKey::PubKey(const Context& _context) :
    context(_context), pubEncrKey(*this), recryptEkey(*this)
{
  recryptKeyID = -1;
}

PubKey::PubKey(const PubKey& other) :
    // copy constructor
    context(other.context),
    pubEncrKey(*this),
    skBounds(other.skBounds),
    keySwitching(other.keySwitching),
    keySwitchMap(other.keySwitchMap),
    KS_strategy(other.KS_strategy),
    recryptKeyID(other.recryptKeyID),
    recryptEkey(*this)
{ // copy pubEncrKey,recryptEkey w/o checking the ref to the public key
  pubEncrKey.privateAssign(other.pubEncrKey);
  recryptEkey.privateAssign(other.recryptEkey);
}

void PubKey::clear()
{
  pubEncrKey.clear();
  skBounds.clear();
  keySwitching.clear();
  keySwitchMap.clear();
  recryptKeyID = -1;
  recryptEkey.clear();
}

void PubKey::setKeySwitchMap(long keyId)
{
  // Sanity-check, do we have such a key?
  assertInRange(keyId,
                0l,
                (long)skBounds.size(),
                "No such key found"); // Sanity-check, do we have such a key?
  long m = context.getM();

  // Initialize an array of "edges" (this is easier than searching through
  // all the matrices for every step). This is a list of all the powers n
  // for which we have a matrix W[s_i(X^n) => s_i(X)], as well as the index
  // of that matrix in the keySwitching array.
  typedef std::pair<long, long> keySwitchingEdge;
  std::vector<keySwitchingEdge> edges;
  for (long i = 0; i < (long)keySwitching.size(); i++) {
    const KeySwitch& mat = keySwitching.at(i);
    if (mat.toKeyID == keyId && mat.fromKey.getPowerOfS() == 1 &&
        mat.fromKey.getSecretKeyID() == keyId)
      edges.push_back(keySwitchingEdge(mat.fromKey.getPowerOfX(), i));
  }
  if (keyId >= (long)keySwitchMap.size()) // allocate more space if needed
    keySwitchMap.resize(keyId + 1);

  // initialize keySwitchMap[keyId] with m empty entries (with -1 in them)
  keySwitchMap.at(keyId).assign(m, -1);

  // A standard BFS implementation using a FIFO queue (complexity O(V+E))

  std::queue<long> bfsQueue;
  bfsQueue.push(1); // Push the target node 1 onto the BFS queue
  while (!bfsQueue.empty()) {
    long currentNode = bfsQueue.front();

    // See what other nodes can reach the current one
    for (long j = 0; j < (long)edges.size(); j++) { // go over the edges
      long n = edges[j].first;
      long matrixIndex = edges[j].second;

      long nextNode = NTL::MulMod(currentNode, n, m);
      if (keySwitchMap.at(keyId).at(nextNode) ==
          -1) { // A new node: mark it now
        // Record the index of the matrix that we use for the first step
        keySwitchMap[keyId][nextNode] = matrixIndex;

        bfsQueue.push(nextNode); // push new node onto BFS queue
      }
    }
    bfsQueue.pop(); // We are done with the current node
  }
}

const KeySwitch& PubKey::getKeySWmatrix(const SKHandle& from, long toIdx) const
{
  // First try to use the keySwitchMap
  if (from.getPowerOfS() == 1 && from.getSecretKeyID() == toIdx &&
      toIdx < (long)keySwitchMap.size()) {
    long matIdx = keySwitchMap.at(toIdx).at(from.getPowerOfX());
    if (matIdx >= 0) {
      const KeySwitch& matrix = keySwitching.at(matIdx);
      if (matrix.fromKey == from)
        return matrix;
    }
  }

  // Otherwise resort to linear search
  for (size_t i = 0; i < keySwitching.size(); i++) {
    if (keySwitching[i].toKeyID == toIdx && keySwitching[i].fromKey == from)
      return keySwitching[i];
  }
  return KeySwitch::dummy(); // return this if nothing is found
}

const KeySwitch& PubKey::getAnyKeySWmatrix(const SKHandle& from) const
{
  // First try to use the keySwitchMap
  if (from.getPowerOfS() == 1 &&
      from.getSecretKeyID() < (long)keySwitchMap.size()) {
    long matIdx = keySwitchMap.at(from.getSecretKeyID()).at(from.getPowerOfX());
    if (matIdx >= 0) {
      const KeySwitch& matrix = keySwitching.at(matIdx);
      if (matrix.fromKey == from)
        return matrix;
    }
  }

  // Otherwise resort to linear search
  for (size_t i = 0; i < keySwitching.size(); i++) {
    if (keySwitching[i].fromKey == from)
      return keySwitching[i];
  }
  return KeySwitch::dummy(); // return this if nothing is found
}

bool PubKey::operator==(const PubKey& other) const
{
  if (this == &other)
    return true;

  if (&context != &other.context)
    return false;
  if (!pubEncrKey.equalsTo(other.pubEncrKey, /*comparePkeys=*/false))
    return false;

  if (skBounds.size() != other.skBounds.size())
    return false;
  for (size_t i = 0; i < skBounds.size(); i++)
    if (fabs(skBounds[i] - other.skBounds[i]) > 0.1)
      return false;

  if (keySwitching.size() != other.keySwitching.size())
    return false;
  for (size_t i = 0; i < keySwitching.size(); i++)
    if (keySwitching[i] != other.keySwitching[i])
      return false;

  if (keySwitchMap.size() != other.keySwitchMap.size())
    return false;
  for (size_t i = 0; i < keySwitchMap.size(); i++) {
    if (keySwitchMap[i].size() != other.keySwitchMap[i].size())
      return false;
    for (size_t j = 0; j < keySwitchMap[i].size(); j++)
      if (keySwitchMap[i][j] != other.keySwitchMap[i][j])
        return false;
  }

  // compare KS_strategy, ignoring trailing HELIB_KSS_UNKNOWN
  long n = KS_strategy.length();
  while (n > 0 && KS_strategy[n - 1] == HELIB_KSS_UNKNOWN)
    n--;
  long n1 = other.KS_strategy.length();
  while (n1 > 0 && other.KS_strategy[n1 - 1] == HELIB_KSS_UNKNOWN)
    n1--;
  if (n != n1)
    return false;
  for (long i : range(n)) {
    if (KS_strategy[i] != other.KS_strategy[i])
      return false;
  }

  if (recryptKeyID != other.recryptKeyID)
    return false;
  if (recryptKeyID >= 0 &&
      !recryptEkey.equalsTo(other.recryptEkey, /*comparePkeys=*/false))
    return false;

  return true;
}

bool PubKey::operator!=(const PubKey& other) const { return !(*this == other); }

const Context& PubKey::getContext() const { return context; }
long PubKey::getPtxtSpace() const { return pubEncrKey.ptxtSpace; }
bool PubKey::keyExists(long keyID) const
{
  return (keyID < (long)skBounds.size());
}

double PubKey::getSKeyBound(long keyID) const { return skBounds.at(keyID); }

const std::vector<KeySwitch>& PubKey::keySWlist() const { return keySwitching; }

const KeySwitch& PubKey::getKeySWmatrix(long fromSPower,
                                        long fromXPower,
                                        long fromID,
                                        long toID) const
{
  return getKeySWmatrix(SKHandle(fromSPower, fromXPower, fromID), toID);
}

bool PubKey::haveKeySWmatrix(const SKHandle& from, long toID) const
{
  return getKeySWmatrix(from, toID).toKeyID >= 0;
}

bool PubKey::haveKeySWmatrix(long fromSPower,
                             long fromXPower,
                             long fromID,
                             long toID) const
{
  return haveKeySWmatrix(SKHandle(fromSPower, fromXPower, fromID), toID);
}

bool PubKey::haveAnyKeySWmatrix(const SKHandle& from) const
{
  return getAnyKeySWmatrix(from).toKeyID >= 0;
}

const KeySwitch& PubKey::getNextKSWmatrix(long fromXPower, long fromID) const
{
  long matIdx = keySwitchMap.at(fromID).at(fromXPower);
  return (matIdx >= 0 ? keySwitching.at(matIdx) : KeySwitch::dummy());
}

bool PubKey::isReachable(long k, long keyID) const
{
  return keyID < long(keySwitchMap.size()) && keySwitchMap.at(keyID).at(k) >= 0;
}

long PubKey::getKSStrategy(long dim) const
{
  long index = dim + 1;
  assertTrue<InvalidArgument>(index >= 0l,
                              "Invalid dimension (dim must be at least -1)");
  if (index >= KS_strategy.length()) {
    return HELIB_KSS_UNKNOWN;
  }
  // HERE
  // std::cout << "*** getKSSStrategy for dim " << dim << " = " <<
  // KS_strategy[index] << "\n";
  return KS_strategy[index];
}

void PubKey::setKSStrategy(long dim, int val)
{
  long index = dim + 1;
  assertTrue<InvalidArgument>(index >= 0l,
                              "Invalid dimension (dim must be at least -1)");
  if (index >= KS_strategy.length())
    KS_strategy.SetLength(index + 1, HELIB_KSS_UNKNOWN);
  KS_strategy[index] = val;
  // HERE
  // std::cout << "*** setKSSStrategy for dim " << dim << " = " << val << "\n";
}

// Encrypts plaintext, result returned in the ciphertext argument. When
// called with highNoise=true, returns a ciphertext with noise level
// approximately q/8. For BGV, ptxtSpace is the intended plaintext
//     space, which cannot be co-prime with pubEncrKey.ptxtSpace.
//     The returned value is the plaintext-space for the resulting
//     ciphertext, which is GCD(ptxtSpace, pubEncrKey.ptxtSpace).
// For CKKS, ptxtSpace is a bound on the size of the complex plaintext
//     elements that are encoded in ptxt (before scaling), it is assumed
//     that they are scaled by eacx.encodeScalingFactor(). The
//     returned value is the same as the argument ptxtSpace.

long PubKey::Encrypt(Ctxt& ctxt,
                     const NTL::ZZX& ptxt,
                     long ptxtSpace,
                     bool highNoise) const
{
  HELIB_TIMER_START;
  // NOTE: isCKKS() checks the tag in the alMod  the context
  if (isCKKS()) {
    double pSize = (ptxtSpace <= 0) ? 1.0 : double(ptxtSpace);
    // For CKKSencrypt, ptxtSpace==1 is the defaults size value
    CKKSencrypt(ctxt, ptxt, pSize); // FIXME: handle highNoise in CKKSencrypt
    return ptxtSpace;
  }

  assertEq(this, &ctxt.pubKey, "Public key and context public key mismatch");
  if (ptxtSpace != pubEncrKey.ptxtSpace) { // plaintext-space mismatch
    ptxtSpace = NTL::GCD(ptxtSpace, pubEncrKey.ptxtSpace);
    if (ptxtSpace <= 1)
      throw RuntimeError("Plaintext-space mismatch on encryption");
  }

  // generate a random encryption of zero from the public encryption key
  ctxt = pubEncrKey; // already an encryption of zero, just not a random one
                     // ctxt with two parts, each with all the ctxtPrimes
  ctxt.noiseBound = 0;

  // choose a random small scalar r and a small random error vector (e0,e1),
  // then set ctxt = r*pk + p*(e0,e1) + (ptxt,0),
  // where pk = pubEncrKey, and p = ptxtSpace.

  // The resulting ciphertext decrypts to
  //   r*<sk,pk> + p*(e0 + sk1*e1) + ptxt,
  // where sk = (1, sk1) is the secret key.
  // This leads to a noise bound of:
  //   r_bound*pubEncrKey.noiseBound
  //     + p*e0_bound + p*e1_bound*getSKeyBound()
  //     + ptxt_bound
  //  Here, r_bound, e0_bound, and e1_bound are values
  //  returned by the corresponding sampling routines.
  //  ptxt_bound is somewhat heuristically set assuming
  //  that the coefficients of the ciphertext are uniformly
  //  and independently chosen from the interval [-p/2, p/2].

  DoubleCRT e(context, context.getCtxtPrimes());
  DoubleCRT r(context, context.getCtxtPrimes());
  double r_bound = r.sampleSmallBounded();

  ctxt.noiseBound += r_bound * pubEncrKey.noiseBound;

  // std::cerr << "*** r_bound*pubEncrKey.noiseBound " << r_bound *
  // pubEncrKey.noiseBound << "\n";

  double stdev = to_double(context.getStdev());
  if (context.getZMStar().getPow2() == 0) // not power of two
    stdev *= sqrt(context.getM());

  for (size_t i = 0; i < ctxt.parts.size(); i++) { // add noise to all the parts
    ctxt.parts[i] *= r;
    NTL::xdouble e_bound;

    if (highNoise && i == 0) {
      // we sample e so that coefficients are uniform over
      // [-Q/(8*ptxtSpace)..Q/(8*ptxtSpace)]

      NTL::ZZ B;
      B = context.productOfPrimes(context.getCtxtPrimes());
      B /= (ptxtSpace * 8);

      e_bound = e.sampleUniform(B);
      // FIXME: should we use a bounded version of this?
      // We haven't implemented this yet
    } else {
      e_bound = e.sampleGaussianBounded(stdev);
    }

    e *= ptxtSpace;
    e_bound *= ptxtSpace;

    if (i == 1) {
      e_bound *= getSKeyBound(ctxt.parts[i].skHandle.getSecretKeyID());
    }

    ctxt.parts[i] += e;
    ctxt.noiseBound += e_bound;

    // std::cerr << "*** e_bound " << e_bound << "\n";
  }

  // add in the plaintext
  // FIXME: we should really randomize ptxt, so that each coefficient
  //    has expected value 0
  // NOTE: This relies on the first part, ctxt[0], to have handle to 1

  // This code sequence could be optimized, but there is no point
  long QmodP = rem(context.productOfPrimes(ctxt.primeSet), ptxtSpace);
  NTL::ZZX ptxt_fixed;
  balanced_MulMod(ptxt_fixed, ptxt, QmodP, ptxtSpace);
  ctxt.parts[0] += ptxt_fixed;

  // NOTE: this is a heuristic, as the ptxt is not really random.
  // although, when ptxtSpace == 2, the balanced_MulMod will
  // randomize it
  double ptxt_bound = context.noiseBoundForMod(ptxtSpace, context.getPhiM());

  // FIXME: for now, we print out a warning, but we can consider
  // implementing a more robust randomization and rejection sampling
  // strategy.
  double ptxt_sz =
      NTL::conv<double>(embeddingLargestCoeff(ptxt_fixed, context.getZMStar()));

  if (ptxt_sz > ptxt_bound) {
    Warning("noise bound exceeded in encryption");
  }

  double ptxt_rat = ptxt_sz / ptxt_bound;
  HELIB_STATS_UPDATE("ptxt_rat", ptxt_rat);

  ctxt.noiseBound += ptxt_bound;

  // std::cerr << "*** ptxt_bound " << ptxt_bound << "\n";

  // fill in the other ciphertext data members
  ctxt.ptxtSpace = ptxtSpace;
  ctxt.intFactor = 1;

  // std::cerr << "*** ctxt.noiseBound " << ctxt.noiseBound << "\n";

  // CheckCtxt(ctxt, "after encryption");

  return ptxtSpace;
}

long PubKey::Encrypt(Ctxt& ciphertxt,
                     const zzX& plaintxt,
                     long ptxtSpace,
                     bool highNoise) const
{
  NTL::ZZX tmp;
  convert(tmp, plaintxt);
  return Encrypt(ciphertxt, tmp, ptxtSpace, highNoise);
}

// FIXME: Some code duplication between here and Encrypt above
void PubKey::CKKSencrypt(Ctxt& ctxt,
                         const NTL::ZZX& ptxt,
                         double ptxtSize,
                         double scaling) const
{
  // VJS-FIXME: this routine has a number of issues and should
  // be deprecated in favor of the new EncodedPtxt-based routines

  assertEq(this, &ctxt.pubKey, "Public key and context public key mismatch");

  if (ptxtSize <= 0)
    ptxtSize = 1.0;
  if (scaling <= 0) // assume the default scaling factor
    scaling = getContext().getEA().getCx().encodeScalingFactor() / ptxtSize;

  long m = context.getM();
  long prec = getContext().getAlMod().getPPowR();

  // generate a random encryption of zero from the public encryption key
  ctxt = pubEncrKey; // already an encryption of zero, just not a random one

  // choose a random small scalar r and a small random error vector
  // (e0,e1), then set ctxt = r*pk + (e0,e1) + (ef*ptxt,0), where
  // pk = pubEncrKey, and ef (the "extra factor") is described below

  // The resulting ciphertext decrypts to
  //   r*<sk,pk> + e0 + sk1*e1 + ef*ptxt,
  // where sk = (1,s) is the secret key. This leads to a noise bound
  // of error_bound = r_bound*pubEncrKey.noiseBound
  //                  + e0_bound + e1_bound*getSKeyBound()
  // Here, r_bound, e0_bound, and e1_bound are values returned by the
  // corresponding sampling routines.
  // We also have ptxt_bound = ef*f*ptxtSize, which is tracked separately.
  //
  // The input ptxt is already scaled by a factor f=scaling, and is being
  // further scaled by the extra factor ef, so ef*f is the new scaling
  // factor. The extra factor ef is set as ceil(error_bound*prec/f),
  // so that we have ef*f >= error_bound*prec.

  DoubleCRT e(context, context.getCtxtPrimes());
  DoubleCRT r(context, context.getCtxtPrimes());

  double r_bound = r.sampleSmallBounded(); // r is a {0,+-1} polynomial

  NTL::xdouble error_bound = r_bound * pubEncrKey.noiseBound;
  // VJS-NOTE: why don't the error bounds include the encoding error?

  double stdev = to_double(context.getStdev());
  if (context.getZMStar().getPow2() == 0) // not power of two
    stdev *= sqrt(m);

  for (size_t i = 0; i < ctxt.parts.size(); i++) {
    // add noise to all the parts

    ctxt.parts[i] *= r;

    double e_bound = e.sampleGaussianBounded(stdev);
    // zero-mean Gaussian, sigma=stdev

    ctxt.parts[i] += e;
    if (i == 1) {
      e_bound *= getSKeyBound(ctxt.parts[i].skHandle.getSecretKeyID());
    }
    error_bound += e_bound;
  }

  // Compute the extra scaling factor, if needed
  long ef = NTL::conv<long>(ceil(error_bound * prec / (scaling * ptxtSize)));
  if (ef > 1) { // scale up some more
    ctxt.parts[0] += ptxt * ef;
    scaling *= ef;
  } else { // no need for extra scaling
    ctxt.parts[0] += ptxt;
  }

  // Round size to next power of two so as not to leak too much
  ctxt.ptxtMag = EncryptedArrayCx::roundedSize(ptxtSize);
  ctxt.ratFactor = scaling;
  ctxt.noiseBound = error_bound;
  ctxt.ptxtSpace = 1;
}

void PubKey::CKKSencrypt(Ctxt& ciphertxt,
                         const zzX& plaintxt,
                         double ptxtSize,
                         double scaling) const
{
  NTL::ZZX tmp;
  convert(tmp, plaintxt);
  CKKSencrypt(ciphertxt, tmp, ptxtSize, scaling);
}

// These methods are overridden by secret-key Encrypt
long PubKey::Encrypt(Ctxt& ciphertxt,
                     const NTL::ZZX& plaintxt,
                     long ptxtSpace) const
{
  return Encrypt(ciphertxt, plaintxt, ptxtSpace, /*highNoise=*/false);
}

long PubKey::Encrypt(Ctxt& ciphertxt, const zzX& plaintxt, long ptxtSpace) const
{
  return Encrypt(ciphertxt, plaintxt, ptxtSpace, /*highNoise=*/false);
}

// These two specialisations are here to avoid a circular dependency on
// EncryptedArray
template <>
void PubKey::Encrypt(Ctxt& ciphertxt, const Ptxt<BGV>& plaintxt) const
{
  EncodedPtxt eptxt;
  plaintxt.encode(eptxt);
  Encrypt(ciphertxt, eptxt);
}

template <>
void PubKey::Encrypt(Ctxt& ciphertxt, const Ptxt<CKKS>& plaintxt) const
{
  EncodedPtxt eptxt;
  plaintxt.encode(eptxt, /*mag=*/NextPow2(Norm(plaintxt.getSlotRepr())));
  // set mag=2^(ceil(log2(max(Norm(pa),1))))
  // This hides the actual magnitude somewhat.
  // Note that Encrypt(Ctxt,EncodedPtxt) does not attempt
  // any hiding: this left up to the caller.
  // This logic mimics the logic in the original CKKSencrypt function.

  // Note also that this API does not allow the user to set precision
  // parameter in encode.

  Encrypt(ciphertxt, eptxt);
}

void PubKey::Encrypt(Ctxt& ctxt, const EncodedPtxt_BGV& eptxt) const
{
  HELIB_TIMER_START;

  assertTrue(!isCKKS(), "Encrypt: mismatched BGV ptxt / CKKS ctxt");
  assertEq(this, &ctxt.pubKey, "Encrypt: public key mismatch");
  assertEq(&context, &eptxt.getContext(), "Encrypt: context mismatch");

  long ptxtSpace = eptxt.getPtxtSpace();
  NTL::ZZX ptxt;

  convert(ptxt, eptxt.getPoly());

  // The rest of the code is copy/pasted from the
  // original Encrypt code, except that for now, highNoise
  // is not implemented.  We can put it back if necessary.
  // We may eventually want to completely deprecate the original
  // Encrypt code, which is why it is copy/pasted for now.
  // We could also just invoke
  //    Encrypt(ctxt, ptxt, ptxtSpace, /*highNoise=*/false);
  // at this point for the same effect.

  // VJS-FIXME: I really should get rid of the unnecessary
  // connversions from zzX to ZZX...I've added a zzX version
  // of balanced_mulMod...but I also need zzX versions
  // of DoubleCRT += and friends.

  if (ptxtSpace != pubEncrKey.ptxtSpace) { // plaintext-space mismatch
    ptxtSpace = NTL::GCD(ptxtSpace, pubEncrKey.ptxtSpace);
    if (ptxtSpace <= 1)
      throw RuntimeError("Plaintext-space mismatch on encryption");
  }

  // generate a random encryption of zero from the public encryption key
  ctxt = pubEncrKey; // already an encryption of zero, just not a random one
                     // ctxt with two parts, each with all the ctxtPrimes
  ctxt.noiseBound = 0;

  // choose a random small scalar r and a small random error vector (e0,e1),
  // then set ctxt = r*pk + p*(e0,e1) + (ptxt,0),
  // where pk = pubEncrKey, and p = ptxtSpace.

  // The resulting ciphertext decrypts to
  //   r*<sk,pk> + p*(e0 + sk1*e1) + ptxt,
  // where sk = (1, sk1) is the secret key.
  // This leads to a noise bound of:
  //   r_bound*pubEncrKey.noiseBound
  //     + p*e0_bound + p*e1_bound*getSKeyBound()
  //     + ptxt_bound
  //  Here, r_bound, e0_bound, and e1_bound are values
  //  returned by the corresponding sampling routines.
  //  ptxt_bound is somewhat heuristically set assuming
  //  that the coefficients of the ciphertext are uniformly
  //  and independently chosen from the interval [-p/2, p/2].

  DoubleCRT e(context, context.getCtxtPrimes());
  DoubleCRT r(context, context.getCtxtPrimes());
  double r_bound = r.sampleSmallBounded();

  ctxt.noiseBound += r_bound * pubEncrKey.noiseBound;

  // std::cerr << "*** r_bound*pubEncrKey.noiseBound " << r_bound *
  // pubEncrKey.noiseBound << "\n";

  double stdev = to_double(context.getStdev());
  if (context.getZMStar().getPow2() == 0) // not power of two
    stdev *= sqrt(context.getM());

  for (size_t i = 0; i < ctxt.parts.size(); i++) { // add noise to all the parts
    ctxt.parts[i] *= r;
    NTL::xdouble e_bound;

    e_bound = e.sampleGaussianBounded(stdev);

    e *= ptxtSpace;
    e_bound *= ptxtSpace;

    if (i == 1) {
      e_bound *= getSKeyBound(ctxt.parts[i].skHandle.getSecretKeyID());
    }

    ctxt.parts[i] += e;
    ctxt.noiseBound += e_bound;

    // std::cerr << "*** e_bound " << e_bound << "\n";
  }

  // add in the plaintext
  // FIXME: we should really randomize ptxt, so that each coefficient
  //    has expected value 0
  // NOTE: This relies on the first part, ctxt[0], to have handle to 1

  // This code sequence could be optimized, but there is no point
  long QmodP = rem(context.productOfPrimes(ctxt.primeSet), ptxtSpace);
  NTL::ZZX ptxt_fixed;
  balanced_MulMod(ptxt_fixed, ptxt, QmodP, ptxtSpace);
  ctxt.parts[0] += ptxt_fixed;

  // NOTE: this is a heuristic, as the ptxt is not really random.
  // although, when ptxtSpace == 2, the balanced_MulMod will
  // randomize it
  double ptxt_bound = context.noiseBoundForMod(ptxtSpace, context.getPhiM());

  // FIXME: for now, we print out a warning, but we can consider
  // implementing a more robust randomization and rejection sampling
  // strategy.
  double ptxt_sz =
      NTL::conv<double>(embeddingLargestCoeff(ptxt_fixed, context.getZMStar()));

  if (ptxt_sz > ptxt_bound) {
    Warning("noise bound exceeded in encryption");
  }

  double ptxt_rat = ptxt_sz / ptxt_bound;
  HELIB_STATS_UPDATE("ptxt_rat", ptxt_rat);

  ctxt.noiseBound += ptxt_bound;

  // std::cerr << "*** ptxt_bound " << ptxt_bound << "\n";

  // fill in the other ciphertext data members
  ctxt.ptxtSpace = ptxtSpace;
  ctxt.intFactor = 1;
  ctxt.ratFactor = ctxt.ptxtMag = 1.0;

  // std::cerr << "*** ctxt.noiseBound " << ctxt.noiseBound << "\n";

  // CheckCtxt(ctxt, "after encryption");
}

void PubKey::Encrypt(Ctxt& ctxt, const EncodedPtxt_CKKS& eptxt) const
{
  assertTrue(isCKKS(), "Encrypt: mismatched CKKS ptxt / BGV ctxt");
  assertEq(this, &ctxt.pubKey, "Public key and context public key mismatch");
  assertEq(&context, &eptxt.getContext(), "Encrypt: context mismatch");

  NTL::ZZX ptxt;
  convert(ptxt, eptxt.getPoly());
  double mag = eptxt.getMag();
  double scale = eptxt.getScale();
  double err = eptxt.getErr();

  assertTrue(mag > 0, "CKKS encryption: mag <= 0");
  assertTrue(scale > 0, "CKKS encryption: scale <= 0");
  assertTrue(err > 0, "CKKS encryption: err <= 0");

  long m = context.getM();

  // generate a random encryption of zero from the public encryption key
  ctxt = pubEncrKey; // already an encryption of zero, just not a random one

  // choose a random small scalar r and a small random error vector
  // (e0,e1), then set ctxt = r*pk + (e0,e1) + (ef*ptxt,0), where
  // pk = pubEncrKey, and ef (the "extra factor") is described below

  // The resulting ciphertext decrypts to
  //   r*<sk,pk> + e0 + sk1*e1 + ef*ptxt,
  // where sk = (1,s) is the secret key.
  // So the noise added to ptxt by the encryption process is
  //    error_bound = r_bound*pubEncrKey.noiseBound
  //                  + e0_bound + e1_bound*getSKeyBound()
  // Here, r_bound, e0_bound, and e1_bound are values returned by the
  // corresponding sampling routines.
  //
  // The input ptxt is already scaled by a factor f=scale, and is being
  // further scaled by the extra factor ef, so ef*f is the new scaling
  // factor. The extra factor ef is set as ceil(error_bound/err),
  // so that we have error_bound/(ef*f) < err/fac....that is,
  // the scaled noise added by encryption is less than the scaled
  // noise already present in the encoded ptxt.

  DoubleCRT e(context, context.getCtxtPrimes());
  DoubleCRT r(context, context.getCtxtPrimes());

  double r_bound = r.sampleSmallBounded(); // r is a {0,+-1} polynomial

  NTL::xdouble error_bound = r_bound * pubEncrKey.noiseBound;

  double stdev = to_double(context.getStdev());

  // VJS-NOTE: this should never happen for CKKS
  if (context.getZMStar().getPow2() == 0) // not power of two
    stdev *= sqrt(m);

  for (size_t i = 0; i < ctxt.parts.size(); i++) {
    // add noise to all the parts

    ctxt.parts[i] *= r;

    double e_bound = e.sampleGaussianBounded(stdev);
    // zero-mean Gaussian, sigma=stdev

    ctxt.parts[i] += e;
    if (i == 1) {
      e_bound *= getSKeyBound(ctxt.parts[i].skHandle.getSecretKeyID());
    }
    error_bound += e_bound;
  }

  // Compute the extra scaling factor, if needed

  // VJS-NOTE: note the new logic for computing ef...
  // see comment above.
  long ef = NTL::conv<long>(ceil(error_bound / err));

  if (ef > 1) { // scale up some more
    ctxt.parts[0] += ptxt * ef;
    scale *= ef;
    err *= ef;
  } else { // no need for extra scaling
    ctxt.parts[0] += ptxt;
  }

  // VJS-NOTE: we no longer round to the next power of two:
  // Then encoding routine should take care of setting mag correctly.
  ctxt.ptxtMag = mag;
  ctxt.ratFactor = scale;
  ctxt.noiseBound = error_bound + err;
  ctxt.ptxtSpace = 1;
  ctxt.intFactor = 1;
}

void PubKey::Encrypt(Ctxt& ctxt, const EncodedPtxt& eptxt) const
{
  if (eptxt.isBGV())
    Encrypt(ctxt, eptxt.getBGV());
  else if (eptxt.isCKKS())
    Encrypt(ctxt, eptxt.getCKKS());
  else
    throw LogicError("Encrypt: bad EncodedPtxt");
}

bool PubKey::isCKKS() const
{
  return (getContext().getAlMod().getTag() == PA_cx_tag);
}
// NOTE: Is taking the alMod from the context the right thing to do?

bool PubKey::isBootstrappable() const { return (recryptKeyID >= 0); }

std::ostream& operator<<(std::ostream& str, const PubKey& pk)
{
  pk.writeToJSON(str);
  return str;
}

std::istream& operator>>(std::istream& str, PubKey& pk)
{
  pk.clear();

  pk.readJSON(str);

  return str;
}

void PubKey::writeTo(std::ostream& str) const
{
  SerializeHeader<PubKey>().writeTo(str);
  writeEyeCatcher(str, EyeCatcher::PK_BEGIN);

  // Write out for PubKey
  //  1. Context Base
  //  2. Ctxt pubEncrKey;
  //  3. vector<long> skBounds;
  //  4. vector<KeySwitch> keySwitching;
  //  5. vector< vector<long> > keySwitchMap;
  //  6. Vec<long> KS_strategy
  //  7. long recryptKeyID;
  //  8. Ctxt recryptEkey;

  this->getContext().writeTo(str);
  this->pubEncrKey.writeTo(str);
  write_raw_vector(str, this->skBounds);

  // Keyswitch Matrices
  write_raw_vector(str, this->keySwitching);

  long sz = this->keySwitchMap.size();
  write_raw_int(str, sz);
  for (auto v : this->keySwitchMap)
    write_raw_vector(str, v);

  write_ntl_vec_long(str, this->KS_strategy);

  write_raw_int(str, this->recryptKeyID);
  this->recryptEkey.writeTo(str);

  writeEyeCatcher(str, EyeCatcher::PK_END);
}

PubKey PubKey::readFrom(std::istream& str, const Context& context)
{
  const auto header = SerializeHeader<PubKey>::readFrom(str);
  assertEq<IOError>(header.version,
                    Binio::VERSION_0_0_1_0,
                    "Header: version " + header.versionString() +
                        " not supported");

  bool eyeCatcherFound = readEyeCatcher(str, EyeCatcher::PK_BEGIN);
  assertTrue<IOError>(eyeCatcherFound,
                      "Could not find pre-public key eyecatcher");

  // TODO code to check context object is what it should be same as the text IO.
  // May be worth putting it in helper func.
  // std::unique_ptr<Context> dummy = buildContextFromBinary(str);

  Context ser_context = Context::readFrom(str);
  assertEq(context, ser_context, "Context mismatch");

  PubKey ret(context);

  // Read in the rest
  ret.pubEncrKey.read(str); // Using in-place ctxt read function for performance
  read_raw_vector(str, ret.skBounds); // Using in-place function for performance

  // Keyswitch Matrices
  ret.keySwitching = read_raw_vector<KeySwitch, Context>(str, context);

  long sz = read_raw_int(str);
  ret.keySwitchMap.clear();
  ret.keySwitchMap.resize(sz);
  for (auto& v : ret.keySwitchMap) {
    read_raw_vector(str, v);
  }

  // TODO: Check with VJS if the following loop is really needed
  for (long i = ret.skBounds.size() - 1; i >= 0; i--) {
    ret.setKeySwitchMap(i);
  }

  read_ntl_vec_long(str, ret.KS_strategy);

  ret.recryptKeyID = read_raw_int(str);
  ret.recryptEkey.read(str); // Using in-place ctxt read function for
                             // performance

  eyeCatcherFound = readEyeCatcher(str, EyeCatcher::PK_END);
  assertTrue<IOError>(eyeCatcherFound,
                      "Could not find post-public key eyecatcher");

  return ret;
}

void PubKey::writeToJSON(std::ostream& str) const
{
  executeRedirectJsonError<void>([&]() { str << writeToJSON(); });
}

JsonWrapper PubKey::writeToJSON() const
{
  auto body = [this]() {
    json j = {{"context", unwrap(this->getContext().writeToJSON())},
              {"pubEncrKey", unwrap(this->pubEncrKey.writeToJSON())},
              {"skBounds", this->skBounds},
              {"keySwitching", writeVectorToJSON(keySwitching)},
              {"keySwitchMap", this->keySwitchMap},
              {"KS_strategy", this->KS_strategy},
              {"recryptKeyID", this->recryptKeyID},
              {"recryptEkey",
               this->recryptKeyID >= 0 ? unwrap(this->recryptEkey.writeToJSON())
                                       : "nullptr"}};
    return wrap(toTypedJson<PubKey>(j));
  };
  return executeRedirectJsonError<JsonWrapper>(body);
}

PubKey PubKey::readFromJSON(std::istream& str, const Context& context)
{
  return executeRedirectJsonError<PubKey>([&]() {
    json j;
    str >> j;
    return readFromJSON(wrap(j), context);
  });
}

PubKey PubKey::readFromJSON(const JsonWrapper& jw, const Context& context)
{
  PubKey pk{context};
  pk.readJSON(jw);
  return pk;
}

void PubKey::readJSON(std::istream& str)
{
  executeRedirectJsonError<void>([&]() {
    json j;
    str >> j;
    readJSON(wrap(j));
  });
}

void PubKey::readJSON(const JsonWrapper& tjw)
{
  auto body = [&, this]() {
    json j = fromTypedJson<PubKey>(unwrap(tjw));
    Context ser_context = Context::readFromJSON(wrap(j.at("context")));
    assertEq(context, ser_context, "Context mismatch");

    this->clear();
    //  std::cerr << "PubKey[";

    // Get the public encryption key itself
    this->pubEncrKey.readJSON(wrap(j.at("pubEncrKey")));

    // Get the vector of secret-key Hamming-weights
    j.at("skBounds").get_to(this->skBounds);

    keySwitching = readVectorFromJSON<KeySwitch>(j.at("keySwitching"), context);

    // Get the key-switching map
    this->keySwitchMap =
        j.at("keySwitchMap").get<std::vector<std::vector<long>>>();

    // TODO: Check with VJS if the following loop is really needed
    // build the key-switching map for all keys
    for (long i = this->skBounds.size() - 1; i >= 0; i--)
      this->setKeySwitchMap(i);

    this->KS_strategy = j.at("KS_strategy");

    // Get the bootstrapping key, if any
    this->recryptKeyID = j.at("recryptKeyID");
    if (this->recryptKeyID >= 0) {
      this->recryptEkey.readJSON(wrap(j.at("recryptEkey")));
    }
  };

  executeRedirectJsonError<void>(body);
}

/******************** SecKey implementation **********************/
/********************************************************************/

SecKey::SecKey(const PubKey& pk) : PubKey(pk) {}

SecKey::SecKey(const Context& _context) : PubKey(_context) {}

bool SecKey::operator==(const SecKey& other) const
{
  if (this == &other)
    return true;

  if (((const PubKey&)*this) != ((const PubKey&)other))
    return false;
  if (sKeys.size() != other.sKeys.size())
    return false;
  for (size_t i = 0; i < sKeys.size(); i++)
    if (sKeys[i] != other.sKeys[i])
      return false;
  return true;
}

bool SecKey::operator!=(const SecKey& other) const { return !(*this == other); }

void SecKey::clear()
{
  PubKey::clear();
  sKeys.clear();
}

// We allow the calling application to choose a secret-key polynomial by
// itself, then insert it into the SecKey object, getting the index of
// that secret key in the sKeys list. If this is the first secret-key for this
// SecKey object, then the procedure below generates a corresponding public
// encryption key.
// It is assumed that the context already contains all parameters.
long SecKey::ImportSecKey(const DoubleCRT& sKey,
                          double bound,
                          long ptxtSpace,
                          long maxDegKswitch)
{
  if (sKeys.empty()) { // 1st secret-key, generate corresponding public key
    if (ptxtSpace < 2)
      ptxtSpace = isCKKS() ? 1 : context.getAlMod().getPPowR();
    // default plaintext space is p^r for BGV, 1 for CKKS

    // allocate space, the parts are DoubleCRTs with all the ctxtPrimes
    pubEncrKey.parts.assign(2, CtxtPart(context, context.getCtxtPrimes()));
    // Choose a new RLWE instance
    pubEncrKey.noiseBound =
        RLWE(pubEncrKey.parts[0], pubEncrKey.parts[1], sKey, ptxtSpace);
    if (isCKKS()) {
      pubEncrKey.ptxtMag = 0.0;
      pubEncrKey.ratFactor = pubEncrKey.noiseBound *
                             getContext().getEA().getCx().encodeScalingFactor();
    }

    // make parts[0],parts[1] point to (1,s)
    pubEncrKey.parts[0].skHandle.setOne();
    pubEncrKey.parts[1].skHandle.setBase();

    // Set the other Ctxt bookeeping parameters in pubEncrKey
    pubEncrKey.primeSet = context.getCtxtPrimes();
    pubEncrKey.ptxtSpace = ptxtSpace;
  }
  skBounds.push_back(bound); // record the size of the new secret-key
  sKeys.push_back(sKey);     // add to the list of secret keys
  long keyID =
      sKeys.size() - 1; // FIXME: not thread-safe, do we need to fix it?

  for (long e = 2; e <= maxDegKswitch; e++)
    GenKeySWmatrix(e, 1, keyID, keyID); // s^e -> s matrix

  return keyID; // return the index where this key is stored
}

long SecKey::GenSecKey(long ptxtSpace, long maxDegKswitch)
{
  long hwt = context.getHwt();

  DoubleCRT newSk(context,
                  context.getCtxtPrimes() | context.getSpecialPrimes());

  if (hwt > 0) {
    // sample a Hamming-weight-hwt polynomial
    double bound = newSk.sampleHWtBounded(hwt);
    return ImportSecKey(newSk, bound, ptxtSpace, maxDegKswitch);
  } else {
    // sample a 0/+-1 polynomial
    double bound = newSk.sampleSmallBounded();
    return ImportSecKey(newSk, bound, ptxtSpace, maxDegKswitch);
  }
}

// Generate a key-switching matrix and store it in the public key.
// The argument p denotes the plaintext space
void SecKey::GenKeySWmatrix(long fromSPower,
                            long fromXPower,
                            long fromIdx,
                            long toIdx,
                            long p)
{
  HELIB_TIMER_START;

  // sanity checks
  if (fromSPower <= 0 || fromXPower <= 0)
    return;
  if (fromSPower == 1 && fromXPower == 1 && fromIdx == toIdx)
    return;

  // See if this key-switching matrix already exists in our list
  if (haveKeySWmatrix(fromSPower, fromXPower, fromIdx, toIdx))
    return; // nothing to do here

  DoubleCRT fromKey = sKeys.at(fromIdx);    // copy object, not a reference
  const DoubleCRT& toKey = sKeys.at(toIdx); // this can be a reference

  if (fromXPower > 1)
    fromKey.automorph(fromXPower); // compute s(X^t)
  if (fromSPower > 1)
    fromKey.Exp(fromSPower); // compute s^r(X^t)
  // SHAI: The above lines compute the automorphism and exponentiation mod q,
  //   turns out this is really what we want (even through usually we think
  //   of the secret key as being mod p^r)

  KeySwitch ksMatrix(fromSPower, fromXPower, fromIdx, toIdx);
  RandomBits(ksMatrix.prgSeed, 256); // a random 256-bit seed

  long n = context.getDigits().size();

  // size-n vector
  ksMatrix.b.resize(
      n,
      DoubleCRT(context, context.getCtxtPrimes() | context.getSpecialPrimes()));

  std::vector<DoubleCRT> a;
  a.resize(
      n,
      DoubleCRT(context, context.getCtxtPrimes() | context.getSpecialPrimes()));

  {
    RandomState state;
    SetSeed(ksMatrix.prgSeed);
    for (long i = 0; i < n; i++)
      a[i].randomize();
  } // restore state upon destruction of state

  // Record the plaintext space for this key-switching matrix
  if (isCKKS())
    p = 1;
  else { // BGV
    if (p < 2) {
      if (context.isBootstrappable()) {
        // use larger bootstrapping plaintext space
        p = context.getRcData().alMod->getPPowR();
      } else {
        p = pubEncrKey.ptxtSpace; // default plaintext space from public key
      }
    }
    // FIXME: We use context.isBootstrappable() rather than
    //   this->isBootstrappable(). So we get the larger bootstrapping
    //   plaintext space even if *this is not currently bootstrappable,
    //   in case the calling application will make it bootstrappable later.

    assertTrue(p >= 2,
               "Invalid p value found generating BGV key-switching matrix");
  }
  ksMatrix.ptxtSpace = p;

  // generate the RLWE instances with pseudorandom ai's

  for (long i = 0; i < n; i++) {
    ksMatrix.noiseBound = RLWE1(ksMatrix.b[i], a[i], toKey, p);
  }
  // Add in the multiples of the fromKey secret key
  fromKey *= context.productOfPrimes(context.getSpecialPrimes());
  for (long i = 0; i < n; i++) {
    ksMatrix.b[i] += fromKey;
    fromKey *= context.productOfPrimes(context.getDigit(i));
  }

  // Push the new matrix onto our list
  keySwitching.push_back(ksMatrix);

#if 0
  // HERE
  std::cout
    << "*** ksMatrix: "
    << fromSPower << " " << fromXPower << " " << fromIdx << " "
    << toIdx << " " << p << " "
    << (log(ksMatrix.noiseBound)/log(2.0)) << "\n";
#endif
}

// Decryption
void SecKey::Decrypt(NTL::ZZX& plaintxt, const Ctxt& ciphertxt) const
{
  NTL::ZZX f;
  Decrypt(plaintxt, ciphertxt, f);
}

// These two specialisations are here to avoid a circular dependency on
// EncryptedArray
template <>
void SecKey::Decrypt<BGV>(Ptxt<BGV>& plaintxt,
                          const Ctxt& ciphertxt,
                          UNUSED OptLong prec) const
{
  NTL::ZZX pp;
  Decrypt(pp, ciphertxt);
  plaintxt.decodeSetData(pp);
}

template <>
void SecKey::Decrypt<CKKS>(Ptxt<CKKS>& plaintxt,
                           const Ctxt& ciphertxt,
                           OptLong prec) const
{
  const Context& context = ciphertxt.getContext();
  assertTrue(&context == &plaintxt.getContext(),
             "Decrypt: inconsistent contexts");

  const View& view = context.getView();
  std::vector<std::complex<double>> ptxt;
  view.decrypt(ciphertxt, *this, ptxt, prec);
  plaintxt.setData(ptxt);
}

// VJS-NOTE: this is duplicated code...moreover, it does
// not implement the mitigation against CKKS vulnerability.
#if 0
{
  std::vector<std::complex<double>> ptxt;
  NTL::ZZX pp;
  Decrypt(pp, ciphertxt);
  const long MAX_BITS = 400;
  long nBits = NTL::MaxBits(pp) - MAX_BITS;
  double factor;
  if (nBits <= 0) { // convert to zzX, double
    CKKS_canonicalEmbedding(ptxt,
                            pp,
                            ciphertxt.getContext().ea->getCx().getPAlgebra());
    factor = NTL::to_double(ciphertxt.getRatFactor());
  } else {
    long dpp = deg(pp);
    std::vector<double> pp_scaled(dpp + 1);
    NTL::ZZ tmp;
    for (long i : range(dpp + 1)) {
      RightShift(tmp, pp.rep[i], nBits);
      pp_scaled[i] = NTL::to_double(tmp);
    }
    CKKS_canonicalEmbedding(ptxt,
                            pp_scaled,
                            ciphertxt.getContext().ea->getCx().getPAlgebra());
    factor =
        NTL::to_double(ciphertxt.getRatFactor() / NTL::power2_xdouble(nBits));
  }
  for (auto& cx : ptxt) // divide by the factor
    cx /= factor;

  plaintxt.setData(ptxt);
}
#endif

void SecKey::Decrypt(NTL::ZZX& plaintxt,
                     const Ctxt& ciphertxt,
                     NTL::ZZX& f) const // plaintext before modular reduction
{
  HELIB_TIMER_START;

  // VJS-NOTE: why are we comapring contexts, rather
  // than addresses, which is what we do everywhere else?
  // I'm changing this for now...
  // assertEq(getContext(), ciphertxt.getContext(), "Context mismatch");
  // To be addressed later
  assertEq(&getContext(), &ciphertxt.getContext(), "Context mismatch");

  // this will trigger a warning if any operations that were
  // previously performed on the polynomial basis were invalid
  // because of excess noise.

  if (!ciphertxt.isCorrect()) {
    std::string message = "Decrypting with too much noise";

// TODO: Turn the following preprocessor logics into a warnOrThrow function
#ifdef HELIB_DEBUG
    Warning(message);
#else
    throw LogicError(message);
#endif
  }

  const IndexSet& ptxtPrimes = ciphertxt.primeSet;

  DoubleCRT ptxt(context, ptxtPrimes); // Set to zero

  // for each ciphertext part, fetch the right key, multiply and add
  for (size_t i = 0; i < ciphertxt.parts.size(); i++) {
    const CtxtPart& part = ciphertxt.parts[i];
    if (part.skHandle.isOne()) { // No need to multiply
      ptxt += part;
      continue;
    }

    long keyIdx = part.skHandle.getSecretKeyID();
    DoubleCRT key = sKeys.at(keyIdx); // copy object, not a reference
    key.setPrimes(ptxtPrimes);
    // need to equalize the prime sets without changing prime set of ciphertext.
    // Note that ciphertext may contain small primes, which are not in key.

    long xPower = part.skHandle.getPowerOfX();
    long sPower = part.skHandle.getPowerOfS();
    if (xPower > 1) {
      key.automorph(xPower); // s(X^t)
    }
    if (sPower > 1) {
      key.Exp(sPower); // s^r(X^t)
    }

    key *= part;
    ptxt += key;
  }
  // convert to coefficient representation & reduce modulo the plaintext space

  if (DECRYPT_ON_PWFL_BASIS && !getContext().getZMStar().getPow2()) {
    const PowerfulDCRT& pwfl_converter = getContext().getPowerfulConverter();
    NTL::Vec<NTL::ZZ> pwfl;

    pwfl_converter.dcrtToPowerful(pwfl, ptxt);
    // convert to powerful basis, reduced mod product of primes in prime chain.
    // the reduction mod Q is done on the powerful basis, as the
    // coefficients tend to be smaller there

    pwfl_converter.powerfulToZZX(plaintxt, pwfl);
    // now convert to polynomial basis, with no modular reduction
  } else {
    ptxt.toPoly(plaintxt);
  }

  f = plaintxt; // f used only for debugging

  if (isCKKS())
    return; // CKKS encryption, nothing else to do
  // NOTE: calling application must still divide by ratFactor after decoding

  PolyRed(plaintxt, ciphertxt.ptxtSpace, true /*reduce to [0,p-1]*/);

  // if p>2, multiply by (intFactor * Q)^{-1} mod p
  if (ciphertxt.getPtxtSpace() > 2) {
    long factor = rem(context.productOfPrimes(ciphertxt.getPrimeSet()),
                      ciphertxt.ptxtSpace);
    factor = NTL::MulMod(factor, ciphertxt.intFactor, ciphertxt.ptxtSpace);
    if (factor != 1) {
      factor = NTL::InvMod(factor, ciphertxt.ptxtSpace);
      MulMod(plaintxt, plaintxt, factor, ciphertxt.ptxtSpace, /*abs=*/true);
    }
  }
}

// Encryption using the secret key, this is useful, e.g., to put an
// encryption of the secret key into the public key.
long SecKey::skEncrypt(Ctxt& ctxt,
                       const NTL::ZZX& ptxt,
                       long ptxtSpace,
                       long skIdx) const
{
  // VJS-FIXME: this routine has a number of issues and should
  // be deprecated in favor of the new EncodedPtxt-based routines

  HELIB_TIMER_START;

  assertEq(((const PubKey*)this),
           &ctxt.pubKey,
           "Key does not match context's public key");

  double ptxtSize = 1.0;
  if (isCKKS()) {
    if (ptxtSpace > 0)
      ptxtSize = ptxtSpace;
    ptxtSpace = 1;
  } else { // BGV
    if (ptxtSpace < 2)
      ptxtSpace = pubEncrKey.ptxtSpace; // default plaintext space is p^r
    assertTrue(ptxtSpace >= 2, "Found invalid p value in BGV encryption");
  }
  ctxt.ptxtSpace = ptxtSpace;

  ctxt.primeSet = context.getCtxtPrimes(); // initialize the primeSet
  {
    CtxtPart tmpPart(context, context.getCtxtPrimes());
    ctxt.parts.assign(2, tmpPart);
  } // allocate space

  // Set Ctxt bookeeping parameters
  ctxt.intFactor = 1;

  // make parts[0],parts[1] point to (1,s)
  ctxt.parts[0].skHandle.setOne();
  ctxt.parts[1].skHandle.setBase(skIdx);

  // Victor says: I reverted the logic here back to an earlier version
  // ac0308715e5ae6bf5e750e8701e736d855550fc8
  // I don't see the reason for the change, and the logic here is
  // very delicate

  const DoubleCRT& sKey = sKeys.at(skIdx); // get key
  // Sample a new RLWE instance
  ctxt.noiseBound = RLWE(ctxt.parts[0], ctxt.parts[1], sKey, ptxtSpace);

  if (isCKKS()) {

    double f = getContext().getEA().getCx().encodeScalingFactor() / ptxtSize;
    long prec = getContext().getAlMod().getPPowR();
    long ef = NTL::conv<long>(ceil(prec * ctxt.noiseBound / (f * ptxtSize)));
    if (ef > 1) { // scale up some more
      ctxt.parts[0] += ptxt * ef;
      f *= ef;
    } else {
      ctxt.parts[0] += ptxt;
    }
    // Round size to next power of two so as not to leak too much
    ctxt.ptxtMag = EncryptedArrayCx::roundedSize(ptxtSize);
    ctxt.ratFactor = f;
    ctxt.noiseBound += ptxtSize * ctxt.ratFactor;
    // VJS-NOTE: the above noise calculation makes no sense to me
    return long(f);

  } else { // BGV

    // The logic here has changed to be identical
    // to that used in public key encryption

    // add in the plaintext
    // FIXME: we should really randomize ptxt, so that each coefficient
    //    has expected value 0
    // NOTE: This relies on the first part, ctxt[0], to have handle to 1

    // This code sequence could be optimized, but there is no point
    long QmodP = rem(context.productOfPrimes(ctxt.primeSet), ptxtSpace);
    NTL::ZZX ptxt_fixed;
    balanced_MulMod(ptxt_fixed, ptxt, QmodP, ptxtSpace);
    ctxt.parts[0] += ptxt_fixed;

    // NOTE: this is a heuristic, as the ptxt is not really random,
    // although, when ptxtSpace == 2, the balanced_MulMod will
    // randomize it
    double ptxt_bound = context.noiseBoundForMod(ptxtSpace, context.getPhiM());

    // FIXME: for now, we print out a warning, but we can consider
    // implementing a more robust randomization and rejection sampling
    // strategy.
    double ptxt_sz = NTL::conv<double>(
        embeddingLargestCoeff(ptxt_fixed, context.getZMStar()));

    if (ptxt_sz > ptxt_bound) {
      Warning("noise bound exceeded in encryption");
    }

    double ptxt_rat = ptxt_sz / ptxt_bound;
    HELIB_STATS_UPDATE("ptxt_rat_sk", ptxt_rat);

    ctxt.noiseBound += ptxt_bound;

    return ctxt.ptxtSpace;
  }
}

long SecKey::skEncrypt(Ctxt& ctxt,
                       const zzX& ptxt,
                       long ptxtSpace,
                       long skIdx) const
{
  NTL::ZZX tmp;
  convert(tmp, ptxt);
  return skEncrypt(ctxt, tmp, ptxtSpace, skIdx);
}
// These methods override the public-key Encrypt methods
long SecKey::Encrypt(Ctxt& ciphertxt,
                     const NTL::ZZX& plaintxt,
                     long ptxtSpace) const
{
  return skEncrypt(ciphertxt, plaintxt, ptxtSpace, /*skIdx=*/0);
}
long SecKey::Encrypt(Ctxt& ciphertxt, const zzX& plaintxt, long ptxtSpace) const
{
  return skEncrypt(ciphertxt, plaintxt, ptxtSpace, /*skIdx=*/0);
}

//=============== new EncodedPtxt interface ==================

void SecKey::Encrypt(Ctxt& ctxt, const EncodedPtxt& eptxt) const
{
  if (eptxt.isBGV())
    Encrypt(ctxt, eptxt.getBGV());
  else if (eptxt.isCKKS())
    Encrypt(ctxt, eptxt.getCKKS());
  else
    throw LogicError("Encrypt: bad EncodedPtxt");
}

void SecKey::Encrypt(Ctxt& ctxt, const EncodedPtxt_BGV& eptxt) const
{
  HELIB_TIMER_START;

  assertTrue(!isCKKS(), "Encrypt: mismatched BGV ptxt / CKKS ctxt");
  assertEq((const PubKey*)this, &ctxt.pubKey, "Encrypt: public key mismatch");
  assertEq(&context, &eptxt.getContext(), "Encrypt: context mismatch");

  long ptxtSpace = eptxt.getPtxtSpace();
  NTL::ZZX ptxt;

  convert(ptxt, eptxt.getPoly());

  long skIdx = 0; // in case we eventually want to generalize

  ctxt.ptxtSpace = ptxtSpace;
  ctxt.primeSet = context.getCtxtPrimes();
  ctxt.intFactor = 1;
  ctxt.ratFactor = ctxt.ptxtMag = 1.0;
  ctxt.parts.assign(2, CtxtPart(context, context.getCtxtPrimes()));

  // make parts[0],parts[1] point to (1,s)
  ctxt.parts[0].skHandle.setOne();
  ctxt.parts[1].skHandle.setBase(skIdx);

  // Sample a new RLWE instance
  const DoubleCRT& sKey = sKeys.at(skIdx);
  ctxt.noiseBound = RLWE(ctxt.parts[0], ctxt.parts[1], sKey, ptxtSpace);

  // The logic here has changed to be identical
  // to that used in public key encryption

  // add in the plaintext
  // FIXME: we should really randomize ptxt, so that each coefficient
  //    has expected value 0
  // NOTE: This relies on the first part, ctxt[0], to have handle to 1

  // This code sequence could be optimized, but there is no point
  long QmodP = rem(context.productOfPrimes(ctxt.primeSet), ptxtSpace);
  NTL::ZZX ptxt_fixed;
  balanced_MulMod(ptxt_fixed, ptxt, QmodP, ptxtSpace);
  ctxt.parts[0] += ptxt_fixed;

  // NOTE: this is a heuristic, as the ptxt is not really random,
  // although, when ptxtSpace == 2, the balanced_MulMod will
  // randomize it
  double ptxt_bound = context.noiseBoundForMod(ptxtSpace, context.getPhiM());

  // FIXME: for now, we print out a warning, but we can consider
  // implementing a more robust randomization and rejection sampling
  // strategy.
  double ptxt_sz =
      NTL::conv<double>(embeddingLargestCoeff(ptxt_fixed, context.getZMStar()));

  if (ptxt_sz > ptxt_bound) {
    Warning("noise bound exceeded in encryption");
  }

  double ptxt_rat = ptxt_sz / ptxt_bound;
  HELIB_STATS_UPDATE("ptxt_rat_sk", ptxt_rat);

  ctxt.noiseBound += ptxt_bound;
}

void SecKey::Encrypt(Ctxt& ctxt, const EncodedPtxt_CKKS& eptxt) const
{
  HELIB_TIMER_START;

  assertTrue(isCKKS(), "Encrypt: mismatched CKKS ptxt / BGV ctxt");
  assertEq((const PubKey*)this, &ctxt.pubKey, "Encrypt: public key mismatch");
  assertEq(&context, &eptxt.getContext(), "Encrypt: context mismatch");

  NTL::ZZX ptxt;
  convert(ptxt, eptxt.getPoly());
  double mag = eptxt.getMag();
  double scale = eptxt.getScale();
  double err = eptxt.getErr();

  long skIdx = 0; // in case we eventually want to generalize

  ctxt.parts.assign(2, CtxtPart(context, context.getCtxtPrimes()));

  // make parts[0],parts[1] point to (1,s)
  ctxt.parts[0].skHandle.setOne();
  ctxt.parts[1].skHandle.setBase(skIdx);

  // Sample a new RLWE instance
  const DoubleCRT& sKey = sKeys.at(skIdx);
  double error_bound = RLWE(ctxt.parts[0], ctxt.parts[1], sKey, 1);

  // This follows the same logic in PubKey::Encrypt(EncodedPtxt_CKKS).
  // See documentation there
  long ef = NTL::conv<long>(ceil(error_bound / err));

  if (ef > 1) { // scale up some more
    ctxt.parts[0] += ptxt * ef;
    scale *= ef;
    err *= ef;
  } else { // no need for extra scaling
    ctxt.parts[0] += ptxt;
  }

  // VJS-NOTE: we no longer round to the next power of two:
  // Then encoding routine should take care of setting mag correctly.
  ctxt.primeSet = context.getCtxtPrimes();
  ctxt.ptxtMag = mag;
  ctxt.ratFactor = scale;
  ctxt.noiseBound = error_bound + err;
  ctxt.ptxtSpace = 1;
  ctxt.intFactor = 1;
}

//============================================================

// Generate bootstrapping data if needed, returns index of key
long SecKey::genRecryptData()
{
  if (recryptKeyID >= 0)
    return recryptKeyID;

  // Make sure that the context has the bootstrapping EA and PAlgMod
  assertTrue(context.isBootstrappable(),
             "Cannot generate recrypt data for non-bootstrappable context");

  long p2ePr = context.getRcData().alMod->getPPowR(); // p^{e-e'+r}
  long p2r = context.getAlMod().getPPowR();           // p^r

  // Generate a new bootstrapping key
  zzX keyPoly;
  long hwt = context.getRcData().skHwt;
  double bound = sampleHWtBounded(keyPoly, context, hwt);

  DoubleCRT newSk(keyPoly,
                  context,
                  context.getCtxtPrimes() | context.getSpecialPrimes());
  // defined relative to all primes

  long keyID = ImportSecKey(newSk, bound, p2r, /*maxDegKswitch=*/1);

  // Generate a key-switching matrix from key 0 to this key
  GenKeySWmatrix(/*fromSPower=*/1,
                 /*fromXPower=*/1,
                 /*fromIdx=*/0,
                 /*toIdx=*/keyID,
                 /*ptxtSpace=*/p2r);

  // Encrypt new key under key #0 and plaintext space p^{e+r}
  Encrypt(recryptEkey, keyPoly, p2ePr);

  return (recryptKeyID = keyID); // return the new key-ID
}

std::ostream& operator<<(std::ostream& str, const SecKey& sk)
{
  sk.writeToJSON(str);
  return str;
}

// FIXME: For consistency we should change this to write in json format too.
std::ostream& SecKey::writeSecKeyDerivedASCII(std::ostream& str) const
{
  str << "[" << sKeys.size() << std::endl;
  for (long i = 0; i < (long)sKeys.size(); i++)
    str << sKeys[i] << std::endl;
  return str << "]";
}

std::istream& operator>>(std::istream& str, SecKey& sk)
{
  sk.readJSON(str);
  return str;
}

void SecKey::writeTo(std::ostream& str, bool sk_only) const
{
  SerializeHeader<SecKey>().writeTo(str);
  writeEyeCatcher(str, EyeCatcher::SK_BEGIN);

  if (!sk_only) {
    // Write out the public key part first.
    this->PubKey::writeTo(str);
  } else {
    // If only writing sk, just write the context.
    this->getContext().writeTo(str);
  }

  // Write out vector<DoubleCRT> sKeys
  write_raw_vector<DoubleCRT>(str, this->sKeys);

  writeEyeCatcher(str, EyeCatcher::SK_END);
}

SecKey SecKey::readFrom(std::istream& str, const Context& context, bool sk_only)
{
  const auto header = SerializeHeader<SecKey>::readFrom(str);
  assertEq<IOError>(header.version,
                    Binio::VERSION_0_0_1_0,
                    "Header: version " + header.versionString() +
                        " not supported");

  bool eyeCatcherFound = readEyeCatcher(str, EyeCatcher::SK_BEGIN);
  assertTrue<IOError>(eyeCatcherFound,
                      "Could not find pre-secret key eyecatcher");
  if (sk_only) {
    // there should be a context written at this point in the file, check it
    // matches provided context
    assertEq(context, Context::readFrom(str), "Context mismatch");
  }
  // now construct a secret key. If public key is written, construct from pk,
  // otherwise, construct from context
  SecKey ret =
      sk_only ? SecKey(context) : SecKey(PubKey::readFrom(str, context));

  // Set the secret part of the secret key.
  ret.sKeys = read_raw_vector<DoubleCRT>(str, context);

  eyeCatcherFound = readEyeCatcher(str, EyeCatcher::SK_END);
  assertTrue<IOError>(eyeCatcherFound,
                      "Could not find post-secret key eyecatcher");

  return ret;
}

void SecKey::writeToJSON(std::ostream& str, bool sk_only) const
{
  executeRedirectJsonError<void>([&]() { str << writeToJSON(sk_only); });
}

JsonWrapper SecKey::writeToJSON(bool sk_only) const
{
  auto body = [this, sk_only]() {
    if (sk_only) {
      // json contains context and secret key(s)
      json j = {{"context", unwrap(this->getContext().writeToJSON())},
                {"sKeys", writeVectorToJSON(this->sKeys)}};
      return wrap(toTypedJson<SecKey>(j));
    } else {
      // json contains public key and secret key(s)
      json j = {{"PubKey", unwrap(this->PubKey::writeToJSON())},
                {"sKeys", writeVectorToJSON(this->sKeys)}};
      return wrap(toTypedJson<SecKey>(j));
    }
  };
  return executeRedirectJsonError<JsonWrapper>(body);
}

SecKey SecKey::readFromJSON(std::istream& str,
                            const Context& context,
                            bool sk_only)
{
  auto body = [&]() {
    json j;
    str >> j;
    return SecKey::readFromJSON(wrap(j), context, sk_only);
  };
  return executeRedirectJsonError<SecKey>(body);
}

SecKey SecKey::readFromJSON(const JsonWrapper& jw,
                            const Context& context,
                            bool sk_only)
{
  SecKey ret{context};
  ret.readJSON(jw, sk_only);
  return ret;
}

void SecKey::readJSON(std::istream& str, bool sk_only)
{
  executeRedirectJsonError<void>([&]() {
    json j;
    str >> j;
    this->readJSON(wrap(j), sk_only);
  });
}

void SecKey::readJSON(const JsonWrapper& tjw, bool sk_only)
{
  executeRedirectJsonError<void>([&]() {
    json j = fromTypedJson<SecKey>(unwrap(tjw));
    this->clear();
    if (sk_only)
      assertEq(this->getContext(),
               Context::readFromJSON(wrap(j.at("context"))),
               "Context mismatch");
    else
      this->PubKey::readJSON(wrap(j.at("PubKey")));

    this->sKeys = readVectorFromJSON<DoubleCRT>(j.at("sKeys"), context);
  });
}

} // namespace helib
