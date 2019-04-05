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
#include <queue> // used in the breadth-first search in setKeySwitchMap
#include "FHE.h"

#include "timing.h"
#include "binio.h"
#include "sample.h"
#include "EncryptedArray.h"

NTL_CLIENT

/******** Utility function to generate RLWE instances *********/

// Assumes that c1 is already chosen by the caller
double RLWE1(DoubleCRT& c0, const DoubleCRT& c1, const DoubleCRT &s, long p)
// Returns a high-probabiliy bound on the L-infty norm
// of the canonical embedding of the decryption of (c0, c1) w/r/to s
{
  //OLD: assert (p>0); // Used with p=1 for CKKS, p>=2 for BGV
  helib::assertTrue<helib::InvalidArgument>(p>0, "Cannot generate RLWE instance with nonpositive p"); // Used with p=1 for CKKS, p>=2 for BGV
  const FHEcontext& context = s.getContext();
  const PAlgebra& palg = context.zMStar;

  // choose a short error e
  double stdev = to_double(context.stdev);
  if (palg.getPow2() == 0) // not power of two
    stdev *= sqrt(palg.getM());
  double bound = c0.sampleGaussianBounded(stdev);

  // Set c0 =  p*e - c1*s.
  // It is assumed that c0,c1 are defined with respect to the same set of
  // primes, but s may be defined relative to a different set. Either way
  // the primes for of c0,c1 are unchanged.
  if (p>1) {
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
double RLWE(DoubleCRT& c0,DoubleCRT& c1, const DoubleCRT &s, long p,
            ZZ* prgSeed)
{
  // choose c1 at random (using prgSeed if not NULL)
  c1.randomize(prgSeed);
  return RLWE1(c0, c1, s, p);
}

/******************** KeySwitch implementation **********************/
/********************************************************************/

bool KeySwitch::operator==(const KeySwitch& other) const
{
  if (this == &other) return true;

  if (fromKey != other.fromKey) return false;
  if (toKeyID != other.toKeyID) return false;
  if (ptxtSpace != other.ptxtSpace) return false;

  if (prgSeed != other.prgSeed) return false;

  if (b.size() != other.b.size()) return false;
  for (size_t i=0; i<b.size(); i++) if (b[i] != other.b[i]) return false;

  return true;
}


void KeySwitch::verify(FHESecKey& sk) 
{
  long fromSPower = fromKey.getPowerOfS();
  long fromXPower = fromKey.getPowerOfX();
  long fromIdx = fromKey.getSecretKeyID(); 
  long toIdx = toKeyID;
  long p = ptxtSpace;
  long n = b.size();

  cout << "KeySwitch::verify\n";
  cout << "fromS = " << fromSPower 
       << " fromX = " << fromXPower 
       << " fromIdx = " << fromIdx 
       << " toIdx = " << toIdx 
       << " p = " << p 
       << " n = " << n 
       << "\n";


  if (fromSPower != 1 || fromXPower != 1 || (fromIdx == toIdx) || n == 0) {
    cout << "KeySwitch::verify: these parameters not checkable\n";
    return;
  }

  const FHEcontext& context = b[0].getContext();

  // we don't store the context in the ks matrix, so let's
  // check that they are consistent

  for (long i = 0; i < n; i++) {
    if (&context != &(b[i].getContext()))
      cout << "KeySwitch::verify: bad context " << i << "\n";
  }

  cout << "context.ctxtPrimes = " << context.ctxtPrimes << "\n";
  cout << "context.specialPrimes = " << context.specialPrimes << "\n";
  IndexSet fullPrimes = context.fullPrimes(); // ctxtPrimes | specialPrimes;

  cout << "digits: ";
  for (long i = 0; i < n; i++) 
    cout << context.digits[i] << " ";
  cout << "\n";

  cout << "IndexSets of b: ";
  for (long i = 0; i < n; i++) 
    cout << b[i].getMap().getIndexSet() << " ";
  cout << "\n";

  // VJS: suspicious shadowing of fromKey, toKey
  const DoubleCRT& _fromKey = sk.sKeys.at(fromIdx);
  const DoubleCRT& _toKey = sk.sKeys.at(toIdx);

  cout << "IndexSet of fromKey: " << _fromKey.getMap().getIndexSet() << "\n";
  cout << "IndexSet of toKey: " << _toKey.getMap().getIndexSet() << "\n";

  vector<DoubleCRT> a;
  a.resize(n, DoubleCRT(context, fullPrimes)); // defined modulo all primes

  { RandomState state;

    SetSeed(prgSeed);
    for (long i = 0; i < n; i++)
      a[i].randomize();

  } // the RandomState destructor "restores the state" (see NumbTh.h)

  vector<ZZX> A, B;

  A.resize(n);
  B.resize(n);

  for (long i = 0; i < n; i++) {
    a[i].toPoly(A[i]);
    b[i].toPoly(B[i]);
  }

  ZZX FromKey, ToKey;
  _fromKey.toPoly(FromKey, fullPrimes);
  _toKey.toPoly(ToKey, fullPrimes);

  ZZ Q = context.productOfPrimes(fullPrimes);
  ZZ prod = context.productOfPrimes(context.specialPrimes);
  ZZX C, D;
  ZZX PhimX = context.zMStar.getPhimX();

  long nb = 0;
  for (long i = 0; i < n; i++) {
    C = (B[i] - FromKey*prod + ToKey*A[i]) % PhimX;
    PolyRed(C, Q);
    if (!divide(D, C, p)) {
      cout << "*** not divisible by p at " << i << "\n";
    }
    else {
      for (long j = 0; j <= deg(D); j++)
         if (NumBits(coeff(D, j)) > nb) nb = NumBits(coeff(D, j));
    }
    prod *= context.productOfPrimes(context.digits[i]);
  }

  cout << "error ratio: " << ((double) nb)/((double) NumBits(Q)) << "\n";
}

const KeySwitch& KeySwitch::dummy()
{
  static const KeySwitch dummy(-1,-1,-1,-1);
  return dummy;
}

ostream& operator<<(ostream& str, const KeySwitch& matrix)
{
  str << "["<<matrix.fromKey  <<" "<<matrix.toKeyID
      << " "<<matrix.ptxtSpace<<" "<<matrix.b.size() << endl;
  for (long i=0; i<(long)matrix.b.size(); i++)
    str << matrix.b[i] << endl;
  str << matrix.prgSeed << " " << matrix.noiseBound << "]";
  return str;
}

// Used in lieu of istream& operator>>(istream& str, KeySwitch& matrix)
void KeySwitch::readMatrix(istream& str, const FHEcontext& context)
{
  seekPastChar(str,'['); // defined in NumbTh.cpp
  str >> fromKey;
  str >> toKeyID;
  str >> ptxtSpace;

  long nDigits;
  str >> nDigits;
  b.resize(nDigits, DoubleCRT(context, IndexSet::emptySet()));
  for (long i=0; i<nDigits; i++)
    str >> b[i];
  str >> prgSeed;
  str >> noiseBound;
  seekPastChar(str,']');
}


void KeySwitch::write(ostream& str) const
{
  writeEyeCatcher(str, BINIO_EYE_SKM_BEGIN);
/*  
    Write out raw
    1. SKHandle fromKey; 
    2. long     toKeyID;
    3. long     ptxtSpace;
    4. vector<DoubleCRT> b;
    5. ZZ prgSeed;
    6. xdouble noiseBound;
*/

  fromKey.write(str);
  write_raw_int(str, toKeyID);
  write_raw_int(str, ptxtSpace);

  write_raw_vector(str, b);
  
  write_raw_ZZ(str, prgSeed);
  write_raw_xdouble(str, noiseBound);

  writeEyeCatcher(str, BINIO_EYE_SKM_END);
}

void KeySwitch::read(istream& str, const FHEcontext& context)
{
  int eyeCatcherFound = readEyeCatcher(str, BINIO_EYE_SKM_BEGIN);
  //OLD: assert(eyeCatcherFound == 0);
  helib::assertEq(eyeCatcherFound, 0, "Could not find pre-secret key eyecatcher");

  fromKey.read(str);
  toKeyID = read_raw_int(str);
  ptxtSpace = read_raw_int(str);
  DoubleCRT blankDCRT(context, IndexSet::emptySet());
  read_raw_vector(str, b, blankDCRT);
  read_raw_ZZ(str, prgSeed);
  noiseBound = read_raw_xdouble(str); 

  eyeCatcherFound = readEyeCatcher(str, BINIO_EYE_SKM_END);
  //OLD: assert(eyeCatcherFound == 0);
  helib::assertEq(eyeCatcherFound, 0, "Could not find post-secret key eyecatcher");
}



/******************** FHEPubKey implementation **********************/
/********************************************************************/
// Computes the keySwitchMap pointers, using breadth-first search (BFS)

void FHEPubKey::setKeySwitchMap(long keyId)
{
  //OLD: assert(keyId>=0 && keyId<(long)skBounds.size()); // Sanity-check, do we have such a key?
  helib::assertInRange(keyId, 0l, (long)skBounds.size(), "No such key found"); // Sanity-check, do we have such a key?
  long m = context.zMStar.getM();

  // Initialize an aray of "edges" (this is easier than searching through
  // all the matrices for every step). This is a list of all the powers n
  // for which we have a matrix W[s_i(X^n) => s_i(X)], as well as the index
  // of that matrix in the keySwitching array.
  typedef pair<long,long> keySwitchingEdge;
  vector<keySwitchingEdge> edges;
  for (long i=0; i<(long)keySwitching.size(); i++) {
    const KeySwitch& mat = keySwitching.at(i);
    if (mat.toKeyID == keyId && mat.fromKey.getPowerOfS()==1
                             && mat.fromKey.getSecretKeyID()== keyId)
      edges.push_back(keySwitchingEdge(mat.fromKey.getPowerOfX(), i));
  }
  if (keyId>=(long)keySwitchMap.size()) // allocate more space if needed
    keySwitchMap.resize(keyId+1);

  // initialize keySwitchMap[keyId] with m empty entries (with -1 in them)
  keySwitchMap.at(keyId).assign(m,-1);

  // A standard BFS implementation using a FIFO queue (complexity O(V+E))

  std::queue<long> bfsQueue; 
  bfsQueue.push(1);          // Push the target node 1 onto the BFS queue
  while (!bfsQueue.empty()) {
    long currentNode = bfsQueue.front();

    // See what other nodes can reach the current one
    for (long j=0; j<(long)edges.size(); j++) { // go over the edges
      long n = edges[j].first;
      long matrixIndex = edges[j].second;

      long nextNode = MulMod(currentNode, n, m);
      if (keySwitchMap.at(keyId).at(nextNode) == -1) {// A new node: mark it now
	// Record the index of the matrix that we use for the first step
	keySwitchMap[keyId][nextNode] = matrixIndex;

	bfsQueue.push(nextNode);    // push new node onto BFS queue
      }
    }
    bfsQueue.pop();                 // We are done with the current node
  }
}

const KeySwitch& FHEPubKey::getKeySWmatrix(const SKHandle& from, 
					   long toIdx) const
{
  // First try to use the keySwitchMap
  if (from.getPowerOfS()==1 && from.getSecretKeyID()==toIdx 
                            && toIdx < (long)keySwitchMap.size()) {
    long matIdx = keySwitchMap.at(toIdx).at(from.getPowerOfX());
    if (matIdx>=0) { 
      const KeySwitch& matrix = keySwitching.at(matIdx);
      if (matrix.fromKey == from) return matrix;
    }
  }

  // Otherwise resort to linear search
  for (size_t i=0; i<keySwitching.size(); i++) {
    if (keySwitching[i].toKeyID==toIdx && keySwitching[i].fromKey==from)
      return keySwitching[i];
  }
  return KeySwitch::dummy(); // return this if nothing is found
}

const KeySwitch& FHEPubKey::getAnyKeySWmatrix(const SKHandle& from) const
{
  // First try to use the keySwitchMap
  if (from.getPowerOfS()==1 && 
      from.getSecretKeyID() < (long)keySwitchMap.size()) {
    long matIdx = keySwitchMap.at(from.getSecretKeyID()).at(from.getPowerOfX());
    if (matIdx>=0) {
      const KeySwitch& matrix = keySwitching.at(matIdx);
      if (matrix.fromKey == from) return matrix;
    }
  }

  // Otherwise resort to linear search
  for (size_t i=0; i<keySwitching.size(); i++) {
    if (keySwitching[i].fromKey==from) return keySwitching[i];
  }
  return KeySwitch::dummy(); // return this if nothing is found
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
long FHEPubKey::Encrypt(Ctxt &ctxt, const ZZX& ptxt, long ptxtSpace,
			bool highNoise) const
{
  FHE_TIMER_START;
  // NOTE: isCKKS() checks the tag in the alMod  the context
  if (isCKKS()) {
    double pSize = (ptxtSpace <= 0)? 1.0 : double(ptxtSpace);
    // For CKKSencrypt, ptxtSpace==1 is the defaults size value
    CKKSencrypt(ctxt, ptxt, pSize); // FIXME: handle highNoise in CKKSencrypt
    return ptxtSpace;
  }

  //OLD: assert(this == &ctxt.pubKey);
  helib::assertEq(this, &ctxt.pubKey, "Public key and context public key mismatch");
  if (ptxtSpace != pubEncrKey.ptxtSpace) { // plaintext-space mistamtch
    ptxtSpace = GCD(ptxtSpace, pubEncrKey.ptxtSpace);
    if (ptxtSpace <= 1) throw helib::RuntimeError("Plaintext-space mismatch on encryption");
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
  
  DoubleCRT e(context, context.ctxtPrimes);
  DoubleCRT r(context, context.ctxtPrimes);
  double r_bound = r.sampleSmall();

  ctxt.noiseBound += r_bound * pubEncrKey.noiseBound;

  //cerr << "*** r_bound*pubEncrKey.noiseBound " << r_bound * pubEncrKey.noiseBound << "\n";

  double stdev = to_double(context.stdev);
  if (context.zMStar.getPow2()==0) // not power of two
    stdev *= sqrt(context.zMStar.getM());

  for (size_t i=0; i<ctxt.parts.size(); i++) {  // add noise to all the parts
    ctxt.parts[i] *= r;
    xdouble e_bound;

    if (highNoise && i == 0) {
      // we sample e so that coefficients are uniform over 
      // [-Q/(8*ptxtSpace)..Q/(8*ptxtSpace)]

      ZZ B;
      B = context.productOfPrimes(context.ctxtPrimes);
      B /= (ptxtSpace*8);
 
      e_bound = e.sampleUniform(B);
    }
    else { 
      e_bound = e.sampleGaussian(stdev);
    }

    e *= ptxtSpace;
    e_bound *= ptxtSpace;

    if (i == 1) {
      e_bound *= getSKeyBound(ctxt.parts[i].skHandle.getSecretKeyID());
    }

    ctxt.parts[i] += e;
    ctxt.noiseBound += e_bound;

    //cerr << "*** e_bound " << e_bound << "\n";
  }

  // add in the plaintext
  // FIXME: we should really randomize ptxt, so that each coefficient
  //    has expected value 0
  // NOTE: This relies on the first part, ctxt[0], to have handle to 1
  
  if (ptxtSpace==2) {
    ctxt.parts[0] += ptxt;
  }
  else { // The general case of ptxtSpace>2: for a ciphertext
         // relative to modulus Q, we add ptxt * Q mod ptxtSpace.
    long QmodP = rem(context.productOfPrimes(ctxt.primeSet), ptxtSpace);
    ctxt.parts[0] += MulMod(ptxt,QmodP,ptxtSpace); // MulMod from module NumbTh
  }

  // NOTE: this is a heuristic
  double ptxt_bound = context.noiseBoundForUniform(double(ptxtSpace)/2.0, context.zMStar.getPhiM());

  ctxt.noiseBound += ptxt_bound;

  //cerr << "*** ptxt_bound " << ptxt_bound << "\n";

  // fill in the other ciphertext data members
  ctxt.ptxtSpace = ptxtSpace;
  ctxt.intFactor = 1;

  //cerr << "*** ctxt.noiseBound " << ctxt.noiseBound << "\n";

  return ptxtSpace;
}

// FIXME: Some code duplication between here and Encrypt above
void FHEPubKey::CKKSencrypt(Ctxt &ctxt, const ZZX& ptxt,
                            double ptxtSize, double scaling) const
{
  //OLD: assert(this == &ctxt.pubKey);
  helib::assertEq(this, &ctxt.pubKey, "Public key and context public key mismatch");

  if (ptxtSize<=0)
    ptxtSize = 1.0;
  if (scaling <= 0) // assume the default scaling factor
    scaling = getContext().ea->getCx().encodeScalingFactor() / ptxtSize;

  long m = context.zMStar.getM();
  long prec = getContext().alMod.getPPowR();

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
  // We also have ptxt_bound = ef*f*ptxtSize, which is tracked sparately.
  //
  // The input ptxt is already scaled by a factof f=scaling, and is being
  // further scaled by the extra factor ef, so ef*f is the new scaling
  // factor. The extra factor ef is set as ceil(error_bound*prec/f),
  // so that we have ef*f >= error_bound*prec.  

  DoubleCRT e(context, context.ctxtPrimes);
  DoubleCRT r(context, context.ctxtPrimes);

  double r_bound = r.sampleSmall(); // r is a {0,+-1} polynomial
  xdouble error_bound = r_bound * pubEncrKey.noiseBound;

  double stdev = to_double(context.stdev);
  if (context.zMStar.getPow2()==0) // not power of two
    stdev *= sqrt(m);

  for (size_t i=0; i<ctxt.parts.size(); i++) {  // add noise to all the parts
    ctxt.parts[i] *= r;
    double e_bound = e.sampleGaussian(stdev);// zero-mean Gaussian, sigma=stdev
    ctxt.parts[i] += e;
    if (i == 1) {
      e_bound *= getSKeyBound(ctxt.parts[i].skHandle.getSecretKeyID());
    }
    error_bound += e_bound;
  }
  // Compute the extra scaling factor, if needed
  long ef = conv<long>(ceil(error_bound*prec/(scaling*ptxtSize)));
  if (ef > 1) { // scale up some more
    ctxt.parts[0] += ptxt*ef;
    scaling *= ef;
  }
  else { // no need for extra scaling
    ctxt.parts[0] += ptxt;
  }
  // Round size to next power of two so as not to leak too much
  ctxt.ptxtMag = EncryptedArrayCx::roundedSize(ptxtSize);
  ctxt.ratFactor = scaling;
  ctxt.noiseBound = error_bound;
  ctxt.ptxtSpace = 1;
}

bool FHEPubKey::operator==(const FHEPubKey& other) const
{
  if (this == &other) return true;

  if (&context != &other.context) return false;
  if (!pubEncrKey.equalsTo(other.pubEncrKey, /*comparePkeys=*/false))
    return false;

  if (skBounds.size() != other.skBounds.size()) return false;
  for (size_t i=0; i<skBounds.size(); i++)
    if (fabs(skBounds[i]-other.skBounds[i])>0.1) return false;

  if (keySwitching.size() != other.keySwitching.size()) return false;
  for (size_t i=0; i<keySwitching.size(); i++)
    if (keySwitching[i] != other.keySwitching[i]) return false;

  if (keySwitchMap.size() != other.keySwitchMap.size()) return false;
  for (size_t i=0; i<keySwitchMap.size(); i++) {
    if (keySwitchMap[i].size() != other.keySwitchMap[i].size()) return false;
    for (size_t j=0; j<keySwitchMap[i].size(); j++)
      if (keySwitchMap[i][j] != other.keySwitchMap[i][j]) return false;
  }

  // compare KS_strategy, ignoring trailing FHE_KSS_UNKNOWN
  long n = KS_strategy.length();
  while (n > 0 && KS_strategy[n-1] == FHE_KSS_UNKNOWN) n--;
  long n1 = other.KS_strategy.length();
  while (n1 > 0 && other.KS_strategy[n1-1] == FHE_KSS_UNKNOWN) n1--;
  if (n != n1) return false;
  for (long i: range(n)) {
    if (KS_strategy[i] != other.KS_strategy[i]) return false;
  }

  if (recryptKeyID!=other.recryptKeyID) return false;
  if (recryptKeyID>=0 && 
      !recryptEkey.equalsTo(other.recryptEkey, /*comparePkeys=*/false))
    return false;

  return true;
}


ostream& operator<<(ostream& str, const FHEPubKey& pk)
{
  str << "[";
  writeContextBase(str, pk.getContext());

  // output the public encryption key itself
  str << pk.pubEncrKey << endl;

  // output skBounds in the same format as vec_double
  str << "[";
  for (long i=0; i<(long)pk.skBounds.size(); i++)
    str << pk.skBounds[i]<<" ";
  str << "]\n";

  // output the key-switching matrices
  str << pk.keySwitching.size() << endl;
  for (long i=0; i<(long)pk.keySwitching.size(); i++)
    str << pk.keySwitching[i] << endl;

  // output keySwitchMap in the same format as vec_vec_long
  str << "[";
  for (long i=0; i<(long)pk.keySwitchMap.size(); i++) {
    str << "[";
    for (long j=0; j<(long)pk.keySwitchMap[i].size(); j++)
      str << pk.keySwitchMap[i][j] << " ";
    str << "]\n ";
  }
  str << "]\n";

  str << pk.KS_strategy << "\n";

  // output the bootstrapping key, if any
  str << pk.recryptKeyID << " ";
  if (pk.recryptKeyID>=0) str << pk.recryptEkey << endl;
  return str << "]";
}

istream& operator>>(istream& str, FHEPubKey& pk)
{
  pk.clear();
  //  cerr << "FHEPubKey[";
  seekPastChar(str, '['); // defined in NumbTh.cpp

  // sanity check, verify that basic context parameters are correct
  unsigned long m, p, r;
  vector<long> gens, ords;
  readContextBase(str, m, p, r, gens, ords);
  //OLD: assert(comparePAlgebra(pk.getContext().zMStar, m, p, r, gens, ords));
  helib::assertTrue(comparePAlgebra(pk.getContext().zMStar, m, p, r, gens, ords), "PAlgebra mismatch");

  // Get the public encryption key itself
  str >> pk.pubEncrKey;

  // Get the vector of secret-key Hamming-weights
  NTL::Vec<double> vl;
  str >> vl;
  pk.skBounds.resize(vl.length());
  for (long i=0; i<(long)pk.skBounds.size(); i++) pk.skBounds[i] = vl[i];

  // Get the key-switching matrices
  long nMatrices;
  str >> nMatrices;
  pk.keySwitching.resize(nMatrices);
  for (long i=0; i<nMatrices; i++)  // read the matrix from input str
    pk.keySwitching[i].readMatrix(str, pk.getContext());

  // Get the key-switching map
  Vec< Vec<long> > vvl;
  str >> vvl;
  pk.keySwitchMap.resize(vvl.length());
  for (long i=0; i<(long)pk.keySwitchMap.size(); i++) {
    pk.keySwitchMap[i].resize(vvl[i].length());
    for (long j=0; j<(long)pk.keySwitchMap[i].size(); j++)
      pk.keySwitchMap[i][j] = vvl[i][j];
  }

  // build the key-switching map for all keys
  for (long i=pk.skBounds.size()-1; i>=0; i--)
    pk.setKeySwitchMap(i);

  str >> pk.KS_strategy; 

  // Get the bootstrapping key, if any
  str >> pk.recryptKeyID;
  if (pk.recryptKeyID>=0) str >> pk.recryptEkey;

  seekPastChar(str, ']');
  return str;
}
      
void writePubKeyBinary(ostream& str, const FHEPubKey& pk) 
{

  writeEyeCatcher(str, BINIO_EYE_PK_BEGIN);  

// Write out for FHEPubKey
//  1. Context Base 
//  2. Ctxt pubEncrKey;
//  3. vector<long> skBounds;
//  4. vector<KeySwitch> keySwitching;
//  5. vector< vector<long> > keySwitchMap;
//  6. Vec<long> KS_strategy
//  7. long recryptKeyID; 
//  8. Ctxt recryptEkey;

  writeContextBaseBinary(str, pk.getContext());
  pk.pubEncrKey.write(str);
  write_raw_vector(str, pk.skBounds);

  // Keyswitch Matrices
  write_raw_vector(str, pk.keySwitching);

  long sz = pk.keySwitchMap.size();
  write_raw_int(str, sz);
  for(auto v: pk.keySwitchMap)
    write_raw_vector(str, v);

  write_ntl_vec_long(str, pk.KS_strategy); 

  write_raw_int(str, pk.recryptKeyID);
  pk.recryptEkey.write(str);

  writeEyeCatcher(str, BINIO_EYE_PK_END);
}

void readPubKeyBinary(istream& str, FHEPubKey& pk)
{
  int eyeCatcherFound = readEyeCatcher(str, BINIO_EYE_PK_BEGIN);
  //OLD: assert(eyeCatcherFound == 0);
  helib::assertEq(eyeCatcherFound, 0, "Could not find pre-public key eyecatcher");
 
  //  // TODO code to check context object is what it should be 
  //  // same as the text IO. May be worth putting it in helper func.
  //  std::unique_ptr<FHEcontext> dummy = buildContextFromBinary(str);
  unsigned long m, p, r;
  vector<long> gens, ords;
  readContextBaseBinary(str, m, p, r, gens, ords);
  //OLD: assert(comparePAlgebra(pk.getContext().zMStar, m, p, r, gens, ords));
  helib::assertTrue(comparePAlgebra(pk.getContext().zMStar, m, p, r, gens, ords), "PAlgebra mismatch");

  // Read in the rest
  pk.pubEncrKey.read(str);
  read_raw_vector(str, pk.skBounds);

  // Keyswitch Matrices
  read_raw_vector(str, pk.keySwitching, pk.getContext());

  long sz = read_raw_int(str);
  pk.keySwitchMap.clear();
  pk.keySwitchMap.resize(sz);
  for(auto& v: pk.keySwitchMap)
    read_raw_vector(str, v);

  read_ntl_vec_long(str, pk.KS_strategy); 

  pk.recryptKeyID = read_raw_int(str);
  pk.recryptEkey.read(str);

  eyeCatcherFound = readEyeCatcher(str, BINIO_EYE_PK_END);
  //OLD: assert(eyeCatcherFound == 0);
  helib::assertEq(eyeCatcherFound, 0, "Could not find post-public key eyecatcher");
}


/******************** FHESecKey implementation **********************/
/********************************************************************/

bool FHESecKey::operator==(const FHESecKey& other) const
{
  if (this == &other) return true;

  if (((const FHEPubKey&)*this)!=((const FHEPubKey&)other)) return false;
  if (sKeys.size() != other.sKeys.size()) return false;
  for (size_t i=0; i<sKeys.size(); i++)
    if (sKeys[i] != other.sKeys[i]) return false;
  return true;
}

// We allow the calling application to choose a secret-key polynomial by
// itself, then insert it into the FHESecKey object, getting the index of
// that secret key in the sKeys list. If this is the first secret-key for this
// FHESecKey object, then the procedure below generates a corresponding public
// encryption key.
// It is assumed that the context already contains all parameters.
long FHESecKey::ImportSecKey(const DoubleCRT& sKey, double bound,
			     long ptxtSpace, long maxDegKswitch)
{
  if (sKeys.empty()) { // 1st secret-key, generate corresponding public key
    if (ptxtSpace<2)
      ptxtSpace = isCKKS()? 1 : context.alMod.getPPowR();
    // default plaintext space is p^r for BGV, 1 for CKKS

    // allocate space, the parts are DoubleCRTs with all the ctxtPrimes
    pubEncrKey.parts.assign(2,CtxtPart(context,context.ctxtPrimes));
    // Choose a new RLWE instance
    pubEncrKey.noiseBound
      = RLWE(pubEncrKey.parts[0], pubEncrKey.parts[1], sKey, ptxtSpace);
    if (isCKKS()) {
      pubEncrKey.ptxtMag = 0.0;
      pubEncrKey.ratFactor = pubEncrKey.noiseBound
                           * getContext().ea->getCx().encodeScalingFactor();
    }

    // make parts[0],parts[1] point to (1,s)
    pubEncrKey.parts[0].skHandle.setOne();
    pubEncrKey.parts[1].skHandle.setBase();

    // Set the other Ctxt bookeeping parameters in pubEncrKey
    pubEncrKey.primeSet = context.ctxtPrimes;
    pubEncrKey.ptxtSpace = ptxtSpace;
  }
  skBounds.push_back(bound); // record the size of the new secret-key
  sKeys.push_back(sKey);     // add to the list of secret keys
  long keyID = sKeys.size()-1; // FIXME: not thread-safe, do we need to fix it?

  for (long e=2; e<=maxDegKswitch; e++)
    GenKeySWmatrix(e,1,keyID,keyID); // s^e -> s matrix

  return keyID; // return the index where this key is stored
}

// Generate a key-switching matrix and store it in the public key.
// The argument p denotes the plaintext space
void FHESecKey::GenKeySWmatrix(long fromSPower, long fromXPower,
			       long fromIdx, long toIdx, long p)
{
  FHE_TIMER_START;

  // sanity checks
  if (fromSPower<=0 || fromXPower<=0) return;  
  if (fromSPower==1 && fromXPower==1 && fromIdx==toIdx) return;

  // See if this key-switching matrix already exists in our list
  if (haveKeySWmatrix(fromSPower, fromXPower, fromIdx, toIdx))
    return; // nothing to do here

  DoubleCRT fromKey = sKeys.at(fromIdx); // copy object, not a reference
  const DoubleCRT& toKey = sKeys.at(toIdx);   // this can be a reference

  if (fromXPower>1) fromKey.automorph(fromXPower); // compute s(X^t)
  if (fromSPower>1) fromKey.Exp(fromSPower);       // compute s^r(X^t)
  // SHAI: The above lines compute the automorphism and exponentiation mod q,
  //   turns out this is really what we want (even through usually we think
  //   of the secret key as being mod p^r)

  KeySwitch ksMatrix(fromSPower,fromXPower,fromIdx,toIdx);
  RandomBits(ksMatrix.prgSeed, 256); // a random 256-bit seed

  long n = context.digits.size();

  // size-n vector
  ksMatrix.b.resize(n, DoubleCRT(context, context.ctxtPrimes | context.specialPrimes)); 

  vector<DoubleCRT> a; 
  a.resize(n, DoubleCRT(context, context.ctxtPrimes | context.specialPrimes));

  { RandomState state;
    SetSeed(ksMatrix.prgSeed);
    for (long i = 0; i < n; i++) 
      a[i].randomize();
  } // restore state upon destruction of state

  // Record the plaintext space for this key-switching matrix
  if (isCKKS()) p = 1;
  else {        // BGV
    if (p<2) {
      if (context.isBootstrappable()) { 
        // use larger bootstrapping plaintext space
        p = context.rcData.alMod->getPPowR();
      }
      else {
        p = pubEncrKey.ptxtSpace; // default plaintext space from public key
      }
    }
    // FIXME: We use context.isBootstrappable() rather than
    //   this->isBootstrappable(). So we get the larger bootstrapping
    //   plaintext space even if *this is not currently bootstrapppable,
    //   in case the calling application will make it bootstrappable later.

    //OLD: assert(p>=2);
    helib::assertTrue(p>=2, "Invalid p value found generating BGV key-switching matrix");
  }
  ksMatrix.ptxtSpace = p;

  // generate the RLWE instances with pseudorandom ai's

  for (long i = 0; i < n; i++) {
    ksMatrix.noiseBound = RLWE1(ksMatrix.b[i], a[i], toKey, p);
  }
  // Add in the multiples of the fromKey secret key
  fromKey *= context.productOfPrimes(context.specialPrimes);
  for (long i = 0; i < n; i++) {
    ksMatrix.b[i] += fromKey;
    fromKey *= context.productOfPrimes(context.digits[i]);
  }

  // Push the new matrix onto our list
  keySwitching.push_back(ksMatrix);
}

// Decryption
void FHESecKey::Decrypt(ZZX& plaintxt, const Ctxt &ciphertxt) const
{
  ZZX f;
  Decrypt(plaintxt, ciphertxt, f);
}

void FHESecKey::Decrypt(ZZX& plaintxt, const Ctxt &ciphertxt,
			ZZX& f) const // plaintext before modular reduction
{
  FHE_TIMER_START;
  //OLD: assert(getContext()==ciphertxt.getContext());
  helib::assertEq(getContext(), ciphertxt.getContext(), "Context mismatch");
  const IndexSet& ptxtPrimes = ciphertxt.primeSet;
  DoubleCRT ptxt(context, ptxtPrimes); // Set to zero

  // for each ciphertext part, fetch the right key, multiply and add
  for (size_t i=0; i<ciphertxt.parts.size(); i++) {
    const CtxtPart& part = ciphertxt.parts[i];
    if (part.skHandle.isOne()) { // No need to multiply
      ptxt += part;
      continue;
    }

    long keyIdx = part.skHandle.getSecretKeyID();
    DoubleCRT key = sKeys.at(keyIdx); // copy object, not a reference
    const IndexSet extraPrimes = key.getIndexSet() / ptxtPrimes;
    key.removePrimes(extraPrimes);    // drop extra primes, for efficiency

    /* Perhaps a slightly more efficient way of doing the same thing is:
       DoubleCRT key(context, ptxtPrimes); // a zero object wrt ptxtPrimes
       key.Add(sKeys.at(keyIdx), false); // add without mathcing primesSet
    */
    long xPower = part.skHandle.getPowerOfX();
    long sPower = part.skHandle.getPowerOfS();
    if (xPower>1) { 
      key.automorph(xPower); // s(X^t)
    }
    if (sPower>1) {
      key.Exp(sPower);       // s^r(X^t)
    }
    key *= part;
    ptxt += key;
  }
  // convert to coefficient representation & reduce modulo the plaintext space
  ptxt.toPoly(plaintxt);
  f = plaintxt; // f used only for debugging

  if (isCKKS()) return; // CKKS encryption, nothing else to do
  // NOTE: calling application must still divide by ratFactor after decoding

  PolyRed(plaintxt, ciphertxt.ptxtSpace, true/*reduce to [0,p-1]*/);

  // if p>2, multiply by (intFactor * Q)^{-1} mod p
  if (ciphertxt.getPtxtSpace()>2) {
    long factor = rem(context.productOfPrimes(ciphertxt.getPrimeSet()),
                     ciphertxt.ptxtSpace);
    factor = MulMod(factor, ciphertxt.intFactor, ciphertxt.ptxtSpace);
    if (factor != 1) {
        factor = InvMod(factor, ciphertxt.ptxtSpace);
        MulMod(plaintxt, plaintxt, factor, ciphertxt.ptxtSpace);
    }
  }
}

// Encryption using the secret key, this is useful, e.g., to put an
// encryption of the secret key into the public key.
long FHESecKey::skEncrypt(Ctxt &ctxt, const ZZX& ptxt,
                          long ptxtSpace, long skIdx) const
{
  FHE_TIMER_START;

  //OLD: assert(((FHEPubKey*)this) == &ctxt.pubKey);
  helib::assertEq(((const FHEPubKey*)this), &ctxt.pubKey, "Key does not match context's public key");

  long m = getContext().zMStar.getM();
  double ptxtSize = 1.0;
  if (isCKKS()) {
    if (ptxtSpace > 0)
      ptxtSize = ptxtSpace;
    ptxtSpace = 1;
  }
  else { // BGV
    if (ptxtSpace<2) 
      ptxtSpace = pubEncrKey.ptxtSpace; // default plaintext space is p^r
    //OLD: assert(ptxtSpace >= 2);
    helib::assertTrue(ptxtSpace >= 2, "Found invalid p value in BGV encryption");
  }
  ctxt.ptxtSpace = ptxtSpace;

  ctxt.primeSet = context.ctxtPrimes; // initialize the primeSet
  {CtxtPart tmpPart(context, context.ctxtPrimes);
  ctxt.parts.assign(2,tmpPart);}      // allocate space

  // Set Ctxt bookeeping parameters
  ctxt.intFactor = 1; // FIXME: is this necessary?

  // make parts[0],parts[1] point to (1,s)
  ctxt.parts[0].skHandle.setOne();
  ctxt.parts[1].skHandle.setBase(skIdx);

  const DoubleCRT& sKey = sKeys.at(skIdx);   // get key
  // Sample a new RLWE instance
  double noiseBound = RLWE(ctxt.parts[0], ctxt.parts[1], sKey, ptxtSpace);

  if (isCKKS()) {
    double f = getContext().ea->getCx().encodeScalingFactor() / ptxtSize;
    long prec = getContext().alMod.getPPowR();
    long ef = conv<long>(ceil(prec*noiseBound/(f*ptxtSize)));
    if (ef>1) { // scale up some more
      ctxt.parts[0] += ptxt * ef;
      f *= ef;
    }
    else {
      ctxt.parts[0] += ptxt;
    }
    // Round size to next power of two so as not to leak too much
    ctxt.ptxtMag = EncryptedArrayCx::roundedSize(ptxtSize);
    ctxt.ratFactor = f;
    ctxt.noiseBound = noiseBound;
    return long(f);
  }
  else { // BGV
    ctxt.addConstant(ptxt);  // add in the plaintext
    double ptxt_bound = context.noiseBoundForUniform(double(ptxtSpace)/2.0, context.zMStar.getPhiM());

    ctxt.noiseBound = noiseBound + ptxt_bound;
    return ctxt.ptxtSpace;
  }
}


// Generate bootstrapping data if needed, returns index of key
long FHESecKey::genRecryptData()
{
  if (recryptKeyID>=0) return recryptKeyID;

  // Make sure that the context has the bootstrapping EA and PAlgMod
  //OLD: assert(context.isBootstrappable());
  helib::assertTrue(context.isBootstrappable(), "Cannot generate recrypt data for non-bootstrappable context");

  long p2ePr = context.rcData.alMod->getPPowR();// p^{e-e'+r}
  long p2r = context.alMod.getPPowR(); // p^r

  // Generate a new bootstrapping key
  zzX keyPoly;
  long hwt = context.rcData.skHwt;
  double bound = sampleHWtBounded(keyPoly, context, hwt);

  DoubleCRT newSk(keyPoly, context, context.ctxtPrimes | context.specialPrimes); 
  // defined relative to all primes

  long keyID = ImportSecKey(newSk, bound, p2r, /*maxDegKswitch=*/1);

  // Generate a key-switching matrix from key 0 to this key
  GenKeySWmatrix(/*fromSPower=*/1,/*fromXPower=*/1,
		 /*fromIdx=*/0,   /*toIdx=*/keyID, /*ptxtSpace=*/p2r);

  // Encrypt new key under key #0 and plaintext space p^{e+r}
  Encrypt(recryptEkey, keyPoly, p2ePr);

  return (recryptKeyID=keyID); // return the new key-ID
}


ostream& operator<<(ostream& str, const FHESecKey& sk)
{
  str << "[" << ((const FHEPubKey&)sk) << endl
      << sk.sKeys.size() << endl;
  for (long i=0; i<(long)sk.sKeys.size(); i++)
    str << sk.sKeys[i] << endl;
  return str << "]";
}

istream& operator>>(istream& str, FHESecKey& sk)
{
  sk.clear();
  //  cerr << "FHESecKey[";
  seekPastChar(str, '['); // defined in NumbTh.cpp
  str >> (FHEPubKey&) sk;

  long nKeys;
  str >> nKeys;
  sk.sKeys.resize(nKeys, DoubleCRT(sk.getContext(),IndexSet::emptySet()));
  for (long i=0; i<nKeys; i++) str >> sk.sKeys[i];
  seekPastChar(str, ']');
  //  cerr << "]\n";
  return str;
}


void writeSecKeyBinary(ostream& str, const FHESecKey& sk)
{
  writeEyeCatcher(str, BINIO_EYE_SK_BEGIN);

  // Write out the public key part first.
  writePubKeyBinary(str, sk);

// Write out 
// 1. vector<DoubleCRT> sKeys  

  write_raw_vector<DoubleCRT>(str, sk.sKeys); 

  writeEyeCatcher(str, BINIO_EYE_SK_END);
}

void readSecKeyBinary(istream& str, FHESecKey& sk)
{
  int eyeCatcherFound = readEyeCatcher(str, BINIO_EYE_SK_BEGIN);
  //OLD: assert(eyeCatcherFound == 0);
  helib::assertEq(eyeCatcherFound, 0, "Could not find pre-secret key eyecatcher");

  // Read in the public key part first.
  readPubKeyBinary(str, sk);

  DoubleCRT blankDCRT(sk.getContext(), IndexSet::emptySet());
  read_raw_vector<DoubleCRT>(str, sk.sKeys, blankDCRT);

  eyeCatcherFound = readEyeCatcher(str, BINIO_EYE_SK_END);
  //OLD: assert(eyeCatcherFound == 0);
  helib::assertEq(eyeCatcherFound, 0, "Could not find post-secret key eyecatcher");
}

