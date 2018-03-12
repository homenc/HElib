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

#include "FHE.h"

#include <queue> // used in the breadth-first search in setKeySwitchMap
#include "timing.h"

/******** Utility function to generate RLWE instances *********/

// Choose random c0,c1 such that c0+s*c1 = p*e for a short e
void RLWE(DoubleCRT& c0,DoubleCRT& c1, const DoubleCRT &s, long p, ZZ* prgSeed)
{
  assert (p>0); // Can be used with p=1, but we always use with p>=2

  // choose c1 at random (using prgSeed if not NULL)
  c1.randomize(prgSeed);

  // choose a short error e, set c0 =  p*e - c1*s
  c0.sampleGaussian();
  c0 *= p;

  // It is assumed that c0,c1 are defined with respect to the same set of
  // primes, but s may be defined relative to a different set. Either way
  // the primes for of c0,c1 are unchanged.
  DoubleCRT tmp(c1);
  tmp.Mul(s, /*matchIndexSets=*/false); // multiply but don't mod-up
  c0 -= tmp;
}

// same as above, but assumes that c1 is already chosen by the caller
void RLWE1(DoubleCRT& c0, const DoubleCRT& c1, const DoubleCRT &s, long p)
{
  assert (p>0); // Can be used with p=1, but we always use with p>=2

  // choose a short error e, set c0 =  p*e - c1*s
  c0.sampleGaussian();
  c0 *= p;

  // It is assumed that c0,c1 are defined with respect to the same set of
  // primes, but s may be defined relative to a different set. Either way
  // the primes for of c0,c1 are unchanged.
  DoubleCRT tmp(c1);
  tmp.Mul(s, /*matchIndexSets=*/false); // multiply but don't mod-up
  c0 -= tmp;
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

  IndexSet allPrimes = context.ctxtPrimes | context.specialPrimes;

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
  a.resize(n, DoubleCRT(context, allPrimes)); // defined modulo all primes

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
  _fromKey.toPoly(FromKey, allPrimes);
  _toKey.toPoly(ToKey, allPrimes);

  ZZ Q = context.productOfPrimes(allPrimes);
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
  str << matrix.prgSeed << "]";
  return str;
}

// Used in lieu of istream& operator>>(istream& str, KeySwitch& matrix)
void KeySwitch::readMatrix(istream& str, const FHEcontext& context)
{
  //  cerr << "KeySwitch[";
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
  seekPastChar(str,']');
  //  cerr << "]";
}

/******************** FHEPubKey implementation **********************/
/********************************************************************/
// Computes the keySwitchMap pointers, using breadth-first search (BFS)

void FHEPubKey::setKeySwitchMap(long keyId)
{
  assert(keyId>=0 && keyId<(long)skHwts.size()); // Sanity-check, do we have such a key?
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


// Encrypts plaintext, result returned in the ciphertext argument. The
// returned value is the plaintext-space for that ciphertext. When called
// with highNoise=true, returns a ciphertext with noise level~q/8.
long FHEPubKey::Encrypt(Ctxt &ctxt, const ZZX& ptxt, long ptxtSpace,
			bool highNoise) const
{
  FHE_TIMER_START;
  assert(this == &ctxt.pubKey);

  if (ptxtSpace != pubEncrKey.ptxtSpace) { // plaintext-space mistamtch
    ptxtSpace = GCD(ptxtSpace, pubEncrKey.ptxtSpace);
    if (ptxtSpace <= 1) Error("Plaintext-space mismatch on encryption");
  }

  // generate a random encryption of zero from the public encryption key
  ctxt = pubEncrKey;  // already an encryption of zero, just not a random one

  // choose a random small scalar r and a small random error vector e,
  // then set ctxt = r*pubEncrKey + ptstSpace*e + (ptxt,0)
  DoubleCRT e(context, context.ctxtPrimes);
  DoubleCRT r(context, context.ctxtPrimes);
  r.sampleSmall();

  for (size_t i=0; i<ctxt.parts.size(); i++) {  // add noise to all the parts
    ctxt.parts[i] *= r;

    if (highNoise && i == 0) {
      // we sample e so that coefficients are uniform over 
      // [-Q/(8*ptxtSpace)..Q/(8*ptxtSpace)]

      ZZ B;
      B = context.productOfPrimes(context.ctxtPrimes);
      B /= ptxtSpace;
      B /= 8;
      e.sampleUniform(B);
    }
    else { 
      e.sampleGaussian();
    }

    e *= ptxtSpace;
    ctxt.parts[i] += e;
  }

  // add in the plaintext
  // FIXME: This relies on the first part, ctxt[0], to have handle to 1
  if (ptxtSpace==2) ctxt.parts[0] += ptxt;

  else { // The general case of ptxtSpace>2: for a ciphertext
         // relative to modulus Q, we add ptxt * Q mod ptxtSpace.
    long QmodP = rem(context.productOfPrimes(ctxt.primeSet), ptxtSpace);
    ctxt.parts[0] += MulMod(ptxt,QmodP,ptxtSpace); // MulMod from module NumbTh
  }

  // fill in the other ciphertext data members
  ctxt.ptxtSpace = ptxtSpace;

  if (highNoise) {
    // hack: we set noiseVar to Q^2/8, which is just below threshold 
    // that will signal an error

    ctxt.noiseVar = xexp(2*context.logOfProduct(context.ctxtPrimes) - log(8.0));

  }
  else {
    // We have <skey,ctxt>= r*<skey,pkey> +p*(e0+e1*s) +m, where VAR(<skey,pkey>)
    // is recorded in pubEncrKey.noiseVar, VAR(ei)=sigma^2*phi(m), and VAR(s) is
    // determined by the secret-key Hamming weight (skHwt). 
    // VAR(r)=phi(m)/2, hence the expected size squared is bounded by:
    // E(X^2) <= pubEncrKey.noiseVar *phi(m) *stdev^2
    //                               + p^2*sigma^2 *phi(m) *(skHwt+1) + p^2
  
    long hwt = skHwts[0];
    xdouble phim = to_xdouble(context.zMStar.getPhiM());
    xdouble sigma2 = context.stdev * context.stdev;
    xdouble p2 = to_xdouble(ptxtSpace) * to_xdouble(ptxtSpace);
    ctxt.noiseVar = pubEncrKey.noiseVar*phim*0.5 
                    + p2*sigma2*phim*(hwt+1)*context.zMStar.get_cM() + p2;
  }
  return ptxtSpace;
}

bool FHEPubKey::operator==(const FHEPubKey& other) const
{
  if (this == &other) return true;

  if (&context != &other.context) return false;
  if (!pubEncrKey.equalsTo(other.pubEncrKey, /*comparePkeys=*/false))
    return false;

  if (skHwts.size() != other.skHwts.size()) return false;
  for (size_t i=0; i<skHwts.size(); i++)
    if (skHwts[i] != other.skHwts[i]) return false;

  if (keySwitching.size() != other.keySwitching.size()) return true;
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
  while (n >= 0 && KS_strategy[n-1] == FHE_KSS_UNKNOWN) n--;
  long n1 = other.KS_strategy.length();
  while (n1 >= 0 && other.KS_strategy[n1-1] == FHE_KSS_UNKNOWN) n1--;
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

  // output skHwts in the same format as vec_long
  str << "[";
  for (long i=0; i<(long)pk.skHwts.size(); i++)
    str << pk.skHwts[i]<<" ";
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
  const PAlgebra& palg = pk.getContext().zMStar;
  assert( m == palg.getM() );
  assert( p == palg.getP() );
  assert( gens.size() == palg.numOfGens() );
  for (long i=0; i<(long)gens.size(); i++) {
    assert(gens[i]==(long)palg.ZmStarGen(i));
    if (palg.SameOrd(i))
      assert(ords[i]==(long) pk.getContext().zMStar.OrderOf(i));
    else
      assert(-ords[i]==(long) pk.getContext().zMStar.OrderOf(i));
  }

  // Get the public encryption key itself
  str >> pk.pubEncrKey;

  // Get the vector of secret-key Hamming-weights
  vec_long vl;
  str >> vl;
  pk.skHwts.resize(vl.length());
  for (long i=0; i<(long)pk.skHwts.size(); i++) pk.skHwts[i] = vl[i];

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
  for (long i=pk.skHwts.size()-1; i>=0; i--)
    pk.setKeySwitchMap(i);

  str >> pk.KS_strategy; 

  // Get the bootstrapping key, if any
  str >> pk.recryptKeyID;
  if (pk.recryptKeyID>=0) str >> pk.recryptEkey;

  seekPastChar(str, ']');
  return str;
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
long FHESecKey::ImportSecKey(const DoubleCRT& sKey, long Hwt,
			     long ptxtSpace, long maxDegKswitch)
{
  if (sKeys.empty()) { // 1st secret-key, generate corresponding public key
    if (ptxtSpace<2)
      ptxtSpace = context.alMod.getPPowR(); // default plaintext space is p^r

    pubEncrKey.parts.assign(2,CtxtPart(context,context.ctxtPrimes));// allocate space
    RLWE(pubEncrKey.parts[0], pubEncrKey.parts[1], sKey, ptxtSpace); // a new RLWE instance

    // make parts[0],parts[1] point to (1,s)
    pubEncrKey.parts[0].skHandle.setOne();
    pubEncrKey.parts[1].skHandle.setBase();

    // Set the other Ctxt bookeeping parameters in pubEncrKey
    pubEncrKey.primeSet = context.ctxtPrimes;
    pubEncrKey.ptxtSpace = ptxtSpace;

    xdouble phim = to_xdouble(context.zMStar.getPhiM());
    pubEncrKey.noiseVar = context.stdev * context.stdev
      * phim * ptxtSpace * ptxtSpace;
  }
  skHwts.push_back(Hwt); // record the Hamming weight of the new secret-key
  sKeys.push_back(sKey); // add to the list of secret keys
  long keyID = sKeys.size()-1; // not thread-safe?

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

  ksMatrix.b.resize(n, DoubleCRT(context)); // size-n vector

  vector<DoubleCRT> a; 
  a.resize(n, DoubleCRT(context));

  { RandomState state;
    SetSeed(ksMatrix.prgSeed);
    for (long i = 0; i < n; i++) 
      a[i].randomize();
  } // restore state upon destruction of state

  // Record the plaintext space for this key-switching matrix
  if (p<2) {
    if (context.isBootstrappable()) // use larger bootstrapping plaintext space
         p = context.rcData.alMod->getPPowR();
    else p = pubEncrKey.ptxtSpace; // default plaintext space from public key
  }
  // FIXME: We use context.isBootstrappable() rather than
  //   this->isBootstrappable(). So we get the larger bootstrapping
  //   plaintext space even if *this is not currently bootstrapppable,
  //   in case the calling application will make it bootstrappable later.

  assert(p>=2);
  ksMatrix.ptxtSpace = p;

  // generate the RLWE instances with pseudorandom ai's

  for (long i = 0; i < n; i++) {
    RLWE1(ksMatrix.b[i], a[i], toKey, p); 
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
#ifdef DEBUG_PRINTOUT
  // The call to findBaseSet is only for the purpose of printing a
  // warning if the noise is large enough so as to risk decryption error
  IndexSet s; ciphertxt.findBaseSet(s);
#endif
  FHE_TIMER_START;
  assert(getContext()==ciphertxt.getContext());
  const IndexSet& ptxtPrimes = ciphertxt.primeSet;
  DoubleCRT ptxt(context, ptxtPrimes); // Set to zero

  // for each ciphertext part, fetch the right key, multiply and add
  for (size_t i=0; i<ciphertxt.parts.size(); i++) {
    const CtxtPart& part = ciphertxt.parts[i];
    //  cout << "decrypt part: "<<part.skHandle<<" "<< part.getIndexSet()<<"\n";
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
  f = plaintxt;

  if (ciphertxt.ptxtSpace>2) { // if p>2, multiply by Q^{-1} mod p
    long qModP = rem(context.productOfPrimes(ciphertxt.getPrimeSet()), 
		     ciphertxt.ptxtSpace);
    if (qModP != 1) {
      qModP = InvMod(qModP, ciphertxt.ptxtSpace);
      MulMod(plaintxt, plaintxt, qModP, ciphertxt.ptxtSpace);
    }
  }
  PolyRed(plaintxt, ciphertxt.ptxtSpace, true/*reduce to [0,p-1]*/);
}
void FHESecKey::Decrypt(DoubleCRT &ptxt, const Ctxt &ciphertxt) const 
{
#ifdef DEBUG_PRINTOUT
  // The call to findBaseSet is only for the purpose of printing a
  // warning if the noise is large enough so as to risk decryption error
  IndexSet s; ciphertxt.findBaseSet(s);
#endif
  FHE_TIMER_START;
  assert(getContext()==ciphertxt.getContext());
  const IndexSet& ptxtPrimes = ciphertxt.primeSet;

  // for each ciphertext part, fetch the right key, multiply and add
  for (size_t i=0; i<ciphertxt.parts.size(); i++) {
    const CtxtPart& part = ciphertxt.parts[i];
    //  cout << "decrypt part: "<<part.skHandle<<" "<< part.getIndexSet()<<"\n";
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
  if (ciphertxt.ptxtSpace>2) { // if p>2, multiply by Q^{-1} mod p
    long qModP = rem(context.productOfPrimes(ciphertxt.getPrimeSet()), 
		     ciphertxt.ptxtSpace);
    if (qModP != 1) {
      qModP = InvMod(qModP, ciphertxt.ptxtSpace);
	  ptxt *= qModP;
    }
  }
}

// Encryption using the secret key, this is useful, e.g., to put an
// encryption of the secret key into the public key.
long FHESecKey::Encrypt(Ctxt &ctxt, const ZZX& ptxt,
			long ptxtSpace, long skIdx) const
{
  FHE_TIMER_START;
  assert(((FHEPubKey*)this) == &ctxt.pubKey);

  if (ptxtSpace<2) 
    ptxtSpace = pubEncrKey.ptxtSpace; // default plaintext space is p^r
  assert(ptxtSpace >= 2);

  ctxt.primeSet = context.ctxtPrimes; // initialize the primeSet
  {CtxtPart tmpPart(context, context.ctxtPrimes);
  ctxt.parts.assign(2,tmpPart);}      // allocate space

  // Set Ctxt bookeeping parameters
  ctxt.ptxtSpace = ptxtSpace;
  xdouble phim = to_xdouble(context.zMStar.getPhiM());
  ctxt.noiseVar = context.stdev*context.stdev * phim * ptxtSpace*ptxtSpace;

  // make parts[0],parts[1] point to (1,s)
  ctxt.parts[0].skHandle.setOne();
  ctxt.parts[1].skHandle.setBase(skIdx);

  const DoubleCRT& sKey = sKeys.at(skIdx);   // get key
  RLWE(ctxt.parts[0], ctxt.parts[1], sKey, ptxtSpace); // a new RLWE instance

  // add in the plaintext
  ctxt.addConstant(ptxt);
  return ptxtSpace;
}


// Generate bootstrapping data if needed, returns index of key
long FHESecKey::genRecryptData()
{
  if (recryptKeyID>=0) return recryptKeyID;

  // Make sure that the context has the bootstrapping EA and PAlgMod
  assert(context.isBootstrappable());

  long p2ePr = context.rcData.alMod->getPPowR();// p^{e-e'+r}
  long p2r = context.alMod.getPPowR(); // p^r

  // Generate a new bootstrapping key
  ZZX keyPoly;
  const long hwt = context.rcData.skHwt;
  sampleHWt(keyPoly, hwt, context.zMStar.getPhiM());
  DoubleCRT newSk(keyPoly, context); // defined relative to all primes
  long keyID = ImportSecKey(newSk, hwt, p2r, /*maxDegKswitch=*/1);

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

