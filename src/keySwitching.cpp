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
/**
 * @file keySwitching.cpp
 * @brief A few strategies for generating key-switching matrices
 *
 * Copyright IBM Corporation 2012 All rights reserved.
 */
#include <unordered_set>
#include <NTL/ZZ.h>
#include <helib/permutations.h>

#include <helib/binio.h>
#include <helib/keySwitching.h>
#include <helib/keys.h>
#include <helib/apiAttributes.h>
#include <helib/log.h>

namespace helib {

/******************** KeySwitch implementation **********************/
/********************************************************************/

KeySwitch::KeySwitch(long sPow, long xPow, long fromID, long toID, long p) :
    fromKey(sPow, xPow, fromID), toKeyID(toID), ptxtSpace(p)
{}

KeySwitch::KeySwitch(const SKHandle& _fromKey,
                     UNUSED long fromID,
                     long toID,
                     long p) :
    fromKey(_fromKey), toKeyID(toID), ptxtSpace(p)
{}

bool KeySwitch::operator==(const KeySwitch& other) const
{
  if (this == &other)
    return true;

  if (fromKey != other.fromKey)
    return false;
  if (toKeyID != other.toKeyID)
    return false;
  if (ptxtSpace != other.ptxtSpace)
    return false;

  if (prgSeed != other.prgSeed)
    return false;

  if (b.size() != other.b.size())
    return false;
  for (size_t i = 0; i < b.size(); i++)
    if (b[i] != other.b[i])
      return false;

  return true;
}
bool KeySwitch::operator!=(const KeySwitch& other) const
{
  return !(*this == other);
}

unsigned long KeySwitch::NumCols() const { return b.size(); }

bool KeySwitch::isDummy() const { return (toKeyID == -1); }

void KeySwitch::verify(SecKey& sk)
{
  long fromSPower = fromKey.getPowerOfS();
  long fromXPower = fromKey.getPowerOfX();
  long fromIdx = fromKey.getSecretKeyID();
  long toIdx = toKeyID;
  long p = ptxtSpace;
  long n = b.size();

  std::cout << "KeySwitch::verify\n";
  std::cout << "fromS = " << fromSPower << " fromX = " << fromXPower
            << " fromIdx = " << fromIdx << " toIdx = " << toIdx << " p = " << p
            << " n = " << n << "\n";

  if (fromSPower != 1 || fromXPower != 1 || (fromIdx == toIdx) || n == 0) {
    std::cout << "KeySwitch::verify: these parameters not checkable\n";
    return;
  }

  const Context& context = b[0].getContext();

  // we don't store the context in the ks matrix, so let's
  // check that they are consistent

  for (long i = 0; i < n; i++) {
    if (&context != &(b[i].getContext()))
      std::cout << "KeySwitch::verify: bad context " << i << "\n";
  }

  std::cout << "context.ctxtPrimes = " << context.ctxtPrimes << "\n";
  std::cout << "context.specialPrimes = " << context.specialPrimes << "\n";
  IndexSet fullPrimes = context.fullPrimes(); // ctxtPrimes | specialPrimes;

  std::cout << "digits: ";
  for (long i = 0; i < n; i++)
    std::cout << context.digits[i] << " ";
  std::cout << "\n";

  std::cout << "IndexSets of b: ";
  for (long i = 0; i < n; i++)
    std::cout << b[i].getMap().getIndexSet() << " ";
  std::cout << "\n";

  // VJS: suspicious shadowing of fromKey, toKey
  const DoubleCRT& _fromKey = sk.sKeys.at(fromIdx);
  const DoubleCRT& _toKey = sk.sKeys.at(toIdx);

  std::cout << "IndexSet of fromKey: " << _fromKey.getMap().getIndexSet()
            << "\n";
  std::cout << "IndexSet of toKey: " << _toKey.getMap().getIndexSet() << "\n";

  std::vector<DoubleCRT> a;
  a.resize(n, DoubleCRT(context, fullPrimes)); // defined modulo all primes

  {
    RandomState state;

    SetSeed(prgSeed);
    for (long i = 0; i < n; i++)
      a[i].randomize();

  } // the RandomState destructor "restores the state" (see NumbTh.h)

  std::vector<NTL::ZZX> A, B;

  A.resize(n);
  B.resize(n);

  for (long i = 0; i < n; i++) {
    a[i].toPoly(A[i]);
    b[i].toPoly(B[i]);
  }

  NTL::ZZX FromKey, ToKey;
  _fromKey.toPoly(FromKey, fullPrimes);
  _toKey.toPoly(ToKey, fullPrimes);

  NTL::ZZ Q = context.productOfPrimes(fullPrimes);
  NTL::ZZ prod = context.productOfPrimes(context.specialPrimes);
  NTL::ZZX C, D;
  NTL::ZZX PhimX = context.zMStar.getPhimX();

  long nb = 0;
  for (long i = 0; i < n; i++) {
    C = (B[i] - FromKey * prod + ToKey * A[i]) % PhimX;
    PolyRed(C, Q);
    if (!divide(D, C, p)) {
      std::cout << "*** not divisible by p at " << i << "\n";
    } else {
      for (long j = 0; j <= deg(D); j++)
        if (NumBits(coeff(D, j)) > nb)
          nb = NumBits(coeff(D, j));
    }
    prod *= context.productOfPrimes(context.digits[i]);
  }

  std::cout << "error ratio: " << ((double)nb) / ((double)NumBits(Q)) << "\n";
}

const KeySwitch& KeySwitch::dummy()
{
  static const KeySwitch dummy(-1, -1, -1, -1);
  return dummy;
}

std::ostream& operator<<(std::ostream& str, const KeySwitch& matrix)
{
  str << "[" << matrix.fromKey << " " << matrix.toKeyID << " "
      << matrix.ptxtSpace << " " << matrix.b.size() << std::endl;
  for (long i = 0; i < (long)matrix.b.size(); i++)
    str << matrix.b[i] << std::endl;
  str << matrix.prgSeed << " " << matrix.noiseBound << "]";
  return str;
}

// Used in lieu of std::istream& operator>>(std::istream& str, KeySwitch&
// matrix)
void KeySwitch::readMatrix(std::istream& str, const Context& context)
{
  seekPastChar(str, '['); // defined in NumbTh.cpp
  str >> fromKey;
  str >> toKeyID;
  str >> ptxtSpace;

  long nDigits;
  str >> nDigits;
  b.resize(nDigits, DoubleCRT(context, IndexSet::emptySet()));
  for (long i = 0; i < nDigits; i++)
    str >> b[i];
  str >> prgSeed;
  str >> noiseBound;
  seekPastChar(str, ']');
}

void KeySwitch::write(std::ostream& str) const
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

void KeySwitch::read(std::istream& str, const Context& context)
{
  int eyeCatcherFound = readEyeCatcher(str, BINIO_EYE_SKM_BEGIN);
  assertEq(eyeCatcherFound, 0, "Could not find pre-secret key eyecatcher");

  fromKey.read(str);
  toKeyID = read_raw_int(str);
  ptxtSpace = read_raw_int(str);
  DoubleCRT blankDCRT(context, IndexSet::emptySet());
  read_raw_vector(str, b, blankDCRT);
  read_raw_ZZ(str, prgSeed);
  noiseBound = read_raw_xdouble(str);

  eyeCatcherFound = readEyeCatcher(str, BINIO_EYE_SKM_END);
  assertEq(eyeCatcherFound, 0, "Could not find post-secret key eyecatcher");
}

long KSGiantStepSize(long D)
{
  assertTrue<InvalidArgument>(D > 0l, "Step size must be positive");
  long g = NTL::SqrRoot(D);
  if (g * g < D)
    g++; // g = ceiling(sqrt(D))
  return g;
}

// A maximalistic approach: generate matrices s(X^e)->s(X) for all e \in Zm*
void addAllMatrices(SecKey& sKey, long keyID)
{
  const Context& context = sKey.getContext();
  long m = context.zMStar.getM();

  // key-switching matrices for the automorphisms
  for (long i = 0; i < m; i++) {
    if (!context.zMStar.inZmStar(i))
      continue;
    sKey.GenKeySWmatrix(1, i, keyID, keyID);
  }
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

// TODO: generate matrices s.t. you can reLinearize each s(X^e) in at most two
// steps void addFewMatrices(SecKey& sKey, long keyID)
// {
//   throw LogicError("addFewMatrices not implemented yet");
// }

// This code block appears at least twice below
#define computeParams(context, m, i)                                           \
  bool native;                                                                 \
  long ord, gi, g2md;                                                          \
  NTL::mulmod_precon_t g2mdminv;                                               \
  if (i == context.zMStar.numOfGens()) { /* Frobenius matrices */              \
    ord = context.zMStar.getOrdP();                                            \
    gi = context.zMStar.getP();                                                \
    native = true;                                                             \
  } else { /* one of the "regular" dimensions */                               \
    ord = context.zMStar.OrderOf(i);                                           \
    gi = context.zMStar.ZmStarGen(i);                                          \
    native = context.zMStar.SameOrd(i);                                        \
    if (!native) {                                                             \
      g2md = PowerMod(gi, -ord, m); /* g^{-ord} mod m */                       \
      g2mdminv = PrepMulModPrecon(g2md, m);                                    \
    }                                                                          \
  }                                                                            \
  NTL::mulmod_precon_t giminv = PrepMulModPrecon(gi, m);

#if 0
static void add1Dmats4dim(SecKey& sKey, long i, long keyID)
{
  const Context &context = sKey.getContext();
  long m = context.zMStar.getM();
  computeParams(context,m,i); // defines vars: native, ord, gi, g2md, giminv, g2mdminv

  /* MAUTO std::vector<long> vals; */
  for (long j=1,val=gi; j < ord; j++) {
    // From s(X^val) to s(X)
    sKey.GenKeySWmatrix(1, val, keyID, keyID);
    if (!native) { // also from s(X^{g^{i-ord}}) to s(X)
      long val2 = MulModPrecon(val,g2md,m,g2mdminv);
      sKey.GenKeySWmatrix(1, val2, keyID, keyID);
      /* MAUTO vals.push_back(val2); */
    }
    /* MAUTO vals.push_back(val); */
    val = MulModPrecon(val, gi, m, giminv); // val *= g mod m (= g^{j+1})
  }

  if (!native) {
    sKey.GenKeySWmatrix(1, context.zMStar.genToPow(i, -ord), keyID, keyID);
  }


/* MAUTO
  sKey.resetTree(i,keyID); // remove existing tree, if any
  sKey.add2tree(i, 1, vals, keyID);
*/
}
#else
// adds all matrices for dim i.
// i == -1 => Frobenius (NOTE: in matmul1D, i ==#gens means something else,
//   so it is best to avoid that).
static void add1Dmats4dim(SecKey& sKey, long i, long keyID)
{
  const PAlgebra& zMStar = sKey.getContext().zMStar;
  long ord;
  bool native;

  if (i != -1) {
    ord = zMStar.OrderOf(i);
    native = zMStar.SameOrd(i);
  } else {
    // Frobenius
    ord = zMStar.getOrdP();
    native = true;
  }

  for (long j = 1; j < ord; j++)
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, j), keyID, keyID);

  if (!native)
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, -ord), keyID, keyID);

  sKey.setKSStrategy(i, HELIB_KSS_FULL);
}

#endif

#if 0
static std::pair<long,long> computeSteps(long ord, long bound, bool native)
{
  long baby,giant;
  if (native) { // using giant+baby matrices
    if (bound*bound >= 4*ord)
      giant = ceil((bound - sqrt((double)bound*bound -4*ord))/2.0);
    else
      giant = sqrt((double)ord);
  }
  else { // using giant+2*baby matrices
    if (bound*bound >= 8*ord)
      giant = ceil((bound - sqrt((double)bound*bound -8*ord))/2.0);
    else
      giant = sqrt((double)ord);
  }
  baby = ord/giant;
  if (baby*giant<ord) baby++;

  //std::cerr << "*** giant steps = " << giant << "\n";
  return std::pair<long,long>(baby,giant);
}

static void addSome1Dmats4dim(SecKey& sKey, long i, long bound, long keyID)
{
  const Context &context = sKey.getContext();
  long m = context.zMStar.getM();
  computeParams(context,m,i); // defines vars: native, ord, gi, g2md, giminv, g2mdminv

  long baby, giant;
  std::tie(baby,giant) = computeSteps(ord, bound, native);

  for (long j=1,val=gi; j<=baby; j++) { // Add matrices for baby steps
    sKey.GenKeySWmatrix(1, val, keyID, keyID);
    if (!native) {
      long val2 = MulModPrecon(val,g2md,m,g2mdminv);
      sKey.GenKeySWmatrix(1, val2, keyID, keyID);
    }
    val = MulModPrecon(val, gi, m, giminv); // val *= g mod m (= g^{j+1})
   }

  long gb = PowerMod(gi,baby,m); // g^baby
  NTL::mulmod_precon_t gbminv = PrepMulModPrecon(gb, m);
  for (long j=2,val=gb; j < giant; j++) { // Add matrices for giant steps
    val = MulModPrecon(val, gb, m, gbminv); // val = g^{(j+1)*baby}
    sKey.GenKeySWmatrix(1, val, keyID, keyID);
  }

  if (!native) {
    sKey.GenKeySWmatrix(1, context.zMStar.genToPow(i, -ord), keyID, keyID);
  }

  // VJS: experimantal feature...because the replication code
  // uses rotations by -1, -2, -4, -8, we add a few
  // of these as well...only the small ones are important,
  // and we only need them if SameOrd(i)...
  // Note: we do indeed get a nontrivial speed-up

  if (native && i<context.zMStar.numOfGens()) {
    for (long k = 1; k < giant; k = 2*k) {
      long j = ord - k;
      long val = PowerMod(gi, j, m); // val = g^j
      sKey.GenKeySWmatrix(1, val, keyID, keyID);
    }
  }

#if 0
MAUTO

  // build the tree for this dimension, the internal nodes are 1 and
  // (subset of) gi^{giant}, gi^{2*giant}, ..., gi^{baby*giant}. We

  MAUTO sKey.resetTree(i,keyID); // remove existing tree, if any

  // keep a list of all the elements that are covered by the tree so far,
  // initialized to only the root (=1).
  std::unordered_set<long> covered({1});

  // Make a list of the automorphisms for this dimension
  std::vector<long> autos;
  for (long j=1,val=gi; j<ord; j++) {
    // Do we have matrices for val and/or val/gi^{di}?
    if (!native) {
      long val2 = MulModPrecon(val, g2md, m, g2mdminv);
      if (sKey.haveKeySWmatrix(1,val2,keyID,keyID)) {
        autos.push_back(val2);
      }
    }
    if (sKey.haveKeySWmatrix(1,val,keyID,keyID)) {
      autos.push_back(val);
    }
    val = MulModPrecon(val, gi, m, giminv); // g^{j+1}
  }

  // Insert internal nodes and their children to tree
  for (long j=0,fromVal=1; j<giant; j++) {
    NTL::mulmod_precon_t fromminv = PrepMulModPrecon(fromVal, m);
    std::vector<long> children;
    for (long k: autos) {
      long toVal = MulModPrecon(k, fromVal, m, fromminv);
      if (covered.count(toVal)==0) { // toVal not covered yet
        covered.insert(toVal);
        children.push_back(toVal);
      }
    }
    if (!children.empty()) { // insert fromVal with its children
      sKey.add2tree(i, fromVal, children, keyID);
    }
    fromVal = MulModPrecon(fromVal, gb, m, gbminv); // g^{(j+1)*baby}
  }

  // Sanity-check, did we cover everything?
  long toCover = native? ord: (2*ord-1);
  if (covered.size()<toCover)
    std::cerr << "**Warning: order-"<<ord<<" dimension, covered "<<covered.size()
         << " of "<<toCover<<std::endl;
#endif
}

#else
// same as above, but uses BS/GS strategy
static void addSome1Dmats4dim(SecKey& sKey,
                              long i,
                              UNUSED long bound,
                              long keyID)
{
  const PAlgebra& zMStar = sKey.getContext().zMStar;
  long ord;
  bool native;

  if (i != -1) {
    ord = zMStar.OrderOf(i);
    native = zMStar.SameOrd(i);
  } else {
    // Frobenius
    ord = zMStar.getOrdP();
    native = true;
  }

  long g = KSGiantStepSize(ord);

  // baby steps
  for (long j = 1; j < g; j++)
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, j), keyID, keyID);

  // giant steps
  for (long j = g; j < ord; j += g)
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, j), keyID, keyID);

  if (!native)
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, -ord), keyID, keyID);

  sKey.setKSStrategy(i, HELIB_KSS_BSGS);

  // NOTE: the old code also added matrices for ord-2^k for small k,
  // in the case of (native && i<context.zMStar.numOfGens()).
  // This supposedly speeds up the replication code, but for now
  // we are leaving this out, for simplicity.   Also, it is a waste
  // of space for applications that don't use replication.
}

#endif

// generate only matrices of the form s(X^{g^i})->s(X), but not all of them.
// For a generator g whose order is larger than bound, generate only enough
// matrices for the giant-step/baby-step procedures (2*sqrt(ord(g))of them).
void addSome1DMatrices(SecKey& sKey, long bound, long keyID)
{
  const Context& context = sKey.getContext();

  // key-switching matrices for the automorphisms
  for (long i : range(context.zMStar.numOfGens())) {
    // For generators of small order, add all the powers
    if (bound >= context.zMStar.OrderOf(i))
      add1Dmats4dim(sKey, i, keyID);
    else // For generators of large order, add only some of the powers
      addSome1Dmats4dim(sKey, i, bound, keyID);
  }
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

void add1DMatrices(SecKey& sKey, long keyID)
{
  addSome1DMatrices(sKey, LONG_MAX, keyID);
}

void addBSGS1DMatrices(SecKey& sKey, long keyID)
{
  addSome1DMatrices(sKey, 0, keyID);
}

// Generate all Frobenius matrices of the form s(X^{p^i})->s(X)
void addSomeFrbMatrices(SecKey& sKey, long bound, long keyID)
{
  const Context& context = sKey.getContext();
  if (bound >= LONG(context.zMStar.getOrdP()))
    add1Dmats4dim(sKey, -1, keyID);
  else // For generators of large order, add only some of the powers
    addSome1Dmats4dim(sKey, -1, bound, keyID);

  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

void addFrbMatrices(SecKey& sKey, long keyID)
{
  addSomeFrbMatrices(sKey, LONG_MAX, keyID);
}

void addBSGSFrbMatrices(SecKey& sKey, long keyID)
{
  addSomeFrbMatrices(sKey, 0, keyID);
}

static void addMinimal1Dmats4dim(SecKey& sKey, long i, long keyID)
{
  const PAlgebra& zMStar = sKey.getContext().zMStar;
  long ord;
  bool native;

  if (i != -1) {
    ord = zMStar.OrderOf(i);
    native = zMStar.SameOrd(i);
  } else {
    // Frobenius
    ord = zMStar.getOrdP();
    native = true;
  }

  sKey.GenKeySWmatrix(1, zMStar.genToPow(i, 1), keyID, keyID);

  if (!native)
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, -ord), keyID, keyID);

  if (ord > HELIB_KEYSWITCH_MIN_THRESH) {
    long g = KSGiantStepSize(ord);
    sKey.GenKeySWmatrix(1, zMStar.genToPow(i, g), keyID, keyID);
  }

  sKey.setKSStrategy(i, HELIB_KSS_MIN);
}

void addMinimal1DMatrices(SecKey& sKey, long keyID)
{
  const Context& context = sKey.getContext();

  // key-switching matrices for the automorphisms
  for (long i : range(context.zMStar.numOfGens())) {
    addMinimal1Dmats4dim(sKey, i, keyID);
  }
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

// Generate all Frobenius matrices of the form s(X^{p^i})->s(X)
void addMinimalFrbMatrices(SecKey& sKey, long keyID)
{
  addMinimal1Dmats4dim(sKey, -1, keyID);
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

// Generate all key-switching matrices for a given permutation network
void addMatrices4Network(SecKey& sKey, const PermNetwork& net, long keyID)
{
  const Context& context = sKey.getContext();
  long m = context.zMStar.getM();

  for (long i = 0; i < net.depth(); i++) {
    long e = net.getLayer(i).getE();
    long gIdx = net.getLayer(i).getGenIdx();
    long g = context.zMStar.ZmStarGen(gIdx);
    long g2e = NTL::PowerMod(g, e, m); // g^e mod m
    const NTL::Vec<long>& shamts = net.getLayer(i).getShifts();
    for (long j = 0; j < shamts.length(); j++) {
      if (shamts[j] == 0)
        continue;
      long val = NTL::PowerMod(g2e, shamts[j], m);
      sKey.GenKeySWmatrix(1, val, keyID, keyID);
    }
  }
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

void addTheseMatrices(SecKey& sKey, const std::set<long>& automVals, long keyID)
{
  std::set<long>::iterator it;
  for (it = automVals.begin(); it != automVals.end(); ++it) {
    long k = *it;
    sKey.GenKeySWmatrix(1, k, keyID, keyID);
  }
  sKey.setKeySwitchMap(); // re-compute the key-switching map
}

} // namespace helib
