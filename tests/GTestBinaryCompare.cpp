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
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <NTL/BasicThreadPool.h>

#include <helib/helib.h>

#include <helib/intraSlot.h>
#include <helib/binaryArith.h>
#include <helib/binaryCompare.h>

#include "gtest/gtest.h"
#include "test_common.h"

#include <helib/debugging.h>

// define flags FLAG_PRINT_ZZX, FLAG_PRINT_POLY, FLAG_PRINT_VEC, functions
//        decryptAndPrint(ostream, ctxt, sk, ea, flags)
//        decryptAndCompare(ctxt, sk, ea, pa);

namespace {

struct Parameters
{
  Parameters(long prm, long bitSize, bool bootstrap, long seed, long nthreads) :
      prm(prm),
      bitSize(bitSize),
      bootstrap(bootstrap),
      seed(seed),
      nthreads(nthreads){};

  long prm;       // parameter size (0-tiny,...,4-huge)
  long bitSize;   // bitSize of input integers (<=32)
  bool bootstrap; // test comparison with bootstrapping
  long seed;      // PRG seed
  long nthreads;  // number of threads

  // Tell google test how to print Parameters
  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "prm=" << params.prm << ","
              << "bitSize=" << params.bitSize << ","
              << "bootstrap=" << params.bootstrap << ","
              << "seed=" << params.seed << ","
              << "nthreads=" << params.nthreads << "}";
  };
};

class GTestBinaryCompare :
    public ::testing::TestWithParam<std::tuple<Parameters, int>>
{
protected:
  static std::vector<helib::zzX> unpackSlotEncoding;
  constexpr static long mValues[][15] = {
      // { p, phi(m),   m,   d, m1, m2, m3,    g1,   g2,   g3, ord1,ord2,ord3,
      // B,c}
      {2, 48, 105, 12, 3, 35, 0, 71, 76, 0, 2, 2, 0, 25, 2},
      {2, 600, 1023, 10, 11, 93, 0, 838, 584, 0, 10, 6, 0, 25, 2},
      {2, 2304, 4641, 24, 7, 3, 221, 3979, 3095, 3760, 6, 2, -8, 25, 3},
      {2, 15004, 15709, 22, 23, 683, 0, 4099, 13663, 0, 22, 31, 0, 25, 3},
      // clang-format off
      {2, 27000, 32767, 15, 31, 7, 151, 11628, 28087, 25824, 30, 6, -10, 28, 4}
      // clang-format on
  };

  static long correctBitSize(long minimum, long oldBitSize)
  {
    long newBitSize;
    if (oldBitSize <= 0)
      newBitSize = minimum;
    else if (oldBitSize > 32)
      newBitSize = 32;
    else
      newBitSize = oldBitSize;
    return newBitSize;
  };

  // Validates the prm value, throwing if invalid
  static long validatePrm(long prm)
  {
    if (prm < 0 || prm >= 5)
      throw std::invalid_argument("prm must be in the interval [0, 4]");
    return prm;
  };

  static NTL::Vec<long> calculateMvec(const long* vals)
  {
    NTL::Vec<long> mvec;
    append(mvec, vals[4]);
    if (vals[5] > 1)
      append(mvec, vals[5]);
    if (vals[6] > 1)
      append(mvec, vals[6]);
    return mvec;
  };

  static std::vector<long> calculateGens(const long* vals)
  {
    std::vector<long> gens;
    gens.push_back(vals[7]);
    if (vals[8] > 1)
      gens.push_back(vals[8]);
    if (vals[9] > 1)
      gens.push_back(vals[9]);
    return gens;
  };

  static std::vector<long> calculateOrds(const long* vals)
  {
    std::vector<long> ords;
    ords.push_back(vals[10]);
    if (abs(vals[11]) > 1)
      ords.push_back(vals[11]);
    if (abs(vals[12]) > 1)
      ords.push_back(vals[12]);
    return ords;
  };

  static long calculateLevels(bool bootstrap, long bitSize)
  {
    return bootstrap
               ? 900
               : 30 * (7 + NTL::NumBits(bitSize + 2)); // that should be enough
  };

  helib::Context& prepareContext(helib::Context& context)
  {
    if (helib_test::verbose) {
      std::cout << "input bitSize=" << bitSize << std::endl;
      if (nthreads > 1)
        std::cout << "  using " << NTL::AvailableThreads() << " threads\n";
      std::cout << "computing key-independent tables..." << std::flush;
    }
    buildModChain(context, L, c, /*willBeBootstrappable=*/bootstrap);
    if (bootstrap) {
      context.makeBootstrappable(mvec, /*t=*/0);
    }
    buildUnpackSlotEncoding(unpackSlotEncoding, *context.ea);

    if (helib_test::verbose) {
      std::cout << " done.\n";
      context.zMStar.printout();
    }

    return context;
  }

  void prepareSecKey(helib::SecKey& secKey)
  {
    if (helib_test::verbose) {
      std::cout << "\ncomputing key-dependent tables..." << std::flush;
    }
    secKey.GenSecKey();
    addSome1DMatrices(secKey); // compute key-switching matrices
    addFrbMatrices(secKey);
    if (bootstrap)
      secKey.genRecryptData();
    if (helib_test::verbose)
      std::cout << " done\n";
  };

  const long prm;
  const long bitSize;
  const bool bootstrap;
  const long seed;
  const long nthreads;
  const long* vals;
  const long p;
  const long m;
  const NTL::Vec<long> mvec;
  const std::vector<long> gens;
  const std::vector<long> ords;
  const long c;
  const long L;
  helib::Context context;
  helib::SecKey secKey;

  GTestBinaryCompare() :
      prm(validatePrm(std::get<0>(GetParam()).prm)),
      bitSize(correctBitSize(5, std::get<0>(GetParam()).bitSize)),
      bootstrap(std::get<0>(GetParam()).bootstrap),
      seed(std::get<0>(GetParam()).seed),
      nthreads(std::get<0>(GetParam()).nthreads),
      vals(mValues[prm]),
      p(vals[0]),
      m(vals[2]),
      mvec(calculateMvec(vals)),
      gens(calculateGens(vals)),
      ords(calculateOrds(vals)),
      c(vals[14]),
      L(calculateLevels(bootstrap, bitSize)),
      context(m, p, /*r=*/1, gens, ords),
      secKey(prepareContext(context)){};

  void SetUp() override
  {
    if (seed)
      NTL::SetSeed(NTL::ZZ(seed));
    if (nthreads > 1)
      NTL::SetNumThreads(nthreads);

    prepareSecKey(secKey);

    helib::activeContext = &context; // make things a little easier sometimes

    helib::setupDebugGlobals(&secKey, context.ea);
  };

  virtual void TearDown() override
  {
#ifdef HELIB_DEBUG
    helib::cleanupDebugGlobals();
#endif
  }

public:
  static void TearDownTestCase()
  {
    if (helib_test::verbose)
      helib::printAllTimers(std::cout);
  };
};

std::vector<helib::zzX> GTestBinaryCompare::unpackSlotEncoding;
constexpr long GTestBinaryCompare::mValues[5][15];

TEST_P(GTestBinaryCompare, comparison)
{
  const helib::EncryptedArray& ea = *context.ea;

  // Choose two random n-bit integers
  long pa = NTL::RandomBits_long(bitSize);
  long pb = NTL::RandomBits_long(bitSize + 1);
  long pMax = std::max(pa, pb);
  long pMin = std::min(pa, pb);
  bool pMu = pa > pb;
  bool pNi = pa < pb;

  // Encrypt the individual bits
  NTL::Vec<helib::Ctxt> eMax, eMin, enca, encb;

  helib::Ctxt mu(secKey), ni(secKey);
  resize(enca, bitSize, mu);
  resize(encb, bitSize + 1, ni);
  for (long i = 0; i <= bitSize; i++) {
    if (i < bitSize)
      secKey.Encrypt(enca[i], NTL::ZZX((pa >> i) & 1));
    secKey.Encrypt(encb[i], NTL::ZZX((pb >> i) & 1));
    if (bootstrap) { // put them at a lower level
      if (i < bitSize)
        enca[i].bringToSet(context.getCtxtPrimes(5));
      encb[i].bringToSet(context.getCtxtPrimes(5));
    }
  }
#ifdef HELIB_DEBUG
  decryptAndPrint((std::cout << " before comparison: "),
                  encb[0],
                  secKey,
                  ea,
                  0);
#endif

  std::vector<long> slotsMin, slotsMax, slotsMu, slotsNi;
  // comparison only
  compareTwoNumbers(mu,
                    ni,
                    helib::CtPtrs_VecCt(enca),
                    helib::CtPtrs_VecCt(encb),
                    false,
                    &unpackSlotEncoding);
  ea.decrypt(mu, secKey, slotsMu);
  ea.decrypt(ni, secKey, slotsNi);
  EXPECT_EQ(std::make_pair(slotsMu[0], slotsNi[0]),
            std::make_pair((long)pMu, (long)pNi))
      << "Comparison (without min max) error: a=" << pa << ", b=" << pb
      << ", mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << std::endl;
  if (helib_test::verbose) {
    std::cout << "Comparison (without min max) succeeded: ";
    std::cout << '(' << pa << ',' << pb << ")=> mu=" << slotsMu[0]
              << ", ni=" << slotsNi[0] << std::endl;
  }

  {
    helib::CtPtrs_VecCt wMin(eMin),
        wMax(eMax); // A wrappers around output vectors
    // comparison with max and min
    compareTwoNumbers(wMax,
                      wMin,
                      mu,
                      ni,
                      helib::CtPtrs_VecCt(enca),
                      helib::CtPtrs_VecCt(encb),
                      false,
                      &unpackSlotEncoding);
    decryptBinaryNums(slotsMax, wMax, secKey, ea);
    decryptBinaryNums(slotsMin, wMin, secKey, ea);
  } // get rid of the wrapper
  ea.decrypt(mu, secKey, slotsMu);
  ea.decrypt(ni, secKey, slotsNi);

  EXPECT_EQ(std::make_tuple(pMax, pMin, pMu, pNi),
            std::make_tuple(slotsMax[0], slotsMin[0], slotsMu[0], slotsNi[0]))
      << "Comparison (with min max) error: a=" << pa << ", b=" << pb
      << ", but min=" << slotsMin[0] << ", max=" << slotsMax[0]
      << ", mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << std::endl;

  if (helib_test::verbose) {
    std::cout << "Comparison (with min max) succeeded: ";
    std::cout << '(' << pa << ',' << pb << ")=>(" << slotsMin[0] << ','
              << slotsMax[0] << "), mu=" << slotsMu[0] << ", ni=" << slotsNi[0]
              << std::endl;
  }

#ifdef HELIB_DEBUG
  const helib::Ctxt* minLvlCtxt = nullptr;
  long minLvl = 1000;
  for (const helib::Ctxt& c : eMax) {
    long lvl = c.logOfPrimeSet();
    if (lvl < minLvl) {
      minLvlCtxt = &c;
      minLvl = lvl;
    }
  }
  decryptAndPrint((std::cout << " after comparison: "),
                  *minLvlCtxt,
                  secKey,
                  ea,
                  0);
  std::cout << std::endl;
#endif
}

TEST_P(GTestBinaryCompare,
       comparingTwoPositiveNumbersInTwosComplementWorksCorrectly)
{
  const helib::EncryptedArray& ea = *context.ea;

  // Choose random n-bit numbers in 2's complement
  long pa = NTL::RandomBits_long(bitSize - 1);
  long pb = NTL::RandomBits_long(bitSize - 1);
  long pMax = std::max(pa, pb);
  long pMin = std::min(pa, pb);
  bool pMu = pa > pb;
  bool pNi = pa < pb;

  // Encrypt the individual bits
  NTL::Vec<helib::Ctxt> eMax, eMin, enca, encb;

  helib::Ctxt mu(secKey), ni(secKey);
  resize(enca, bitSize, mu);
  resize(encb, bitSize, ni);
  for (long i = 0; i < bitSize; i++) {
    secKey.Encrypt(enca[i], NTL::ZZX((pa >> i) & 1));
    secKey.Encrypt(encb[i], NTL::ZZX((pb >> i) & 1));
    if (bootstrap) { // put them at a lower level
      if (i < bitSize)
        enca[i].bringToSet(context.getCtxtPrimes(5));
      encb[i].bringToSet(context.getCtxtPrimes(5));
    }
  }

  std::vector<long> slotsMin, slotsMax, slotsMu, slotsNi;
  // comparison only
  compareTwoNumbers(mu,
                    ni,
                    helib::CtPtrs_VecCt(enca),
                    helib::CtPtrs_VecCt(encb),
                    true,
                    &unpackSlotEncoding);
  ea.decrypt(mu, secKey, slotsMu);
  ea.decrypt(ni, secKey, slotsNi);
  EXPECT_EQ(std::make_pair(slotsMu[0], slotsNi[0]),
            std::make_pair((long)pMu, (long)pNi))
      << "Comparison (without min max) error: a=" << pa << ", b=" << pb
      << ", mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << std::endl;
  if (helib_test::verbose) {
    std::cout << "Comparison (without min max) succeeded: ";
    std::cout << '(' << pa << ',' << pb << ")=> mu=" << slotsMu[0]
              << ", ni=" << slotsNi[0] << std::endl;
  }

  {
    helib::CtPtrs_VecCt wMin(eMin),
        wMax(eMax); // A wrappers around output vectors
    // comparison with max and min
    compareTwoNumbers(wMax,
                      wMin,
                      mu,
                      ni,
                      helib::CtPtrs_VecCt(enca),
                      helib::CtPtrs_VecCt(encb),
                      true,
                      &unpackSlotEncoding);
    decryptBinaryNums(slotsMax, wMax, secKey, ea, true);
    decryptBinaryNums(slotsMin, wMin, secKey, ea, true);
  } // get rid of the wrapper
  ea.decrypt(mu, secKey, slotsMu);
  ea.decrypt(ni, secKey, slotsNi);

  EXPECT_EQ(std::make_tuple(pMax, pMin, pMu, pNi),
            std::make_tuple(slotsMax[0], slotsMin[0], slotsMu[0], slotsNi[0]))
      << "Comparison (with min max) error: a=" << pa << ", b=" << pb
      << ", but min=" << slotsMin[0] << ", max=" << slotsMax[0]
      << ", mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << std::endl;
}

TEST_P(GTestBinaryCompare,
       comparingTwoNegativeNumbersInTwosComplementWorksCorrectly)
{
  const helib::EncryptedArray& ea = *context.ea;

  // Choose random n-bit numbers in 2's complement
  long pa_data = NTL::RandomBits_long(bitSize - 1);
  long pb_data = NTL::RandomBits_long(bitSize - 1);
  pa_data |= (1 << (bitSize - 1));
  pb_data |= (1 << (bitSize - 1));
  long pa = helib::bitSetToLong(pa_data, bitSize);
  long pb = helib::bitSetToLong(pb_data, bitSize);
  long pMax = std::max(pa, pb);
  long pMin = std::min(pa, pb);
  bool pMu = pa > pb;
  bool pNi = pa < pb;

  // Encrypt the individual bits
  NTL::Vec<helib::Ctxt> eMax, eMin, enca, encb;

  helib::Ctxt mu(secKey), ni(secKey);
  resize(enca, bitSize, mu);
  resize(encb, bitSize, ni);
  for (long i = 0; i < bitSize; i++) {
    secKey.Encrypt(enca[i], NTL::ZZX((pa >> i) & 1));
    secKey.Encrypt(encb[i], NTL::ZZX((pb >> i) & 1));
    if (bootstrap) { // put them at a lower level
      if (i < bitSize)
        enca[i].bringToSet(context.getCtxtPrimes(5));
      encb[i].bringToSet(context.getCtxtPrimes(5));
    }
  }

  std::vector<long> slotsMin, slotsMax, slotsMu, slotsNi;
  // comparison only
  compareTwoNumbers(mu,
                    ni,
                    helib::CtPtrs_VecCt(enca),
                    helib::CtPtrs_VecCt(encb),
                    true,
                    &unpackSlotEncoding);
  ea.decrypt(mu, secKey, slotsMu);
  ea.decrypt(ni, secKey, slotsNi);
  EXPECT_EQ(std::make_pair(slotsMu[0], slotsNi[0]),
            std::make_pair((long)pMu, (long)pNi))
      << "Comparison (without min max) error: a=" << pa << ", b=" << pb
      << ", mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << std::endl;
  if (helib_test::verbose) {
    std::cout << "Comparison (without min max) succeeded: ";
    std::cout << '(' << pa << ',' << pb << ")=> mu=" << slotsMu[0]
              << ", ni=" << slotsNi[0] << std::endl;
  }

  {
    helib::CtPtrs_VecCt wMin(eMin),
        wMax(eMax); // A wrappers around output vectors
    // comparison with max and min
    compareTwoNumbers(wMax,
                      wMin,
                      mu,
                      ni,
                      helib::CtPtrs_VecCt(enca),
                      helib::CtPtrs_VecCt(encb),
                      true,
                      &unpackSlotEncoding);
    decryptBinaryNums(slotsMax, wMax, secKey, ea, true);
    decryptBinaryNums(slotsMin, wMin, secKey, ea, true);
  } // get rid of the wrapper
  ea.decrypt(mu, secKey, slotsMu);
  ea.decrypt(ni, secKey, slotsNi);

  EXPECT_EQ(std::make_tuple(pMax, pMin, pMu, pNi),
            std::make_tuple(slotsMax[0], slotsMin[0], slotsMu[0], slotsNi[0]))
      << "Comparison (with min max) error: a=" << pa << ", b=" << pb
      << ", but min=" << slotsMin[0] << ", max=" << slotsMax[0]
      << ", mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << std::endl;
}

TEST_P(GTestBinaryCompare,
       comparingNegativeAndPositiveNumbersInTwosComplementWorksCorrectly)
{
  const helib::EncryptedArray& ea = *context.ea;

  // Choose random n-bit numbers in 2's complement
  long pa_data = NTL::RandomBits_long(bitSize - 1);
  long pb_data = NTL::RandomBits_long(bitSize - 1);
  // pa is non-negative, pb is negative
  pb_data |= (1 << (bitSize - 1));
  long pa = helib::bitSetToLong(pa_data, bitSize);
  long pb = helib::bitSetToLong(pb_data, bitSize);
  long pMax = pa;
  long pMin = pb;
  bool pMu = pa > pb; // Always true
  bool pNi = pa < pb; // Always false

  // Encrypt the individual bits
  NTL::Vec<helib::Ctxt> eMax, eMin, enca, encb;

  helib::Ctxt mu(secKey), ni(secKey);
  resize(enca, bitSize, mu);
  resize(encb, bitSize, ni);
  for (long i = 0; i < bitSize; i++) {
    secKey.Encrypt(enca[i], NTL::ZZX((pa >> i) & 1));
    secKey.Encrypt(encb[i], NTL::ZZX((pb >> i) & 1));
    if (bootstrap) { // put them at a lower level
      if (i < bitSize)
        enca[i].bringToSet(context.getCtxtPrimes(5));
      encb[i].bringToSet(context.getCtxtPrimes(5));
    }
  }

  std::vector<long> slotsMin, slotsMax, slotsMu, slotsNi;
  // comparison only
  compareTwoNumbers(mu,
                    ni,
                    helib::CtPtrs_VecCt(enca),
                    helib::CtPtrs_VecCt(encb),
                    true,
                    &unpackSlotEncoding);
  ea.decrypt(mu, secKey, slotsMu);
  ea.decrypt(ni, secKey, slotsNi);
  EXPECT_EQ(std::make_pair(slotsMu[0], slotsNi[0]),
            std::make_pair((long)pMu, (long)pNi))
      << "Comparison (without min max) error: a=" << pa << ", b=" << pb
      << ", mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << std::endl;
  if (helib_test::verbose) {
    std::cout << "Comparison (without min max) succeeded: ";
    std::cout << '(' << pa << ',' << pb << ")=> mu=" << slotsMu[0]
              << ", ni=" << slotsNi[0] << std::endl;
  }

  {
    helib::CtPtrs_VecCt wMin(eMin),
        wMax(eMax); // A wrappers around output vectors
    // comparison with max and min
    compareTwoNumbers(wMax,
                      wMin,
                      mu,
                      ni,
                      helib::CtPtrs_VecCt(enca),
                      helib::CtPtrs_VecCt(encb),
                      true,
                      &unpackSlotEncoding);
    decryptBinaryNums(slotsMax, wMax, secKey, ea, true);
    decryptBinaryNums(slotsMin, wMin, secKey, ea, true);
  } // get rid of the wrapper
  ea.decrypt(mu, secKey, slotsMu);
  ea.decrypt(ni, secKey, slotsNi);

  EXPECT_EQ(std::make_tuple(pMax, pMin, pMu, pNi),
            std::make_tuple(slotsMax[0], slotsMin[0], slotsMu[0], slotsNi[0]))
      << "Comparison (with min max) error: a=" << pa << ", b=" << pb
      << ", but min=" << slotsMin[0] << ", max=" << slotsMax[0]
      << ", mu=" << slotsMu[0] << ", ni=" << slotsNi[0] << std::endl;
}

INSTANTIATE_TEST_SUITE_P(
    smallParamaterSizesRepeated,
    GTestBinaryCompare,
    ::testing::Combine(::testing::Values(
                           // SLOW
                           Parameters(1, 5, false, 0, 1)
                           // FAST
                           // Parameters(0, 5, false, 0, 1)
                           ),
                       ::testing::Range(0, 3)) // The range is for repeats
);

} // namespace

#if 0
e2=a2+b2+1
e1=a1+b1+1
e0=a0+b0+1, ne0 = a0+b0

e*2 = e2
e*1 = e2 e1
e*0 = e*1 e0, ne*0 = e*1 ne0

a*2 = a2    , b*2 = b2
a*1 = e2 a1 , b*1 = e*1 -e2 - a*1
a*0 = e*1 a0, b*0 = ne*0 -a*0

c2 = a2 b2
c1 = a1 b1
c0 = a0 b0

c*2 = c2
c*1 = e2 c1
c*0 = e*1 c0

A2 = a2,      B2 = b2,      C2 = c2
A1 = a2 +a*1, B1 = b2 +b*1, C1 = c2 + c*1
A0 = A1 +a*0, B0 = B1 +b*0, C0 = C1 + c*0

(a>b) = A0+C0
(a<b) = B0+C0

mx2 = a2+b2+c2
mx1 = a*1+b*1+c1 + a1(a2+c2) + b1(b2+c)

X ab01 = b1 a0
X b10 = b1 b0

mx0 = a0+b0+c*0
     + a0(A1+C1) = a0 A1 + a0(c*2 +a*1 b1) = a0(A1+c*2) + a*1 ab01
     + b0(B1+C1) = b0 B1 + b0(c*2 +a*1 b1) = b0(B1+c*2) + a*1 b10
    = a0(1+A1+c*2) + b0(1+B1+c*2) + a*1(b10+ab01) + c*0

mn0 = c0
     + a0(B1+C1) = a0 B1 + a0(c*2 +a*1 b1) = a0(B1+c*2) + a*1 ab01
     + b0(A1+C1) = b0 A1 + b0(c*2 +a*1 b1) = b0(A1+c*2) + a*1 b10
    = a0(B1+c*2) + b0(A1+C*2) + a*1(b10 + ab01) + c*0

e20 = e2 ne0

      a*1(b10+ab01)=(a2+b2+1)a1 b1(a0+b0) = e20 c1

mx0 = a0(1+A1+c*2) + b0(1+B1+c*2) + e20 c1 + c*0
mn0 = a0(B1+c*2) + b0(A1+c*2) + e20 c1 + c*0

#endif
