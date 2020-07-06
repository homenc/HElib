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
#include <NTL/BasicThreadPool.h>
#include <helib/intraSlot.h>
#include <helib/tableLookup.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

struct Parameters
{
  Parameters(long prm,
             long bitSize,
             long outSize,
             long nTests,
             bool bootstrap,
             long seed,
             long nthreads) :
      prm(prm),
      bitSize(bitSize),
      outSize(outSize),
      nTests(nTests),
      bootstrap(bootstrap),
      seed(seed),
      nthreads(nthreads){};

  long prm;       // parameter size (0-tiny,...,4-huge)
  long bitSize;   // bitSize of input integers (<=32)
  long outSize;   // bitSize of output integers
  long nTests;    // number of tests to run
  bool bootstrap; // test multiplication with bootstrapping
  long seed;      // PRG seed
  long nthreads;  // number of threads

  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "prm=" << params.prm << ","
              << "bitSize=" << params.bitSize << ","
              << "outSize=" << params.outSize << ","
              << "nTests=" << params.nTests << ","
              << "bootstrap=" << params.bootstrap << ","
              << "seed=" << params.seed << ","
              << "nthreads=" << params.nthreads << "}";
  };
};

class GTestTableLookup : public ::testing::TestWithParam<Parameters>
{
protected:
  // clang-format off
  static constexpr long mValues[][15] = {
  //  {p,phi(m),    m,  d, m1, m2, m3,   g1,   g2,   g3,ord1,ord2,ord3,  B, c}
      {2,    48,  105, 12,  3, 35,  0,   71,   76,    0,   2,   2,   0, 25, 2},
      {2,   600, 1023, 10, 11, 93,  0,  838,  584,    0,  10,   6,   0, 25, 2},
      {2,  2304, 4641, 24,  7,  3,221, 3979, 3095, 3760,   6,   2,  -8, 25, 3},
      {2, 15004,15709, 22, 23,683,  0, 4099,13663,    0,  22,  31,   0, 25, 3},
      {2, 27000,32767, 15, 31,  7,151,11628,28087,25824,  30,   6, -10, 28, 4}
  };
  // clang-format on

  // Utility encryption/decryption methods
  static void encryptIndex(std::vector<helib::Ctxt>& ei,
                           long index,
                           const helib::SecKey& sKey)
  {
    for (long i = 0; i < helib::lsize(ei); i++)
      sKey.Encrypt(ei[i], NTL::to_ZZX((index >> i) & 1)); // i'th bit of index
  }

  static long decryptIndex(std::vector<helib::Ctxt>& ei,
                           const helib::SecKey& sKey)
  {
    long num = 0;
    for (long i = 0; i < helib::lsize(ei); i++) {
      NTL::ZZX poly;
      sKey.Decrypt(poly, ei[i]);
      num += to_long(NTL::ConstTerm(poly)) << i;
    }
    return num;
  }

  static long validatePrm(const long prm)
  {
    if (prm < 0 || prm >= 5)
      throw std::invalid_argument("Invalid prm value");
    return prm;
  };

  static long validateBitSize(const long bitSize)
  {
    if (bitSize > 7)
      throw std::invalid_argument("Invalid bitSize value: must be <=7");
    else if (bitSize <= 0)
      throw std::invalid_argument("Invalid bitSize value: must be >0");
    return bitSize;
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

  static long calculateLevels(const bool bootstrap, const long bitSize)
  {
    long L;
    if (bootstrap)
      L = 900; // that should be enough
    else
      L = 30 * (5 + bitSize);
    return L;
  };

  static void printPreContextPrepDiagnostics(const long bitSize,
                                             const long outSize,
                                             const long nTests,
                                             const long nthreads)
  {
    if (helib_test::verbose) {
      std::cout << "input bitSize=" << bitSize
                << ", output size bound=" << outSize << ", running " << nTests
                << " tests for each function\n";
      if (nthreads > 1)
        std::cout << "  using " << NTL::AvailableThreads() << " threads\n";
      std::cout << "computing key-independent tables..." << std::flush;
    }
  };

  static void printPostContextPrepDiagnostics(const helib::Context& context,
                                              const long L)
  {
    if (helib_test::verbose) {
      std::cout << " done.\n";
      context.zMStar.printout();
      std::cout << " L=" << L << std::endl;
    };
  }

  // Not static as many instance variables are required.
  helib::Context& prepareContext(helib::Context& context)
  {
    printPreContextPrepDiagnostics(bitSize, outSize, nTests, nthreads);
    helib::buildModChain(context, L, c, /*willBeBootstrappable*/ bootstrap);
    if (bootstrap) {
      context.makeBootstrappable(mvec, /*t=*/0);
    }
    helib::buildUnpackSlotEncoding(unpackSlotEncoding, *context.ea);
    printPostContextPrepDiagnostics(context, L);
    return context;
  };

  static void prepareSecKey(helib::SecKey& secretKey, const bool bootstrap)
  {
    if (helib_test::verbose)
      std::cout << "\ncomputing key-dependent tables..." << std::flush;
    secretKey.GenSecKey();
    helib::addSome1DMatrices(secretKey); // compute key-switching matrices
    helib::addFrbMatrices(secretKey);
    if (bootstrap)
      secretKey.genRecryptData();
    if (helib_test::verbose)
      std::cout << " done\n";
  };

  static void setSeedIfNeeded(const long seed)
  {
    if (seed)
      NTL::SetSeed(NTL::ZZ(seed));
    ;
  };

  static void setThreadsIfNeeded(const long nthreads)
  {
    if (nthreads > 1)
      NTL::SetNumThreads(nthreads);
  };

  GTestTableLookup() :
      prm(validatePrm(GetParam().prm)),
      bitSize(validateBitSize(GetParam().bitSize)),
      outSize(GetParam().outSize),
      nTests(GetParam().nTests),
      bootstrap(GetParam().bootstrap),
      seed((setSeedIfNeeded(GetParam().seed), GetParam().seed)),
      nthreads((setThreadsIfNeeded(GetParam().nthreads), GetParam().nthreads)),
      vals(mValues[prm]),
      p(vals[0]),
      m(vals[2]),
      mvec(calculateMvec(vals)),
      gens(calculateGens(vals)),
      ords(calculateOrds(vals)),
      c(vals[14]),
      L(calculateLevels(bootstrap, bitSize)),
      context(m, p, /*r=*/1, gens, ords),
      secretKey(prepareContext(context))
  {
    prepareSecKey(secretKey, bootstrap);
  };

  std::vector<helib::zzX> unpackSlotEncoding;
  const long prm;
  const long bitSize;
  const long outSize;
  const long nTests;
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
  helib::SecKey secretKey;

  void SetUp() override
  {
    helib::activeContext = &context; // make things a little easier sometimes
    helib::setupDebugGlobals(&secretKey, context.ea);
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
    if (helib_test::verbose) {
      helib::printAllTimers(std::cout);
    }
  };
};

constexpr long GTestTableLookup::mValues[][15];

TEST_P(GTestTableLookup, lookupFunctionsCorrectly)
{
  // Build a table s.t. T[i] = 2^{outSize -1}/(i+1), i=0,...,2^bitSize -1
  std::vector<helib::zzX> T;
  helib::buildLookupTable(
      T,
      [](double x) { return 1 / (x + 1.0); },
      bitSize,
      /*scale_in=*/0,
      /*sign_in=*/0,
      outSize,
      /*scale_out=*/1 - outSize,
      /*sign_out=*/0,
      *(secretKey.getContext().ea));

  ASSERT_EQ(helib::lsize(T), 1L << bitSize);
  for (long i = 0; i < helib::lsize(T); i++) {
    helib::Ctxt c(secretKey);
    std::vector<helib::Ctxt> ei(bitSize, c);
    encryptIndex(ei, i, secretKey); // encrypt the index
    helib::tableLookup(c,
                       T,
                       helib::CtPtrs_vectorCt(ei)); // get the encrypted entry
    // decrypt and compare
    NTL::ZZX poly;
    secretKey.Decrypt(poly, c); // decrypt
    helib::zzX poly2;
    helib::convert(poly2, poly); // convert to zzX
    EXPECT_EQ(poly2, T[i]) << "testLookup error: decrypted T[" << i << "]\n";
  }
}

TEST_P(GTestTableLookup, writeinFunctionsCorrectly)
{
  long tSize = 1L << bitSize; // table size

  // encrypt a random table
  std::vector<long> pT(tSize, 0);                            // plaintext table
  std::vector<helib::Ctxt> T(tSize, helib::Ctxt(secretKey)); // encrypted table
  for (long i = 0; i < bitSize; i++) {
    long bit = NTL::RandomBits_long(1); // a random bit
    secretKey.Encrypt(T[i], NTL::to_ZZX(bit));
    pT[i] = bit;
  }

  // Add 1 to 20 random entries in the table
  for (long count = 0; count < nTests; count++) {
    // encrypt a random index into the table
    long index = NTL::RandomBnd(tSize); // 0 <= index < tSize
    std::vector<helib::Ctxt> I(bitSize, helib::Ctxt(secretKey));
    encryptIndex(I, index, secretKey);

    // do the table write-in
    tableWriteIn(helib::CtPtrs_vectorCt(T),
                 helib::CtPtrs_vectorCt(I),
                 &unpackSlotEncoding);
    pT[index]++; // add 1 to entry 'index' in the plaintext table
  }

  // Check that the ciphertext and plaintext tables still match
  for (int i = 0; i < tSize; i++) {
    NTL::ZZX poly;
    secretKey.Decrypt(poly, T[i]);
    long decrypted = to_long(NTL::ConstTerm(poly));
    long p = T[i].getPtxtSpace();
    ASSERT_EQ((pT[i] - decrypted) % p, 0) // should be equal mod p
        << "testWritein error: decrypted T[" << i << "]=" << decrypted
        << " but should be " << pT[i] << " (mod " << p << ")\n";
  }
}

INSTANTIATE_TEST_SUITE_P(typicalParameters,
                         GTestTableLookup,
                         ::testing::Values(
                             // SLOW
                             Parameters(1, 5, 0, 3, false, 0, 1)
                             // FAST
                             // Parameters(0, 5, 0, 3, false, 0, 1)
                             ));

} // namespace
