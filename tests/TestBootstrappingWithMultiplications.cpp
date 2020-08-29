/* Copyright (C) 2019-2020 IBM Corp.
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
#include <helib/helib.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"
#include <random>

namespace {
struct Parameters
{
  Parameters(long m,
             long p,
             long r,
             long c,
             long bits,
             long t,
             int c_m,
             long n,
             std::vector<long> mvec,
             std::vector<long> gens,
             std::vector<long> ords) :
      m(m),
      p(p),
      r(r),
      c(c),
      bits(bits),
      t(t),
      c_m(c_m),
      n(n),
      mvec(helib::convert<NTL::Vec<long>>(mvec)),
      gens(gens),
      ords(ords)
  {
    if (mvec.empty(), gens.empty() || ords.empty())
      throw helib::LogicError("mvec, gens, and ords must be non-empty");
  };

  const long m;
  const long p;
  const long r;
  const long c;
  const long bits;
  const long t;
  const int c_m;
  const long n;
  const NTL::Vec<long> mvec;
  const std::vector<long> gens;
  const std::vector<long> ords;

  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "m=" << params.m << ","
              << "p=" << params.p << ","
              << "r=" << params.r << ","
              << "c=" << params.c << ","
              << "bits=" << params.bits << ","
              << "skHwt=" << params.t << ","
              << "gens=" << helib::vecToStr(params.gens) << ","
              << "ords=" << helib::vecToStr(params.ords) << ","
              << "mvec=" << params.mvec << ","
              << "c_m=" << params.c_m << ","
              << "computation depth=" << params.n << "}";
  }
};

class TestFatBootstrappingWithMultiplications :
    public ::testing::TestWithParam<Parameters>
{
protected:
  const long n; // Number of multiplications to perform

  helib::Context context;
  helib::SecKey secretKey;
  helib::PubKey publicKey;
  const helib::EncryptedArray& ea;

  TestFatBootstrappingWithMultiplications() :
      n(GetParam().n),
      context(GetParam().m,
              GetParam().p,
              GetParam().r,
              GetParam().gens,
              GetParam().ords),
      secretKey(postContextSetup(context,
                                 GetParam().c_m,
                                 GetParam().bits,
                                 GetParam().c,
                                 GetParam().t,
                                 GetParam().mvec)),
      publicKey(keySetup(secretKey)),
      ea(*(context.ea))
  {}

  static helib::Context& postContextSetup(helib::Context& context,
                                          int c_m,
                                          long bits,
                                          long c,
                                          long t,
                                          NTL::Vec<long> mvec)
  {
    context.zMStar.set_cM(c_m / 100);
    helib::buildModChain(context, bits, c, true, t);
    context.makeBootstrappable(mvec, t, 0);
    return context;
  }

  static helib::SecKey& keySetup(helib::SecKey& secretKey)
  {
    secretKey.GenSecKey();
    helib::addSome1DMatrices(secretKey);
    helib::addFrbMatrices(secretKey);
    secretKey.genRecryptData();
    return secretKey;
  }

  virtual void SetUp() override
  {
    if (helib_test::verbose) {
      std::cout << "m=" << GetParam().m << ", p=" << GetParam().p
                << ", r=" << GetParam().r << ", bits=" << GetParam().bits
                << ", c=" << GetParam().c << ", skHwt=" << GetParam().t
                << ", c_m=" << GetParam().c_m
                << ", depth to compute=" << GetParam().n
                << ", mvec=" << GetParam().mvec
                << ", gens=" << helib::vecToStr(GetParam().gens)
                << ", ords=" << helib::vecToStr(GetParam().ords) << std::endl;
      ea.getPAlgebra().printout();
      std::cout << "ctxtPrimes=" << context.ctxtPrimes
                << ", specialPrimes=" << context.specialPrimes << std::endl
                << std::endl;
    }
    helib::setupDebugGlobals(&secretKey, context.ea);
  }

  virtual void TearDown() override
  {
    if (helib_test::verbose) {
      helib::printAllTimers();
    }
    helib::cleanupDebugGlobals();
  }
};

std::vector<long> generateRandomBinaryVector(long nslots)
{
  std::vector<long> ptxt(nslots);
  std::mt19937 gen(helib_test::random_seed);
  std::uniform_int_distribution<int> coinFlipDist(0, 1);
  for (auto& num : ptxt)
    num = coinFlipDist(gen);
  return ptxt;
}

long multiplyWithRecryption(helib::Ctxt& ctxt,
                            helib::Ctxt& tmp_ctxt,
                            helib::PubKey& publicKey,
                            long depth,
                            long n,
                            bool thin)
{
  long count = 0; // number of multiplications
  HELIB_NTIMER_START(Multiplications);
  while (ctxt.bitCapacity() >= 40 && depth + count < n) {
    if (helib_test::verbose) {
      std::cout << "multiplication " << count + 1 << std::endl;
    }
    ctxt.multiplyBy(tmp_ctxt);
    count += 1;
  }
  HELIB_NTIMER_STOP(Multiplications);
  if (!helib_test::noPrint) {
    helib::CheckCtxt(ctxt, "Before recryption");
  }
  if (depth + count < n) {
    // Time the recryption step
    HELIB_NTIMER_START(Bootstrap);
    // Recrypt/Bootstrap the ctxt
    if (thin) {
      publicKey.thinReCrypt(ctxt);
    } else {
      publicKey.reCrypt(ctxt);
    }
    HELIB_NTIMER_STOP(Bootstrap);
  }
  if (!helib_test::noPrint) {
    helib::CheckCtxt(ctxt, "After recryption");
  }
  return count;
}

void multiplyPtxt(std::vector<long> ptxt, long count, long nslots, long p2r)
{
  std::vector<long> tmp_ptxt(ptxt);
  for (int j = 0; j < count; j++) {
    for (int i = 0; i < nslots; i++) {
      ptxt[i] *= tmp_ptxt[i];
      ptxt[i] = ptxt[i] % p2r;
    }
  }
}

TEST_P(TestFatBootstrappingWithMultiplications,
       correctlyPerformsFatBootstrappingWithNoMultiplications)
{
  const long nslots = ea.size();
  std::vector<long> ptxt(
      generateRandomBinaryVector(nslots)); // Random 0s and 1s
  helib::Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, ptxt);
  if (!helib_test::noPrint) {
    helib::CheckCtxt(ctxt, "Before recryption");
  }
  // Time the recryption step
  HELIB_NTIMER_START(Bootstrap);
  // Recrypt/Bootstrap the ctxt
  publicKey.reCrypt(ctxt);
  HELIB_NTIMER_STOP(Bootstrap);
  if (!helib_test::noPrint) {
    helib::CheckCtxt(ctxt, "After recryption");
  }
  std::vector<long> decrypted(nslots);
  ea.decrypt(ctxt, secretKey, decrypted);

  EXPECT_EQ(decrypted, ptxt);
}

TEST_P(TestFatBootstrappingWithMultiplications,
       correctlyPerformsFatBootstrappingWithMultiplications)
{
  const long nslots = ea.size();
  const long p2r = context.alMod.getPPowR();
  std::vector<long> ptxt(
      generateRandomBinaryVector(nslots)); // Random 0s and 1s
  helib::Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, ptxt);
  helib::Ctxt tmp_ctxt(ctxt);

  long depth = 0; // count to keep track of number of multiplications
  long round = 0; // count for number of bootstraps
  while (depth < n) {
    if (helib_test::verbose) {
      std::cout << "Round " << round << std::endl;
    }
    if (!helib_test::noPrint) {
      helib::CheckCtxt(ctxt, "Before multiplication");
    }
    long count = 0; // count for number of multiplications this round
    // Multiply the ciphertext with itself n times
    // until number of bits falls below threshold
    HELIB_NTIMER_START(RoundTotal);
    count = multiplyWithRecryption(ctxt, tmp_ctxt, publicKey, depth, n, false);
    HELIB_NTIMER_STOP(RoundTotal);
    // Plaintext operation
    // Multiply with itself n times
    multiplyPtxt(ptxt, count, nslots, p2r);
    depth += count; // current depth of circuit computed
    round += 1;     // end of round
  }
  std::vector<long> decrypted(nslots);
  ea.decrypt(ctxt, secretKey, decrypted);

  EXPECT_EQ(decrypted, ptxt);

  if (!helib_test::noPrint) {
    for (const auto& timerName :
         {"Setup", "Multiplications", "Bootstrap", "RoundTotal"}) {
      helib::printNamedTimer(std::cout << std::endl, timerName);
    }
  }
}

// clang-format off
INSTANTIATE_TEST_SUITE_P(typicalParameters, TestFatBootstrappingWithMultiplications, ::testing::Values(
    // Using fast parameters as slow ones takes too long to converge.
    //FAST
    Parameters( 31*41,  2, 1, 2, 380, 64, 100, 30, {  31, 41}, { 1026,  249}, {30, -2}),
    Parameters(7*5*37, 17, 1, 3, 600, 64, 100, 40, {7, 5, 37}, {  556, 1037}, { 6,  4})
    //SLOW
    //Parameters( 31775,  2, 1, 2, 580, 64, 100, 50, {  41,775}, { 6976,24806}, {40, 30}),
    //Parameters( 35113,  2, 1, 2, 580, 64, 100, 50, {  37,949}, {16134, 8548}, {36, 24})
));
// clang-format on

class TestThinBootstrappingWithMultiplications :
    public ::testing::TestWithParam<Parameters>
{
protected:
  const long n; // Number of multiplications to perform

  helib::Context context;
  helib::SecKey secretKey;
  helib::PubKey publicKey;
  const helib::EncryptedArray& ea;

  TestThinBootstrappingWithMultiplications() :
      n(GetParam().n),
      context(GetParam().m,
              GetParam().p,
              GetParam().r,
              GetParam().gens,
              GetParam().ords),
      secretKey(postContextSetup(context,
                                 GetParam().c_m,
                                 GetParam().bits,
                                 GetParam().c,
                                 GetParam().t,
                                 GetParam().mvec)),
      publicKey(keySetup(secretKey)),
      ea(*(context.ea))
  {}

  static helib::Context& postContextSetup(helib::Context& context,
                                          int c_m,
                                          long bits,
                                          long c,
                                          long t,
                                          NTL::Vec<long> mvec)
  {
    context.zMStar.set_cM(c_m / 100);
    helib::buildModChain(context, bits, c, true, t);
    context.makeBootstrappable(mvec, t, 0, false);
    return context;
  }

  static helib::SecKey& keySetup(helib::SecKey& secretKey)
  {
    secretKey.GenSecKey();
    helib::addSome1DMatrices(secretKey);
    helib::addFrbMatrices(secretKey);
    secretKey.genRecryptData();
    return secretKey;
  }

  virtual void SetUp() override
  {
    if (helib_test::verbose) {
      std::cout << "m=" << GetParam().m << ", p=" << GetParam().p
                << ", r=" << GetParam().r << ", bits=" << GetParam().bits
                << ", c=" << GetParam().c << ", skHwt=" << GetParam().t
                << ", c_m=" << GetParam().c_m
                << ", depth to compute=" << GetParam().n
                << ", mvec=" << GetParam().mvec
                << ", gens=" << helib::vecToStr(GetParam().gens)
                << ", ords=" << helib::vecToStr(GetParam().ords) << std::endl;
      ea.getPAlgebra().printout();
      std::cout << "ctxtPrimes=" << context.ctxtPrimes
                << ", specialPrimes=" << context.specialPrimes << std::endl
                << std::endl;
    }

    helib::setupDebugGlobals(&secretKey, context.ea);
  }

  virtual void TearDown() override
  {
    if (helib_test::verbose) {
      helib::printAllTimers();
    }
    helib::cleanupDebugGlobals();
  }
};

TEST_P(TestThinBootstrappingWithMultiplications,
       correctlyPerformsThinBootstrappingWithNoMultiplications)
{
  const long nslots = ea.size();
  std::vector<long> ptxt(
      generateRandomBinaryVector(nslots)); // Random 0s and 1s
  helib::Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, ptxt);
  if (!helib_test::noPrint) {
    helib::CheckCtxt(ctxt, "Before recryption");
  }
  // Time the recryption step
  HELIB_NTIMER_START(Bootstrap);
  // Recrypt/Bootstrap the ctxt
  publicKey.thinReCrypt(ctxt);
  HELIB_NTIMER_STOP(Bootstrap);
  if (!helib_test::noPrint) {
    helib::CheckCtxt(ctxt, "After recryption");
  }
  std::vector<long> decrypted(nslots);
  ea.decrypt(ctxt, secretKey, decrypted);

  EXPECT_EQ(decrypted, ptxt);
}

TEST_P(TestThinBootstrappingWithMultiplications,
       correctlyPerformsThinBootstrappingWithMultiplications)
{
  const long nslots = ea.size();
  const long p2r = context.alMod.getPPowR();
  std::vector<long> ptxt(
      generateRandomBinaryVector(nslots)); // Random 0s and 1s
  helib::Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, ptxt);
  helib::Ctxt tmp_ctxt(ctxt);

  long depth = 0; // count to keep track of number of multiplications
  long round = 0; // count for number of bootstraps
  while (depth < n) {
    if (helib_test::verbose) {
      std::cout << "Round " << round << std::endl;
    }
    if (!helib_test::noPrint) {
      helib::CheckCtxt(ctxt, "Before multiplication");
    }
    long count = 0; // count for number of multiplications this round
    // Multiply the ciphertext with itself n times
    // until number of bits falls below threshold
    HELIB_NTIMER_START(RoundTotal);
    count = multiplyWithRecryption(ctxt, tmp_ctxt, publicKey, depth, n, true);
    HELIB_NTIMER_STOP(RoundTotal);
    // Plaintext operation
    // Multiply with itself n times
    multiplyPtxt(ptxt, count, nslots, p2r);
    depth += count; // current depth of circuit computed
    round += 1;     // end of round
  }
  std::vector<long> decrypted(nslots);
  ea.decrypt(ctxt, secretKey, decrypted);

  EXPECT_EQ(decrypted, ptxt);

  if (!helib_test::noPrint) {
    for (const auto& timerName :
         {"Setup", "Multiplications", "Bootstrap", "RoundTotal"}) {
      helib::printNamedTimer(std::cout << std::endl, timerName);
    }
  }
}

// clang-format off
INSTANTIATE_TEST_SUITE_P(typicalParameters, TestThinBootstrappingWithMultiplications, ::testing::Values(
    // Using fast parameters as slow ones takes too long to converge.
    //FAST
    Parameters( 31*41,  2, 1, 2, 320, 64, 100, 25,   {31, 41}, { 1026,  249}, {30, -2}),
    Parameters(7*5*37, 17, 1, 3, 500, 64, 100, 30, {7, 5, 37}, {  556, 1037}, { 6,  4})
    //SLOW
    //Parameters( 31775,  2, 1, 2, 580, 64, 100, 50,   {41,775}, { 6976,24806}, {40, 30}),
    //Parameters( 35113,  2, 1, 2, 580, 64, 100, 50,   {37,949}, {16134, 8548}, {36, 24})
));
// clang-format on

} // namespace
