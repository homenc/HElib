/* Copyright (C) 2020 IBM Corp.
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
#include <helib/set.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

struct BGVParameters
{
  BGVParameters(unsigned m, unsigned p, unsigned r, unsigned bits) :
      m(m), p(p), r(r), bits(bits){};

  const unsigned m;
  const unsigned p;
  const unsigned r;
  const unsigned bits;

  friend std::ostream& operator<<(std::ostream& os, const BGVParameters& params)
  {
    return os << "{"
              << "m = " << params.m << ", "
              << "p = " << params.p << ", "
              << "r = " << params.r << ", "
              << "bits = " << params.bits << "}";
  }
};

class TestSet : public ::testing::TestWithParam<BGVParameters>
{
protected:
  const unsigned long m;
  const unsigned long p;
  const unsigned long r;
  const unsigned long bits;
  helib::Context context;
  helib::SecKey secretKey;
  helib::PubKey publicKey;
  const helib::EncryptedArray& ea;

  TestSet() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      bits(GetParam().bits),
      context(m, p, r),
      secretKey((buildModChain(context, bits), context)),
      publicKey((secretKey.GenSecKey(),
                 addSome1DMatrices(secretKey),
                 addFrbMatrices(secretKey),
                 secretKey)),
      ea(*(context.ea))
  {}

  virtual void SetUp() override
  {
    if (helib_test::verbose) {
      ea.getPAlgebra().printout();
      std::cout << "r = " << context.alMod.getR() << std::endl;
      std::cout << "ctxtPrimes=" << context.ctxtPrimes
                << ", specialPrimes=" << context.specialPrimes << std::endl
                << std::endl;
    }

    helib::setupDebugGlobals(&secretKey, context.ea);
  }

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }
};

static inline NTL::ZZX polyFromBinary(long x)
{
  constexpr long N = 16;
  std::bitset<N> b(x);
  NTL::ZZX poly;

  for (long i = 0; i < N; ++i) {
    SetCoeff(poly, i, b[i]);
  }

  return poly;
}

TEST_P(TestSet, calculateSetIntersectionPtxtWorksCorrectly)
{
  std::cout << "Available Threads: " << NTL::AvailableThreads() << std::endl;

  constexpr long N = 1 << 9;

  // Create server set - simple array
  long cnt = 0;
  std::vector<NTL::ZZX> server_set(N);
  std::generate(server_set.begin(), server_set.end(), [&cnt]() {
    return polyFromBinary(++cnt);
  });

  // create client set - single ctxt
  std::vector<NTL::ZZX> query_numbers = {
      polyFromBinary(1),
      polyFromBinary(5),
      polyFromBinary(501),
      polyFromBinary(2048),
  };

  helib::Ptxt<helib::BGV> client_set(context, query_numbers);

  // Now for call the psi op
  auto result = calculateSetIntersection(client_set, server_set);

  // create expected result set - single ptxt
  std::vector<NTL::ZZX> expected_matches = {
      polyFromBinary(1),
      polyFromBinary(5),
      polyFromBinary(501),
  };

  helib::Ptxt<helib::BGV> expected_result(context, expected_matches);

  // test
  EXPECT_EQ(result, expected_result);
}

TEST_P(TestSet, calculateSetIntersectionCtxtWorksCorrectly)
{
  std::cout << "Available Threads: " << NTL::AvailableThreads() << std::endl;

  constexpr long N = 1 << 4;

  // Create server set - simple array
  long cnt = 0;
  std::vector<NTL::ZZX> server_set(N);
  std::generate(server_set.begin(), server_set.end(), [&cnt]() {
    return polyFromBinary(++cnt);
  });

  // create client set - single ctxt
  std::vector<NTL::ZZX> query_numbers = {
      polyFromBinary(1),
      polyFromBinary(5),
      polyFromBinary(15),
      polyFromBinary(2048),
  };

  helib::Ptxt<helib::BGV> client_set(context, query_numbers);

  helib::Ctxt encrypted_client_set(publicKey);
  publicKey.Encrypt(encrypted_client_set, client_set);

  // Now for call the psi op
  auto encrypted_result =
      calculateSetIntersection(encrypted_client_set, server_set);

  // Decrypt result
  helib::Ptxt<helib::BGV> result(context);
  secretKey.Decrypt(result, encrypted_result);

  // create expected result set - single ptxt
  std::vector<NTL::ZZX> expected_matches = {
      polyFromBinary(1),
      polyFromBinary(5),
      polyFromBinary(15),
  };

  helib::Ptxt<helib::BGV> expected_result(context, expected_matches);

  // test
  EXPECT_EQ(result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(variousParameters,
                         TestSet,
                         ::testing::Values(BGVParameters(771, 2, 1, 700)));

} // namespace
