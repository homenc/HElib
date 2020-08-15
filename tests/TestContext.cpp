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
#include <cmath> // isinf
#include <helib/helib.h>

#include "test_common.h"
#include "gtest/gtest.h"

// TODO - Currently does not cover well the Context object.

// Needed for testing
namespace helib {
extern NTL::ZZX getG(const EncryptedArray& ea);
}

namespace {

struct BGVParameters
{
  BGVParameters(unsigned m, unsigned p, unsigned r) : m(m), p(p), r(r){};

  const unsigned m;
  const unsigned p;
  const unsigned r;

  friend std::ostream& operator<<(std::ostream& os, const BGVParameters& params)
  {
    return os << "{"
              << "m = " << params.m << ", "
              << "p = " << params.p << ", "
              << "r = " << params.r << "}";
  }
};

class TestContext : public ::testing::TestWithParam<BGVParameters>
{
protected:
  TestContext() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      context(std::make_shared<helib::Context>(m, p, r))
  {}

  const unsigned long m;
  const unsigned long p;
  const unsigned long r;

  const std::shared_ptr<helib::Context> context;
};

TEST_P(TestContext,
       ContextThrowExceptionWhenCalculatingSecurityBeforeModchainBuilt)
{
  EXPECT_THROW(context->securityLevel(), helib::LogicError);
}

TEST_P(TestContext, ContextCalculatingSecurityAfterModchainBuilt)
{
  buildModChain(*context, /*bits=*/100, /*c=*/2);
  double result = context->securityLevel();
  EXPECT_FALSE(std::isinf(result));
}

TEST_P(TestContext, hasCorrectSlotRingWhenConstructed)
{
  EXPECT_EQ(context->slotRing->p, p);
  EXPECT_EQ(context->slotRing->r, r);
  EXPECT_EQ(context->slotRing->p2r, pow(p, r));
  EXPECT_EQ(context->slotRing->G, helib::getG(*(context->ea)));
}

TEST_P(TestContext, buildModChainThrowsWhenBitsIsZero)
{
  EXPECT_THROW(helib::buildModChain(*context, /*bits=*/0, /*c=*/2),
               helib::InvalidArgument);
}

TEST_P(TestContext, calculateBitSizeOfQ)
{

  long bits = 1016;
  buildModChain(*context, bits, /*c=*/2);
  long bitsize = context->bitSizeOfQ();

  // Get the primes used by HElib.
  helib::IndexSet fullPrimes = context->fullPrimes();
  long calcFullPrimesBitSize =
      ceil(context->logOfProduct(fullPrimes) / log(2.0));

  long calcCtxtPrimesBitSize =
      ceil(context->logOfProduct(context->ctxtPrimes) / log(2.0));

  // Check if the ctxtPrimes are the bits we asked for.
  // Will be close but not exact.
  EXPECT_NEAR(calcCtxtPrimesBitSize, bits, 0.04 * bits);

  // Check if total primes bitsize is same as calulated here.
  EXPECT_EQ(calcFullPrimesBitSize, bitsize);
}

TEST(TestContext, securityHasLowerBoundOfZero)
{
  // Security = -109
  helib::Context large_negative_sec_context(/*m=*/17, /*p=*/2, /*r=*/1);
  buildModChain(large_negative_sec_context, /*bits=*/100, /*c=*/2);
  double large_negative_result = large_negative_sec_context.securityLevel();
  EXPECT_DOUBLE_EQ(large_negative_result, 0.0);

  // Security = -6
  helib::Context small_negative_sec_context(/*m=*/2501, /*p=*/2, /*r=*/1);
  buildModChain(small_negative_sec_context, /*bits=*/100, /*c=*/2);
  double small_negative_result = small_negative_sec_context.securityLevel();
  EXPECT_DOUBLE_EQ(small_negative_result, 0.0);

  // Security 0.6
  helib::Context negligible_sec_context(/*m=*/3205, /*p=*/2, /*r=*/1);
  buildModChain(negligible_sec_context, /*bits=*/100, /*c=*/2);
  double negl_result = negligible_sec_context.securityLevel();
  EXPECT_NEAR(negl_result, 0.656656, 0.000001);

  // Security 43
  helib::Context small_sec_context(/*m=*/5215, /*p=*/2, /*r=*/1);
  buildModChain(small_sec_context, /*bits=*/100, /*c=*/2);
  double small_result = small_sec_context.securityLevel();
  EXPECT_NEAR(small_result, 43.4318, 0.0001);
}

INSTANTIATE_TEST_SUITE_P(variousParameters,
                         TestContext,
                         ::testing::Values(BGVParameters(17, 2, 1)));

} // namespace
