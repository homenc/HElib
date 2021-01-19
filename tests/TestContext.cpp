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

struct CKKSParameters
{
  CKKSParameters(unsigned m, unsigned r) : m(m), r(r){};

  const unsigned m;
  const unsigned r;

  friend std::ostream& operator<<(std::ostream& os,
                                  const CKKSParameters& params)
  {
    return os << "{"
              << "m = " << params.m << ", "
              << "r = " << params.r << "}";
  }
};

class TestContextBGV : public ::testing::TestWithParam<BGVParameters>
{
protected:
  TestContextBGV() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      context(helib::ContextBuilder<helib::BGV>()
                  .m(m)
                  .p(p)
                  .r(r)
                  .buildModChain(false)
                  .buildPtr())
  {}

  const unsigned long m;
  const unsigned long p;
  const unsigned long r;

  const std::shared_ptr<helib::Context> context;
};

class TestContextCKKS : public ::testing::TestWithParam<CKKSParameters>
{
protected:
  TestContextCKKS() :
      m(GetParam().m),
      r(GetParam().r),
      context(helib::ContextBuilder<helib::CKKS>()
                  .m(m)
                  .precision(r)
                  .buildModChain(false)
                  .buildPtr())
  {}

  const unsigned long m;
  const unsigned long r;
  const std::shared_ptr<helib::Context> context;
};

TEST_P(TestContextBGV,
       ContextThrowExceptionWhenCalculatingSecurityBeforeModchainBuilt)
{
  EXPECT_THROW(context->securityLevel(), helib::LogicError);
}

TEST_P(TestContextBGV, contextEquals)
{
  helib::Context someOtherContext = helib::ContextBuilder<helib::BGV>()
                                        .m(17)
                                        .p(2)
                                        .r(1)
                                        .bits(100)
                                        .c(2)
                                        .build();

  context->buildModChain(/*bits=*/100, /*c=*/2);

  EXPECT_EQ(*context, *context);
  EXPECT_EQ(*context, someOtherContext);

  EXPECT_EQ(context->getM(), 17);            // m
  EXPECT_EQ(context->getP(), 2);             // p
  EXPECT_EQ(context->getAlMod().getR(), 1);  // r
  EXPECT_EQ(context->getDigits().size(), 2); // c
  EXPECT_GT(context->numPrimes(), 0);
  EXPECT_EQ(context->numPrimes(), someOtherContext.numPrimes());
  for (long i = 0; i < context->numPrimes(); i++) {
    const helib::Cmodulus& m1 = context->ithModulus(i);
    const helib::Cmodulus& m2 = someOtherContext.ithModulus(i);
    EXPECT_EQ(m1.getQ(), m2.getQ()) << " index: " << i;
  }
  EXPECT_EQ(context->getSmallPrimes(), someOtherContext.getSmallPrimes());
  EXPECT_EQ(context->getCtxtPrimes(), someOtherContext.getCtxtPrimes());
  EXPECT_EQ(context->getSpecialPrimes(), someOtherContext.getSpecialPrimes());
  for (std::size_t i = 0; i < context->getDigits().size(); ++i) {
    EXPECT_EQ(context->getDigit(i), someOtherContext.getDigit(i))
        << " index: " << i;
  }
  EXPECT_EQ(context->getStdev(), someOtherContext.getStdev());
  EXPECT_EQ(context->getScale(), someOtherContext.getScale());
  EXPECT_EQ(context->getRcData(), someOtherContext.getRcData());
}

TEST_P(TestContextBGV, contextNotEquals)
{
  helib::Context someOtherContext = helib::ContextBuilder<helib::BGV>()
                                        .m(13)
                                        .p(3)
                                        .r(2)
                                        .scale(6)
                                        .stdev(3.0)
                                        .bits(200)
                                        .c(2)
                                        .bootstrappable()
                                        .mvec({13})
                                        .build();
  context->buildModChain(/*bits=*/100, /*c=*/3);

  EXPECT_NE(*context, someOtherContext);
  EXPECT_NE(context->numPrimes(), someOtherContext.numPrimes());
  bool atLeastOneDoesNotEqual = false;
  for (long i = 0; i < context->numPrimes(); i++) {
    const helib::Cmodulus& m1 = context->ithModulus(i);
    const helib::Cmodulus& m2 = someOtherContext.ithModulus(i);
    if (m1.getQ() != m2.getQ()) {
      atLeastOneDoesNotEqual = true;
      break;
    }
  }
  EXPECT_TRUE(atLeastOneDoesNotEqual);
  // This only checks the handles, not the primes themselves.
  EXPECT_EQ(context->getSmallPrimes(), someOtherContext.getSmallPrimes());
  EXPECT_NE(context->getCtxtPrimes(), someOtherContext.getCtxtPrimes());
  EXPECT_NE(context->getSpecialPrimes(), someOtherContext.getSpecialPrimes());
  for (std::size_t i = 0; i < context->getDigits().size(); ++i) {
    EXPECT_NE(context->getDigit(i), someOtherContext.getDigit(i))
        << " index: " << i;
  }
  EXPECT_NE(context->getStdev(), someOtherContext.getStdev());
  EXPECT_NE(context->getScale(), someOtherContext.getScale());
  EXPECT_NE(context->getRcData(), someOtherContext.getRcData());
}

TEST_P(TestContextBGV, ContextCalculatingSecurityAfterModchainBuilt)
{
  context->buildModChain(/*bits=*/100, /*c=*/2);
  double result = context->securityLevel();
  EXPECT_FALSE(std::isinf(result));
}

TEST_P(TestContextBGV, hasCorrectSlotRingWhenConstructed)
{
  EXPECT_EQ(context->getSlotRing()->p, p);
  EXPECT_EQ(context->getSlotRing()->r, r);
  EXPECT_EQ(context->getSlotRing()->p2r, pow(p, r));
  EXPECT_EQ(context->getSlotRing()->G, helib::getG(context->getEA()));
}

TEST_P(TestContextBGV, buildModChainThrowsWhenBitsIsZero)
{
  EXPECT_THROW(context->buildModChain(/*bits=*/0, /*c=*/2),
               helib::InvalidArgument);
}

TEST_P(TestContextBGV, calculateBitSizeOfQ)
{
  long bits = 1016;
  context->buildModChain(bits, /*c=*/2);
  long bitsize = context->bitSizeOfQ();

  // Get the primes used by HElib.
  helib::IndexSet fullPrimes = context->fullPrimes();
  long calcFullPrimesBitSize =
      ceil(context->logOfProduct(fullPrimes) / log(2.0));

  long calcCtxtPrimesBitSize =
      ceil(context->logOfProduct(context->getCtxtPrimes()) / log(2.0));

  // Check if the ctxtPrimes are the bits we asked for.
  // Will be close but not exact.
  EXPECT_NEAR(calcCtxtPrimesBitSize, bits, 0.04 * bits);

  // Check if total primes bitsize is same as calulated here.
  EXPECT_EQ(calcFullPrimesBitSize, bitsize);
}

TEST(TestContextBGV, securityHasLowerBoundOfZero)
{
// VJS-FIXME: this kind of test is not very good, as
// security level calculations are subject to change...
// I got rid of this for now
#if 0
  // Security = -14
  helib::Context small_negative_sec_context(/*m=*/17, /*p=*/2, /*r=*/1);
  buildModChain(small_negative_sec_context, /*bits=*/100, /*c=*/2);
  double small_negative_result = small_negative_sec_context.securityLevel();
  EXPECT_DOUBLE_EQ(small_negative_result, 0.0);

  // Security = 13
  helib::Context small_sec_context(/*m=*/2501, /*p=*/2, /*r=*/1);
  buildModChain(small_sec_context, /*bits=*/200, /*c=*/2);
  double small_result = small_sec_context.securityLevel();
  EXPECT_NEAR(small_result, 8.31563, 0.0001);
#endif
}

/* The following tests are for the ContextBuilder class */

TEST(TestContextBGV, contextBuilderWithDefaultArguments)
{
  helib::Context context_built = helib::ContextBuilder<helib::BGV>().build();

  helib::Context expected_default_context =
      helib::ContextBuilder<helib::BGV>().m(3).p(2).r(1).bits(300).c(3).build();

  // Making sure the number of columns is not clipped.
  EXPECT_GT(expected_default_context.getCtxtPrimes().card(),
            expected_default_context.getDigits().size());
  EXPECT_GT(context_built.getCtxtPrimes().card(),
            context_built.getDigits().size());
  EXPECT_EQ(context_built.getDigits().size(), 3);

  EXPECT_FALSE(context_built.isBootstrappable());
  EXPECT_GT(context_built.numPrimes(), 0);
  EXPECT_EQ(context_built, expected_default_context);
}

TEST(TestContextBGV, contextBuilderBuildsPointer)
{
  std::unique_ptr<helib::Context> context_built{
      helib::ContextBuilder<helib::BGV>().buildPtr()};

  helib::Context expected_default_context =
      helib::ContextBuilder<helib::BGV>().m(3).p(2).r(1).bits(300).c(3).build();

  // Making sure the number of columns is not clipped.
  EXPECT_GT(expected_default_context.getCtxtPrimes().card(),
            expected_default_context.getDigits().size());
  EXPECT_GT(context_built->getCtxtPrimes().card(),
            context_built->getDigits().size());
  EXPECT_EQ(context_built->getDigits().size(), 3);

  EXPECT_FALSE(context_built->isBootstrappable());
  EXPECT_GT(context_built->numPrimes(), 0);
  EXPECT_EQ(*context_built, expected_default_context);
}

TEST(TestContextBGV, contextBuilderClipsDigitsSizeWithSmallBits)
{
  long c = 4; // Columns of SKMs
  helib::Context context_built =
      helib::ContextBuilder<helib::BGV>().bits(100).c(c).build();

  // Because bits is small, c gets clipped automatically.
  EXPECT_LT(context_built.getDigits().size(), c);
}

TEST(TestContextBGV, contextBuilderDoesNotClipDigitsSize)
{
  long c = 5; // Columns of SKMs
  helib::Context context_built =
      helib::ContextBuilder<helib::BGV>().bits(500).c(c).build();

  // Should have sufficient number of bits to have c columns.
  EXPECT_EQ(context_built.getDigits().size(), c);
}

TEST_P(TestContextBGV, contextBuilderWithBasicParams)
{
  // clang-format off
  helib::Context context_built = helib::ContextBuilder<helib::BGV>()
                               .m(m)
                               .p(p)
                               .r(r)
                               .build();   
  // clang-format off

  context->buildModChain( /*bits*/300, /*c=*/3);
  EXPECT_EQ(context_built, *context);
  EXPECT_GT(context_built.numPrimes(), 0);
}

TEST_P(TestContextBGV, contextBuilderWithGensOrdsToo)
{ 
  // clang-format off
  helib::Context context_built = helib::ContextBuilder<helib::BGV>()
                                 .m(m)
                                 .p(p)
                                 .r(r)
                                 .gens({3})
                                 .ords({-2})
                                 .build();   
  // clang-format off

  context->buildModChain( /*bits*/300, /*c=*/3);
  EXPECT_EQ(context_built, *context);
  EXPECT_GT(context_built.numPrimes(), 0);
}

TEST_P(TestContextBGV, contextBuilderNoModChain)
{
  // clang-format off
  helib::Context context_built = helib::ContextBuilder<helib::BGV>()
                                 .m(m)
                                 .p(p)
                                 .r(r)
                                 .buildModChain(false)
                                 .build();   
  // clang-format off
  EXPECT_EQ(context_built.numPrimes(), 0);
}

TEST_P(TestContextBGV, contextBuilderBootstrappableContext)
{
  NTL::Vec<long> mvec;
  mvec.SetLength(1);
  mvec[0] = 3;
  // clang-format off
  helib::Context context_built = helib::ContextBuilder<helib::BGV>()
                                 .mvec(mvec)
                                 .bootstrappable(true)
                                 .build();
  // clang-format on

  EXPECT_TRUE(context_built.isBootstrappable());
}

TEST(TestContextCKKS, contextBuilderWithDefaultArguments)
{
  helib::Context context_built = helib::ContextBuilder<helib::CKKS>().build();

  helib::Context expected_default_context = helib::ContextBuilder<helib::CKKS>()
                                                .m(4)
                                                .bits(300)
                                                .c(3)
                                                .precision(20)
                                                .build();

  // Making sure columns is not clipped.
  EXPECT_GT(expected_default_context.getCtxtPrimes().card(),
            expected_default_context.getDigits().size());
  EXPECT_GT(context_built.getCtxtPrimes().card(),
            context_built.getDigits().size());
  EXPECT_EQ(context_built.getDigits().size(), 3);

  EXPECT_FALSE(context_built.isBootstrappable());
  EXPECT_GT(context_built.numPrimes(), 0);
  EXPECT_EQ(context_built, expected_default_context);
}

TEST(TestContextCKKS, contextBuilderClipsDigitsSizeWithSmallBits)
{
  long c = 4; // Columns of SKMs
  helib::Context context_built =
      helib::ContextBuilder<helib::CKKS>().bits(100).c(c).build();

  // Because bits is small, c gets clipped automatically.
  EXPECT_LT(context_built.getDigits().size(), c);
}

TEST(TestContextCKKS, contextBuilderDoesNotClipDigitsSize)
{
  long c = 5; // Columns of SKMs
  helib::Context context_built =
      helib::ContextBuilder<helib::CKKS>().bits(500).c(c).build();

  // Should have sufficient number of bits to have c columns.
  EXPECT_EQ(context_built.getDigits().size(), c);
}

TEST_P(TestContextCKKS, contextBuilderWithBasicParams)
{
  // clang-format off
  helib::Context context_built = helib::ContextBuilder<helib::CKKS>()
                               .m(m)
                               .precision(r)
                               .build();   
  // clang-format off

  context->buildModChain( /*bits*/300, /*c=*/3);
  EXPECT_EQ(context_built, *context);
  EXPECT_GT(context_built.numPrimes(), 0);
}

TEST_P(TestContextCKKS, contextBuilderWithGensOrdsToo)
{ 
  // clang-format off
  helib::Context context_built = helib::ContextBuilder<helib::CKKS>()
                               .m(m)
                               .precision(r)
                               .gens({3})
                               .ords({64})
                               .build();   
  // clang-format off

  context->buildModChain( /*bits*/300, /*c=*/3);
  EXPECT_EQ(context_built, *context);
  EXPECT_GT(context_built.numPrimes(), 0);
}

TEST_P(TestContextCKKS, contextBuilderNoModChain)
{
  // clang-format off
  helib::Context context_built = helib::ContextBuilder<helib::CKKS>()
                               .m(m)
                               .precision(r)
                               .buildModChain(false)
                               .build();   
  // clang-format off
  EXPECT_EQ(context_built.numPrimes(), 0);
}

// Just checking manually if the printout works as expected
// TEST_P(TestContextBGV, contextBuilderPrintoutWorksCorrectly)
// {
  // auto bgv_context_builder = helib::ContextBuilder<helib::BGV>();
  // auto ckks_context_builder = helib::ContextBuilder<helib::CKKS>();
  // std::cout << bgv_context_builder << std::endl;
  // std::cout << ckks_context_builder << std::endl;
// }

INSTANTIATE_TEST_SUITE_P(variousParameters,
                         TestContextBGV,
                         ::testing::Values(BGVParameters(17, 2, 1)));

INSTANTIATE_TEST_SUITE_P(variousParameters,
                         TestContextCKKS,
                         ::testing::Values(CKKSParameters(256, 20)));

} // namespace
