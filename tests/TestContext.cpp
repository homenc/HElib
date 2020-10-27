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
      context(std::make_shared<helib::Context>(m, p, r))
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
      context(std::make_shared<helib::Context>(m, /*p=*/-1, r))
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
  helib::Context someOtherContext(/*m=*/17, /*p=*/2, /*r=*/1);
  buildModChain(*context, /*bits=*/100, /*c=*/2);
  buildModChain(someOtherContext, /*bits=*/100, /*c=*/2);

  EXPECT_EQ(*context, *context);
  EXPECT_EQ(*context, someOtherContext);

  EXPECT_EQ(context->zMStar.getM(), 17); // m
  EXPECT_EQ(context->zMStar.getP(), 2);  // p
  EXPECT_EQ(context->alMod.getR(), 1);   // r
  EXPECT_EQ(context->digits.size(), 2);  // c
  EXPECT_GT(context->numPrimes(), 0);
  EXPECT_EQ(context->numPrimes(), someOtherContext.numPrimes());
  for (long i = 0; i < context->numPrimes(); i++) {
    const helib::Cmodulus& m1 = context->ithModulus(i);
    const helib::Cmodulus& m2 = someOtherContext.ithModulus(i);
    EXPECT_EQ(m1.getQ(), m2.getQ()) << " index: " << i;
  }
  EXPECT_EQ(context->smallPrimes, someOtherContext.smallPrimes);
  EXPECT_EQ(context->ctxtPrimes, someOtherContext.ctxtPrimes);
  EXPECT_EQ(context->specialPrimes, someOtherContext.specialPrimes);
  for (std::size_t i = 0; i < context->digits.size(); ++i) {
    EXPECT_EQ(context->digits[i], someOtherContext.digits[i])
        << " index: " << i;
  }
  EXPECT_EQ(context->stdev, someOtherContext.stdev);
  EXPECT_EQ(context->scale, someOtherContext.scale);
  EXPECT_EQ(context->rcData, someOtherContext.rcData);
}

TEST_P(TestContextBGV, contextNotEquals)
{
  helib::Context someOtherContext(/*m=*/13, /*p=*/3, /*r=*/2);
  buildModChain(*context, /*bits=*/100, /*c=*/3);
  buildModChain(someOtherContext, /*bits=*/200, /*c=*/2, 
                /*willBeBootstrappable =*/true);
  someOtherContext.scale = 6;
  someOtherContext.stdev = 3.0;
  NTL::Vec<long> mvec;
  mvec.SetLength(1);
  mvec[0] = 13;
  someOtherContext.enableBootStrapping(mvec);

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
  EXPECT_EQ(context->smallPrimes, someOtherContext.smallPrimes);
  EXPECT_NE(context->ctxtPrimes, someOtherContext.ctxtPrimes);
  EXPECT_NE(context->specialPrimes, someOtherContext.specialPrimes);
  for (std::size_t i = 0; i < context->digits.size(); ++i) {
    EXPECT_NE(context->digits[i], someOtherContext.digits[i])
        << " index: " << i;
  }
  EXPECT_NE(context->stdev, someOtherContext.stdev);
  EXPECT_NE(context->scale, someOtherContext.scale);
  EXPECT_NE(context->rcData, someOtherContext.rcData);
}

TEST_P(TestContextBGV, ContextCalculatingSecurityAfterModchainBuilt)
{
  buildModChain(*context, /*bits=*/100, /*c=*/2);
  double result = context->securityLevel();
  EXPECT_FALSE(std::isinf(result));
}

TEST_P(TestContextBGV, hasCorrectSlotRingWhenConstructed)
{
  EXPECT_EQ(context->slotRing->p, p);
  EXPECT_EQ(context->slotRing->r, r);
  EXPECT_EQ(context->slotRing->p2r, pow(p, r));
  EXPECT_EQ(context->slotRing->G, helib::getG(*(context->ea)));
}

TEST_P(TestContextBGV, buildModChainThrowsWhenBitsIsZero)
{
  EXPECT_THROW(helib::buildModChain(*context, /*bits=*/0, /*c=*/2),
               helib::InvalidArgument);
}

TEST_P(TestContextBGV, calculateBitSizeOfQ)
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
  helib::Context context_built { helib::ContextBuilder<helib::BGV>().build() };

  helib::Context expected_default_context(/*m=*/3, /*p=*/2, /*r=*/1);
  buildModChain(expected_default_context, /*bits=*/300, /*c=*/3);

  // Making sure the number of columns is not clipped.
  EXPECT_GT(expected_default_context.ctxtPrimes.card(),
            expected_default_context.digits.size());
  EXPECT_GT(context_built.ctxtPrimes.card(), context_built.digits.size());
  EXPECT_EQ(context_built.digits.size(), 3);

  EXPECT_FALSE(context_built.isBootstrappable());
  EXPECT_GT(context_built.numPrimes(), 0);
  EXPECT_EQ(context_built, expected_default_context);
}

TEST(TestContextBGV, contextBuilderClipsDigitsSizeWithSmallBits)
{
  long c = 4; // Columns of SKMs
  helib::Context context_built 
      { helib::ContextBuilder<helib::BGV>().bits(100).c(c).build() };

  // Because bits is small, c gets clipped automatically.
  EXPECT_LT(context_built.digits.size(), c);
}

TEST(TestContextBGV, contextBuilderDoesNotClipDigitsSize)
{
  long c = 5; // Columns of SKMs
  helib::Context context_built 
    { helib::ContextBuilder<helib::BGV>().bits(500).c(c).build() };

  // Should have sufficient number of bits to have c columns.
  EXPECT_EQ(context_built.digits.size(), c);
}

TEST_P(TestContextBGV, contextBuilderWithBasicParams)
{
  // clang-format off
  helib::Context context_built { helib::ContextBuilder<helib::BGV>()
                               .m(m)
                               .p(p)
                               .r(r)
                               .build() };   
  // clang-format off

  buildModChain(*context, /*bits*/300, /*c=*/3);
  EXPECT_EQ(context_built, *context);
  EXPECT_GT(context_built.numPrimes(), 0);
}

TEST_P(TestContextBGV, contextBuilderWithGensOrdsToo)
{ 
  // clang-format off
  helib::Context context_built { helib::ContextBuilder<helib::BGV>()
                               .m(m)
                               .p(p)
                               .r(r)
                               .gens({3})
                               .ords({-2})
                               .build() };   
  // clang-format off

  buildModChain(*context, /*bits*/300, /*c=*/3);
  EXPECT_EQ(context_built, *context);
  EXPECT_GT(context_built.numPrimes(), 0);
}

TEST_P(TestContextBGV, contextBuilderNoModChain)
{
  // clang-format off
  helib::Context context_built { helib::ContextBuilder<helib::BGV>()
                               .m(m)
                               .p(p)
                               .r(r)
                               .buildModChain(false)
                               .build() };   
  // clang-format off
  EXPECT_EQ(context_built.numPrimes(), 0);
}

TEST_P(TestContextBGV, contextBuilderBootstrappableContext)
{
  NTL::Vec<long> mvec;
  mvec.SetLength(1);
  mvec[0] = 3;
  // clang-format off
  helib::Context context_built { helib::ContextBuilder<helib::BGV>()
                               .mvec(mvec)
                               .bootstrappable(true)
                               .build() };
  // clang-format on

  EXPECT_TRUE(context_built.isBootstrappable());
}

TEST_P(TestContextBGV, contextBuilderLogsCorrectly)
{
  NTL::Vec<long> mvec;
  mvec.SetLength(1);
  mvec[0] = 3;

  long c = 3;
  std::vector<long> gens = {3};
  std::vector<long> ords = {-2};
  bool buildModChainFlag = false;
  long bits = 2;
  long skHwt = 64;
  long resolution = 1;
  long bitsInSpecialPrimes = 15;
  bool bootstrappableFlag = true;
  bool buildCacheFlag = true;
  bool thickFlag = true;

  // clang-format off
  auto cb = helib::ContextBuilder<helib::BGV>()
                          .m(m)
                          .p(p)
                          .r(r)
                          .c(c)
                          .gens(gens)
                          .ords(ords)
                          .buildModChain(buildModChainFlag)
                          .bits(bits)
                          .skHwt(skHwt)
                          .resolution(resolution)
                          .bitsInSpecialPrimes(bitsInSpecialPrimes)
                          .bootstrappable(bootstrappableFlag)
                          .mvec(mvec)
                          .buildCache(buildCacheFlag)
                          .thickboot();
  // clang-format off

  std::stringstream expected_ss;
  expected_ss << "{\n"
              << "  scheme: BGV\n"
              << "  m: " << m << "\n"
              << "  p: " << p << "\n"
              << "  r: " << r << "\n"
              << "  c: " << c << "\n"
              << "  gens: " << helib::vecToStr(gens) << "\n"
              << "  ords: " << helib::vecToStr(ords) << "\n"
              << "  buildModChainFlag: " << buildModChainFlag << "\n"
              << "  bits: " << bits << "\n"
              << "  skHwt: " << skHwt << "\n"
              << "  resolution: " << resolution << "\n"
              << "  bitsInSpecialPrimes: " << bitsInSpecialPrimes << "\n"
              << "  bootstrappableFlag: " << bootstrappableFlag << "\n"
              << "  mvec: " << mvec << "\n"
              << "  buildCacheFlag: " << buildCacheFlag << "\n"
              << "  thickFlag: " << thickFlag << "\n"
              << "}" << std::endl;
  
  std::stringstream actual_ss;
  actual_ss << cb;

  EXPECT_EQ(expected_ss.str(), actual_ss.str());
}

TEST(TestContextCKKS, contextBuilderWithDefaultArguments)
{
  helib::Context context_built { helib::ContextBuilder<helib::CKKS>().build() };

  helib::Context expected_default_context(/*m=*/4, /*p=*/-1, /*r=*/20);
  buildModChain(expected_default_context, /*bits=*/300, /*c=*/3);

  // Making sure columns is not clipped.
  EXPECT_GT(expected_default_context.ctxtPrimes.card(),
            expected_default_context.digits.size());
  EXPECT_GT(context_built.ctxtPrimes.card(), context_built.digits.size());
  EXPECT_EQ(context_built.digits.size(), 3);

  EXPECT_FALSE(context_built.isBootstrappable());
  EXPECT_GT(context_built.numPrimes(), 0);
  EXPECT_EQ(context_built, expected_default_context);
}

TEST(TestContextCKKS, contextBuilderClipsDigitsSizeWithSmallBits)
{
  long c = 4; // Columns of SKMs
  helib::Context context_built {
      helib::ContextBuilder<helib::CKKS>().bits(100).c(c).build() };

  // Because bits is small, c gets clipped automatically.
  EXPECT_LT(context_built.digits.size(), c);
}

TEST(TestContextCKKS, contextBuilderDoesNotClipDigitsSize)
{
  long c = 5; // Columns of SKMs
  helib::Context context_built {
      helib::ContextBuilder<helib::CKKS>().bits(500).c(c).build() };

  // Should have sufficient number of bits to have c columns.
  EXPECT_EQ(context_built.digits.size(), c);
}

TEST_P(TestContextCKKS, contextBuilderWithBasicParams)
{
  // clang-format off
  helib::Context context_built { helib::ContextBuilder<helib::CKKS>()
                               .m(m)
                               .precision(r)
                               .build() };   
  // clang-format off

  buildModChain(*context, /*bits*/300, /*c=*/3);
  EXPECT_EQ(context_built, *context);
  EXPECT_GT(context_built.numPrimes(), 0);
}

TEST_P(TestContextCKKS, contextBuilderWithGensOrdsToo)
{ 
  // clang-format off
  helib::Context context_built { helib::ContextBuilder<helib::CKKS>()
                               .m(m)
                               .precision(r)
                               .gens({3})
                               .ords({64})
                               .build() };   
  // clang-format off

  buildModChain(*context, /*bits*/300, /*c=*/3);
  EXPECT_EQ(context_built, *context);
  EXPECT_GT(context_built.numPrimes(), 0);
}

TEST_P(TestContextCKKS, contextBuilderNoModChain)
{
  // clang-format off
  helib::Context context_built { helib::ContextBuilder<helib::CKKS>()
                               .m(m)
                               .precision(r)
                               .buildModChain(false)
                               .build() };   
  // clang-format off
  EXPECT_EQ(context_built.numPrimes(), 0);
}

TEST_P(TestContextCKKS, contextBuilderLogsCorrectly)
{
  NTL::Vec<long> mvec;
  mvec.SetLength(1);
  mvec[0] = 3;

  long precision = r;
  long c = 3;
  std::vector<long> gens = {3};
  std::vector<long> ords = {-2};
  bool buildModChainFlag = false;
  long bits = 2;
  long skHwt = 64;
  long resolution = 1;
  long bitsInSpecialPrimes = 15;

  // clang-format off
  auto cb = helib::ContextBuilder<helib::CKKS>()
                          .m(m)
                          .precision(precision)
                          .c(c)
                          .gens(gens)
                          .ords(ords)
                          .buildModChain(buildModChainFlag)
                          .bits(bits)
                          .skHwt(skHwt)
                          .resolution(resolution)
                          .bitsInSpecialPrimes(bitsInSpecialPrimes);
  // clang-format off

  std::stringstream expected_ss;
  expected_ss << "{\n"
              << "  scheme: CKKS\n"
              << "  m: " << m << ",\n"
              << "  precision: " << precision << ",\n"
              << "  c: " << c << ",\n"
              << "  gens: " << helib::vecToStr(gens) << ",\n"
              << "  ords: " << helib::vecToStr(ords) << ",\n"
              << "  buildModChainFlag: " << buildModChainFlag << ",\n"
              << "  bits: " << bits << ",\n"
              << "  skHwt: " << skHwt << ",\n"
              << "  resolution: " << resolution << ",\n"
              << "  bitsInSpecialPrimes: " << bitsInSpecialPrimes << "\n"
              << "}" << std::endl;
  
  std::stringstream actual_ss;
  actual_ss << cb;

  EXPECT_EQ(expected_ss.str(), actual_ss.str());
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
