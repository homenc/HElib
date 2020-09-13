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

// This test file does not fully cover Ctxt, and is intended as a starting
// point for testing new functionality of Ctxt going forward.
// The older tests with more extensive coverage can be found in the files
// with names matching "GTest*".

#include <helib/helib.h>
#include <helib/debugging.h>

#include "test_common.h"
#include "gtest/gtest.h"

namespace {

struct BGVParameters
{
  BGVParameters(unsigned m,
                unsigned p,
                unsigned r,
                unsigned bits,
                const std::vector<long>& gens = {},
                const std::vector<long>& ords = {}) :
      m(m), p(p), r(r), bits(bits), gens(gens), ords(ords){};

  const unsigned m;
  const unsigned p;
  const unsigned r;
  const unsigned bits;
  const std::vector<long> gens;
  const std::vector<long> ords;

  friend std::ostream& operator<<(std::ostream& os, const BGVParameters& params)
  {
    return os << "{"
              << "m = " << params.m << ", "
              << "p = " << params.p << ", "
              << "r = " << params.r << ", "
              << "gens = " << helib::vecToStr(params.gens) << ", "
              << "ords = " << helib::vecToStr(params.ords) << ", "
              << "bits = " << params.bits << "}";
  }
};

class TestCtxt : public ::testing::TestWithParam<BGVParameters>
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

  TestCtxt() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      bits(GetParam().bits),
      context(m, p, r),
      secretKey((buildModChain(context, bits), context)),
      publicKey((secretKey.GenSecKey(),
                 addFrbMatrices(secretKey),
                 addSome1DMatrices(secretKey),
                 secretKey)),
      ea(*(context.ea))
  {}

  virtual void SetUp() override
  {
    helib::setupDebugGlobals(&secretKey, context.ea);
  };

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }

  virtual ~TestCtxt() = default;
};

class TestCtxtWithBadDimensions : public TestCtxt
{
protected:
  TestCtxtWithBadDimensions() : TestCtxt()
  {
    for (long i = 0; i < context.zMStar.numOfGens(); ++i) {
      if (!ea.nativeDimension(i)) {
        return;
      }
    }
    throw std::logic_error("Algebra provided does not have a bad dimension");
  }
};

TEST_P(TestCtxt, timesEqualsWithLongWorks)
{
  helib::Ptxt<helib::BGV> ptxt(context, std::vector<long>(ea.size(), 5));
  helib::Ctxt ctxt(publicKey);
  publicKey.Encrypt(ctxt, ptxt);
  ctxt *= 2l;

  helib::Ptxt<helib::BGV> expected_result(context,
                                          std::vector<long>(ea.size(), 10));
  helib::Ptxt<helib::BGV> decrypted_result(context);
  secretKey.Decrypt(decrypted_result, ctxt);

  EXPECT_EQ(decrypted_result, expected_result);
}

TEST_P(TestCtxt, timesEqualsWithZZXWorks)
{
  helib::Ptxt<helib::BGV> ptxt(context, std::vector<long>(ea.size(), 5));
  helib::Ctxt ctxt(publicKey);
  publicKey.Encrypt(ctxt, ptxt);
  ctxt *= NTL::ZZX(2l);

  helib::Ptxt<helib::BGV> expected_result(context,
                                          std::vector<long>(ea.size(), 10));
  helib::Ptxt<helib::BGV> decrypted_result(context);
  secretKey.Decrypt(decrypted_result, ctxt);

  EXPECT_EQ(decrypted_result, expected_result);
}

TEST_P(TestCtxt, mapTo01WorksCorrectlyForConstantInputs)
{
  std::vector<long> data(ea.size());
  std::iota(data.begin(), data.end(), 0);
  for (auto& num : data)
    num %= p;
  helib::Ptxt<helib::BGV> ptxt(context, data);
  helib::Ctxt ctxt(publicKey);
  publicKey.Encrypt(ctxt, ptxt);
  mapTo01(ea, ctxt);
  mapTo01(ea, ptxt);

  std::vector<long> expected(data);
  for (auto& num : expected)
    num = num ? 1 : 0;

  helib::Ptxt<helib::BGV> expected_result(context, expected);
  helib::Ptxt<helib::BGV> result(context);
  secretKey.Decrypt(result, ctxt);

  EXPECT_EQ(expected_result, ptxt);
  EXPECT_EQ(expected_result, result);
  EXPECT_EQ(ptxt, result);
}

TEST_P(TestCtxt, mapTo01WorksCorrectlyForNonConstantInputs)
{
  NTL::ZZX poly;
  NTL::SetCoeff(poly, 0, 1001);
  NTL::SetCoeff(poly, 1, 1001);

  helib::Ptxt<helib::BGV> ptxt(context, std::vector<NTL::ZZX>(ea.size(), poly));
  helib::Ctxt ctxt(publicKey);
  publicKey.Encrypt(ctxt, ptxt);
  mapTo01(ea, ctxt);
  mapTo01(ea, ptxt);

  helib::Ptxt<helib::BGV> expected_result(context,
                                          std::vector<long>(ea.size(), 1));
  helib::Ptxt<helib::BGV> result(context);
  secretKey.Decrypt(result, ctxt);

  EXPECT_EQ(expected_result, ptxt);
  EXPECT_EQ(expected_result, result);
  EXPECT_EQ(ptxt, result);
}

TEST_P(TestCtxtWithBadDimensions,
       frobeniusAutomorphWorksCorrectlyWithBadDimensions)
{
  std::vector<long> data(ea.size());
  std::iota(data.begin(), data.end(), 0);
  helib::Ptxt<helib::BGV> ptxt(context, data);
  helib::Ctxt ctxt(publicKey);
  publicKey.Encrypt(ctxt, ptxt);
  helib::Ptxt<helib::BGV> expected_result(ptxt);
  for (long i = 0; i <= ea.getDegree(); ++i) {
    ctxt.frobeniusAutomorph(i);
    for (long j = 0; j < i; ++j) {
      expected_result.power(p);
    }

    helib::Ptxt<helib::BGV> result(context);
    secretKey.Decrypt(result, ctxt);

    EXPECT_EQ(expected_result, result)
        << "Frobenius automorph failed with i=" << i << std::endl;
  }
}

TEST_P(TestCtxt, frobeniusAutomorphWorksCorrectly)
{
  std::vector<long> data(ea.size());
  std::iota(data.begin(), data.end(), 0);
  helib::Ptxt<helib::BGV> ptxt(context, data);
  helib::Ctxt ctxt(publicKey);
  publicKey.Encrypt(ctxt, ptxt);
  helib::Ptxt<helib::BGV> expected_result(ptxt);
  for (long i = 0; i <= ea.getDegree(); ++i) {
    ctxt.frobeniusAutomorph(i);
    for (long j = 0; j < i; ++j) {
      expected_result.power(p);
    }

    helib::Ptxt<helib::BGV> result(context);
    secretKey.Decrypt(result, ctxt);

    EXPECT_EQ(ptxt, result)
        << "Frobenius automorph failed with i=" << i << std::endl;
  }
}

TEST_P(TestCtxtWithBadDimensions, rotate1DRotatesCorrectlyWithBadDimensions)
{
  std::vector<long> data(ea.size());
  std::iota(data.begin(), data.end(), 0);
  helib::Ptxt<helib::BGV> ptxt(context, data);
  helib::Ctxt ctxt(publicKey);
  publicKey.Encrypt(ctxt, ptxt);

  for (long i = 0; i < context.zMStar.numOfGens(); ++i) {
    helib::Ctxt tmp(ctxt);
    ea.rotate1D(tmp, i, 3);
    helib::Ptxt<helib::BGV> expected_result(ptxt);
    expected_result.rotate1D(i, 3);
    helib::Ptxt<helib::BGV> result(context);
    secretKey.Decrypt(result, tmp);

    EXPECT_EQ(expected_result, result);
  }
}

// Use this when thoroughly exploring an (m, p) grid of parameters.
// std::vector<BGVParameters> getParameters(bool good)
// {
//   std::vector<BGVParameters> parameterSets;
//
//   const long r = 1;
//   const long bits = 500;
//   const long min_p = 257;
//   const long max_p = 2003;
//   const long min_m = 100;
//   const long max_m = 3000;
//
//   std::vector<BGVParameters> params;
//   auto getParamFunc = good ? helib_test::getGoodDimensionParams
//                            : helib_test::getBadDimensionParams;
//   auto m_p_pairs = getParamFunc(min_m, max_m, min_p, max_p, 10, 10);
//   std::transform(m_p_pairs.begin(),
//                  m_p_pairs.end(),
//                  std::back_inserter(params),
//                  [](const auto& pair) {
//                    return BGVParameters(pair.first, pair.second, r, bits);
//                  });
//   return params;
// }

// INSTANTIATE_TEST_SUITE_P(variousParameters, TestCtxt,
// ::testing::ValuesIn(getParameters(true)));
// INSTANTIATE_TEST_SUITE_P(variousParameters, TestCtxtWithBadDimensions,
// ::testing::ValuesIn(getParameters(false)));

INSTANTIATE_TEST_SUITE_P(
    variousParameters,
    TestCtxt,
    ::testing::Values(BGVParameters(2049, 2, 1, 300),
                      BGVParameters(45, 1009, 1, 500),
                      BGVParameters(45, 317, 1, 500),
                      BGVParameters(45, 353, 1, 500),
                      BGVParameters(45, 367, 1, 500),
                      BGVParameters(45, 397, 1, 500),
                      BGVParameters(45, 419, 1, 500),
                      BGVParameters(45, 443, 1, 500),
                      BGVParameters(45, 971, 1, 500, {11}, {12}),
                      // This is here to test an algebra with 3 good dimensions
                      // BGVParameters(10005, 37, 1, 500),
                      BGVParameters(45, 19, 1, 300)));

INSTANTIATE_TEST_SUITE_P(variousParameters,
                         TestCtxtWithBadDimensions,
                         ::testing::Values(BGVParameters(45, 1009, 1, 500),
                                           BGVParameters(45, 349, 1, 300),
                                           BGVParameters(45, 379, 1, 300),
                                           BGVParameters(45, 499, 1, 300),
                                           BGVParameters(45, 619, 1, 300),
                                           BGVParameters(45, 709, 1, 300),
                                           BGVParameters(45, 769, 1, 300),
                                           BGVParameters(45, 829, 1, 300),
                                           BGVParameters(45, 919, 1, 300),
                                           BGVParameters(45, 19, 1, 300)));

} // namespace
