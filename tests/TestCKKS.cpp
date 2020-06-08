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
#include <NTL/ZZ.h>
#include <algorithm>
#include <complex>

#include <helib/norms.h>
#include <helib/helib.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {
struct Parameters
{
  Parameters(long m, long r, long L, double epsilon) :
      m(m), r(r), L(L), epsilon(epsilon){};

  const long m;
  const long r;
  const long L;
  const double epsilon;

  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "m=" << params.m << ","
              << "r=" << params.r << ","
              << "L=" << params.L << ","
              << "epsilon=" << params.epsilon << "}";
  }
};

// Utility functions for the tests

// Compute the L-infinity distance between two vectors
double calcMaxDiff(const std::vector<std::complex<double>>& v1,
                   const std::vector<std::complex<double>>& v2)
{
  if (helib::lsize(v1) != helib::lsize(v2)) {
    throw std::runtime_error("Vector sizes differ.");
  }

  double maxDiff = 0.0;
  for (long i = 0; i < helib::lsize(v1); i++) {
    double diffAbs = std::abs(v1[i] - v2[i]);
    if (diffAbs > maxDiff)
      maxDiff = diffAbs;
  }
  return maxDiff;
}
// Compute the max relative difference between two vectors
double calcMaxRelDiff(const std::vector<std::complex<double>>& v1,
                      const std::vector<std::complex<double>>& v2)
{
  if (helib::lsize(v1) != helib::lsize(v2)) {
    throw std::runtime_error("Vector sizes differ.");
  }

  // Compute the largest-magnitude value in the vector
  double maxAbs = 0.0;
  for (auto& x : v1) {
    if (std::abs(x) > maxAbs)
      maxAbs = std::abs(x);
  }
  if (maxAbs < 1e-10)
    maxAbs = 1e-10;

  double maxDiff = 0.0;
  for (long i = 0; i < helib::lsize(v1); i++) {
    double relDiff = std::abs(v1[i] - v2[i]) / maxAbs;
    if (relDiff > maxDiff)
      maxDiff = relDiff;
  }

  return maxDiff;
}

inline bool cx_equals(const std::vector<std::complex<double>>& v1,
                      const std::vector<std::complex<double>>& v2,
                      double epsilon)
{
  return (calcMaxDiff(v1, v2) < epsilon);
}

void negateVec(std::vector<std::complex<double>>& p1)
{
  for (auto& x : p1)
    x = -x;
}
void conjVec(std::vector<std::complex<double>>& p1)
{
  for (auto& x : p1)
    x = conj(x);
}
void add(std::vector<std::complex<double>>& to,
         const std::vector<std::complex<double>>& from)
{
  if (to.size() < from.size())
    to.resize(from.size(), 0);
  for (std::size_t i = 0; i < from.size(); i++)
    to[i] += from[i];
}
void add(std::vector<std::complex<double>>& to, double from)
{
  for (std::size_t i = 0; i < to.size(); i++)
    to[i] += from;
}
void sub(std::vector<std::complex<double>>& to,
         const std::vector<std::complex<double>>& from)
{
  if (to.size() < from.size())
    to.resize(from.size(), 0);
  for (std::size_t i = 0; i < from.size(); i++)
    to[i] -= from[i];
}
void mul(std::vector<std::complex<double>>& to,
         const std::vector<std::complex<double>>& from)
{
  if (to.size() < from.size())
    to.resize(from.size(), 0);
  for (std::size_t i = 0; i < from.size(); i++)
    to[i] *= from[i];
}
void mul(std::vector<std::complex<double>>& to, double from)
{
  for (std::size_t i = 0; i < to.size(); i++)
    to[i] *= from;
}
void rotate(std::vector<std::complex<double>>& p, long amt)
{
  long sz = p.size();
  std::vector<std::complex<double>> tmp(sz);
  for (long i = 0; i < sz; i++)
    tmp[((i + amt) % sz + sz) % sz] = p[i];
  p = tmp;
}

class TestCKKS : public ::testing::TestWithParam<Parameters>
{
protected:
  const long m;         // Zm*
  const long r;         // bit precision
  const long L;         // Number of bits
  const double epsilon; // error threshold

  helib::Context context;
  helib::SecKey secretKey;
  const helib::PubKey publicKey;
  const helib::EncryptedArrayCx& ea;

  TestCKKS() :
      m(GetParam().m),
      r(GetParam().r),
      L(GetParam().L),
      epsilon(GetParam().epsilon),
      context(m, /*p=*/-1, r),
      secretKey((context.scale = 4,
                 helib::buildModChain(context, L, /*c=*/2),
                 context)),
      publicKey((secretKey.GenSecKey(),
                 helib::addSome1DMatrices(secretKey),
                 secretKey)),
      ea(context.ea->getCx())
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

TEST_P(TestCKKS, negatingCiphertextWorks)
{
  helib::Ctxt c1(publicKey);
  std::vector<std::complex<double>> vd1, vd2;
  NTL::ZZX poly;
  NTL::xdouble rf, pm;

  ea.random(vd1);
  ea.encrypt(c1, publicKey, vd1);
  rf = c1.getRatFactor();
  pm = c1.getPtxtMag();
  ea.encode(poly, vd1, /*size=*/1.0);
  c1.negate();
  ea.decrypt(c1, secretKey, vd2);

  negateVec(vd1);

  EXPECT_TRUE(cx_equals(vd2, vd1, NTL::conv<double>(epsilon)))
      << "  max(vd1)=" << helib::largestCoeff(vd1)
      << ", max(vd2)=" << helib::largestCoeff(vd2)
      << ", maxDiff=" << calcMaxDiff(vd1, vd2) << std::endl
      << std::endl;
  EXPECT_EQ(rf, c1.getRatFactor());
  EXPECT_EQ(pm, c1.getPtxtMag());
}

TEST_P(TestCKKS, addingPolyConstantToCiphertextWorks)
{
  helib::Ctxt c1(publicKey);
  std::vector<std::complex<double>> vd1, vd2, vd3;
  NTL::ZZX poly;
  NTL::xdouble rf, pm;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  rf = c1.getRatFactor();
  pm = c1.getPtxtMag();
  ea.encode(poly, vd2, /*size=*/1.0);
  c1.addConstantCKKS(poly);
  ea.decrypt(c1, secretKey, vd3);

  add(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
      << "  max(vd1)=" << helib::largestCoeff(vd1)
      << ", max(vd3)=" << helib::largestCoeff(vd3)
      << ", maxDiff=" << calcMaxDiff(vd1, vd3) << std::endl
      << std::endl;
}

TEST_P(TestCKKS, addingNegatedPolyConstantToCiphertextWorks)
{
  helib::Ctxt c1(publicKey);
  std::vector<std::complex<double>> vd1, vd2, vd3;
  NTL::ZZX poly;
  NTL::xdouble rf, pm;

  ea.random(vd1);
  ea.random(vd2);
  negateVec(vd2);
  ea.encrypt(c1, publicKey, vd1);
  rf = c1.getRatFactor();
  pm = c1.getPtxtMag();
  ea.encode(poly, vd2, /*size=*/1.0);
  c1.addConstantCKKS(poly);
  ea.decrypt(c1, secretKey, vd3);

  add(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
      << "  max(vd1)=" << helib::largestCoeff(vd1)
      << ", max(vd3)=" << helib::largestCoeff(vd3) << std::endl
      << ", maxDiff=" << calcMaxDiff(vd1, vd3)
      << ", maxRelDiff=" << calcMaxRelDiff(vd1, vd3) << std::endl
      << std::endl;
}

TEST_P(TestCKKS, multiplyingPolyConstantToCiphertextWorks)
{
  helib::Ctxt c1(publicKey);
  std::vector<std::complex<double>> vd1, vd2, vd3;
  NTL::ZZX poly;
  NTL::xdouble rf, pm;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  rf = c1.getRatFactor();
  pm = c1.getPtxtMag();
  ea.encode(poly, vd2, /*size=*/1.0);
  c1.multByConstantCKKS(poly);
  ea.decrypt(c1, secretKey, vd3);

  mul(vd1, vd2);
  rf *= ea.encodeScalingFactor();

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
      << "  max(vd1)=" << helib::largestCoeff(vd1)
      << ", max(vd3)=" << helib::largestCoeff(vd3)
      << ", maxDiff=" << calcMaxDiff(vd1, vd3) << std::endl
      << std::endl;
  EXPECT_EQ(rf, c1.getRatFactor());
  EXPECT_EQ(pm, c1.getPtxtMag());
}

TEST_P(TestCKKS, addingDoubleToCiphertextWorks)
{
  helib::Ctxt c1(publicKey);
  std::vector<std::complex<double>> vd1, vd2;
  NTL::xdouble rf, pm;
  std::vector<double> vd(1);
  ea.random(vd);

  ea.random(vd1);
  ea.encrypt(c1, publicKey, vd1);
  rf = c1.getRatFactor();
  pm = c1.getPtxtMag();
  c1.addConstantCKKS(vd[0]);
  ea.decrypt(c1, secretKey, vd2);

  add(vd1, vd[0]);

  EXPECT_TRUE(cx_equals(vd2, vd1, NTL::conv<double>(epsilon)))
      << "  max(vd1)=" << helib::largestCoeff(vd1)
      << ", max(vd2)=" << helib::largestCoeff(vd2)
      << ", maxDiff=" << calcMaxDiff(vd1, vd2) << std::endl
      << std::endl;
}

TEST_P(TestCKKS, multiplyingDoubleToCiphertextWorks)
{
  helib::Ctxt c1(publicKey);
  std::vector<std::complex<double>> vd1, vd2, vd0;
  NTL::xdouble rf, pm;
  std::vector<double> vd(1);
  ea.random(vd);

  ea.random(vd1);
  vd0 = vd1;
  ea.encrypt(c1, publicKey, vd1);
  rf = c1.getRatFactor();
  pm = c1.getPtxtMag();
  c1.multByConstantCKKS(vd[0]);
  ea.decrypt(c1, secretKey, vd2);

  mul(vd1, vd[0]);
  rf /= vd[0];
  pm *= std::abs(vd[0]);

  EXPECT_TRUE(cx_equals(vd2, vd1, NTL::conv<double>(epsilon)))
      << "  max(vd1)=" << helib::largestCoeff(vd1)
      << ", max(vd2)=" << helib::largestCoeff(vd2)
      << ", maxDiff=" << calcMaxDiff(vd1, vd2) << std::endl
      << ", ptxtMag=" << c1.getPtxtMag() << std::endl
      << ", vd[0]=" << vd[0] << std::endl;
  EXPECT_EQ(rf, c1.getRatFactor());
  EXPECT_EQ(pm, c1.getPtxtMag());
}

TEST_P(TestCKKS, gettingTheComplexConjugateWorks)
{
  helib::Ctxt c1(publicKey);
  std::vector<std::complex<double>> vd1, vd2;
  NTL::xdouble rf, pm;
  NTL::ZZX poly;

  ea.random(vd1);
  ea.encrypt(c1, publicKey, vd1);
  rf = c1.getRatFactor();
  pm = c1.getPtxtMag();
  ea.encode(poly, vd1, /*size=*/1.0);
  c1.complexConj();
  ea.decrypt(c1, secretKey, vd2);

  conjVec(vd1);

  EXPECT_TRUE(cx_equals(vd2, vd1, NTL::conv<double>(epsilon)))
      << "  max(vd1)=" << helib::largestCoeff(vd1)
      << ", max(vd2)=" << helib::largestCoeff(vd2)
      << ", maxDiff=" << calcMaxDiff(vd1, vd2) << std::endl
      << std::endl;
  EXPECT_EQ(rf, c1.getRatFactor());
  EXPECT_EQ(pm, c1.getPtxtMag());
}

TEST_P(TestCKKS, rotatingCiphertextWorks)
{
  helib::Ctxt c1(publicKey);
  std::vector<std::complex<double>> vd1, vd2;
  NTL::xdouble rf, pm;
  NTL::ZZX poly;

  ea.random(vd1);
  ea.encrypt(c1, publicKey, vd1);
  pm = c1.getPtxtMag();
  ea.encode(poly, vd1, /*size=*/1.0);
  ea.rotate(c1, 3);
  ea.decrypt(c1, secretKey, vd2);
  // TODO: understand how this affects ratfactor and add appropriate expectation
  rotate(vd1, 3);
  // vd1 is now the expected result

  EXPECT_TRUE(cx_equals(vd2, vd1, NTL::conv<double>(epsilon)))
      << "  max(vd1)=" << helib::largestCoeff(vd1)
      << ", max(vd2)=" << helib::largestCoeff(vd2)
      << ", maxDiff=" << calcMaxDiff(vd1, vd2) << std::endl
      << std::endl;
  EXPECT_EQ(pm, c1.getPtxtMag());
}

TEST_P(TestCKKS, addingCiphertextsWorks)
{
  helib::Ctxt c1(publicKey), c2(publicKey);
  std::vector<std::complex<double>> vd1, vd2, vd3;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  ea.encrypt(c2, publicKey, vd2);
  c1 += c2;
  ea.decrypt(c1, secretKey, vd3);

  add(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
      << "  max(vd1)=" << helib::largestCoeff(vd1)
      << ", max(vd3)=" << helib::largestCoeff(vd3)
      << ", maxDiff=" << calcMaxDiff(vd1, vd3) << std::endl
      << std::endl;
}

TEST_P(TestCKKS, subtractingCiphertextsWorks)
{
  helib::Ctxt c1(publicKey), c2(publicKey);
  std::vector<std::complex<double>> vd1, vd2, vd3;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  ea.encrypt(c2, publicKey, vd2);
  c1 -= c2;
  ea.decrypt(c1, secretKey, vd3);

  sub(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
      << "  max(vd1)=" << helib::largestCoeff(vd1)
      << ", max(vd3)=" << helib::largestCoeff(vd3)
      << ", maxDiff=" << calcMaxDiff(vd1, vd3) << std::endl
      << std::endl;
}

TEST_P(TestCKKS, rawMultiplicationOfCiphertextsWorks)
{
  helib::Ctxt c1(publicKey), c2(publicKey);
  std::vector<std::complex<double>> vd1, vd2, vd3;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  ea.encrypt(c2, publicKey, vd2);
  NTL::xdouble expectedPtxtMag = c1.getPtxtMag() * c2.getPtxtMag();
  c1 *= c2;
  ea.decrypt(c1, secretKey, vd3);

  mul(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
      << "  max(vd1)=" << helib::largestCoeff(vd1)
      << ", max(vd3)=" << helib::largestCoeff(vd3)
      << ", maxDiff=" << calcMaxDiff(vd1, vd3) << std::endl
      << std::endl;
  EXPECT_EQ(expectedPtxtMag, c1.getPtxtMag());
}

TEST_P(TestCKKS, highLevelMultiplicationOfCiphertextsWorks)
{
  helib::Ctxt c1(publicKey), c2(publicKey);
  std::vector<std::complex<double>> vd1, vd2, vd3;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  ea.encrypt(c2, publicKey, vd2);
  NTL::xdouble expectedPtxtMag = c1.getPtxtMag() * c2.getPtxtMag();
  c1.multiplyBy(c2);
  ea.decrypt(c1, secretKey, vd3);

  mul(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
      << "  max(vd1)=" << helib::largestCoeff(vd1)
      << ", max(vd3)=" << helib::largestCoeff(vd3)
      << ", maxDiff=" << calcMaxDiff(vd1, vd3) << std::endl
      << std::endl;
  EXPECT_EQ(expectedPtxtMag, c1.getPtxtMag());
}

TEST_P(
    TestCKKS,
    multiplyingCiphertextByNegativeConstantAndThenAddingToOtherCiphertextWorks)
{
  helib::Ctxt c1(publicKey), c2(publicKey);
  std::vector<std::complex<double>> vd1, vd2, vd3, vd4;
  NTL::ZZX poly;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  ea.encrypt(c2, publicKey, vd2);
  c1.multByConstantCKKS(std::make_pair<long, long>(-5, 2));
  c1 += c2;
  ea.decrypt(c1, secretKey, vd4);

  mul(vd1, std::vector<std::complex<double>>(ea.size(), -5 / 2.0));
  add(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd4, vd1, epsilon))
      << "  max(vd1)=" << helib::largestCoeff(vd1)
      << ", max(vd4)=" << helib::largestCoeff(vd4) << std::endl
      << ", maxDiff=" << calcMaxDiff(vd1, vd4) << std::endl
      << std::endl;
}

TEST(TestCKKS, buildingCKKSContextWithMAsNotAPowerOfTwoThrows)
{
  EXPECT_THROW(helib::Context context(99, -1, 20), helib::InvalidArgument);
}

INSTANTIATE_TEST_SUITE_P(typicalParameters,
                         TestCKKS,
                         ::testing::Values(
                             // SLOW
                             Parameters(1024, 20, 150, 0.01)
                             // FAST
                             // Parameters(128, 20, 150, 0.01)
                             ));

} // namespace
