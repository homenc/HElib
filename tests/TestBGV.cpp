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

#include <NTL/ZZ.h>
#include <algorithm>

#include <helib/helib.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {
struct Parameters
{
  Parameters(unsigned m,
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

  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
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

class TestBGV : public ::testing::TestWithParam<Parameters>
{
protected:
  const unsigned long m;
  const unsigned long p;
  const unsigned long r;
  const unsigned long bits;
  helib::Context context;
  helib::SecKey secretKey;
  const helib::PubKey publicKey;
  const helib::EncryptedArray& ea;

  TestBGV() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      bits(GetParam().bits),
      context(helib::ContextBuilder<helib::BGV>()
                  .m(m)
                  .p(p)
                  .r(r)
                  .bits(bits)
                  .build()),
      secretKey(context),
      publicKey((secretKey.GenSecKey(),
                 helib::addSome1DMatrices(secretKey),
                 secretKey)),
      ea(context.getEA())
  {}

  virtual void SetUp() override
  {
    if (helib_test::verbose) {
      ea.getPAlgebra().printout();
      std::cout << "r = " << context.getAlMod().getR() << std::endl;
      std::cout << "ctxtPrimes=" << context.getCtxtPrimes()
                << ", specialPrimes=" << context.getSpecialPrimes() << "\n"
                << std::endl;
    }

    helib::setupDebugGlobals(&secretKey, context.shareEA());
  }

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }
};

TEST_P(TestBGV, negatingCiphertextWorks)
{
  helib::PtxtArray p1(ea), p2(ea);
  p1.random();

  helib::Ctxt c1(publicKey);
  p1.encrypt(c1);

  c1.negate();
  p1.negate();

  p2.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p2);
}

TEST_P(TestBGV, addingPolyConstantToCiphertextWorks)
{
  helib::PtxtArray p1(ea), p2(ea), const1(ea);
  p1.random();
  const1.random();

  helib::Ctxt c1(publicKey);
  p1.encrypt(c1);

  c1.addConstant(const1);
  p1 += const1;

  p2.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p2);
}

TEST_P(TestBGV, addingNegatedPolyConstantToCiphertextWorks)
{
  helib::PtxtArray p1(ea), p2(ea), const1(ea);
  p1.random();
  const1.random();

  helib::Ctxt c1(publicKey);
  p1.encrypt(c1);

  p1 -= const1;
  const1.negate();
  c1.addConstant(const1);

  p2.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p2);
}

TEST_P(TestBGV, multiplyingPolyConstantToCiphertextWorks)
{
  helib::PtxtArray p1(ea), p2(ea), const1(ea);
  p1.random();
  const1.random();

  helib::Ctxt c1(publicKey);
  p1.encrypt(c1);

  c1.multByConstant(const1);
  p1 *= const1;

  p2.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p2);
}

TEST_P(TestBGV, addingLongToCiphertextWorks)
{
  long const1 = 1;
  helib::PtxtArray p1(ea), p2(ea);
  p1.random();

  helib::Ctxt c1(publicKey);
  p1.encrypt(c1);

  c1.addConstant(const1);
  p1 += const1;

  p2.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p2);
}

TEST_P(TestBGV, multiplyingLongToCiphertextWorks)
{
  long const1 = 2;
  helib::PtxtArray p1(ea), p2(ea);
  p1.random();

  helib::Ctxt c1(publicKey);
  p1.encrypt(c1);

  c1.multByConstant(const1);
  p1 *= const1;

  p2.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p2);
}

TEST_P(TestBGV, rotatingCiphertextWorks)
{
  helib::PtxtArray p1(ea), p2(ea);
  p1.random();

  helib::Ctxt c1(publicKey);
  p1.encrypt(c1);

  ea.rotate(c1, 3);
  rotate(p1, 3);

  p2.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p2);
}

TEST_P(TestBGV, addingCiphertextsWorks)
{
  helib::PtxtArray p1(ea), p2(ea), p3(ea);
  p1.random();
  p2.random();

  helib::Ctxt c1(publicKey), c2(publicKey);
  p1.encrypt(c1);
  p2.encrypt(c2);

  c1 += c2;
  p1 += p2;

  p3.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p3);
}

TEST_P(TestBGV, subtractingCiphertextsWorks)
{
  helib::PtxtArray p1(ea), p2(ea), p3(ea);
  p1.random();
  p2.random();

  helib::Ctxt c1(publicKey), c2(publicKey);
  p1.encrypt(c1);
  p2.encrypt(c2);

  c1 -= c2;
  p1 -= p2;

  p3.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p3);
}

TEST_P(TestBGV, timesEqualsOfCiphertextsWorks)
{
  helib::PtxtArray p1(ea), p2(ea), p3(ea);
  p1.random();
  p2.random();

  helib::Ctxt c1(publicKey), c2(publicKey);
  p1.encrypt(c1);
  p2.encrypt(c2);

  c1 *= c2;
  p1 *= p2;

  p3.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p3);
  // Check that relinearization has occurred
  EXPECT_TRUE(c1.inCanonicalForm());
}

TEST_P(TestBGV, rawMultiplicationOfCiphertextsWorks)
{
  helib::PtxtArray p1(ea), p2(ea), p3(ea);
  p1.random();
  p2.random();

  helib::Ctxt c1(publicKey), c2(publicKey);
  p1.encrypt(c1);
  p2.encrypt(c2);

  c1.multLowLvl(c2);
  p1 *= p2;

  p3.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p3);
  // Check that relinearization has not occurred
  EXPECT_FALSE(c1.inCanonicalForm());
}

TEST_P(TestBGV, highLevelMultiplicationOfCiphertextsWorks)
{
  helib::PtxtArray p1(ea), p2(ea), p3(ea);
  p1.random();
  p2.random();

  helib::Ctxt c1(publicKey), c2(publicKey);
  p1.encrypt(c1);
  p2.encrypt(c2);

  c1.multiplyBy(c2);
  p1 *= p2;

  p3.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p3);
  // Check that relinearization has not occurred
  EXPECT_TRUE(c1.inCanonicalForm());
}

TEST_P(TestBGV, squaringCiphertextWorks)
{
  helib::PtxtArray p1(ea), p2(ea);
  p1.random();

  helib::Ctxt c1(publicKey);
  p1.encrypt(c1);

  c1.square();
  p1 *= p1;

  p2.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p2);
  // Check that relinearization has not occurred
  EXPECT_TRUE(c1.inCanonicalForm());
}

TEST_P(
    TestBGV,
    multiplyingCiphertextByNegativeConstantAndThenAddingToOtherCiphertextWorks)
{
  helib::PtxtArray p1(ea), p2(ea), p3(ea), const1(ea, -1);
  p1.random();
  p2.random();

  helib::Ctxt c1(publicKey), c2(publicKey);
  p1.encrypt(c1);
  p2.encrypt(c2);

  c1.multByConstant(const1);
  c1 += c2;
  p1 *= const1;
  p1 += p2;

  p3.decrypt(c1, secretKey);

  EXPECT_EQ(p1, p3);
}

INSTANTIATE_TEST_SUITE_P(typicalParameters,
                         TestBGV,
                         ::testing::Values(
                             // FAST
                             Parameters(257, 2, 1, 150)
                             // SLOW
                             // Parameters(2049, 2, 1, 150)
                             ));

} // namespace
