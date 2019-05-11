/* Copyright (C) 2012-2017 IBM Corp.
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

#include "norms.h"
#include "EncryptedArray.h"
#include "FHE.h"
#include "debugging.h"

#include "gtest/gtest.h"
#include "test_common.h"

namespace {
struct Parameters {
    Parameters(long m, long r, long L, double epsilon) :
        m(m),
        r(r),
        L(L),
        epsilon(epsilon)
    {};

    const long m;
    const long r;
    const long L;
    const double epsilon;

    friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
    {
        return os << "{" <<
            "m=" << params.m << "," <<
            "r=" << params.r << "," <<
            "L=" << params.L << "," <<
            "epsilon=" << params.epsilon <<
            "}";
    }
};

// Utility functions for the tests

  // Compute the L-infinity distance between two vectors
  double calcMaxDiff(const std::vector<cx_double>& v1,
                     const std::vector<cx_double>& v2){
    if(lsize(v1) != lsize(v2)) {
      throw std::runtime_error("Vector sizes differ.");
    }
    
    double maxDiff = 0.0;
    for (long i=0; i<lsize(v1); i++) {
      double diffAbs = std::abs(v1[i]-v2[i]);
      if (diffAbs > maxDiff)
        maxDiff = diffAbs;
    }
    return maxDiff;
  }
  // Compute the max relative difference between two vectors
  double calcMaxRelDiff(const std::vector<cx_double>& v1,
                        const std::vector<cx_double>& v2)
  {
    if(lsize(v1)!=lsize(v2)) {
      throw std::runtime_error("Vector sizes differ.");
      
    }
    
    // Compute the largest-magnitude value in the vector
    double maxAbs = 0.0;
    for (auto& x : v1) {
      if (std::abs(x) > maxAbs)
        maxAbs = std::abs(x);
    }
    if (maxAbs<1e-10)
      maxAbs = 1e-10;
    
    double maxDiff = 0.0;
    for (long i=0; i<lsize(v1); i++) {
      double relDiff = std::abs(v1[i]-v2[i]) / maxAbs;
      if (relDiff > maxDiff)
        maxDiff = relDiff;
    }
    
    return maxDiff;
  }
  
  inline bool cx_equals(const std::vector<cx_double>& v1,
                        const std::vector<cx_double>& v2,
                        double epsilon)
  {
    return (calcMaxDiff(v1,v2) < epsilon);
  }
  
  ::testing::AssertionResult ciphertextMatches(const EncryptedArrayCx &ea, const FHESecKey &sk,
                                               const std::vector<cx_double> &p, const Ctxt &c, double epsilon)
  {
    std::vector<cx_double> pp;
    ea.decrypt(c, sk, pp);
    if (helib_test::verbose) {
      std::cout << "    relative-error=" << calcMaxRelDiff(p, pp)
                << ", absolute-error=" << calcMaxRelDiff(p, pp) << std::endl;
    }
    
    if(cx_equals(pp, p, epsilon)) {
      return ::testing::AssertionSuccess();
    } else {
      return ::testing::AssertionFailure()
          << "Ciphertext does not match plaintext:" << std::endl
          << "p = " << p << std::endl
          << "pp = " << pp << std::endl;
    }
  }
  
  void negateVec(std::vector<cx_double>& p1)
  {
    for (auto& x: p1) x = -x;
  }
  void conjVec(std::vector<cx_double>& p1)
  {
    for (auto& x: p1) x = conj(x);
  }
  void add(std::vector<cx_double>& to, const std::vector<cx_double>& from)
  {
    if (to.size() < from.size())
      to.resize(from.size(), 0);
    for (long i=0; i<from.size(); i++) to[i] += from[i];
  }
  void add(std::vector<cx_double>& to, double from)
  {
    for (long i=0; i<to.size(); i++) to[i] += from;
  }
  void sub(std::vector<cx_double>& to, const std::vector<cx_double>& from)
  {
    if (to.size() < from.size())
      to.resize(from.size(), 0);
    for (long i=0; i<from.size(); i++) to[i] -= from[i];
  }
  void mul(std::vector<cx_double>& to, const std::vector<cx_double>& from)
  {
    if (to.size() < from.size())
      to.resize(from.size(), 0);
    for (long i=0; i<from.size(); i++) to[i] *= from[i];
  }
  void mul(std::vector<cx_double>& to, double from)
  {
    for (long i=0; i<to.size(); i++) to[i] *= from;
  }
  void rotate(std::vector<cx_double>& p, long amt)
  {
    long sz = p.size();
    std::vector<cx_double> tmp(sz);
    for (long i=0; i<sz; i++)
      tmp[((i+amt)%sz +sz)%sz] = p[i];
    p = tmp;
  }

#ifndef FFT_ARMA
class Test_CKKS : public ::testing::TestWithParam<Parameters> {};
TEST_P(Test_CKKS, error)
{
    FAIL() << "Cannot run CKKS test when not using armadillo.";
}
#else

class Test_CKKS : public ::testing::TestWithParam<Parameters> {
    protected:
        const long m;         // Zm*
        const long r;         // bit precision
        const long L;         // Number of bits
        const double epsilon; // error threshold

        FHEcontext context;
        FHESecKey secretKey;
        const FHEPubKey publicKey;
        const EncryptedArrayCx& ea;

        Test_CKKS () :
            m(GetParam().m),
            r(GetParam().r),
            L(GetParam().L),
            epsilon(GetParam().epsilon),
            context(m, /*p=*/-1, r),
            secretKey((context.scale=4, buildModChain(context, L, /*c=*/2), context)),
            publicKey((secretKey.GenSecKey(), addSome1DMatrices(secretKey), secretKey)),
            ea(context.ea->getCx())
        {}

        virtual void SetUp() override
        {
            if (helib_test::verbose) {
                ea.getPAlgebra().printout();
                std::cout << "r = " << context.alMod.getR() << std::endl;
                std::cout << "ctxtPrimes="<<context.ctxtPrimes
                    << ", specialPrimes="<<context.specialPrimes<<std::endl<<std::endl;
            }
        }

        virtual void TearDown() override
        {
            cleanupGlobals();
        }

};
 
TEST_P(Test_CKKS, negating_ciphertext_works)
{
  Ctxt c1(publicKey);
  std::vector<cx_double> vd1, vd2;
  NTL::ZZX poly;
  NTL::xdouble rf, pm;

  ea.random(vd1);
  ea.encrypt(c1, publicKey, vd1);
  rf = c1.getRatFactor();
  pm = c1.getPtxtMag();
  ea.encode(poly, vd1, /*size=*/ 1.0);
  c1.negate();
  ea.decrypt(c1, secretKey, vd2);

  negateVec(vd1);

  EXPECT_TRUE(cx_equals(vd2, vd1, NTL::conv<double>(epsilon)))
                  << "  max(vd1)=" << largestCoeff(vd1)
                  << ", max(vd2)=" << largestCoeff(vd2)
                  << ", maxDiff=" << calcMaxDiff(vd1,vd2) << std::endl << std::endl;
  EXPECT_EQ(rf, c1.getRatFactor());
  EXPECT_EQ(pm, c1.getPtxtMag());
}

TEST_P(Test_CKKS, adding_poly_constant_to_ciphertext_works)
{
  Ctxt c1(publicKey);
  std::vector<cx_double> vd1, vd2, vd3;
  NTL::ZZX poly;
  NTL::xdouble rf, pm;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  rf = c1.getRatFactor();
  pm = c1.getPtxtMag();
  ea.encode(poly, vd2, /*size=*/ 1.0);
  c1.addConstantCKKS(poly);
  ea.decrypt(c1, secretKey, vd3);

  add(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
                  << "  max(vd1)=" << largestCoeff(vd1)
                  << ", max(vd3)=" << largestCoeff(vd3)
                  << ", maxDiff=" << calcMaxDiff(vd1,vd3) << std::endl << std::endl;
}

TEST_P(Test_CKKS, adding_negated_poly_constant_to_ciphertext_works)
{
  Ctxt c1(publicKey);
  std::vector<cx_double> vd1, vd2, vd3;
  NTL::ZZX poly;
  NTL::xdouble rf, pm;

  ea.random(vd1);
  ea.random(vd2);
  negateVec(vd2);
  ea.encrypt(c1, publicKey, vd1);
  rf = c1.getRatFactor();
  pm = c1.getPtxtMag();
  ea.encode(poly, vd2, /*size=*/ 1.0);
  c1.addConstantCKKS(poly);
  ea.decrypt(c1, secretKey, vd3);

  add(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
                  << "  max(vd1)=" << largestCoeff(vd1)
                  << ", max(vd3)=" << largestCoeff(vd3) << std::endl
                  << ", maxDiff=" << calcMaxDiff(vd1,vd3) 
                  << ", maxRelDiff=" << calcMaxRelDiff(vd1,vd3) << std::endl << std::endl;
}

TEST_P(Test_CKKS, multiplying_poly_constant_to_ciphertext_works)
{
  Ctxt c1(publicKey);
  std::vector<cx_double> vd1, vd2, vd3;
  NTL::ZZX poly;
  NTL::xdouble rf, pm;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  rf = c1.getRatFactor();
  pm = c1.getPtxtMag();
  ea.encode(poly, vd2, /*size=*/ 1.0);
  c1.multByConstantCKKS(poly);
  ea.decrypt(c1, secretKey, vd3);

  mul(vd1, vd2);
  rf *= ea.encodeScalingFactor();

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
                  << "  max(vd1)=" << largestCoeff(vd1)
                  << ", max(vd3)=" << largestCoeff(vd3)
                  << ", maxDiff=" << calcMaxDiff(vd1,vd3) << std::endl << std::endl;
  EXPECT_EQ(rf, c1.getRatFactor());
  EXPECT_EQ(pm, c1.getPtxtMag());
}

TEST_P(Test_CKKS, adding_double_to_ciphertext_works)
{
  Ctxt c1(publicKey);
  std::vector<cx_double> vd1, vd2;
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
                  << "  max(vd1)=" << largestCoeff(vd1)
                  << ", max(vd2)=" << largestCoeff(vd2)
                  << ", maxDiff=" << calcMaxDiff(vd1,vd2) << std::endl << std::endl;
}

TEST_P(Test_CKKS, multiplying_double_to_ciphertext_works)
{
  Ctxt c1(publicKey);
  std::vector<cx_double> vd1, vd2, vd0;
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
                  << "  max(vd1)=" << largestCoeff(vd1)
                  << ", max(vd2)=" << largestCoeff(vd2)
                  << ", maxDiff=" << calcMaxDiff(vd1,vd2) << std::endl
		  << ", ptxtMag=" << c1.getPtxtMag() << std::endl
		  << ", vd[0]=" << vd[0] << std::endl;
  EXPECT_EQ(rf, c1.getRatFactor());
  EXPECT_EQ(pm, c1.getPtxtMag());
}

TEST_P(Test_CKKS, getting_the_complex_conjugate_works)
{
  Ctxt c1(publicKey);
  std::vector<cx_double> vd1, vd2;
  NTL::xdouble rf, pm;
  NTL::ZZX poly;

  ea.random(vd1);
  ea.encrypt(c1, publicKey, vd1);
  rf = c1.getRatFactor();
  pm = c1.getPtxtMag();
  ea.encode(poly, vd1, /*size=*/ 1.0);
  c1.complexConj();
  ea.decrypt(c1, secretKey, vd2);

  conjVec(vd1);

  EXPECT_TRUE(cx_equals(vd2, vd1, NTL::conv<double>(epsilon)))
                  << "  max(vd1)=" << largestCoeff(vd1)
                  << ", max(vd2)=" << largestCoeff(vd2)
                  << ", maxDiff=" << calcMaxDiff(vd1,vd2) << std::endl << std::endl;
  EXPECT_EQ(rf, c1.getRatFactor());
  EXPECT_EQ(pm, c1.getPtxtMag());
}

TEST_P(Test_CKKS, adding_ciphertexts_works)
{
  Ctxt c1(publicKey), c2(publicKey);
  std::vector<cx_double> vd1, vd2, vd3;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  ea.encrypt(c2, publicKey, vd2);
  c1 += c2;
  ea.decrypt(c1, secretKey, vd3);

  add(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
                  << "  max(vd1)=" << largestCoeff(vd1)
                  << ", max(vd3)=" << largestCoeff(vd3)
                  << ", maxDiff=" << calcMaxDiff(vd1,vd3) << std::endl << std::endl;
}

TEST_P(Test_CKKS, subtracting_ciphertexts_works)
{
  Ctxt c1(publicKey), c2(publicKey);
  std::vector<cx_double> vd1, vd2, vd3;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  ea.encrypt(c2, publicKey, vd2);
  c1 -= c2;
  ea.decrypt(c1, secretKey, vd3);

  sub(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
                  << "  max(vd1)=" << largestCoeff(vd1)
                  << ", max(vd3)=" << largestCoeff(vd3)
                  << ", maxDiff=" << calcMaxDiff(vd1,vd3) << std::endl << std::endl;
}

TEST_P(Test_CKKS, raw_multiplication_of_ciphertexts_works)
{
  Ctxt c1(publicKey), c2(publicKey);
  std::vector<cx_double> vd1, vd2, vd3;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  ea.encrypt(c2, publicKey, vd2);
  NTL::xdouble expectedPtxtMag = c1.getPtxtMag() * c2.getPtxtMag();
  c1 *= c2;
  ea.decrypt(c1, secretKey, vd3);

  mul(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
                  << "  max(vd1)=" << largestCoeff(vd1)
                  << ", max(vd3)=" << largestCoeff(vd3)
                  << ", maxDiff=" << calcMaxDiff(vd1,vd3) << std::endl << std::endl;
  EXPECT_EQ(expectedPtxtMag, c1.getPtxtMag());
}

TEST_P(Test_CKKS, high_level_multiplication_of_ciphertexts_works)
{
  Ctxt c1(publicKey), c2(publicKey);
  std::vector<cx_double> vd1, vd2, vd3;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  ea.encrypt(c2, publicKey, vd2);
  NTL::xdouble expectedPtxtMag = c1.getPtxtMag() * c2.getPtxtMag();
  c1.multiplyBy(c2);
  ea.decrypt(c1, secretKey, vd3);

  mul(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd3, vd1, NTL::conv<double>(epsilon)))
                  << "  max(vd1)=" << largestCoeff(vd1)
                  << ", max(vd3)=" << largestCoeff(vd3)
                  << ", maxDiff=" << calcMaxDiff(vd1,vd3) << std::endl << std::endl;
  EXPECT_EQ(expectedPtxtMag, c1.getPtxtMag());
}

TEST_P(Test_CKKS, multiplying_ciphertext_by_negative_constant_and_then_adding_to_other_ciphertext_works)
{
  Ctxt c1(publicKey), c2(publicKey);
  std::vector<cx_double> vd1, vd2, vd3, vd4;
  NTL::ZZX poly;

  ea.random(vd1);
  ea.random(vd2);
  ea.encrypt(c1, publicKey, vd1);
  ea.encrypt(c2, publicKey, vd2);
  c1.multByConstantCKKS(std::make_pair<long, long>(-5, 2));
  c1 += c2;
  ea.decrypt(c1, secretKey, vd4);

  mul(vd1,std::vector<cx_double>(ea.size(),-5/2.0));
  add(vd1, vd2);

  EXPECT_TRUE(cx_equals(vd4, vd1, epsilon))
                  << "  max(vd1)=" << largestCoeff(vd1)
                  << ", max(vd4)=" << largestCoeff(vd4) << std::endl
                  << ", maxDiff=" << calcMaxDiff(vd1,vd4) << std::endl << std::endl;
}

#endif // FFT_ARMA

INSTANTIATE_TEST_SUITE_P(typical_parameters, Test_CKKS, ::testing::Values(
            //SLOW
            Parameters(1024, 20, 150, 0.01)
            //FAST
            //Parameters(128, 20, 150, 0.01)
            ));

} // namespace
