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
    Parameters(long R, long m, long r, long L, double epsilon) :
        R(R),
        m(m),
        r(r),
        L(L),
        epsilon(epsilon)
    {};

    const long R; // number of rounds
    const long m;
    const long r;
    const long L;
    const double epsilon;

    friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
    {
        return os << "{" <<
            "R=" << params.R << "," <<
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
    return (calcMaxRelDiff(v1,v2) < epsilon);
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
  void add(std::vector<cx_double>& to, const std::vector<cx_double>& from)
  {
    if (to.size() < from.size())
      to.resize(from.size(), 0);
    for (long i=0; i<from.size(); i++) to[i] += from[i];
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
  void rotate(std::vector<cx_double>& p, long amt)
  {
    long sz = p.size();
    std::vector<cx_double> tmp(sz);
    for (long i=0; i<sz; i++)
      tmp[((i+amt)%sz +sz)%sz] = p[i];
    p = tmp;
  }

#ifndef FFT_ARMA
class GTest_approxNums : public ::testing::TestWithParam<Parameters> {};
TEST_P(GTest_approxNums, error)
{
    FAIL() << "Cannot run approxNums test when not using armadillo.";
}
#else

class GTest_approxNums : public ::testing::TestWithParam<Parameters> {
    protected:
        const long R;
        const long m;
        const long r;
        const long L;
        const double epsilon;

        FHEcontext context;
        FHESecKey secretKey;
        const FHEPubKey publicKey;
        const EncryptedArrayCx& ea;

        GTest_approxNums () :
            R(GetParam().R),
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
//          if (helib_test::debug) {
//            dbgKey = &secretKey;
//            dbgEa = const_cast<EncryptedArray*>(context.ea);
//          }
        }

        virtual void TearDown() override
        {
            cleanupGlobals();
        }

};

TEST_P(GTest_approxNums, basic_arithmetic_works)
{
  if (helib_test::verbose)  std::cout << "Test Arithmetic ";
  // Test objects

  Ctxt c1(publicKey), c2(publicKey), c3(publicKey);
  
  std::vector<cx_double> vd;
  std::vector<cx_double> vd1, vd2, vd3;
  ea.random(vd1);
  ea.random(vd2);

  // test encoding of shorter vectors
  vd1.resize(vd1.size()-2);
  ea.encrypt(c1, publicKey, vd1, /*size=*/1.0);
  vd1.resize(vd1.size()+2, 0.0);

  ea.encrypt(c2, publicKey, vd2, /*size=*/1.0);


  // Test - Multiplication  
  c1 *= c2;
  for (long i=0; i<lsize(vd1); i++) vd1[i] *= vd2[i];

  NTL::ZZX poly;
  ea.random(vd3);
  ea.encode(poly, vd3, /*size=*/1.0);
  c1.addConstant(poly); // vd1*vd2 + vd3
  for (long i=0; i<lsize(vd1); i++) vd1[i] += vd3[i];

  // Test encoding, encryption of a single number
  double xx = NTL::RandomLen_long(16)/double(1L<<16); // random in [0,1]
  ea.encryptOneNum(c2, publicKey, xx);
  c1 += c2;
  for (auto& x : vd1) x += xx;

  // Test - Multiply by a mask
  std::vector<long> mask(lsize(vd1), 1);
  for (long i=0; i*(i+1)<lsize(mask); i++) {
    mask[i*i] = 0;
    mask[i*(i+1)] = -1;
  }

  ea.encode(poly,mask, /*size=*/1.0);
  c1.multByConstant(poly); // mask*(vd1*vd2 + vd3)
  for (long i=0; i<lsize(vd1); i++) vd1[i] *= mask[i];

  // Test - Addition
  ea.random(vd3);
  ea.encrypt(c3, publicKey, vd3, /*size=*/1.0);
  c1 += c3;
  for (long i=0; i<lsize(vd1); i++) vd1[i] += vd3[i];

  c1.negate();
  c1.addConstant(NTL::to_ZZ(1));
  for (long i=0; i<lsize(vd1); i++) vd1[i] = 1.0 - vd1[i];

  // Diff between approxNums HE scheme and plaintext floating  
  ea.decrypt(c1, secretKey, vd);
#ifdef DEBUG_PRINTOUT
  printVec(std::cout<<"res=", vd, 10)<<std::endl;
  printVec(std::cout<<"vec=", vd1, 10)<<std::endl;
#endif
  if (helib_test::verbose)
    std::cout << "(max |res-vec|_{infty}="<< calcMaxDiff(vd, vd1) << "): ";

  EXPECT_TRUE(cx_equals(vd, vd1, NTL::conv<double>(epsilon * c1.getPtxtMag())))
                  << "  max(vd)=" << largestCoeff(vd)
                  << ", max(vd1)=" << largestCoeff(vd1)
                  << ", maxDiff=" << calcMaxDiff(vd,vd1) << std::endl << std::endl;
}

TEST_P(GTest_approxNums, complex_arithmetic_works)
{
  // Test complex conjugate
  Ctxt c1(publicKey), c2(publicKey);

  std::vector<cx_double> vd;
  std::vector<cx_double> vd1, vd2;
  ea.random(vd1);
  ea.random(vd2);
   
  ea.encrypt(c1, publicKey, vd1, /*size=*/1.0);
  ea.encrypt(c2, publicKey, vd2, /*size=*/1.0);

  if (helib_test::verbose)
    std::cout << "Test Conjugate: ";
  for_each(vd1.begin(), vd1.end(), [](cx_double& d){d=std::conj(d);});
  c1.complexConj();  
  ea.decrypt(c1, secretKey, vd);
#ifdef DEBUG_PRINTOUT
  printVec(cout<<"vd1=", vd1, 10) << std::endl;
  printVec(cout<<"res=", vd, 10) << std::endl;
#endif
  EXPECT_TRUE(cx_equals(vd, vd1, NTL::conv<double>(epsilon * c1.getPtxtMag())))
                    << "  max(vd)="  << largestCoeff(vd)
                    << ", max(vd1)=" << largestCoeff(vd1)
                    << ", maxDiff="  << calcMaxDiff(vd,vd1) << std::endl << std::endl;;

  // Test that real and imaginary parts are actually extracted.
  Ctxt realCtxt(c2), imCtxt(c2);
  std::vector<cx_double> realParts(vd2), real_dec;
  std::vector<cx_double> imParts(vd2), im_dec;

  if (helib_test::verbose)
    std::cout << "Test Real and Im parts: ";
  for_each(realParts.begin(), realParts.end(), [](cx_double& d){d=std::real(d);});
  for_each(imParts.begin(), imParts.end(), [](cx_double& d){d=std::imag(d);});

  ea.extractRealPart(realCtxt);
  ea.decrypt(realCtxt, secretKey, real_dec);

  ea.extractImPart(imCtxt);
  ea.decrypt(imCtxt, secretKey, im_dec);

#ifdef DEBUG_PRINTOUT
  printVec(std::cout << "vd2=", vd2, 10) << std::endl;
  printVec(std::cout << "real=", realParts, 10) << std::endl;
  printVec(std::cout << "res=", real_dec, 10) << std::endl;
  printVec(std::cout << "im=", imParts, 10) << std::endl;
  printVec(std::cout << "res=", im_dec, 10) << std::endl;
#endif
  EXPECT_TRUE(cx_equals(realParts, real_dec, NTL::conv<double>(epsilon * realCtxt.getPtxtMag())))
                    << "  max(re)="  << largestCoeff(realParts)
                    << ", max(re1)=" << largestCoeff(real_dec)
                    << ", maxDiff="  << calcMaxDiff(realParts,real_dec) << std::endl;
  EXPECT_TRUE(cx_equals(imParts, im_dec, NTL::conv<double>(epsilon * imCtxt.getPtxtMag())))
                    << "  max(im)="  << largestCoeff(imParts)
                    << ", max(im1)=" << largestCoeff(im_dec)
                    << ", maxDiff="  << calcMaxDiff(imParts,im_dec) << std::endl << std::endl;
}

TEST_P(GTest_approxNums, rotates_and_shifts_work)
{
  std::srand(std::time(0)); // set seed, current time.
  int nplaces = rand() % static_cast<int>(ea.size()/2.0) + 1;

  if (helib_test::verbose)
    std::cout << "Test Rotation of " << nplaces << ": ";  

  Ctxt c1(publicKey);
  std::vector<cx_double> vd1;
  std::vector<cx_double> vd_dec;
  ea.random(vd1);
  ea.encrypt(c1, publicKey, vd1, /*size=*/1.0);

#ifdef DEBUG_PRINTOUT
  printVec(std::cout << "vd1=", vd1, 10) << std::endl;
#endif
  std::rotate(vd1.begin(), vd1.end() - nplaces, vd1.end());
  ea.rotate(c1, nplaces);
  ea.decrypt(c1, secretKey, vd_dec);
#ifdef DEBUG_PRINTOUT
  printVec(std::cout << "vd1(rot)=", vd1, 10) << std::endl;
  printVec(std::cout << "res: ", vd_dec, 10) << std::endl;
#endif

  EXPECT_TRUE(cx_equals(vd1, vd_dec, NTL::conv<double>(epsilon * c1.getPtxtMag())))
                  << "  max(vd)=" << largestCoeff(vd_dec)
                  << ", max(vd1)=" << largestCoeff(vd1)
                  << ", maxDiff=" << calcMaxDiff(vd_dec,vd1) << std::endl << std::endl;
}
  
TEST_P(GTest_approxNums, general_ops_works) {
  /************** Each round consists of the following:
   1. c1.multiplyBy(c0)
   2. c0 += random constant
   3. c2 *= random constant
   4. tmp = c1
   5. ea.rotate(tmp, random amount in [-nSlots/2, nSlots/2])
   6. c2 += tmp
   7. ea.rotate(c2, random amount in [1-nSlots, nSlots-1])
   8. c1.negate()
   9. c3.multiplyBy(c2)
   10. c0 -= c3
   **************/
  long nslots = ea.size();
  char buffer[32];
    
  std::vector<cx_double> p0, p1, p2, p3;
  ea.random(p0);
  ea.random(p1);
  ea.random(p2);
  ea.random(p3);
    
  Ctxt c0(publicKey), c1(publicKey), c2(publicKey), c3(publicKey);
  ea.encrypt(c0, publicKey, p0, /*size=*/1.0);
  ea.encrypt(c1, publicKey, p1, /*size=*/1.0);
  ea.encrypt(c2, publicKey, p2, /*size=*/1.0);
  ea.encrypt(c3, publicKey, p3, /*size=*/1.0);
    
  resetAllTimers();
  FHE_NTIMER_START(Circuit);
    
    for (long i = 0; i < R; i++) {
      
      if (helib_test::verbose)
        std::cout << "*** round " << i << "..." << std::endl;
      
      long shamt = NTL::RandomBnd(2 * (nslots / 2) + 1) - (nslots / 2);
      // random number in [-nslots/2..nslots/2]
      long rotamt = NTL::RandomBnd(2 * nslots - 1) - (nslots - 1);
      // random number in [-(nslots-1)..nslots-1]
      
      // two random constants
      std::vector<cx_double> const1, const2;
      ea.random(const1);
      ea.random(const2);
      
      NTL::ZZX const1_poly, const2_poly;
      ea.encode(const1_poly, const1, /*size=*/1.0);
      ea.encode(const2_poly, const2, /*size=*/1.0);
      
      mul(p1, p0);     // c1.multiplyBy(c0)
      c1.multiplyBy(c0);
      if (helib_test::verbose) {
        CheckCtxt(c1, "c1*=c0");
      }
      EXPECT_TRUE(ciphertextMatches(ea, secretKey, p1, c1, epsilon));
      
      add(p0, const1); // c0 += random constant
      c0.addConstant(const1_poly);
      if (helib_test::verbose) {
        CheckCtxt(c0, "c0+=k1");
      }
      EXPECT_TRUE(ciphertextMatches(ea, secretKey, p0, c0, epsilon));
      
      mul(p2, const2); // c2 *= random constant
      c2.multByConstant(const2_poly);
      if (helib_test::verbose) {
        CheckCtxt(c2, "c2*=k2");
      }
      EXPECT_TRUE(ciphertextMatches(ea, secretKey, p2, c2, epsilon));
      
      std::vector<cx_double> tmp_p(p1); // tmp = c1
      Ctxt tmp(c1);
      sprintf(buffer, "tmp=c1>>=%d", (int)shamt);
      rotate(tmp_p, shamt); // ea.shift(tmp, random amount in [-nSlots/2,nSlots/2])
      ea.rotate(tmp, shamt);
      if (helib_test::verbose) {
        CheckCtxt(tmp, buffer);
      }
      EXPECT_TRUE(ciphertextMatches(ea, secretKey, tmp_p, tmp, epsilon));
      
      add(p2, tmp_p);  // c2 += tmp
      c2 += tmp;
      if (helib_test::verbose) {
        CheckCtxt(c2, "c2+=tmp");
      }
      EXPECT_TRUE(ciphertextMatches(ea, secretKey, p2, c2, epsilon));
      
      sprintf(buffer, "c2>>>=%d", (int)rotamt);
      rotate(p2, rotamt); // ea.rotate(c2, random amount in [1-nSlots, nSlots-1])
      ea.rotate(c2, rotamt);
      if (helib_test::verbose) {
        CheckCtxt(c2, buffer);
      }
      EXPECT_TRUE(ciphertextMatches(ea, secretKey, p2, c2, epsilon));
      
      negateVec(p1); // c1.negate()
      c1.negate();
      if (helib_test::verbose) {
        CheckCtxt(c1, "c1=-c1");
      }
      EXPECT_TRUE(ciphertextMatches(ea, secretKey, p1, c1, epsilon));
      
      mul(p3, p2); // c3.multiplyBy(c2)
      c3.multiplyBy(c2);
      if (helib_test::verbose) {
        CheckCtxt(c3, "c3*=c2");
      }
      EXPECT_TRUE(ciphertextMatches(ea, secretKey, p3, c3, epsilon));
      
      sub(p0, p3); // c0 -= c3
      c0 -= c3;
      if (helib_test::verbose) {
        CheckCtxt(c0, "c0=-c3");
      }
      EXPECT_TRUE(ciphertextMatches(ea, secretKey, p0, c0, epsilon));
      
    }
    
    c0.cleanUp();
    c1.cleanUp();
    c2.cleanUp();
    c3.cleanUp();
    
    FHE_NTIMER_STOP(Circuit);
    
    std::vector<cx_double> pp0, pp1, pp2, pp3;
    
    ea.decrypt(c0, secretKey, pp0);
    ea.decrypt(c1, secretKey, pp1);
    ea.decrypt(c2, secretKey, pp2);
    ea.decrypt(c3, secretKey, pp3);
  
    if (helib_test::verbose) {
      std::cout << "Test " << R << " rounds of mixed operations, ";
    }
    EXPECT_TRUE(cx_equals(pp0, p0, NTL::conv<double>(epsilon * c0.getPtxtMag()))
                && cx_equals(pp1, p1, NTL::conv<double>(epsilon * c1.getPtxtMag()))
                && cx_equals(pp2, p2, NTL::conv<double>(epsilon * c2.getPtxtMag()))
                && cx_equals(pp3, p3, NTL::conv<double>(epsilon * c3.getPtxtMag())))
                               << "  max(p0)="  << largestCoeff(p0)
                               << ", max(pp0)=" << largestCoeff(pp0)
                               << ", maxDiff="  << calcMaxDiff(p0,pp0)
                               << std::endl
                               << "  max(p1)="  << largestCoeff(p1)
                               << ", max(pp1)=" << largestCoeff(pp1)
                               << ", maxDiff="  << calcMaxDiff(p1,pp1)
                               << std::endl
                               << "  max(p2)="  << largestCoeff(p2)
                               << ", max(pp2)=" << largestCoeff(pp2)
                               << ", maxDiff="  << calcMaxDiff(p2,pp2)
                               << std::endl
                               << "  max(p3)="  << largestCoeff(p3)
                               << ", max(pp3)=" << largestCoeff(pp3)
                               << ", maxDiff="  << calcMaxDiff(p3,pp3)
                               << std::endl << std::endl;
  
    if (helib_test::verbose) {
        std::cout << std::endl;
        printAllTimers();
        std::cout << std::endl;
      }
      resetAllTimers();
  }
  
#endif // FFT_ARMA

INSTANTIATE_TEST_SUITE_P(typical_parameters, GTest_approxNums, ::testing::Values(
            //SLOW
            Parameters(1, 1024, 8, 150, 0.01)
            //FAST
            //Parameters(1, 128, 8, 150, 0.01)
            ));
//if (R<=0) R=1;
//if (R<=2)
//  L = 100*R;
//else
//  L = 220*(R-1);

} // namespace
