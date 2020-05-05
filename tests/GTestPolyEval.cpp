/* Copyright (C) 2012-2020 IBM Corp.
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
#include <helib/polyEval.h>
#include <helib/EncryptedArray.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {
struct Parameters
{
  const long p; //   p is the plaintext base
  const long r; //   r is the lifting
  const long m; //   m is a specific cyclotomic ring
  const long d; //   d is the polynomial degree
  const long k; //   k is the baby-step parameter
  const long max_d;
  const long L;
  const bool isMonic;

  Parameters(long p,
             long r,
             long m,
             long d,
             long k,
             long max_d,
             long L,
             bool isMonic) :
      p(p), r(r), m(m), d(d), k(k), max_d(max_d), L(L), isMonic(isMonic){};

  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "p=" << params.p << ","
              << "r=" << params.r << ","
              << "m=" << params.m << ","
              << "d=" << params.d << ","
              << "k=" << params.k << ","
              << "max_d=" << params.max_d << ","
              << "L=" << params.L << ","
              << "isMonic=" << params.isMonic << "}";
  };
};

class GTestPolyEval : public ::testing::TestWithParam<Parameters>
{
protected:
  long p;
  long r;
  long d;
  long max_d;
  long L;
  bool isMonic;
  long m;
  long k;
  helib::Context context;
  long p2r;
  std::shared_ptr<helib::EncryptedArray> ea;
  helib::SecKey secretKey;
  const helib::PubKey& publicKey;

  GTestPolyEval() :
      p(GetParam().p),
      r(GetParam().r),
      d(GetParam().d),
      max_d(GetParam().max_d),
      L(GetParam().L),
      isMonic(GetParam().isMonic),
      m(GetParam().m),
      k(GetParam().k),
      context((helib::setDryRun(helib_test::dry), m), p, r),
      p2r(context.alMod.getPPowR()),
      ea(std::make_shared<helib::EncryptedArray>(
          (helib::buildModChain(context, L, /*c=*/3), context))),
      secretKey(context),
      publicKey((secretKey.GenSecKey(), secretKey))
      //  addSome1DMatrices(secretKey); // compute key-switching matrices
      {};

  virtual void SetUp() override
  {
    helib::setupDebugGlobals(&secretKey, ea);

    if (!helib_test::noPrint)
      std::cout << (helib::isDryRun() ? "* dry run, " : "* ") << "degree-" << d
                << ", m=" << m << ", L=" << L << ", p^r=" << p2r << std::endl;
  };

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }
};

TEST_P(GTestPolyEval, encryptedPolynomialsEvaluateAtEncryptedPointCorrectly)
{
  NTL::zz_pBak bak;
  bak.save();
  NTL::zz_p::init(p);
  NTL::zz_pXModulus phimX = NTL::conv<NTL::zz_pX>(ea->getPAlgebra().getPhimX());

  // Choose random plaintext polynomials
  NTL::zz_pX pX = NTL::random_zz_pX(deg(phimX) - 1);
  NTL::Vec<NTL::zz_pX> ppoly(NTL::INIT_SIZE, d);
  for (long i = 0; i < ppoly.length(); i++)
    random(ppoly[i], deg(phimX) - 1);

  // Evaluate the non-encrypted polynomial
  NTL::zz_pX pres =
      (ppoly.length() > 0) ? ppoly[ppoly.length() - 1] : NTL::zz_pX::zero();
  for (long i = ppoly.length() - 2; i >= 0; i--) {
    MulMod(pres, pres, pX, phimX);
    pres += ppoly[i];
  }

  // Encrypt the random polynomials
  helib::Ctxt cX(publicKey);
  NTL::Vec<helib::Ctxt> cpoly(NTL::INIT_SIZE, d, cX);

  secretKey.Encrypt(cX, NTL::conv<NTL::ZZX>(pX));

  for (long i = 0; i < ppoly.length(); i++)
    secretKey.Encrypt(cpoly[i], NTL::conv<NTL::ZZX>(ppoly[i]));

  // Evaluate the encrypted polynomial
  helib::polyEval(cX, cpoly, cX);

  // Compare the results
  NTL::ZZX ret;
  secretKey.Decrypt(ret, cX);
  NTL::zz_pX cres = NTL::conv<NTL::zz_pX>(ret);
  EXPECT_EQ(cres, pres) << "encrypted poly MISMATCH";
}

TEST_P(GTestPolyEval, evaluatePolynomialOnCiphertext)
{
  // evaluate at random points (at least one co-prime with p)
  std::vector<long> x;
  ea->random(x);
  while (NTL::GCD(x[0], p) != 1) {
    x[0] = NTL::RandomBnd(p2r);
  }
  helib::Ctxt inCtxt(publicKey), outCtxt(publicKey);
  ea->encrypt(inCtxt, publicKey, x);

  NTL::ZZX poly;
  for (long i = d; i >= 0; i--)
    SetCoeff(poly, i, NTL::RandomBnd(p2r)); // coefficients are random
  if (isMonic)
    SetCoeff(poly, d); // set top coefficient to 1

  // Evaluate poly on the ciphertext
  helib::polyEval(outCtxt, poly, inCtxt, k);

  // Check the result
  std::vector<long> y;
  ea->decrypt(outCtxt, secretKey, y);
  for (long i = 0; i < ea->size(); i++) {
    EXPECT_EQ(helib::polyEvalMod(poly, x[i], p2r), y[i])
        << "plaintext poly MISMATCH\n";
  }
}

std::vector<Parameters> getParameters()
{
  std::vector<Parameters> allParams;

  // SLOW
  const long p = 7;
  const long r = 2;
  long m = 0;
  const long k = 0;
  long d = 34;

  // FAST
  // const long p = 3;
  // const long r = 2;
  //      long m = 91;
  // const long k = 0;
  // long d = -1;

  const long max_d = (d <= 0) ? 35 : d;
  const long L = (7 + NTL::NextPowerOfTwo(max_d)) * 30;

  if (m < 2) {
    m = helib::FindM(
        /*secprm=*/80,
        L,
        /*c=*/3,
        p,
        1,
        0,
        m,
        !helib_test::noPrint);
  }

  // Test both monic and non-monic polynomials of this degree
  if (d >= 0) {
    allParams.emplace_back(Parameters{p, r, m, d, k, max_d, L, false});
    allParams.emplace_back(Parameters{p, r, m, d, k, max_d, L, true});
  } else {
    // Test degrees 1 to 3 and 25 through 35
    for (d = 1; d <= 3; d += 2)
      allParams.emplace_back(Parameters{p, r, m, d, k, max_d, L, true});
    for (d = 25; d <= 33; d += 2)
      allParams.emplace_back(Parameters{p, r, m, d, k, max_d, L, true});
  }
  return allParams;
}

INSTANTIATE_TEST_SUITE_P(manyDegrees,
                         GTestPolyEval,
                         ::testing::ValuesIn(getParameters()));

} // namespace
