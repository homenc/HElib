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

/**
 * @file matmul.h
 * @brief some matrix / linear algebra stuff
 */

#include <numeric>

#include <NTL/BasicThreadPool.h>

#include <helib/matmul.h>
#include <helib/fhe_stats.h>
#include <helib/debugging.h>

#if (defined(__unix__) || defined(__unix) || defined(unix))
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

struct Parameters
{
  Parameters(long m,
             long r,
             long bits,
             long nt,
             int force_bsgs,
             int force_hoist,
             int ks_strategy) :
      m(m),
      r(r),
      bits(bits),
      nt(nt),
      force_bsgs(force_bsgs),
      force_hoist(force_hoist),
      ks_strategy(ks_strategy){};

  const long m;          // defines the cyclotomic polynomial Phi_m(X)
  const long r;          // Bit precision
  const long bits;       // # bits in the modulus chain
  const long nt;         // # threads
  const int force_bsgs;  // 1 to force on, -1 to force off
  const int force_hoist; // -1 to force off
  const int ks_strategy; // 0: default, 1: full, 2: bsgs, 3: minimal

  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "m=" << params.m << ","
              << "r=" << params.r << ","
              << "bits=" << params.bits << ","
              << "nt=" << params.nt << ","
              << "force_bsgs=" << params.force_bsgs << ","
              << "force_hoist=" << params.force_hoist << ","
              << "ks_strategy=" << params.ks_strategy << "}";
  }
};

class TestMatmulCKKS : public ::testing::TestWithParam<Parameters>
{
protected:
  static void setGlobals(long force_bsgs, long force_hoist)
  {
    helib::fhe_test_force_bsgs = force_bsgs;
    helib::fhe_test_force_hoist = force_hoist;
  }

  const long m;
  const long r;
  const long bits;
  const long nt;

  helib::Context context;
  helib::SecKey secretKey;
  const helib::PubKey publicKey;
  const helib::EncryptedArray& ea;

  TestMatmulCKKS() :
      m(GetParam().m),
      r(GetParam().r),
      bits(GetParam().bits),
      nt(GetParam().nt),
      context(helib::ContextBuilder<helib::CKKS>()
                  .m(m)
                  .precision(r)
                  .bits(bits)
                  .build()),
      secretKey(context),
      publicKey(keySetup(secretKey, GetParam().ks_strategy)),
      ea(context.getEA())
  {}

  static helib::SecKey& keySetup(helib::SecKey& secretKey, int ks_strategy)
  {
    secretKey.GenSecKey();
    // We call addSomeFrbMatrices for all strategies except minimal
    switch (ks_strategy) {
    case 0:
      addSome1DMatrices(secretKey);
      addSomeFrbMatrices(secretKey);
      break;
    case 1:
      add1DMatrices(secretKey);
      addSomeFrbMatrices(secretKey);
      break;
    case 2:
      addBSGS1DMatrices(secretKey);
      addSomeFrbMatrices(secretKey);
      break;
    case 3:
      addMinimal1DMatrices(secretKey);
      addMinimalFrbMatrices(secretKey);
      break;

    default:
      NTL::Error("bad ks_strategy");
    }
    return secretKey;
  }

  virtual void SetUp() override
  {
    if (helib_test::verbose) {
      context.getZMStar().printout();
      std::cout
          << "# small primes = " << context.getSmallPrimes().card() << "\n"
          << "# ctxt primes = " << context.getCtxtPrimes().card() << "\n"
          << "# bits in ctxt primes = "
          << long(context.logOfProduct(context.getCtxtPrimes()) / log(2.0) +
                  0.5)
          << "\n"
          << "# special primes = " << context.getSpecialPrimes().card() << "\n"
          << "# bits in special primes = "
          << long(context.logOfProduct(context.getSpecialPrimes()) / log(2.0) +
                  0.5)
          << "\n";

      helib::fhe_stats = true;
    }
    helib::setupDebugGlobals(&secretKey, context.shareEA());
  }

  virtual void TearDown() override
  {
    if (helib_test::verbose) {
      helib::printAllTimers();
#if (defined(__unix__) || defined(__unix) || defined(unix))
      struct rusage rusage;
      getrusage(RUSAGE_SELF, &rusage);
      std::cout << "  rusage.ru_maxrss=" << rusage.ru_maxrss << std::endl;
#endif
      helib::print_stats(std::cout);
    }
    helib::cleanupDebugGlobals();
  }
};

TEST_P(TestMatmulCKKS, vectorToMatrixMultiplication)
{
  std::vector<double> v(ea.size());
  std::iota(v.begin(), v.end(), 1);

  helib::MatMul_CKKS_Complex mat(context, [&v](long i, long j) {
    return ((i + j) % v.size()) / double(v.size());
  });
  // Note the use of a "lambda": this allows for quite
  // general ways to describe a matrix with minimal fuss

  helib::Ctxt ctxt(publicKey);
  // Initialize ptxt with the vector v
  helib::PtxtArray ptxt(context, v);

  // Encrypt ptxt. We have to supply *some* upper bound on the magnitude
  // of the slots.
  ptxt.encrypt(ctxt);

  // Perform the linear transformation on the encrypted data
  helib::EncodedMatMul_CKKS emat(mat);
  // emat.upgrade();
  ctxt *= emat;

  // We can also do it this way:
  // ctxt *= mat;

  // Perform the linear transformation on the plaintext data
  ptxt *= mat;

  helib::PtxtArray ptxt1(context);
  ptxt1.decrypt(ctxt, secretKey);

  // w1 is the result of performing the transformation on
  // the encrypted data
  std::vector<double> w1;
  ptxt1.store(w1);

  // w is the result of performing the transformation on
  // the encrypted data
  std::vector<double> w;
  ptxt.store(w);

  for (long i = 0; i < ea.size(); ++i) {
    EXPECT_NEAR(w[i], w1[i], 0.015);
  }
}

// clang-format off
INSTANTIATE_TEST_SUITE_P(typicalParameters, TestMatmulCKKS, ::testing::Values(
      Parameters(/*m=*/16, /*r=*/10, /*bits=*/200, /*nt=*/1, /*force_bsgs=*/0, /*force_hoist=*/0, /*ks_strategy=*/0)
      ));
// clang-format on

} // namespace
