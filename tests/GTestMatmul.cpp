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
/**
 * @file matmul.h
 * @brief some matrix / linear algebra stuff
 */
#include <helib/matmul.h>
#include <helib/fhe_stats.h>
#include <NTL/BasicThreadPool.h>

#if (defined(__unix__) || defined(__unix) || defined(unix))
#include <sys/time.h>
#include <sys/resource.h>
#endif

// Implementation of the various random matrices is found here
#include <helib/randomMatrices.h>
/*
 * Defined in this file are the following class templates:
 *
 *   class RandomMatrix: public MatMul1D_derived<type>
 *   class RandomMultiMatrix: public MatMul1D_derived<type>
 *   class RandomBlockMatrix: public BlockMatMul1D_derived<type>
 *   class RandomMultiBlockMatrix: public BlockMatMul1D_derived<type>
 *   class RandomFullMatrix: public MatMulFull_derived<type>
 *   class RandomFullBlockMatrix : public BlockMatMulFull_derived<type>
 *
 * Each of them has a corresponding build function, namely:
 *
 *   MatMul1D* buildRandomMatrix(const EncryptedArray& ea, long dim);
 *   MatMul1D* buildRandomMultiMatrix(const EncryptedArray& ea, long dim);
 test_common.h   BlockMatMul1D* buildRandomBlockMatrix(const EncryptedArray& ea,
 long dim);
 *   BlockMatMul1D* buildRandomMultiBlockMatrix(const EncryptedArray& ea, long
 dim);
 *   MatMulFull* buildRandomFullMatrix(const EncryptedArray& ea);
 *   BlockMatMulFull* buildRandomFullBlockMatrix(const EncryptedArray& ea);
 */

#include <helib/debugging.h>
#include "gtest/gtest.h"
#include "test_common.h"

namespace {

struct Parameters
{
  Parameters(long m,
             long p,
             long r,
             long L,
             long dim,
             long nt,
             long full,
             long block,
             std::vector<long> gens,
             std::vector<long> ords,
             int ks_strategy,
             int force_bsgs,
             int force_hoist) :
      m(m),
      p(p),
      r(r),
      L(L),
      dim(dim),
      nt(nt),
      full(full),
      block(block),
      gens(gens),
      ords(ords),
      ks_strategy(ks_strategy),
      force_bsgs(force_bsgs),
      force_hoist(force_hoist){};

  long m;                 // defines the cyclotomic polynomial Phi_m(X)
  long p;                 // plaintext base
  long r;                 // lifting
  long L;                 // # of levels in the modulus chain
  long dim;               // dimension along which to multiply
  long nt;                // # threads
  long full;              // 0: 1D, 1: full
  long block;             // 0: normal, 1: block
  std::vector<long> gens; // use specified vector of generators
                          // e.g., gens='[420 1105 1425]'
  std::vector<long> ords; // use specified vector of orders
                          // e.g., ords='[11 8 2]', negative means 'bad'
  int ks_strategy;        // 0: default, 1:full, 2:bsgs, 3:minimal
  int force_bsgs;         // 1 to force on, -1 to force off
  int force_hoist;        // -1 to force off

  // tell google test how to print Parameters
  // Removed as not used by templated tests
  // friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  //{
  //    return os << "{" <<
  //        "m=" << params.m << "," <<
  //        "p=" << params.p << "," <<
  //        "r=" << params.r << "," <<
  //        "L=" << params.L << "," <<
  //        "dim=" << params.dim << "," <<
  //        "nt=" << params.nt << "," <<
  //        "full=" << params.full <<  "," <<
  //        "block=" << params.block <<  "," <<
  //        "gens=" << helib::vecToStr(params.gens) << "," <<
  //        "ords=" << helib::vecToStr(params.ords) << "," <<
  //        "ks_strategy=" << params.ks_strategy << "," <<
  //        "force_bsgs=" << params.force_bsgs << "," <<
  //        "force_hoist=" << params.force_hoist <<
  //        "}";
  //};
};

template <typename Mat>
std::unique_ptr<Mat> buildMat(const helib::EncryptedArray& ea, long dim);

template <>
std::unique_ptr<helib::MatMul1D> buildMat(const helib::EncryptedArray& ea,
                                          long dim)
{
  return std::unique_ptr<helib::MatMul1D>{buildRandomMatrix(ea, dim)};
}
// template<> std::unique_ptr<helib::BlockMatMul1D> buildMat(const
// helib::EncryptedArray &ea, long dim)
//{
//    return std::unique_ptr<helib::BlockMatMul1D>{buildRandomBlockMatrix(ea,
//    dim)};
//};
// template<> std::unique_ptr<helib::MatMulFull> buildMat(const
// helib::EncryptedArray &ea, long dim)
//{
//    return std::unique_ptr<helib::MatMulFull>{buildRandomFullMatrix(ea)};
//};
// template<> std::unique_ptr<helib::BlockMatMulFull> buildMat(const
// helib::EncryptedArray &ea, long dim)
//{
//    return
//    std::unique_ptr<helib::BlockMatMulFull>{buildRandomFullBlockMatrix(ea)};
//};

template <typename T>
class GTestMatmul : public ::testing::Test
{
protected:
  static void setGlobals(long force_bsgs, long force_hoist)
  {
    helib::fhe_test_force_bsgs = force_bsgs;
    helib::fhe_test_force_hoist = force_hoist;
  };

  bool minimal;
  long m;
  long p;
  long r;
  long L;
  long dim;
  long nt;

  long full;
  long block;
  std::vector<long> gens;
  std::vector<long> ords;

  long ks_strategy;

  helib::Context context;
  helib::SecKey secretKey;
  const helib::PubKey& publicKey;
  helib::EncryptedArray ea;
  std::unique_ptr<typename T::MatrixType> matrixPtr;

  static helib::Context& setupContext(helib::Context& context, long L)
  {
    buildModChain(context, L, /*c=*/3);
    if (helib_test::verbose) {
      context.zMStar.printout();
      std::cout << "# small primes = " << context.smallPrimes.card() << "\n";
      std::cout << "# ctxt primes = " << context.ctxtPrimes.card() << "\n";
      std::cout << "# bits in ctxt primes = "
                << long(context.logOfProduct(context.ctxtPrimes) / log(2.0) +
                        0.5)
                << "\n";
      std::cout << "# special primes = " << context.specialPrimes.card()
                << "\n";
      std::cout << "# bits in special primes = "
                << long(context.logOfProduct(context.specialPrimes) / log(2.0) +
                        0.5)
                << "\n";
      helib::fhe_stats = true;
    }
    return context;
  }

  GTestMatmul() :
      minimal((setGlobals(T::parameters.force_bsgs, T::parameters.force_hoist),
               T::parameters.ks_strategy == 3)),
      m(T::parameters.m),
      p(T::parameters.p),
      r(T::parameters.r),
      L(T::parameters.L),
      dim(T::parameters.dim),
      nt((NTL::SetNumThreads(T::parameters.nt), T::parameters.nt)),
      full(T::parameters.full),
      block(T::parameters.block),
      gens(T::parameters.gens),
      ords(T::parameters.ords),
      ks_strategy(T::parameters.ks_strategy),
      context((helib::setTimersOn(), m), p, r, gens, ords),
      secretKey(setupContext(context, L)),
      publicKey((secretKey.GenSecKey(), secretKey)),
      ea(context, context.alMod), // encrypted array with "full slots"
      matrixPtr(buildMat<typename T::MatrixType>(ea, dim)){};

  virtual void SetUp()
  {
    // we call addSomeFrbMatrices for all strategies except minimal
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
      throw std::invalid_argument("bad ks_strategy");
    }

    if (helib_test::verbose) {
      {
        std::stringstream ss;
        ss << publicKey;
        std::cout << "\n  |pubKey|=" << ss.tellp() << " bytes";
      }
      std::cout << ", security=" << context.securityLevel() << std::endl;
      std::cout << "  #threads=" << NTL::AvailableThreads();
      if (block)
        std::cout << "block-size=" << (ea.getPAlgebra().getOrdP());
      std::cout << ", vector-dimension="
                << (full ? ea.size() : ea.sizeOfDimension(dim)) << ", ";
      if (!full)
        std::cout << (ea.size() / ea.sizeOfDimension(dim))
                  << " products in parallel, ";
    }
  };

  virtual void TearDown()
  {
    if (helib_test::verbose) {
      helib::printAllTimers(std::cout);
#if (defined(__unix__) || defined(__unix) || defined(unix))
      struct rusage rusage;
      getrusage(RUSAGE_SELF, &rusage);
      std::cout << "  rusage.ru_maxrss=" << rusage.ru_maxrss << std::endl;
#endif
      helib::print_stats(std::cout);
    }
    helib::cleanupDebugGlobals();
  };
};

// This is how we pack a matrix type and a reference to a Parameters object into
// a single type for gtest
template <typename T, const Parameters& params>
struct MatrixTypeAndParams
{
  typedef T MatrixType;
  constexpr const static Parameters& parameters = params;
};

// These need to be at namespace scope, and cannot be anonymous since one cannot
// provide rvalue references as template non-type parameters.
// clang-format off
//SLOW
Parameters oneDimensionalMatrixParams     (18631, 2, 1, 300, 0, 1, 0, 0, std::vector<long>{}            , std::vector<long>{}      , 0, 0, 0);
Parameters oneDimensionalBlockMatrixParams(24295, 2, 1, 300, 0, 1, 0, 1, std::vector<long>{16386, 16427}, std::vector<long>{42, 16}, 0, 0, 0);

// The above commented-out parameters were used before the new noise was brought in - it is too slow now.
//FAST
// Parameters oneDimensionalMatrixParams(91, 2, 1, 300, 0, 1, 0, 0, std::vector<long>{9, 3}, std::vector<long>{3, -2}, 0, 0, 0);
// clang-format on

using TypesToTest = ::testing::Types<
    MatrixTypeAndParams<helib::MatMul1D, oneDimensionalMatrixParams>,
    MatrixTypeAndParams<helib::MatMul1D, oneDimensionalBlockMatrixParams>>;

// Currently gtest does not intend on supporting the -Wall and -Wextra flags so
// this does not conform to the C++ standard. Until gtest changes, we need a
// pragma to ignore this warning.
#pragma GCC diagnostic push
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wgnu-zero-variadic-macro-arguments"
#endif
TYPED_TEST_SUITE(GTestMatmul, TypesToTest);
#pragma GCC diagnostic pop

TYPED_TEST(GTestMatmul, multipliesWithoutErrors)
{
  HELIB_NTIMER_START(EncodeMartix_MatMul);
  // Templated class so explicit 'this' necessary
  const typename TypeParam::MatrixType& mat = *(this->matrixPtr);
  typename TypeParam::MatrixType::ExecType mat_exec(mat, (this->minimal));
  mat_exec.upgrade();
  HELIB_NTIMER_STOP(EncodeMartix_MatMul);

  // choose a random plaintext vector and encrypt it
  helib::PlaintextArray v(this->ea);
  random(this->ea, v);

  // encrypt the random vector
  helib::Ctxt ctxt(this->secretKey);
  this->ea.encrypt(ctxt, this->secretKey, v);
  helib::Ctxt ctxt2 = ctxt;

  mat_exec.mul(ctxt);

  mul(v, mat); // multiply the plaintext vector

  helib::PlaintextArray v1(this->ea);
  this->ea.decrypt(ctxt, this->secretKey, v1); // decrypt the ciphertext vector

  EXPECT_TRUE(equals(this->ea, v, v1)); // check that we've got the right answer
}

} // namespace
