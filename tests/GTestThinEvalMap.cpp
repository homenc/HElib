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

/* Test_ThinEvalMap.cpp - Testing the evaluation map for thin bootstrapping
 */
#include <helib/helib.h>
#include <helib/EvalMap.h>
#include <NTL/BasicThreadPool.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

struct Parameters
{

  const long p;
  const long r;
  const long c;
  const long k;
  const long L;
  const long s;
  const long seed;
  const NTL::Vec<long> mvec;
  const std::vector<long> gens;
  const std::vector<long> ords;
  const long nthreads;
  const long useCache;

  Parameters(long p,
             long r,
             long c,
             long k,
             long L,
             long s,
             long seed,
             NTL::Vec<long> mvec,
             std::vector<long> gens,
             std::vector<long> ords,
             long nthreads,
             long useCache) :
      p(p),
      r(r),
      c(c),
      k(k),
      L(L),
      s(s),
      seed(seed),
      mvec(mvec),
      gens(gens),
      ords(ords),
      nthreads(nthreads),
      useCache(useCache){};

  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "p=" << params.p << ","
              << "r=" << params.r << ","
              << "c=" << params.c << ","
              << "k=" << params.k << ","
              << "L=" << params.L << ","
              << "s=" << params.s << ","
              << "seed=" << params.seed << ","
              << "mvec=" << params.mvec << ","
              << "gens=" << helib::vecToStr(params.gens) << ","
              << "ords=" << helib::vecToStr(params.ords) << ","
              << "nthreads=" << params.nthreads << ","
              << "useCache=" << params.useCache << "}";
  };
};

class GTestThinEvalMap : public ::testing::TestWithParam<Parameters>
{
protected:
  const long p;
  const long r;
  const long c;
  const long k;
  const long w;
  const long L;
  const long s;
  const long seed;
  const NTL::Vec<long> mvec;
  const std::vector<long> gens;
  const std::vector<long> ords;
  const long m;
  const long nthreads;
  const long useCache;
  helib::Context context;
  const long d;
  const long phim;
  const long nslots;
  helib::SecKey secretKey;
  const helib::PubKey& publicKey;

  static std::vector<long> getDefaultGens()
  {
    return std::vector<long>{3979, 3095, 3760};
  };

  static std::vector<long> getDefaultOrds()
  {
    return std::vector<long>{6, 2, -8};
  };

  static NTL::Vec<long> getDefaultMvec()
  {
    NTL::Vec<long> defaultMvec;
    defaultMvec.SetLength(3);
    defaultMvec[0] = 7;
    defaultMvec[1] = 3;
    defaultMvec[2] = 221;
    return defaultMvec;
  };

  static void validateMvec(const NTL::Vec<long>& mvec)
  {
    long nfactors = mvec.length();
    for (long i = 0; i < nfactors; i++)
      for (long j = i + 1; j < nfactors; j++)
        if (NTL::GCD(mvec[i], mvec[j]) != 1)
          throw std::invalid_argument(
              "mvec must include the prime-power factorization of m");
  };

  static NTL::Vec<long> calculateMvec(const NTL::Vec<long>& userMvec)
  {
    NTL::Vec<long> mvec;
    if (helib::lsize(userMvec) >= 1)
      mvec = userMvec;
    else
      mvec = getDefaultMvec();
    validateMvec(mvec);
    return mvec;
  };

  static helib::Context& prepareContext(helib::Context& context,
                                        const long L,
                                        const long c)
  {
    helib::buildModChain(context, L, c);
    if (!helib_test::noPrint) {
      context.zMStar.printout(); // print structure of Zm* /(p) to std::cout
      std::cout << std::endl;
    }
    return context;
  };

  static helib::SecKey& prepareSecretKey(helib::SecKey& secretKey, const long w)
  {
    secretKey.GenSecKey(w); // A Hamming-weight-w secret key
    helib::addSome1DMatrices(
        secretKey); // compute key-switching matrices that we need
    helib::addFrbMatrices(
        secretKey); // compute key-switching matrices that we need
    return secretKey;
  }

  GTestThinEvalMap() :
      p(GetParam().p),
      r(GetParam().r),
      c(GetParam().c),
      k(GetParam().k),
      w(64),
      L(GetParam().L),
      s(GetParam().s),
      seed(
          (NTL::SetSeed(NTL::conv<NTL::ZZ>(GetParam().seed)), GetParam().seed)),
      mvec(calculateMvec(GetParam().mvec)),
      gens(helib::lsize(GetParam().mvec) >= 1 ? GetParam().gens
                                              : getDefaultGens()),
      ords(helib::lsize(GetParam().mvec) >= 1 ? GetParam().ords
                                              : getDefaultOrds()),
      m(helib::computeProd(mvec)),
      nthreads((NTL::SetNumThreads(GetParam().nthreads), GetParam().nthreads)),
      useCache(GetParam().useCache),
      context((helib::setTimersOn(),
               helib::setDryRun(
                   false), // Need to get a "real context" to test ThinEvalMap
               m),
              p,
              r,
              gens,
              ords),
      d(prepareContext(context, L, c).zMStar.getOrdP()),
      phim(context.zMStar.getPhiM()),
      nslots(phim / d),
      secretKey(
          (helib::setDryRun(
               helib_test::dry), // Now we can set the dry-run flag if desired
           context)),
      publicKey(prepareSecretKey(secretKey, w)){};

  virtual void TearDown() override
  {
    if (!helib_test::noPrint) {
      std::cout << "\n*********\n";
      helib::printAllTimers();
      std::cout << std::endl;
    }
    helib::cleanupDebugGlobals();
  }
};

TEST_P(GTestThinEvalMap, thinEvalMapIsCorrect)
{
  // GG defines the plaintext space Z_p[X]/GG(X)
  NTL::ZZX GG;
  GG = context.alMod.getFactorsOverZZ()[0];
  helib::EncryptedArray ea(context, GG);

  NTL::zz_p::init(context.alMod.getPPowR());

  NTL::Vec<NTL::zz_p> val0(NTL::INIT_SIZE, nslots);
  for (auto& x : val0)
    random(x);

  std::vector<NTL::ZZX> val1;
  val1.resize(nslots);
  for (long i = 0; i < nslots; i++) {
    val1[i] = NTL::conv<NTL::ZZX>(NTL::conv<NTL::ZZ>(rep(val0[i])));
  }

  helib::Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, val1);

  helib::resetAllTimers();
  HELIB_NTIMER_START(ALL);

  // Compute homomorphically the transformation that takes the
  // coefficients packed in the slots and produces the polynomial
  // corresponding to cube

  if (!helib_test::noPrint)
    helib::CheckCtxt(ctxt, "init");

  if (!helib_test::noPrint)
    std::cout << "build ThinEvalMap\n";

  helib::ThinEvalMap map(ea,
                         /*minimal=*/false,
                         mvec,
                         /*invert=*/false,
                         /*build_cache=*/false);
  // compute the transformation to apply

  if (!helib_test::noPrint)
    std::cout << "apply ThinEvalMap\n";
  if (useCache)
    map.upgrade();
  map.apply(ctxt); // apply the transformation to ctxt
  if (!helib_test::noPrint)
    helib::CheckCtxt(ctxt, "ThinEvalMap");
  if (!helib_test::noPrint)
    std::cout << "check results\n";

  if (!helib_test::noPrint)
    std::cout << "build ThinEvalMap\n";
  helib::ThinEvalMap imap(ea,
                          /*minimal=*/false,
                          mvec,
                          /*invert=*/true,
                          /*build_cache=*/false);
  // compute the transformation to apply
  if (!helib_test::noPrint)
    std::cout << "apply ThinEvalMap\n";
  if (useCache)
    imap.upgrade();
  imap.apply(ctxt); // apply the transformation to ctxt
  if (!helib_test::noPrint) {
    helib::CheckCtxt(ctxt, "ThinEvalMap");
    std::cout << "check results\n";
  }

#if 1

  /* create dirty version of ctxt */
  NTL::Vec<NTL::zz_pX> dirty_val0;
  dirty_val0.SetLength(nslots);
  for (long i = 0; i < nslots; i++) {
    random(dirty_val0[i], d);
    SetCoeff(dirty_val0[i], 0, val0[i]);
  }

  std::vector<NTL::ZZX> dirty_val1;
  dirty_val1.resize(nslots);
  for (long i = 0; i < nslots; i++) {
    dirty_val1[i] = NTL::conv<NTL::ZZX>(dirty_val0[i]);
  }

  helib::Ctxt dirty_ctxt(publicKey);
  ea.encrypt(dirty_ctxt, publicKey, dirty_val1);

  helib::EvalMap dirty_map(ea,
                           /*minimal=*/false,
                           mvec,
                           /*invert=*/false,
                           /*build_cache=*/false);

  dirty_map.apply(dirty_ctxt);
  imap.apply(dirty_ctxt);
#endif

  std::vector<NTL::ZZX> val2;
  ea.decrypt(ctxt, secretKey, val2);

  EXPECT_EQ(val1, val2);

#if 1
  std::vector<NTL::ZZX> dirty_val2;
  ea.decrypt(dirty_ctxt, secretKey, dirty_val2);

  EXPECT_EQ(val1, dirty_val2);
#endif

  HELIB_NTIMER_STOP(ALL);
}

// clang-format off
INSTANTIATE_TEST_SUITE_P(variousParameters, GTestThinEvalMap, ::testing::Values(
    //SLOW
    Parameters(2, 1, 2, 80, 300, 0, 0,
               helib::convert<NTL::Vec<long>, std::vector<long>>(std::vector<long>{7, 3, 221}),
               std::vector<long>{3979, 3095, 3760},
               std::vector<long>{6, 2, -8},
               1, 0)
    //FAST
    //Parameters(2, 1, 2, 80, 300, 0, 0,
    //           helib::convert<NTL::Vec<long>, std::vector<long>>(std::vector<long>{3, 35}),
    //           std::vector<long>{71, 76},
    //           std::vector<long>{2, 2},
    //           1, 0)
    ));
// clang-format on
} // namespace
