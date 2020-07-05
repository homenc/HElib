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

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <NTL/BasicThreadPool.h>
#include <helib/FHE.h>
#include <helib/EncryptedArray.h>
#include <helib/matmul.h>
#include <helib/debugging.h>
#include <helib/fhe_stats.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

/********************************************************************
 ********************************************************************/
struct Parameters
{
  const long p; // plaintext base
  const long r; // exponent
  // p^r is the plaintext-space modulus
  const long c;        // number of columns in the key-switching matrices
  const long bits;     // # of bits in the modulus chain
  const long skHwt;    // Hamming weight of recryption secret key
  const long nthreads; // number of threads
  const long seed;     // random number seed
  const int useCache;  // 0: zzX cache, 1: DCRT cache
  const int force_bsgs;
  const int force_hoist;
  // int disable_intFactor // fhe_disable_intFactor (disabled by Victor)
  const int chen_han;
  const bool debug; // generate debugging output
  const int scale;  // scale parameter
  const NTL::Vec<long> global_gens;
  const NTL::Vec<long> global_ords;
  const NTL::Vec<long> global_mvec;
  const int c_m; // = 100;
  const long outer_rep;
  const long inner_rep;

  Parameters(long p,
             long r,
             long c,
             long bits,
             long skHwt,
             long nthreads,
             long seed,
             long useCache,
             int c_m,
             int force_bsgs,
             int force_hoist,
             int chen_han,
             bool debug,
             int scale,
             const std::vector<long>& global_gens,
             const std::vector<long>& global_ords,
             const std::vector<long>& global_mvec,
             long outer_rep,
             long inner_rep) :
      p(p),
      r(r),
      c(c),
      bits(bits),
      skHwt(skHwt),
      nthreads(nthreads),
      seed(seed),
      useCache(useCache),
      force_bsgs(force_bsgs),
      force_hoist(force_hoist),
      chen_han(chen_han),
      debug(debug),
      scale(scale),
      global_gens(helib::convert<NTL::Vec<long>>(global_gens)),
      global_ords(helib::convert<NTL::Vec<long>>(global_ords)),
      global_mvec(helib::convert<NTL::Vec<long>>(global_mvec)),
      c_m(c_m),
      outer_rep(outer_rep),
      inner_rep(inner_rep)
  {
    if (global_gens.empty() || global_ords.empty() || global_mvec.empty())
      throw helib::LogicError("gens, ords, and mvec must be non-empty");
  };

  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "p" << params.p << ","
              << "r" << params.r << ","
              << "c" << params.c << ","
              << "bits" << params.bits << ","
              << "skHwt" << params.skHwt << ","
              << "nthreads" << params.nthreads << ","
              << "seed" << params.seed << ","
              << "useCache" << params.useCache << ","
              << "force_bsgs" << params.force_bsgs << ","
              << "force_hoist" << params.force_hoist << ","
              << "chen_han" << params.chen_han << ","
              << "debug" << params.debug << ","
              << "scale" << params.scale << ","
              << "global_gens" << params.global_gens << ","
              << "global_ords" << params.global_ords << ","
              << "global_mvec" << params.global_mvec << ","
              << "c_m" << params.c_m << "}";
  }
};

class GTestFatboot : public ::testing::TestWithParam<Parameters>
{
private:
  void preContextSetup()
  {

    if (!helib_test::noPrint)
      helib::fhe_stats = true;

    if (!helib_test::noPrint) {
      std::cout << "*** GTestFatboot";
      if (helib::isDryRun())
        std::cout << " (dry run)";
      std::cout << ": p=" << p << ", r=" << r << ", bits=" << bits
                << ", t=" << skHwt << ", c=" << c << ", m=" << m
                << " mvec=" << mvec << ", gens=" << helib::vecToStr(gens)
                << ", ords=" << helib::vecToStr(ords) << std::endl;
      std::cout << "Computing key-independent tables..." << std::flush;
    }
    helib::setTimersOn();
    helib::setDryRun(
        false); // Need to get a "real context" to test bootstrapping
    time = -NTL::GetTime();
  }

  void postContextSetup()
  {
    if (scale) {
      context.scale = scale;
    }

    context.zMStar.set_cM(c_m / 100.0);
  }

  static void setGlobals(int force_bsgs, int force_hoist, int chen_han)
  {
    helib::fhe_test_force_bsgs = force_bsgs;
    helib::fhe_test_force_hoist = force_hoist;
    helib::fhe_force_chen_han = chen_han;
  };

  void cleanupBootstrappingGlobals()
  {
    helib::fhe_test_force_bsgs = old_fhe_test_force_bsgs;
    helib::fhe_test_force_hoist = old_fhe_test_force_hoist;
    helib::fhe_force_chen_han = old_fhe_force_chen_han;
  }

  static void setSeedIfNeeded(const long seed)
  {
    if (seed)
      SetSeed(NTL::ZZ(seed));
  };

  static void checkPM(const long p, const long m)
  {
    helib::assertTrue(NTL::GCD(p, m) == 1, "GCD(p, m) == 1");
  }

protected:
  const int old_fhe_test_force_bsgs;
  const int old_fhe_test_force_hoist;
  const int old_fhe_force_chen_han;
  const long p;
  const long r;
  const long c;
  const long bits;
  const long skHwt;
  const long nthreads;
  const long seed;
  const long useCache;
  const int force_bsgs;
  const int force_hoist;
  const int chen_han;
  const bool debug;
  const int scale;
  const NTL::Vec<long> mvec;
  const std::vector<long> gens;
  const std::vector<long> ords;
  const int c_m; // = 100;
  const long outer_rep;
  const long inner_rep;

  const long m, phim;
  double time;
  helib::Context context;

  GTestFatboot() :
      old_fhe_test_force_bsgs(helib::fhe_test_force_bsgs),
      old_fhe_test_force_hoist(helib::fhe_test_force_hoist),
      old_fhe_force_chen_han(helib::fhe_force_chen_han),
      p((setGlobals(GetParam().force_bsgs,
                    GetParam().force_hoist,
                    GetParam().chen_han),
         GetParam().p)),
      r(GetParam().r),
      c(GetParam().c),
      bits(GetParam().bits),
      skHwt(GetParam().skHwt),
      nthreads((NTL::SetNumThreads(GetParam().nthreads), GetParam().nthreads)),
      seed((setSeedIfNeeded(GetParam().seed), GetParam().seed)),
      useCache(GetParam().useCache),
      force_bsgs(GetParam().force_bsgs),
      force_hoist(GetParam().force_hoist),
      chen_han(GetParam().chen_han),
      debug(GetParam().debug),
      scale(GetParam().scale),
      mvec(GetParam().global_mvec),
      gens(helib::convert<std::vector<long>>(GetParam().global_gens)),
      ords(helib::convert<std::vector<long>>(GetParam().global_ords)),
      c_m(GetParam().c_m),
      outer_rep(GetParam().outer_rep),
      inner_rep(GetParam().inner_rep),
      m(helib::computeProd(mvec)),
      phim((checkPM(p, m), helib::phi_N(m))),
      time(0),
      context((preContextSetup(), m), p, r, gens, ords)
  {
    postContextSetup();
  }

  void TearDown() override
  {
    if (!helib_test::noPrint) {
      helib::printAllTimers();
    }
    if (helib::fhe_stats)
      helib::print_stats(std::cout);

    cleanupBootstrappingGlobals();
    helib::cleanupDebugGlobals();
  }
};

TEST_P(GTestFatboot, correctlyPerformsFatboot)
{
  helib::buildModChain(context,
                       bits,
                       c,
                       /*willBeBootstrappable=*/true,
                       /*t=*/skHwt);

  if (!helib_test::noPrint) {
    std::cout << "security=" << context.securityLevel() << std::endl;
    std::cout << "# small primes = " << context.smallPrimes.card() << std::endl;
    std::cout << "# ctxt primes = " << context.ctxtPrimes.card() << std::endl;
    std::cout << "# bits in ctxt primes = "
              << long(context.logOfProduct(context.ctxtPrimes) / log(2.0) + 0.5)
              << std::endl;
    std::cout << "# special primes = " << context.specialPrimes.card()
              << std::endl;
    std::cout << "# bits in special primes = "
              << long(context.logOfProduct(context.specialPrimes) / log(2.0) +
                      0.5)
              << std::endl;
    std::cout << "scale=" << context.scale << std::endl;
  }

  context.makeBootstrappable(mvec, /*t=*/skHwt, useCache);
  time += NTL::GetTime();

  if (!helib_test::noPrint) {
    std::cout << " done in " << time << " seconds" << std::endl;
    std::cout << "  e=" << context.rcData.e << ", e'=" << context.rcData.ePrime
              << ", t=" << context.rcData.skHwt << std::endl
              << "  ";
    context.zMStar.printout();
  }
  helib::setDryRun(
      helib_test::dry); // Now we can set the dry-run flag if desired

  long p2r = context.alMod.getPPowR();

  for (long numkey = 0; numkey < outer_rep; numkey++) { // test with 3 keys

    time = -NTL::GetTime();
    if (!helib_test::noPrint)
      std::cout << "Generating keys, " << std::flush;
    helib::SecKey secretKey(context);
    helib::PubKey& publicKey = secretKey;
    secretKey.GenSecKey(skHwt); // A +-1/0 secret key
    helib::addSome1DMatrices(
        secretKey); // compute key-switching matrices that we need
    helib::addFrbMatrices(secretKey);
    if (!helib_test::noPrint)
      std::cout << "computing key-dependent tables..." << std::flush;
    secretKey.genRecryptData();
    time += NTL::GetTime();
    if (!helib_test::noPrint)
      std::cout << " done in " << time << " seconds\n";

    NTL::zz_p::init(p2r);
    NTL::zz_pX poly_p = NTL::random_zz_pX(context.zMStar.getPhiM());
    helib::zzX poly_p1 = helib::balanced_zzX(poly_p);
    NTL::ZZX ptxt_poly = helib::convert<NTL::ZZX>(poly_p1);
    NTL::ZZX ptxt_poly1;
    helib::PolyRed(ptxt_poly1, ptxt_poly, p2r, true);
    // this is the format produced by decryption

    helib::setupDebugGlobals(&secretKey, context.ea);

    NTL::ZZX poly2;
    helib::Ctxt c1(publicKey);

    secretKey.Encrypt(c1, ptxt_poly, p2r);

    helib::resetAllTimers();
    for (long num = 0; num < inner_rep; num++) { // multiple tests with same key
      publicKey.reCrypt(c1);
      secretKey.Decrypt(poly2, c1);

      EXPECT_EQ(ptxt_poly1, poly2);
    }
  }
  if (!helib_test::noPrint)
    helib::printAllTimers();
#if (defined(__unix__) || defined(__unix) || defined(unix))
  struct rusage rusage;
  getrusage(RUSAGE_SELF, &rusage);
  if (!helib_test::noPrint)
    std::cout << "  rusage.ru_maxrss=" << rusage.ru_maxrss << std::endl;
#endif
  if (helib::fhe_stats)
    helib::print_stats(std::cout);
}

// LEGACY TEST DEFAULT PARAMETERS:
// long p=2;
// long r=1;
// long c=3;
// long bits=600;
// long t=64;
// long nthreads=1;
// long seed=0;
// long useCache=1;
// int c_m = 100;

// clang-format off
INSTANTIATE_TEST_SUITE_P(typicalParameters, GTestFatboot, ::testing::Values(
    //SLOW
    Parameters( 2, 1, 3, 600, 64, 1, 0, 1, 100, 0, 0, 0, 0, 0, {1026,  249}, {30, -2}, {  31, 41}, 1, 1),
    Parameters(17, 1, 3, 600, 64, 1, 0, 1, 100, 0, 0, 0, 0, 0, { 556, 1037}, { 6,  4}, {7, 5, 37}, 1, 1)
    //FAST
    // Parameters( 2, 1, 3, 600, 64, 1, 0, 1, 100, 0, 0, 0, 0, 0, {1026,  249}, {30, -2}, {  31, 41}, 1, 1)
    ));
// clang-format on
} // namespace
