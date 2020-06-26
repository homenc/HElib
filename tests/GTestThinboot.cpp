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
#include <helib/helib.h>
#include <helib/matmul.h>
#include <helib/debugging.h>
#include <helib/fhe_stats.h>

#include "gtest/gtest.h"
#include "test_common.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <NTL/ZZ.h>

namespace {

/********************************************************************
 ********************************************************************/
struct Parameters
{
  long p; // plaintext base
  long r; // exponent
  // p^r is the plaintext-space modulus
  long c;        // number of columns in the key-switching matrices
  long bits;     // # of bits in the modulus chain
  long skHwt;    // Hamming weight of recryption secret key
  long nthreads; // number of threads
  long seed;     // random number seed
  int useCache;  // 0: zzX cache, 1: DCRT cache
  int force_bsgs;
  int force_hoist;
  // int disable_intFactor // fhe_disable_intFactor (disabled by Victor)
  int chen_han;
  bool debug; // generate debugging output
  int scale;  // scale parameter
  NTL::Vec<long> global_gens;
  NTL::Vec<long> global_ords;
  NTL::Vec<long> global_mvec;
  const int c_m; // = 100;
  const long iter;
  const std::string v_values_name;

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
             long iter,
             const std::string& v_values_name) :
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
      iter(iter),
      v_values_name(v_values_name)
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
              << "c_m" << params.c_m << ","
              << "iter" << params.iter << ","
              << "v_values_name" << params.v_values_name << "}";
  }

  //
  //
  //  if (seed)
  //    SetSeed(ZZ(seed));
  //
  //  SetNumThreads(nthreads);
  //
  //  TestIt(p,r,bits,c,t,useCache);
};

class GTestThinboot : public ::testing::TestWithParam<Parameters>
{
private:
  void preContextSetup()
  {

    if (!helib_test::noPrint) {
      std::cout << "*** GTestThinboot";
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
    if (!helib_test::noPrint) {
      helib::fhe_stats = true;
    }
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
  const long iter;

  const long m, phim;
  double time;
  helib::Context context;

  std::string v_values_name;

  GTestThinboot() :
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
      iter(GetParam().iter),
      m(helib::computeProd(mvec)),
      phim((checkPM(p, m), helib::phi_N(m))),
      time(0),
      context((preContextSetup(), m), p, r, gens, ords),
      v_values_name(GetParam().v_values_name)

  {
    postContextSetup();
  }

  void TearDown() override
  {
    if (helib_test::verbose) {
      helib::printAllTimers();

#if (defined(__unix__) || defined(__unix) || defined(unix))
      struct rusage rusage;
      getrusage(RUSAGE_SELF, &rusage);
      if (!helib_test::noPrint)
        std::cout << "  rusage.ru_maxrss=" << rusage.ru_maxrss << std::endl;
#endif

      if (helib::fhe_stats) {
        helib::print_stats(std::cout);
        std::cout << "mvec=" << mvec << ", ";
        print_anderson_darling();
        dump_v_values();
        std::cout << std::endl;
      }
    }

    cleanupBootstrappingGlobals();
    helib::cleanupDebugGlobals();
  }

  void dump_v_values()
  {
    const std::vector<double>* v_values = helib::fetch_saved_values("v_values");
    if (v_values && v_values_name != "") {
      // write v_values to a file

      std::cerr << "writing v_values to " << v_values_name << std::endl;

      std::ofstream F;
      F.open(v_values_name.c_str());
      for (long i : helib::range(v_values->size()))
        F << (*v_values)[i] << std::endl;
    }
  }

  void anderson_darling(const std::vector<double>& X, double& AD, double& p_val)
  {
    long N = X.size();

    if (N < 2) {
      AD = 0;
      p_val = 1;
      return;
    }

    std::vector<double> Y(X);
    sort(Y.begin(), Y.end());

    // compute the sample mean
    double SM = 0;
    for (long i : helib::range(N))
      SM += Y[i];
    SM /= N;

    // compute the sample variance
    double SV = 0;
    for (long i : helib::range(N))
      SV += (Y[i] - SM) * (Y[i] - SM);
    SV /= (N - 1);

    // replace Y[i] by CDF of Y[i]
    for (long i : helib::range(N))
      Y[i] = 0.5 * (1 + erf((Y[i] - SM) / sqrt(2 * SV)));

    double S = 0;
    for (long i : helib::range(N)) {
      S += (2 * i + 1) * (log(Y[i]) + log1p(-Y[N - 1 - i]));
    }
    AD = -N - S / N;

    AD *= (1 + 0.75 / N + 2.25 / N / N);
    // This adjustment and the p-values below come from:
    // R.B. D'Augostino and M.A. Stephens, Eds., 1986,
    // Goodness-of-Fit Techniques, Marcel Dekker.

    if (AD >= 0.6)
      p_val = exp(1.2937 - 5.709 * (AD) + 0.0186 * helib::fsquare(AD));
    else if (AD > 0.34)
      p_val = exp(0.9177 - 4.279 * (AD)-1.38 * helib::fsquare(AD));
    else if (AD > 0.2)
      p_val = 1 - exp(-8.318 + 42.796 * (AD)-59.938 * helib::fsquare(AD));
    else
      p_val = 1 - exp(-13.436 + 101.14 * (AD)-223.73 * helib::fsquare(AD));
  }

  void print_anderson_darling()
  {
    const std::vector<double>* v_values = helib::fetch_saved_values("v_values");
    if (v_values) {
      double AD, p_val;
      anderson_darling(*v_values, AD, p_val);
      std::cout << "AD=" << AD << ", p_val=" << p_val << std::endl;
    }
  }
};

TEST_P(GTestThinboot, correctlyPerformsThinboot)
{
  helib::buildModChain(context,
                       bits,
                       c,
                       /*willBeBootstrappable=*/true,
                       /*skHwt=*/skHwt,
                       /*resolution=*/3,
                       /*bitsInSpecialPrimes=*/helib_test::special_bits);

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

  context.makeBootstrappable(mvec, /*t=*/skHwt, useCache, /*alsoThick=*/false);
  // save time...disable some fat boot precomputation

  time += NTL::GetTime();

  // if (skHwt>0) context.rcData.skHwt = skHwt;
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

  for (long numkey = 0; numkey < iter; numkey++) { // test with 3 keys
    if (helib::fhe_stats && numkey > 0 && numkey % 100 == 0) {
      helib::print_stats(std::cout);

      std::cout << "mvec=" << mvec << ", ";
      print_anderson_darling();
      dump_v_values();
    }

    if (!helib_test::noPrint)
      std::cerr << "*********** iter=" << (numkey + 1) << std::endl;

    time = -NTL::GetTime();
    if (!helib_test::noPrint)
      std::cout << "Generating keys, " << std::flush;
    helib::SecKey secretKey(context);
    secretKey.GenSecKey(skHwt); // A Hamming-weight-64 secret key
    helib::addSome1DMatrices(
        secretKey); // compute key-switching matrices that we need
    helib::addFrbMatrices(secretKey);
    if (!helib_test::noPrint)
      std::cout << "computing key-dependent tables..." << std::flush;
    secretKey.genRecryptData();
    time += NTL::GetTime();
    if (!helib_test::noPrint)
      std::cout << " done in " << time << " seconds\n";

    const helib::PubKey publicKey = secretKey;

    long d = context.zMStar.getOrdP();
    long phim = context.zMStar.getPhiM();
    long nslots = phim / d;

    // GG defines the plaintext space Z_p[X]/GG(X)
    NTL::ZZX GG;
    GG = context.alMod.getFactorsOverZZ()[0];
    std::shared_ptr<helib::EncryptedArray> ea(
        std::make_shared<helib::EncryptedArray>(context, GG));

    helib::setupDebugGlobals(&secretKey, ea);

    NTL::zz_p::init(p2r);
    NTL::Vec<NTL::zz_p> val0(NTL::INIT_SIZE, nslots);
    for (auto& x : val0)
      random(x);

    std::vector<NTL::ZZX> val1;
    val1.resize(nslots);
    for (long i = 0; i < nslots; i++) {
      val1[i] = NTL::conv<NTL::ZZX>(NTL::conv<NTL::ZZ>(rep(val0[i])));
    }

    helib::Ctxt c1(publicKey);
    ea->encrypt(c1, publicKey, val1);

    // Make some noise!
    // This ensures that we do not start the sequence of squarings
    // from a "fresh" ciphertext (which may not be representative).
    ea->rotate(c1, 1);
    ea->rotate(c1, -1);

    helib::Ctxt c2(c1);

    if (!helib_test::noPrint)
      helib::CheckCtxt(c2, "before squarings");

    long sqr_count = -1;
    helib::Ctxt next_c2(c2);
    do {
      c2 = next_c2;
      next_c2.multiplyBy(next_c2);
      sqr_count++;
    } while (next_c2.bitCapacity() >= 100);

    if (!helib_test::noPrint) {
      // compute minimal capacity before bootstrapping (rawModSwitch)
      long e = context.rcData.e;
      long q = NTL::power_long(p, e) + 1;
      double Bnd = context.boundForRecryption();
      double mfac = context.zMStar.getNormBnd();
      double min_bit_cap =
          log(mfac * q / (p2r * Bnd * HELIB_MIN_CAP_FRAC)) / log(2.0);

      std::cout << "min_bit_cap=" << min_bit_cap << std::endl;

      std::cout << "log2(modSwitchAddedNoiseBound)="
                << log(c1.modSwitchAddedNoiseBound()) / log(2.0) << std::endl;

      std::cout << "sqr_count=" << sqr_count << std::endl;
      if (sqr_count > 0) {
        std::cout << "BITS-PER-LEVEL: "
                  << ((c1.bitCapacity() - c2.bitCapacity()) / double(sqr_count))
                  << std::endl;
      }
      helib::CheckCtxt(c2, "before recryption");
    }

    helib::resetAllTimers();

    {
      HELIB_NTIMER_START(AAA_thinRecrypt);

      publicKey.thinReCrypt(c2);
    }

    if (!helib_test::noPrint)
      helib::CheckCtxt(c2, "after recryption");

    for (auto& x : val0) {
      for (long i = 0; i < sqr_count; ++i)
        x = x * x;
    }

    for (long i = 0; i < nslots; i++) {
      val1[i] = NTL::conv<NTL::ZZX>(NTL::conv<NTL::ZZ>(rep(val0[i])));
    }

    std::vector<NTL::ZZX> val2;
    ea->decrypt(c2, secretKey, val2);

    EXPECT_EQ(val1, val2);
  }
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
INSTANTIATE_TEST_SUITE_P(typicalParameters, GTestThinboot, ::testing::Values(
    //SLOW
    Parameters( 2, 1, 3, 600, 64, 1, 0, 1, 100, 0, 0, 0, 0, 0, {1026,  249}, {30, -2}, {  31, 41}, 1, ""),
    Parameters(17, 1, 3, 600, 64, 1, 0, 1, 100, 0, 0, 0, 0, 0, { 556, 1037}, { 6,  4}, {7, 5, 37}, 1, "")
    //FAST
    //Parameters( 2, 1, 3, 600, 64, 1, 0, 1, 100, 0, 0, 0, 0, 0, {1026,  249}, {30, -2}, {  31, 41}, 1, "")
    ));
// clang-format on
} // namespace
