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

/* Test_bootstrapping.cpp - Testing the recryption procedure */

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <NTL/ZZ.h>
#include <NTL/fileio.h>
#include <NTL/BasicThreadPool.h>
NTL_CLIENT
#include <cassert>
#include <helib/EncryptedArray.h>
#include <helib/EvalMap.h>
#include <helib/powerful.h>
#include <helib/matmul.h>
#include <helib/debugging.h>
#include <helib/fhe_stats.h>
#include <helib/ArgMap.h>

using namespace helib;

static bool noPrint = false;
static bool dry = false; // a dry-run flag
static bool debug = 0;   // a debug flag
static int scale = 0;

extern long printFlag;


#define OUTER_REP (1)
#define INNER_REP (1)


static Vec<long> global_mvec, global_gens, global_ords;
static int c_m = 100;


void TestIt(long p, long r, long L, long c, long skHwt, int build_cache=0)
{
  Vec<long> mvec;
  vector<long> gens;
  vector<long> ords;

  long m, phim;

    mvec = global_mvec;
    convert(gens, global_gens);
    convert(ords, global_ords);

    m = computeProd(mvec);
    phim = phi_N(m);
    assertTrue(GCD(p, m) == 1, "GCD(p, m) == 1");

  if (!noPrint) fhe_stats = true;

  if (!noPrint) {
    cout << "*** TestIt";
    if (isDryRun()) cout << " (dry run)";
    cout << ": p=" << p
	 << ", r=" << r
	 << ", L=" << L
	 << ", t=" << skHwt
	 << ", c=" << c
	 << ", m=" << m
	 << " (=" << mvec << "), gens="<<gens<<", ords="<<ords
	 << endl;
    cout << "Computing key-independent tables..." << std::flush;
  }
  setTimersOn();
  setDryRun(false); // Need to get a "real context" to test bootstrapping

  double t = -GetTime();
  Context context(m, p, r, gens, ords);
  if (scale) {
    context.scale = scale;
  }


  context.zMStar.set_cM(c_m/100.0);
  buildModChain(context, L, c, /*willBeBootstrappable=*/true, /*t=*/skHwt);

  if (!noPrint) {
    std::cout << "security=" << context.securityLevel()<<endl;
    std::cout << "# small primes = " << context.smallPrimes.card() << "\n";
    std::cout << "# ctxt primes = " << context.ctxtPrimes.card() << "\n";
    std::cout << "# bits in ctxt primes = "
         << long(context.logOfProduct(context.ctxtPrimes)/log(2.0) + 0.5) << "\n";
    std::cout << "# special primes = " << context.specialPrimes.card() << "\n";
    std::cout << "# bits in special primes = "
         << long(context.logOfProduct(context.specialPrimes)/log(2.0) + 0.5) << "\n";
    std::cout << "scale=" << context.scale<<endl;
  }



  context.makeBootstrappable(mvec, /*t=*/skHwt, build_cache);
  t += GetTime();

  if (!noPrint) {
    cout << " done in "<<t<<" seconds\n";
    cout << "  e="    << context.rcData.e
	 << ", e'="   << context.rcData.ePrime
	 << ", t="    << context.rcData.skHwt
	 << "\n  ";
    context.zMStar.printout();
  }
  setDryRun(dry); // Now we can set the dry-run flag if desired

  long p2r = context.alMod.getPPowR();

  for (long numkey=0; numkey<OUTER_REP; numkey++) { // test with 3 keys

  t = -GetTime();
  if (!noPrint) cout << "Generating keys, " << std::flush;
  SecKey secretKey(context);
  PubKey& publicKey = secretKey;
  secretKey.GenSecKey(skHwt);      // A +-1/0 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey);
  if (!noPrint) cout << "computing key-dependent tables..." << std::flush;
  secretKey.genRecryptData();
  t += GetTime();
  if (!noPrint) cout << " done in "<<t<<" seconds\n";

  zz_p::init(p2r);
  zz_pX poly_p = random_zz_pX(context.zMStar.getPhiM());
  zzX poly_p1 = balanced_zzX(poly_p);
  ZZX ptxt_poly = convert<ZZX>(poly_p1);
  ZZX ptxt_poly1;
  PolyRed(ptxt_poly1, ptxt_poly, p2r, true);
  // this is the format produced by decryption

#ifdef HELIB_DEBUG
      dbgEa = context.ea;
      dbgKey = &secretKey;
#endif

  if (debug) {
    dbgKey = &secretKey; // debugging key
  }

  ZZX poly2;
  Ctxt c1(publicKey);

  secretKey.Encrypt(c1,ptxt_poly,p2r);


  resetAllTimers();
  for (long num=0; num<INNER_REP; num++) { // multiple tests with same key
    publicKey.reCrypt(c1);
    secretKey.Decrypt(poly2,c1);

    if (ptxt_poly1 == poly2)
      cout << "GOOD\n";
    else
      cout << "BAD\n";
  }
  }
  if (!noPrint) printAllTimers();
#if (defined(__unix__) || defined(__unix) || defined(unix))
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
    if (!noPrint) cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif
  if (fhe_stats) print_stats(cout);
}

/********************************************************************
 ********************************************************************/

//extern long fhe_disable_intFactor;
// extern long fhe_force_chen_han;

int main(int argc, char *argv[])
{
  ArgMap amap;

  long p=2;
  long r=1;
  long c=3;
  long L=600;
  long N=0;
  long t=64;
  long nthreads=1;

  long seed=0;
  long useCache=1;

  amap.arg("p", p, "plaintext base");

  amap.arg("r", r,  "exponent");
  amap.note("p^r is the plaintext-space modulus");

  amap.arg("c", c, "number of columns in the key-switching matrices");
  amap.arg("L", L, "# of bits in the modulus chain");
  amap.arg("N", N, "lower-bound on phi(m)");
  amap.arg("t", t, "Hamming weight of recryption secret key", "heuristic");
  amap.arg("dry", dry, "dry=1 for a dry-run");
  amap.arg("nthreads", nthreads, "number of threads");
  amap.arg("seed", seed, "random number seed");
  amap.arg("noPrint", noPrint, "suppress printouts");
  amap.arg("useCache", useCache, "0: zzX cache, 1: DCRT cache");


  amap.arg("force_bsgs", fhe_test_force_bsgs);
  amap.arg("force_hoist", fhe_test_force_hoist);


  //  amap.arg("disable_intFactor", fhe_disable_intFactor);
  amap.arg("chen_han", fhe_force_chen_han);

  amap.arg("debug", debug, "generate debugging output");
  amap.arg("scale", scale, "scale parameter");


  amap.arg("gens", global_gens);
  amap.arg("ords", global_ords);
  amap.arg("mvec", global_mvec);
  amap.arg("c_m", c_m);

  amap.parse(argc, argv);

  if (global_gens.length() == 0 || global_ords.length() == 0 || global_mvec.length() == 0)
    Error("gens, ords, and mvec must be initialized");

  if (seed)
    SetSeed(ZZ(seed));

  SetNumThreads(nthreads);


  TestIt(p,r,L,c,t,useCache);

  return 0;
}
