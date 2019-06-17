/* Copyright (C) 2012-2018 IBM Corp.
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
#include <NTL/BasicThreadPool.h>
#include "FHE.h"
#include "EncryptedArray.h"
#include "matmul.h"
#include "debugging.h"
#include "fhe_stats.h"

NTL_CLIENT

static bool noPrint = false;
static bool dry = false; // a dry-run flag
static bool debug = 0;   // a debug flag
static int scale = 0;


static long OUTER_REP = 1;


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
    helib::assertTrue(GCD(p, m) == 1, "GCD(p, m) == 1");

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
  FHEcontext context(m, p, r, gens, ords);
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

  context.makeBootstrappable(mvec,/*t=*/skHwt,build_cache,/*alsoThick=*/false);
  // save time...disable some fat boot precomputation

  t += GetTime();

  //if (skHwt>0) context.rcData.skHwt = skHwt;
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

  if (!noPrint) fhe_stats = true;

  for (long numkey=0; numkey<OUTER_REP; numkey++) { // test with 3 keys
  if (fhe_stats && numkey > 0 && numkey%100 == 0) 
    print_stats(cout);

  if (!noPrint) cerr << "*********** iter=" << (numkey+1) << "\n";

  t = -GetTime();
  if (!noPrint) cout << "Generating keys, " << std::flush;
  FHESecKey secretKey(context);
  secretKey.GenSecKey(skHwt);      // A Hamming-weight-64 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey);
  if (!noPrint) cout << "computing key-dependent tables..." << std::flush;
  secretKey.genRecryptData();
  t += GetTime();
  if (!noPrint) cout << " done in "<<t<<" seconds\n";

  FHEPubKey publicKey = secretKey;

  long d = context.zMStar.getOrdP();
  long phim = context.zMStar.getPhiM();
  long nslots = phim/d;

  // GG defines the plaintext space Z_p[X]/GG(X)
  ZZX GG;
  GG = context.alMod.getFactorsOverZZ()[0];
  EncryptedArray ea(context, GG);

  if (debug) {
    dbgKey = &secretKey;
    dbgEa = &ea;
  }

  zz_p::init(p2r);
  Vec<zz_p> val0(INIT_SIZE, nslots);
  for (auto& x: val0)
    random(x);

  vector<ZZX> val1;
  val1.resize(nslots);
  for (long i = 0; i < nslots; i++) {
    val1[i] = conv<ZZX>(conv<ZZ>(rep(val0[i])));
  }


  Ctxt c1(publicKey);
  ea.encrypt(c1, publicKey, val1);

  Ctxt c2(c1);

  if (!noPrint) CheckCtxt(c2, "before squarings");

  long sqr_count = -1;
  Ctxt next_c2(c2);
  do {
    c2 = next_c2;
    next_c2.multiplyBy(next_c2);
    sqr_count++;
  } 
  while (next_c2.bitCapacity() >= 100);

  if (!noPrint) {
    cout << "sqr_count=" << sqr_count << "\n";
    if (sqr_count > 0) {
      cout << "BITS-PER-LEVEL: " 
           << ((c1.bitCapacity()-c2.bitCapacity())/double(sqr_count)) << "\n";
    }
    CheckCtxt(c2, "before recryption");
  }

  resetAllTimers();

  { FHE_NTIMER_START(AAA_thinRecrypt);

  publicKey.thinReCrypt(c2);
 
  }


  if (!noPrint) CheckCtxt(c2, "after recryption");

  for (auto& x: val0) {
    for (long i: range(sqr_count)) x = x*x;
  }

  for (long i = 0; i < nslots; i++) {
    val1[i] = conv<ZZX>(conv<ZZ>(rep(val0[i])));
  }

  vector<ZZX> val2;
  ea.decrypt(c2, secretKey, val2);

  if (val1 == val2) 
    cout << "GOOD\n";
  else
    cout << "BAD\n";
  }

 
  if (!noPrint) printAllTimers();

  if (fhe_stats) print_stats(cout);

}


//extern long fhe_disable_intFactor;
extern long fhe_disable_fat_boot;
extern long fhe_force_chen_han;

/********************************************************************
 ********************************************************************/
int main(int argc, char *argv[]) 
{
  ArgMapping amap;

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
  amap.arg("L", L, "# of levels in the modulus chain");
  amap.arg("t", t, "Hamming weight of recryption secret key");
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

  amap.arg("iter", OUTER_REP);

  amap.parse(argc, argv);

  if (global_gens.length() == 0 || global_ords.length() == 0 || global_mvec.length() == 0)
    Error("gens, ords, and mvec must be initialized");

  if (seed) 
    SetSeed(ZZ(seed));

  SetNumThreads(nthreads);

  TestIt(p,r,L,c,t,useCache);

  return 0;
}
