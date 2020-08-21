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
#include <helib/ArgMap.h>

#include <algorithm>
#include <cmath>
#include <string>

NTL_CLIENT
using namespace helib;

static bool noPrint = false;
static bool dry = false; // a dry-run flag
static bool debug = 0;   // a debug flag
static int scale = 0;
static long special_bits;

static string v_values_name = "";


static long OUTER_REP = 1;


static Vec<long> global_mvec, global_gens, global_ords;
static int c_m = 100;


static void
dump_v_values()
{
  const vector<double> *v_values = fetch_saved_values("v_values");
  if (v_values && v_values_name != "") {
    // write v_values to a file

    cerr << "writing v_values to " << v_values_name << "\n";

    ofstream F;
    F.open(v_values_name.c_str());
    for (long i: range(v_values->size()))
      F << (*v_values)[i] << "\n";
  }
}

static void
anderson_darling(const vector<double>& X, double& AD, double& p_val)
{
  long N = X.size();

  if (N < 2) {
    AD = 0;
    p_val = 1;
    return;
  }

  vector<double> Y(X);
  sort(Y.begin(), Y.end());

  // compute the sample mean
  double SM = 0;
  for (long i: range(N)) SM += Y[i];
  SM /= N;


  // compute the sample variance
  double SV = 0;
  for (long i: range(N)) SV += (Y[i]-SM)*(Y[i]-SM);
  SV /= (N-1);

  // replace Y[i] by CDF of Y[i]
  for (long i: range(N))
    Y[i] = 0.5*(1 + erf((Y[i]-SM)/sqrt(2*SV)));

  double S = 0;
  for (long i: range(N)) {
    S += (2*i+1)*(log(Y[i]) + log1p(-Y[N-1-i]));
  }
  AD = -N - S/N;

  AD *= (1 + 0.75/N + 2.25/N/N);
  // This adjustment and the p-values below come from:
  // R.B. D'Augostino and M.A. Stephens, Eds., 1986,
  // Goodness-of-Fit Techniques, Marcel Dekker.

  if (AD >= 0.6)
    p_val = exp(1.2937 - 5.709*(AD)+ 0.0186*fsquare(AD));
  else if (AD > 0.34)
    p_val = exp(0.9177 - 4.279*(AD) - 1.38*fsquare(AD));
  else if (AD > 0.2)
    p_val = 1 - exp(-8.318 + 42.796*(AD)- 59.938*fsquare(AD));
  else
    p_val = 1 - exp(-13.436 + 101.14*(AD)- 223.73*fsquare(AD));
}

static void
print_anderson_darling()
{
  const vector<double> *v_values = fetch_saved_values("v_values");
  if (v_values) {
    double AD, p_val;
    anderson_darling(*v_values, AD, p_val);
    cout << "AD=" << AD << ", p_val=" << p_val << "\n";
  }
  else
    cout << "\n";
}




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

  if (!noPrint) fhe_stats = true;

  double t = -GetTime();
  Context context(m, p, r, gens, ords);
  if (scale) {
    context.scale = scale;
  }

  context.zMStar.set_cM(c_m/100.0);
  buildModChain(context, L, c,
    /*willBeBootstrappable=*/true,
    /*t=*/skHwt,
    /*resolution=*/3,
    /*bitsInSpecialPrimes=*/special_bits);

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


  for (long numkey=0; numkey<OUTER_REP; numkey++) { // test with 3 keys
  if (1 && fhe_stats && numkey > 0 && numkey%100 == 0) {
    print_stats(cout);
    cout << "mvec=" << mvec << ", ";
    print_anderson_darling();
    dump_v_values();
  }

  if (!noPrint) cerr << "*********** iter=" << (numkey+1) << "\n";

  t = -GetTime();
  if (!noPrint) cout << "Generating keys, " << std::flush;
  SecKey secretKey(context);
  secretKey.GenSecKey(skHwt);      // A Hamming-weight-64 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey);
  if (!noPrint) cout << "computing key-dependent tables..." << std::flush;
  secretKey.genRecryptData();
  t += GetTime();
  if (!noPrint) cout << " done in "<<t<<" seconds\n";

#ifdef HELIB_DEBUG
      dbgEa = context.ea;
      dbgKey = &secretKey;
#endif

  const PubKey publicKey = secretKey;

  long d = context.zMStar.getOrdP();
  long phim = context.zMStar.getPhiM();
  long nslots = phim/d;

  // GG defines the plaintext space Z_p[X]/GG(X)
  ZZX GG;
  GG = context.alMod.getFactorsOverZZ()[0];
  std::shared_ptr<EncryptedArray> ea_ptr(std::make_shared<EncryptedArray>(context, GG));
  // Alias to avoid issues with previous code
  EncryptedArray& ea(*ea_ptr);

  if (debug) {
    dbgKey = &secretKey;
    dbgEa = ea_ptr;
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

  // Make some noise!
  // This ensures that we do not start the sequence of squarings
  // from a "fresh" ciphertext (which may not be representive).
  ea.rotate(c1, 1);
  ea.rotate(c1, -1);

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
    // compute minimal capacity before bootstrapping (rawModSwitch)
    long e = context.rcData.e;
    long q = power_long(p, e) + 1;
    double Bnd = context.boundForRecryption();
    double mfac = context.zMStar.getNormBnd();
    double min_bit_cap = log(  mfac*q / (p2r*Bnd*HELIB_MIN_CAP_FRAC) )/log(2.0);

    cout << "min_bit_cap=" << min_bit_cap << "\n";

    cout << "log2(modSwitchAddedNoiseBound)="
         << log(c1.modSwitchAddedNoiseBound())/log(2.0) << "\n";
    cout << "sqr_count=" << sqr_count << "\n";
    if (sqr_count > 0) {
      cout << "BITS-PER-LEVEL: "
           << ((c1.bitCapacity()-c2.bitCapacity())/double(sqr_count)) << "\n";
    }
    CheckCtxt(c2, "before recryption");
  }

  resetAllTimers();

  { HELIB_NTIMER_START(AAA_thinRecrypt);

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

#if (defined(__unix__) || defined(__unix) || defined(unix))
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
    if (!noPrint) cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif

  if (fhe_stats) {
    print_stats(cout);
    cout << "mvec=" << mvec << ", ";
    print_anderson_darling();
    dump_v_values();
  }

}


//extern long fhe_disable_intFactor;
// extern long fhe_disable_fat_boot;
// extern long fhe_force_chen_han;

/********************************************************************
 ********************************************************************/
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
  amap.arg("special_bits", special_bits, "# of bits in special primes");

  amap.arg("gens", global_gens);
  amap.arg("ords", global_ords);
  amap.arg("mvec", global_mvec);
  amap.arg("c_m", c_m);

  amap.arg("iter", OUTER_REP);

  amap.arg("v_values", v_values_name);

  amap.parse(argc, argv);

  if (global_gens.length() == 0 || global_ords.length() == 0 || global_mvec.length() == 0)
    Error("gens, ords, and mvec must be initialized");

  if (seed)
    SetSeed(ZZ(seed));

  SetNumThreads(nthreads);

  TestIt(p,r,L,c,t,useCache);

  return 0;
}
