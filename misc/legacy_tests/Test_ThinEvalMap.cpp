/* Copyright (C) 2012-2019 IBM Corp.
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

/* Test_ThinEvalMap.cpp - Testing the evalution map for thin bootstrapping
 */
#include <cassert>
#include <helib/helib.h>
#include <helib/EvalMap.h>
#include <NTL/BasicThreadPool.h>
#include <helib/ArgMap.h>

NTL_CLIENT
using namespace helib;

static bool dry = false; // a dry-run flag
static bool noPrint = true;

void  TestIt(long p, long r, long c, long _k, long w,
             long L, Vec<long>& mvec,
             Vec<long>& gens, Vec<long>& ords, long useCache)
{
  if (lsize(mvec)<1) { // use default values
    mvec.SetLength(3); gens.SetLength(3); ords.SetLength(3);
    mvec[0] = 7;    mvec[1] = 3;    mvec[2] = 221;
    gens[0] = 3979; gens[1] = 3095; gens[2] = 3760;
    ords[0] = 6;    ords[1] = 2;    ords[2] = -8;
  }
  if (!noPrint)
    cout << "*** TestIt"
       << (dry? " (dry run):" : ":")
       << " p=" << p
       << ", r=" << r
       << ", c=" << c
       << ", k=" << _k
       << ", w=" << w
       << ", L=" << L
       << ", mvec=" << mvec << ", "
       << ", useCache = " << useCache
       << endl;

  setTimersOn();
  setDryRun(false); // Need to get a "real context" to test ThinEvalMap

  // mvec is supposed to include the prime-power factorization of m
  long nfactors = mvec.length();
  for (long i = 0; i < nfactors; i++)
    for (long j = i+1; j < nfactors; j++)
      assert(GCD(mvec[i], mvec[j]) == 1);

  // multiply all the prime powers to get m itself
  long m = computeProd(mvec);
  assert(GCD(p, m) == 1);

  // build a context with these generators and orders
  vector<long> gens1, ords1;
  convert(gens1, gens);
  convert(ords1, ords);
  Context context(m, p, r, gens1, ords1);
  buildModChain(context, L, c);

  if (!noPrint) {
    context.zMStar.printout(); // print structure of Zm* /(p) to cout
    cout << endl;
  }
  long d = context.zMStar.getOrdP();
  long phim = context.zMStar.getPhiM();
  long nslots = phim/d;

  setDryRun(dry); // Now we can set the dry-run flag if desired

  SecKey secretKey(context);
  const PubKey& publicKey = secretKey;
  secretKey.GenSecKey(w); // A Hamming-weight-w secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need

  // GG defines the plaintext space Z_p[X]/GG(X)
  ZZX GG;
  GG = context.alMod.getFactorsOverZZ()[0];
  EncryptedArray ea(context, GG);

  zz_p::init(context.alMod.getPPowR());

  Vec<zz_p> val0(INIT_SIZE, nslots);
  for (auto& x: val0)
    random(x);

  vector<ZZX> val1;
  val1.resize(nslots);
  for (long i = 0; i < nslots; i++) {
    val1[i] = conv<ZZX>(conv<ZZ>(rep(val0[i])));
  }

  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, val1);

  resetAllTimers();
  HELIB_NTIMER_START(ALL);

  // Compute homomorphically the transformation that takes the
  // coefficients packed in the slots and produces the polynomial
  // corresponding to cube

  if (!noPrint) CheckCtxt(ctxt, "init");

  if (!noPrint) cout << "build ThinEvalMap\n";
  ThinEvalMap map(ea, /*minimal=*/false, mvec,
    /*invert=*/false, /*build_cache=*/false);
  // compute the transformation to apply

  if (!noPrint) cout << "apply ThinEvalMap\n";
  if (useCache) map.upgrade();
  map.apply(ctxt); // apply the transformation to ctxt
  if (!noPrint) CheckCtxt(ctxt, "ThinEvalMap");
  if (!noPrint) cout << "check results\n";

  if (!noPrint) cout << "build ThinEvalMap\n";
  ThinEvalMap imap(ea, /*minimal=*/false, mvec,
    /*invert=*/true, /*build_cache=*/false);
  // compute the transformation to apply
  if (!noPrint) cout << "apply ThinEvalMap\n";
  if (useCache) imap.upgrade();
  imap.apply(ctxt); // apply the transformation to ctxt
  if (!noPrint) {
    CheckCtxt(ctxt, "ThinEvalMap");
    cout << "check results\n";
  }

#if 1

  /* create dirty version of ctxt */
  Vec<zz_pX> dirty_val0;
  dirty_val0.SetLength(nslots);
  for (long i = 0; i < nslots; i++) {
    random(dirty_val0[i], d);
    SetCoeff(dirty_val0[i], 0, val0[i]);
  }

  vector<ZZX> dirty_val1;
  dirty_val1.resize(nslots);
  for (long i = 0; i < nslots; i++) {
    dirty_val1[i] = conv<ZZX>(dirty_val0[i]);
  }

  Ctxt dirty_ctxt(publicKey);
  ea.encrypt(dirty_ctxt, publicKey, dirty_val1);


  EvalMap dirty_map(ea, /*minimal=*/false, mvec,
    /*invert=*/false, /*build_cache=*/false);

  dirty_map.apply(dirty_ctxt);
  imap.apply(dirty_ctxt);
#endif


  vector<ZZX> val2;
  ea.decrypt(ctxt, secretKey, val2);
  cout << ((val1 == val2)? "GOOD\n" : "BAD\n");

  vector<ZZX> dirty_val2;
  ea.decrypt(dirty_ctxt, secretKey, dirty_val2);
  cout << ((val1 == dirty_val2)? "GOOD\n" : "BAD\n");


  HELIB_NTIMER_STOP(ALL);

  if (!noPrint) {
    cout << "\n*********\n";
    printAllTimers();
    cout << endl;
  }
}


/* Usage: Test_EvalMap_x.exe [ name=value ]...
 *  p       plaintext base  [ default=2 ]
 *  r       lifting  [ default=1 ]
 *  c       number of columns in the key-switching matrices  [ default=2 ]
 *  k       security parameter  [ default=80 ]
 *  L       # of bits in the modulus chain
 *  s       minimum number of slots  [ default=0 ]
 *  seed    PRG seed  [ default=0 ]
 *  mvec    use specified factorization of m
 *             e.g., mvec='[5 3 187]'
 *  gens    use specified vector of generators
 *             e.g., gens='[562 1871 751]'
 *  ords    use specified vector of orders
 *             e.g., ords='[4 2 -4]', negative means 'bad'
 */
int main(int argc, char *argv[])
{
  ArgMap amap;

  long p=2;
  amap.arg("p", p, "plaintext base");

  long r=1;
  amap.arg("r", r,  "lifting");

  long c=2;
  amap.arg("c", c, "number of columns in the key-switching matrices");

  long k=80;
  amap.arg("k", k, "security parameter");

  long L=300;
  amap.arg("L", L, "# of levels in the modulus chain");

  long s=0;
  amap.arg("s", s, "minimum number of slots");

  long seed=0;
  amap.arg("seed", seed, "PRG seed");

  Vec<long> mvec;
  amap.arg("mvec", mvec, "use specified factorization of m", nullptr);
  amap.note("e.g., mvec='[7 3 221]'");

  Vec<long> gens;
  amap.arg("gens", gens, "use specified vector of generators", nullptr);
  amap.note("e.g., gens='[3979 3095 3760]'");

  Vec<long> ords;
  amap.arg("ords", ords, "use specified vector of orders", nullptr);
  amap.note("e.g., ords='[6 2 -8]', negative means 'bad'");

  amap.arg("dry", dry, "a dry-run flag to check the noise");

  long nthreads=1;
  amap.arg("nthreads", nthreads, "number of threads");

  amap.arg("noPrint", noPrint, "suppress printouts");

  long useCache=0;
  amap.arg("useCache", useCache, "0: zzX cache, 2: DCRT cache");

  amap.parse(argc, argv);

  SetNumThreads(nthreads);

  SetSeed(conv<ZZ>(seed));
  TestIt(p, r, c, k, /*Key Hamming weight=*/64, L, mvec, gens, ords, useCache);
}

// ./Test_ThinEvalMap_x mvec="[73 433]" gens="[18620 12995]" ords="[72 -6]"
