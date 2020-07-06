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
namespace std {} using namespace std;
namespace NTL {} using namespace NTL;

#include <NTL/BasicThreadPool.h>

#include <cassert>

#include <helib/EvalMap.h>
#include <helib/hypercube.h>
#include <helib/powerful.h>
#include <helib/ArgMap.h>

NTL_CLIENT
using namespace helib;

static bool dry = false; // a dry-run flag
static bool noPrint = true;

void  TestIt(long p, long r, long c, long _k,
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
       << ", L=" << L
       << ", mvec=" << mvec << ", "
       << ", useCache = " << useCache
       << endl;

  setTimersOn();
  setDryRun(false); // Need to get a "real context" to test EvalMap

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
  secretKey.GenSecKey(); // A Hamming-weight-w secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need

  // GG defines the plaintext space Z_p[X]/GG(X)
  ZZX GG;
  GG = context.alMod.getFactorsOverZZ()[0];
  EncryptedArray ea(context, GG);

  zz_p::init(context.alMod.getPPowR());
  zz_pX F;
  random(F, phim); // a random polynomial of degree phi(m)-1 modulo p

  // convert F to powerful representation: cube represents a multi-variate
  // polynomial with as many variables Xi as factors mi in mvec. cube has
  // degree phi(mi) in the variable Xi, and the coefficients are given
  // in lexicographic order.

  // compute tables for converting between powerful and zz_pX
  PowerfulTranslationIndexes ind(mvec); // indpendent of p
  PowerfulConversion pConv(ind);        // depends on p

  HyperCube<zz_p> cube(pConv.getShortSig());
  pConv.polyToPowerful(cube, F);

  // Sanity check: convert back and compare
  zz_pX F2;
  pConv.powerfulToPoly(F2, cube);
  if (F != F2) {
    cout << "BAD\n";
    if (!noPrint) cout << " @@@ conversion error ):\n";
  }
  // pack the coefficients from cube in the plaintext slots: the j'th
  // slot contains the polynomial pj(X) = \sum_{t=0}^{d-1} cube[jd+t] X^t
  vector<ZZX> val1;
  val1.resize(nslots);
  for (long i = 0; i < phim; i++) {
    val1[i/d] += conv<ZZX>(conv<ZZ>(cube[i])) << (i % d);
  }
  PlaintextArray pa1(ea);
  encode(ea, pa1, val1);

  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, pa1);

  resetAllTimers();
  HELIB_NTIMER_START(ALL);

  // Compute homomorphically the transformation that takes the
  // coefficients packed in the slots and produces the polynomial
  // corresponding to cube

  if (!noPrint) CheckCtxt(ctxt, "init");

  if (!noPrint) cout << "build EvalMap\n";
  EvalMap map(ea, /*minimal=*/false, mvec,
    /*invert=*/false, /*build_cache=*/false, /*normal_basis=*/false);
  // compute the transformation to apply

  if (!noPrint) cout << "apply EvalMap\n";
  if (useCache) map.upgrade();
  map.apply(ctxt); // apply the transformation to ctxt
  if (!noPrint) CheckCtxt(ctxt, "EvalMap");
  if (!noPrint) cout << "check results\n";

  ZZX FF1;
  secretKey.Decrypt(FF1, ctxt);
  zz_pX F1 = conv<zz_pX>(FF1);

  if (F1 == F)
    cout << "GOOD\n";
  else
    cout << "BAD\n";

  publicKey.Encrypt(ctxt, balanced_zzX(F1));
  if (!noPrint) CheckCtxt(ctxt, "init");

  // Compute homomorphically the inverse transformation that takes the
  // polynomial corresponding to cube and produces the coefficients
  // packed in the slots

  if (!noPrint) cout << "build EvalMap\n";
  EvalMap imap(ea, /*minimal=*/false, mvec,
    /*invert=*/true, /*build_cache=*/false, /*normal_basis=*/false);
  // compute the transformation to apply
  if (!noPrint) cout << "apply EvalMap\n";
  if (useCache) imap.upgrade();
  imap.apply(ctxt); // apply the transformation to ctxt
  if (!noPrint) {
    CheckCtxt(ctxt, "EvalMap");
    cout << "check results\n";
  }
  PlaintextArray pa2(ea);
  ea.decrypt(ctxt, secretKey, pa2);

  if (equals(ea, pa1, pa2))
    cout << "GOOD\n";
  else
    cout << "BAD\n";
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
  amap.arg("L", L, "# of bits in the modulus chain");

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
  TestIt(p, r, c, k, L, mvec, gens, ords, useCache);
}

// ./Test_EvalMap_x mvec="[73 433]" gens="[18620 12995]" ords="[72 -6]"
