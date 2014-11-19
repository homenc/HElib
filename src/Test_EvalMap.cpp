namespace std {} using namespace std;
namespace NTL {} using namespace NTL;

#include "EvalMap.h"
#include "hypercube.h"
#include "powerful.h"

void  TestIt(long p, long r, long c, long _k, long w,
             long L, const Vec<long>& mvec, 
             const Vec<long>& gens, const Vec<long>& ords )
{
  cerr << "*** TestIt: "
       << "  p=" << p
       << ", r=" << r
       << ", c=" << c
       << ", k=" << _k
       << ", w=" << w
       << ", L=" << L
       << ", mvec=" << mvec
       << endl;

  setTimersOn();

  // mvec is supposed to include the prime-power factorization of m
  long nfactors = mvec.length();
  for (long i = 0; i < nfactors; i++)
    for (long j = i+1; j < nfactors; j++)
      assert(GCD(mvec[i], mvec[j]) == 1);

  // multiply all the prime powers to get m itself
  long m = computeProd(mvec);
  assert(GCD(p, m) == 1);

  vector<long> gens1, ords1;
  convert(gens1, gens);
  convert(ords1, ords);


  FHEcontext context(m, p, r, gens1, ords1);
  buildModChain(context, L, c);
  context.zMStar.printout();
  cerr << endl;

  long d = context.zMStar.getOrdP();
  long phim = context.zMStar.getPhiM();
  long nslots = phim/d;


  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(w); // A Hamming-weight-w secret key

  cerr << "generating key-switching matrices... ";
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need
  cerr << "done\n";


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
  if (F != F2) cerr << " @@@ conversion error ):\n";

  // pack the coefficients from cube in the plaintext slots: the j'th
  // slot contains the polynomial pj(X) = \sum_{t=0}^{d-1} cube[jd+t] X^t
  vector<ZZX> val1;
  val1.resize(nslots);
  for (long i = 0; i < phim; i++) {
    val1[i/d] += conv<ZZX>(conv<ZZ>(cube[i])) << (i % d);
  }
  PlaintextArray pa1(ea);
  pa1.encode(val1);

  PlaintextArray check_val(ea);
  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, pa1);

  resetAllTimers();

  FHE_NTIMER_START(ALL);

  // Compute homomorphically the transformation that takes the
  // coefficients packed in the slots and produces the polynomial
  // corresponding to cube

  CheckCtxt(ctxt, "init");

  cout << "build EvalMap\n";
  EvalMap map(ea, mvec, false); // compute the transformation to apply
  cout << "apply EvalMap\n";
  map.apply(ctxt);                     // apply the transformation to ctxt
  CheckCtxt(ctxt, "EvalMap");
  cout << "check results\n";

  ZZX FF1;
  secretKey.Decrypt(FF1, ctxt);
  zz_pX F1 = conv<zz_pX>(FF1);

  if (F1 == F)
    cout << "EvalMap: GOOD\n";
  else
    cout << "EvalMap: BAD\n";

#if 1
  // This test no longer works, because the inverse map 
  // now applies the normal basis transformation

  publicKey.Encrypt(ctxt, FF1);
  CheckCtxt(ctxt, "init");

  // Compute homomorphically the inverse transformation that takes the
  // polynomial corresponding to cube and produces the coefficients
  // packed in the slots

  cout << "build EvalMap\n";
  EvalMap imap(ea, mvec, true, false); // compute the transformation to apply
  cout << "apply EvalMap\n";
  imap.apply(ctxt);                    // apply the transformation to ctxt
  CheckCtxt(ctxt, "EvalMap");
  cout << "check results\n";

  PlaintextArray pa2(ea);
  ea.decrypt(ctxt, secretKey, pa2);

  if (pa1.equals(pa2))
    cout << "EvalMap: GOOD\n";
  else
    cout << "EvalMap: BAD\n";
#endif

  FHE_NTIMER_STOP(ALL);

  cerr << "\n*********\n";
  printAllTimers();
  cerr << endl;
}


int main(int argc, char *argv[])
{
  ArgMapping amap;

  long p=2;
  amap.arg("p", p, "plaintext base");

  long r=1;
  amap.arg("r", r,  "lifting");

  long c=2;
  amap.arg("c", c, "number of columns in the key-switching matrices");
  
  long k=80;
  amap.arg("k", k, "security parameter");

  long L=6;
  amap.arg("L", L, "# of levels in the modulus chain");

  long s=0;
  amap.arg("s", s, "minimum number of slots");

  long seed=0;
  amap.arg("seed", seed, "PRG seed");

  Vec<long> gens;
  amap.arg("gens", gens, "use specified vector of generators", NULL);
  amap.note("e.g., 'gens=[17 63 42]'");

  Vec<long> ords;
  amap.arg("ords", ords, "use specified vector of orders", NULL);
  amap.note("e.g., 'ords=[100 20 -4]', negative means 'bad'");

  Vec<long> mvec;
  amap.arg("mvec", mvec, "use specified factorization of m", NULL);

  amap.parse(argc, argv);

  long w = 64; // Hamming weight of secret key

  if (seed) SetSeed(conv<ZZ>(seed));

  TestIt(p, r, c, k, w, L, mvec, gens, ords);
}

