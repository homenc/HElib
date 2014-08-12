

namespace std {} using namespace std;
namespace NTL {} using namespace NTL;

#include "EvalMap.h"
#include "hypercube.h"
#include "powerful.h"

void convertToPowerful(Vec<zz_p>& v, const zz_pX& F, const Vec<long>& mvec)
{ 
  long nfactors = mvec.length();

  long m = computeProd(mvec);
  
  Vec<long> phivec;
  phivec.SetLength(nfactors);
  for (long i = 0; i < nfactors; i++) phivec[i] = phi_N(mvec[i]);

  long phim = computeProd(phivec);

  Vec<long> divvec;
  computeDivVec(divvec, m, mvec);

  Vec<long> invvec;
  computeInvVec(invvec, divvec, mvec);

  CubeSignature shortsig(phivec);
  CubeSignature longsig(mvec);

  Vec<long> polyToCubeMap;
  Vec<long> cubeToPolyMap;
  computePowerToCubeMap(polyToCubeMap, cubeToPolyMap, m, mvec, invvec, longsig);

  Vec<long> shortToLongMap;
  computeShortToLongMap(shortToLongMap, shortsig, longsig);


  Vec<zz_pX> cycvec;
  computeCycVec(cycvec, mvec);


  ZZX PhimX = Cyclotomic(m);
  zz_pX phimX = conv<zz_pX>(PhimX);

  HyperCube<zz_p> cube(shortsig);
  HyperCube<zz_p> tmpCube(longsig);

  convertPolyToPowerful(cube, tmpCube, F, cycvec, 
                        polyToCubeMap, shortToLongMap);

  zz_pX poly1;

  convertPowerfulToPoly(poly1, cube, m, shortToLongMap, cubeToPolyMap, phimX);

  if (F == poly1)
    cout << "*********** :-)\n";
  else {
    cout << "*********** :-(\n";
    cout << F << "\n";
    cout << poly1 << "\n";
  }

  v.SetLength(phim);
  for (long i = 0; i < phim; i++) v[i] = cube[i];
}


void  TestIt(long R, long p, long r, long c, long _k, long w, 
               long L, const Vec<long>& mvec, long width)
{
  cerr << "*** TestIt: R=" << R 
       << ", p=" << p
       << ", r=" << r
       << ", c=" << c
       << ", k=" << _k
       << ", w=" << w
       << ", L=" << L
       << ", mvec=" << mvec
       << ", width=" << width 
       << endl;

  setTimersOn();

  long nfactors = mvec.length();
  for (long i = 0; i < nfactors; i++)
    for (long j = i+1; j < nfactors; j++)
      assert(GCD(mvec[i], mvec[j]) == 1);


  long m = computeProd(mvec);
  assert(GCD(p, m) == 1); 



  FHEcontext context(m, p, r); 

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


  ZZX GG;
  GG = context.alMod.getFactorsOverZZ()[0];

  EncryptedArray ea(context, GG);

  zz_p::init(context.alMod.getPPowR());


  zz_pX F;
  random(F, phim);


  Vec<zz_p> cube;
  convertToPowerful(cube, F, mvec);
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


  CheckCtxt(ctxt, "init");

  cout << "build EvalMap\n";
  EvalMap map(ea, mvec, width, false);
  cout << "apply EvalMap\n";
  map.apply(ctxt);
  CheckCtxt(ctxt, "EvalMap");
  cout << "check results\n";

  ZZX FF1;
  secretKey.Decrypt(FF1, ctxt);
  zz_pX F1 = conv<zz_pX>(FF1);

  if (F1 == F) 
    cout << "EvalMap: good\n";
  else
    cout << "EvalMap: bad\n";

  publicKey.Encrypt(ctxt, FF1);
  CheckCtxt(ctxt, "init");

  cout << "build EvalMap\n";
  EvalMap imap(ea, mvec, width, true);
  cout << "apply EvalMap\n";
  imap.apply(ctxt);
  CheckCtxt(ctxt, "EvalMap");
  cout << "check results\n";

  PlaintextArray pa2(ea);
  ea.decrypt(ctxt, secretKey, pa2);

  if (pa1.equals(pa2))
    cout << "EvalMap: good\n";
  else
    cout << "EvalMap: bad\n";

  FHE_NTIMER_STOP(ALL);

  cerr << "\n*********\n";
  printAllTimers();
  cerr << endl;
}


void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'R=1 p=2 k=80'\n\n";
  cerr << "  R is the number of rounds\n";
  cerr << "  p is the plaintext base [default=2]" << endl;
  cerr << "  r is the lifting [default=1]" << endl;
  cerr << "  d is the degree of the field extension [default==0]\n";
  cerr << "    (d == 0 => factors[0] defined the extension)\n";
  cerr << "  c is number of columns in the key-switching matrices [default=2]\n";
  cerr << "  k is the security parameter [default=80]\n";
  cerr << "  L is the # of primes in the modulus chai [default=4*R]\n";
  cerr << "  s is the minimum number of slots [default=4]\n";
  cerr << "  m defined the cyclotomic polynomial Phi_m(X)\n";
  cerr << "  seed is the PRG seed\n";
  exit(0);
}


int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["R"] = "1";
  argmap["p"] = "2";
  argmap["r"] = "1";
  argmap["c"] = "2";
  argmap["k"] = "80";
  argmap["L"] = "0";
  argmap["s"] = "0";
  argmap["m1"] = "0";
  argmap["m2"] = "0";
  argmap["m3"] = "0";
  argmap["m4"] = "0";
  argmap["width"] = "5";
  argmap["seed"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long R = atoi(argmap["R"]);
  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long c = atoi(argmap["c"]);
  long k = atoi(argmap["k"]);
  //  long z = atoi(argmap["z"]);
  long L = atoi(argmap["L"]);
  if (L==0) { // determine L based on R,r
    if (r==1) L = 2*R+2;
    else      L = 4*R;
  }
  long s = atoi(argmap["s"]);

  long m1 = atoi(argmap["m1"]);
  long m2 = atoi(argmap["m2"]);
  long m3 = atoi(argmap["m3"]);
  long m4 = atoi(argmap["m4"]);
  long width = atoi(argmap["width"]);
  long seed = atoi(argmap["seed"]);

  long w = 64; // Hamming weight of secret key
  //  long L = z*R; // number of levels

  Vec<long> mvec;
  if (m1 != 0) append(mvec, m1);
  if (m2 != 0) append(mvec, m2);
  if (m3 != 0) append(mvec, m3);
  if (m4 != 0) append(mvec, m4);
  

  if (seed) SetSeed(conv<ZZ>(seed));

  TestIt(R, p, r, c, k, w, L, mvec, width);


}

//   [1 1 3 8] Test_EvalMap_x p=2 m1=3 m2=5 m3=7 m4=17
//   Test_EvalMap_x p=2 m1=11 m2=41 m3=31 (phim1=6, phim2=40, d2=4)
//   [1 1 20]  Test_EvalMap_x p=2 m1=3 m2=11 m3=25
//   [1 1 20]  Test_EvalMap_x p=2 m1=3 m2=11 m3=41
//   Test_EvalMap_x p=2 m1=3 m2=11 m3=17
//   Test_EvalMap_x p=2 m1=3 m2=5 m3=43 (phim1 == 3)
//   Test_EvalMap_x p=2 m1=7 m2=13 m3=73 (phim1=8, phim2=12, d2=4)
//   Test_EvalMap_x p=2 m1=7 m2=33 m3=73 (phim1=8, phim2=20, d2=10)


    




