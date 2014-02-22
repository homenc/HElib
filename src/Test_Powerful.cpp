#include "hypercube.h"
#include "powerful.h"

void usage()
{
  cerr << "bad args\n";
  exit(0);
}

int main(int argc, char *argv[])
{

  argmap_t argmap;

  argmap["m"] = "225"; // 225 = 5^2*3^2
  argmap["iter"] = "10";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage();

  long m = atoi(argmap["m"]);
  long iter = atoi(argmap["iter"]);

  cout << "m=" << m << "\n";

  long q; // find least prime q s/t q = 2*k*m + 1 for some k

  for (long k = 1; q = 2*k*m + 1, !ProbPrime(q, 20); k++);
  cout << "q=" << q << "\n";
  

  Vec< Pair<long, long> > factors;

  factorize(factors, m);

  cout << factors << "\n";

  Vec<long> phiVec;
  computePhiVec(phiVec, factors);
  cout << phiVec << "\n";

  long phim = computeProd(phiVec);
  cout << phim << "\n";

  Vec<long> powVec;
  computePowVec(powVec, factors);
  cout << powVec << "\n";

  Vec<long> divVec;
  computeDivVec(divVec, m, powVec);

  Vec<long> invVec;
  computeInvVec(invVec, divVec, powVec);


  CubeSignature shortSig(phiVec);
  CubeSignature longSig(powVec);

  Vec<long> polyToCubeMap;
  Vec<long> cubeToPolyMap;
  computePowerToCubeMap(polyToCubeMap, cubeToPolyMap, m, powVec, invVec, longSig);

  Vec<long> shortToLongMap;
  computeShortToLongMap(shortToLongMap, shortSig, longSig);

  Vec<long> longToShortMap;
  computeLongToShortMap(longToShortMap, m, shortToLongMap);
  

  zz_p::init(q);

  Vec<zz_pX> cycVec;
  computeCycVec(cycVec, powVec);


  ZZX PhimX = Cyclotomic(m);
  zz_pX phimX = conv<zz_pX>(PhimX);


  zz_pX poly;
  random(poly, phim);

  HyperCube<zz_p> cube(shortSig);
  HyperCube<zz_p> tmpCube(longSig);

  zz_pX tmp1, tmp2;

  convertPolyToPowerful(cube, tmpCube, poly, cycVec, 
                        polyToCubeMap, shortToLongMap);

  zz_pX poly1;

  convertPowerfulToPoly(poly1, cube, m, shortToLongMap, cubeToPolyMap, phimX);

  if (poly == poly1)
    cout << ":-)\n";
  else {
    cout << ":-(\n";
    cout << poly << "\n";
    cout << poly1 << "\n";
  }

  long lbase = 1;
  while (lbase < q && multOrd(lbase, q) != m) lbase++;

  assert(lbase < q);

  zz_p base = conv<zz_p>(lbase);

  // Vec< Vec<zz_p> > multiEvalPoints;
  Vec< copied_ptr<FFTHelper> > multiEvalPoints;

  computeMultiEvalPoints(multiEvalPoints, base, m, powVec, phiVec);

  Vec<zz_p> linearEvalPoints;
  computeLinearEvalPoints(linearEvalPoints, base, m, phim);

  Vec< Vec<long> > compressedIndex;
  computeCompressedIndex(compressedIndex, powVec);

  Vec<long> powToCompressedIndexMap;
  computePowToCompressedIndexMap(powToCompressedIndexMap, m, powVec, 
                                 compressedIndex, shortSig);


  
  HyperCube<zz_p> cube1 = cube;

  eval(cube1, multiEvalPoints);

  Vec<zz_p> res;
  eval(res, poly, linearEvalPoints);

  bool eval_ok = true;

  for (long i = 0, j = 0; i < m; i++) {
    if (GCD(i, m) == 1) {
      if (cube1[powToCompressedIndexMap[i]] != res[j++]) eval_ok = false;
    }
  } 


  if (eval_ok)
    cout << "eval ok\n";
  else
    cout << "eval not ok\n";


  interp(cube1, multiEvalPoints);

  if (cube1 == cube)
    cout << "interp OK\n";
  else
    cout << "interp NOT OK\n";
  

  FFTHelper fftHelper(m, base);
 
  Vec<zz_p> res1;

  fftHelper.FFT(poly, res1);

  if (res1 == res) 
    cout << "fast eval OK\n";
  else
    cout << "fast eval NOT OK\n";

  zz_pX poly2;

  fftHelper.iFFT(poly2, res);

  if (poly2 == poly)
    cout << "fast interp OK\n";
  else
    cout << "fast interp NOT OK\n";



  double t;

  cout << "time for cube eval...";

  t = GetTime();

  HyperCube<zz_p> cube2(shortSig);

  eval(cube1, multiEvalPoints);
  t = GetTime();
  for (long i = 0; i < iter; i++) {
    cube2 = cube;
    eval(cube2, multiEvalPoints);
  }
  t = GetTime()-t;
  cout << t << "\n";


  cout << "time for linear eval...";
  
  t = GetTime();
  for (long i = 0; i < iter; i++) {
    fftHelper.FFT(poly, res1);
  }
  t = GetTime()-t;
  cout << t << "\n";
  
}
