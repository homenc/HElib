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
#include <helib/hypercube.h>
#include <helib/powerful.h>
#include <helib/Context.h>
#include <helib/ArgMap.h>

NTL_CLIENT
using namespace helib;

void testSimpleConversion(const Vec<long>& mvec)
{
  PowerfulTranslationIndexes ind(mvec);
  PowerfulConversion pConv;
  long q = NextPrime(ind.m);
  zz_p::init(q);
  pConv.initPConv(ind); // uses ind and also the current modulus q
  zz_pX poly, poly2;
  random(poly,ind.phim);

  HyperCube<zz_p> cube(pConv.getShortSig());
  pConv.polyToPowerful(cube, poly);
  pConv.powerfulToPoly(poly2, cube);
  cout << ((poly == poly2)? "GOOD" : "BAD") << endl;
}

void testHighLvlConversion(const Context& context, const Vec<long>& mvec)
{
  PowerfulDCRT p2d(context, mvec);
  DoubleCRT dcrt(context, context.fullPrimes());
  ZZX poly1, poly2;
  Vec<ZZ> pwfl;

  dcrt.randomize();
  dcrt.toPoly(poly1);

  p2d.ZZXtoPowerful(pwfl, poly1);
  p2d.powerfulToZZX(poly2, pwfl);

  cout << ((poly1 != poly2)? "BAD" : "GOOD") << endl;
}


int main(int argc, char *argv[])
{

  long m1 = 7;
  long m2 = 13;
  long m3 = 17;
  long p = 2;
  long r = 1;

  ArgMap amap;
  amap.arg("m1", m1, "Factor of m");
  amap.arg("m2", m2, "Factor of m");
  amap.arg("m3", m3, "Factor of m");
  amap.arg("p",  p,  "Plaintext base");
  amap.arg("r",  r,  "Lifting");
  amap.parse(argc, argv);

  // get parameters from the command line
//  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  if (m1<2 || m2<2) {
    cerr << "m1,m2 are mandatory\n";
    exit(0);
  }
  Vec<long> mvec(INIT_SIZE,2);
  mvec[0] = m1;
  mvec[1] = m2;
  if (m3>1) append(mvec,m3);
  long m = computeProd(mvec);

  cout << "m="<<m<<" "<<mvec<<", p="<<p<<", r="<<r<<endl;

  // Test conversion between zz_pX abd HyperCube<zz_p>
  testSimpleConversion(mvec);

  Context context(m,p,r);
  buildModChain(context, /*L=*/100, /*c=*/3);

  testHighLvlConversion(context, mvec);
  return 0;
  /****************** UNUSED OLD CODE, COMMENTED OUT *****************/
#if 0
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
#endif
}
