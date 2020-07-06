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
#include <iostream>
#include <cassert>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <NTL/BasicThreadPool.h>
NTL_CLIENT

#include <helib/helib.h>

#include <helib/intraSlot.h>
#include <helib/binaryArith.h>
#include <helib/binaryCompare.h>
#include <helib/ArgMap.h>

#ifdef HELIB_DEBUG
#include <helib/debugging.h>
#endif

using namespace helib;

// define flags FLAG_PRINT_ZZX, FLAG_PRINT_POLY, FLAG_PRINT_VEC, functions
//        decryptAndPrint(ostream, ctxt, sk, ea, flags)
//        decryptAndCompare(ctxt, sk, ea, pa);

static std::vector<zzX> unpackSlotEncoding; // a global variable
static bool verbose=false;

static long mValues[][15] = {
// { p, phi(m),   m,   d, m1, m2, m3,    g1,   g2,   g3, ord1,ord2,ord3, B,c}
  {  2,    48,   105, 12,  3, 35,  0,    71,    76,    0,   2,  2,   0, 25, 2},
  {  2 ,  600,  1023, 10, 11, 93,  0,   838,   584,    0,  10,  6,   0, 25, 2},
  {  2,  2304,  4641, 24,  7,  3,221,  3979,  3095, 3760,   6,  2,  -8, 25, 3},
  {  2, 15004, 15709, 22, 23,683,  0,  4099, 13663,    0,  22, 31,   0, 25, 3},
  {  2, 27000, 32767, 15, 31,  7, 151, 11628, 28087,25824, 30,  6, -10, 28, 4}
};

void testCompare(SecKey& secKey, long bitSize, bool bootstrap=false);


int main(int argc, char *argv[])
{
  ArgMap amap;
  long prm=1;
  amap.arg("prm", prm, "parameter size (0-tiny,...,4-huge)");
  long bitSize = 5;
  amap.arg("bitSize", bitSize, "bitSize of input integers (<=32)");
  long nTests = 3;
  amap.arg("nTests", nTests, "number of tests to run");
  bool bootstrap = false;
  amap.arg("bootstrap", bootstrap, "test comparison with bootstrapping");
  long seed=0;
  amap.arg("seed", seed, "PRG seed");
  long nthreads=1;
  amap.arg("nthreads", nthreads, "number of threads");
  amap.arg("verbose", verbose, "print more information");

  amap.parse(argc, argv);
  assert(prm >= 0 && prm < 5);
  if (seed) NTL::SetSeed(ZZ(seed));
  if (nthreads>1) NTL::SetNumThreads(nthreads);

  if (bitSize<=0) bitSize=5;
  else if (bitSize>32) bitSize=32;

  long* vals = mValues[prm];
  long p = vals[0];
  //  long phim = vals[1];
  long m = vals[2];

  NTL::Vec<long> mvec;
  append(mvec, vals[4]);
  if (vals[5]>1) append(mvec, vals[5]);
  if (vals[6]>1) append(mvec, vals[6]);

  std::vector<long> gens;
  gens.push_back(vals[7]);
  if (vals[8]>1) gens.push_back(vals[8]);
  if (vals[9]>1) gens.push_back(vals[9]);

  std::vector<long> ords;
  ords.push_back(vals[10]);
  if (abs(vals[11])>1) ords.push_back(vals[11]);
  if (abs(vals[12])>1) ords.push_back(vals[12]);

  long B = vals[13];
  long c = vals[14];

  // Compute the number of levels
  long L;
  if (bootstrap) L = 900; // that should be enough
  else           L = 30*(7+ NTL::NumBits(bitSize+2));

  if (verbose) {
    cout <<"input bitSize="<<bitSize
         <<", running "<<nTests<<" tests for each function\n";
    if (nthreads>1) cout << "  using "<<NTL::AvailableThreads()<<" threads\n";
    cout << "computing key-independent tables..." << std::flush;
  }
  Context context(m, p, /*r=*/1, gens, ords);
  buildModChain(context, L, c,/*willBeBootstrappable=*/bootstrap);
  if (bootstrap) {
    context.makeBootstrappable(mvec, /*t=*/0);
  }
  buildUnpackSlotEncoding(unpackSlotEncoding, *context.ea);
  if (verbose) {
    cout << " done.\n";
    context.zMStar.printout();
    cout << " L="<<L<<", B="<<B<<endl;
    cout << "\ncomputing key-dependent tables..." << std::flush;
  }
  SecKey secKey(context);
  secKey.GenSecKey();
  addSome1DMatrices(secKey); // compute key-switching matrices
  addFrbMatrices(secKey);
  if (bootstrap) secKey.genRecryptData();
  if (verbose) cout << " done\n";

  activeContext = &context; // make things a little easier sometimes
#ifdef HELIB_DEBUG
  dbgEa = context.ea;
  dbgKey = &secKey;
#endif

  for (long i=0; i<nTests; i++)
    testCompare(secKey, bitSize, bootstrap);
  cout << "GOOD\n";

  if (verbose) printAllTimers(cout);
  return 0;
}


void testCompare(SecKey& secKey, long bitSize, bool bootstrap)
{
  const Context& context = secKey.getContext();
  const EncryptedArray& ea = *(context.ea);

  // Choose two random n-bit integers
  long pa = RandomBits_long(bitSize);
  long pb = RandomBits_long(bitSize+1);
  long pMax = std::max(pa,pb);
  long pMin = std::min(pa,pb);
  bool pMu = pa>pb;
  bool pNi = pa<pb;

  // Encrypt the individual bits
  NTL::Vec<Ctxt> eMax, eMin, enca, encb;

  Ctxt mu(secKey), ni(secKey);
  resize(enca, bitSize, mu);
  resize(encb, bitSize+1, ni);
  for (long i=0; i<=bitSize; i++) {
    if (i<bitSize) secKey.Encrypt(enca[i], ZZX((pa>>i)&1));
    secKey.Encrypt(encb[i], ZZX((pb>>i)&1));
    if (bootstrap) { // put them at a lower level
      if (i<bitSize) enca[i].bringToSet(context.getCtxtPrimes(5));
      encb[i].bringToSet(context.getCtxtPrimes(5));
    }
  }
#ifdef HELIB_DEBUG
  decryptAndPrint((cout<<" before comparison: "), encb[0], secKey, ea,0);
#endif

  vector<long> slotsMin, slotsMax, slotsMu, slotsNi;

  //cmp only (as unsigned integer)
  compareTwoNumbers(mu, ni, CtPtrs_VecCt(enca), CtPtrs_VecCt(encb), false,
                    &unpackSlotEncoding);
  ea.decrypt(mu, secKey, slotsMu);
  ea.decrypt(ni, secKey, slotsNi);
  if (slotsMu[0]!=pMu || slotsNi[0]!=pNi) {
    cout << "BAD\n";
    if (verbose)
      cout << "Comparison (without min max) error: a="<<pa<<", b="<<pb
           << ", mu="<<slotsMu[0]<<", ni="<<slotsNi[0]<<endl;
    exit(0);
  }
  else if (verbose) {
    cout << "Comparison (without min max) succeeded: ";
    cout << '('<<pa<<','<<pb<<")=> mu="<<slotsMu[0]<<", ni="<<slotsNi[0]<<endl;
  }

  {CtPtrs_VecCt wMin(eMin), wMax(eMax); // A wrappers around output vectors

  //cmp with max and min (as unsigned integer)
  compareTwoNumbers(wMax, wMin, mu, ni,
                    CtPtrs_VecCt(enca), CtPtrs_VecCt(encb), false,
                    &unpackSlotEncoding);
  decryptBinaryNums(slotsMax, wMax, secKey, ea);
  decryptBinaryNums(slotsMin, wMin, secKey, ea);
  } // get rid of the wrapper
  ea.decrypt(mu, secKey, slotsMu);
  ea.decrypt(ni, secKey, slotsNi);

  if (slotsMax[0]!=pMax || slotsMin[0]!=pMin
      || slotsMu[0]!=pMu || slotsNi[0]!=pNi) {
    cout << "BAD\n";
    if (verbose)
      cout << "Comparison (with min max) error: a="<<pa<<", b="<<pb
           << ", but min="<<slotsMin[0]<<", max="<<slotsMax[0]
           << ", mu="<<slotsMu[0]<<", ni="<<slotsNi[0]<<endl;
    exit(0);
  }
  else if (verbose) {
    cout << "Comparison (with min max) succeeded: ";
    cout << '('<<pa<<','<<pb<<")=>("<<slotsMin[0]<<','<<slotsMax[0]
         <<"), mu="<<slotsMu[0]<<", ni="<<slotsNi[0]<<endl;
  }

#ifdef HELIB_DEBUG
  const Ctxt* minLvlCtxt = nullptr;
  long minLvl=1000;
  for (const Ctxt& c: eMax) {
    long lvl = c.logOfPrimeSet();
    if (lvl < minLvl) {
      minLvlCtxt = &c;
      minLvl = lvl;
    }
  }
  decryptAndPrint((cout<<" after comparison: "), *minLvlCtxt, secKey, ea,0);
  cout << endl;
#endif
}


#if 0
e2=a2+b2+1
e1=a1+b1+1
e0=a0+b0+1, ne0 = a0+b0

e*2 = e2
e*1 = e2 e1
e*0 = e*1 e0, ne*0 = e*1 ne0

a*2 = a2    , b*2 = b2
a*1 = e2 a1 , b*1 = e*1 -e2 - a*1
a*0 = e*1 a0, b*0 = ne*0 -a*0

c2 = a2 b2
c1 = a1 b1
c0 = a0 b0

c*2 = c2
c*1 = e2 c1
c*0 = e*1 c0

A2 = a2,      B2 = b2,      C2 = c2
A1 = a2 +a*1, B1 = b2 +b*1, C1 = c2 + c*1
A0 = A1 +a*0, B0 = B1 +b*0, C0 = C1 + c*0

(a>b) = A0+C0
(a<b) = B0+C0

mx2 = a2+b2+c2
mx1 = a*1+b*1+c1 + a1(a2+c2) + b1(b2+c)

X ab01 = b1 a0
X b10 = b1 b0

mx0 = a0+b0+c*0
     + a0(A1+C1) = a0 A1 + a0(c*2 +a*1 b1) = a0(A1+c*2) + a*1 ab01
     + b0(B1+C1) = b0 B1 + b0(c*2 +a*1 b1) = b0(B1+c*2) + a*1 b10
    = a0(1+A1+c*2) + b0(1+B1+c*2) + a*1(b10+ab01) + c*0

mn0 = c0
     + a0(B1+C1) = a0 B1 + a0(c*2 +a*1 b1) = a0(B1+c*2) + a*1 ab01
     + b0(A1+C1) = b0 A1 + b0(c*2 +a*1 b1) = b0(A1+c*2) + a*1 b10
    = a0(B1+c*2) + b0(A1+C*2) + a*1(b10 + ab01) + c*0

e20 = e2 ne0

      a*1(b10+ab01)=(a2+b2+1)a1 b1(a0+b0) = e20 c1

mx0 = a0(1+A1+c*2) + b0(1+B1+c*2) + e20 c1 + c*0
mn0 = a0(B1+c*2) + b0(A1+c*2) + e20 c1 + c*0

#endif
