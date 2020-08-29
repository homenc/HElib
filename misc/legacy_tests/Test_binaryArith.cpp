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
  {  2,    48,   105, 12,   3, 35,  0,    71,    76,    0,   2,  2,   0, 25, 2},
  {  2 ,  600,  1023, 10,  11, 93,  0,   838,   584,    0,  10,  6,   0, 25, 2},
  {  2,  2304,  4641, 24,   7,  3,221,  3979,  3095, 3760,   6,  2,  -8, 25, 3},
  {  2,  5460,  8193, 26,8193,  0,  0,    46,     0,    0, 210,  0,   0, 25, 3},
  {  2,  8190,  8191, 13,8191,  0,  0,    39,     0,    0, 630,  0,   0, 25, 3},
  {  2, 10752, 11441, 48,  17,673,  0,  4712,  2024,    0,  16,-14,   0, 25, 3},
  {  2, 15004, 15709, 22,  23,683,  0,  4099, 13663,    0,  22, 31,   0, 25, 3},
  {  2, 27000, 32767, 15,  31,  7,151, 11628, 28087,25824,  30,  6, -10, 28, 4}
};

void test15for4(SecKey& secKey);
void testProduct(SecKey& secKey, long bitSize1, long bitSize2,
                 long outSize, bool bootstrap = false);
void testAdd(SecKey& secKey, long bitSize1, long bitSize2,
             long outSize, bool bootstrap = false);

int main(int argc, char *argv[])
{
  ArgMap amap;
  long prm=1;
  amap.arg("prm", prm, "parameter size (0-tiny,...,7-huge)");
  long bitSize = 5;
  amap.arg("bitSize", bitSize, "bitSize of input integers (<=32)");
  long bitSize2 = 0;
  amap.arg("bitSize2", bitSize2, "bitSize of 2nd input integer (<=32)",
           "same as bitSize");
  long outSize = 0;
  amap.arg("outSize", outSize, "bitSize of output integers", "as many as needed");
  long nTests = 2;
  amap.arg("nTests", nTests, "number of tests to run");
  bool bootstrap = false;
  amap.arg("bootstrap", bootstrap, "test multiplication with bootstrapping");
  long seed=0;
  amap.arg("seed", seed, "PRG seed");
  long nthreads=1;
  amap.arg("nthreads", nthreads, "number of threads");
  amap.arg("verbose", verbose, "print more information");

  long tests2avoid = 1;
  amap.arg("tests2avoid", tests2avoid, "bitmap of tests to disable (1-15for4, 2-add, 4-multiply");

  amap.parse(argc, argv);
  assert(prm >= 0 && prm < 5);
  if (seed) NTL::SetSeed(ZZ(seed));
  if (nthreads>1) NTL::SetNumThreads(nthreads);

  if (bitSize<=0) bitSize=5;
  else if (bitSize>32) bitSize=32;
  if (bitSize2<=0) bitSize2=bitSize;
  else if (bitSize2>32) bitSize2=32;

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
  if (bootstrap) L=900; // that should be enough
  else {
    double nBits =
      (outSize>0 && outSize<2*bitSize)? outSize : (2*bitSize);
    double three4twoLvls = log(nBits/2) / log(1.5);
    double add2NumsLvls = log(nBits) / log(2.0);
    L = (5 + ceil(three4twoLvls + add2NumsLvls))*30;
  }

  if (verbose) {
    cout <<"input bitSizes="<<bitSize<<','<<bitSize2
         <<", output size bound="<<outSize
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

  if (!(tests2avoid & 1)) {
    for (long i=0; i<nTests; i++)
      test15for4(secKey);
    cout << "GOOD\n";
  }
  if (!(tests2avoid & 2)) {
    for (long i=0; i<nTests; i++)
      testAdd(secKey, bitSize, bitSize2, outSize, bootstrap);
    cout << "GOOD\n";
  }
  if (!(tests2avoid & 4)) {
    for (long i=0; i<nTests; i++)
      testProduct(secKey, bitSize, bitSize2, outSize, bootstrap);
    cout << "GOOD\n";
  }
  if (verbose) printAllTimers(cout);
  return 0;
}

void test15for4(SecKey& secKey)
{
  std::vector<Ctxt> inBuf(15, Ctxt(secKey));
  std::vector<Ctxt*> inPtrs(15, nullptr);

  std::vector<Ctxt> outBuf(5, Ctxt(secKey));

  long sum=0;
  std::string inputBits = "(";
  for (int i=0; i<15; i++) {
    if (NTL::RandomBnd(10)>0) { // leave empty with small probability
      inPtrs[i] = &(inBuf[i]);
      long bit = NTL::RandomBnd(2);  // a random bit
      secKey.Encrypt(inBuf[i], ZZX(bit));
      inputBits += std::to_string(bit) + ",";
      sum += bit;
    }
    else inputBits += "-,";
  }
  inputBits += ")";

  // Add these bits
  if (verbose) {
    cout << endl;
    CheckCtxt(inBuf[lsize(inBuf)-1], "b4 15for4");
  }
  long numOutputs
    = fifteenOrLess4Four(CtPtrs_vectorCt(outBuf), CtPtrs_vectorPt(inPtrs));
  if (verbose)
    CheckCtxt(outBuf[lsize(outBuf)-1], "after 15for4");

  // Check the result
  long sum2=0;
  for (int i=0; i<numOutputs; i++) {
    ZZX poly;
    secKey.Decrypt(poly, outBuf[i]);
    sum2 += to_long(ConstTerm(poly)) << i;
  }
  if (sum != sum2) {
    cout << "BAD\n";
    if (verbose) {
      cout << "  15to4: inputs="<<inputBits<<", sum="<<sum
           << " but sum2="<<sum2<<endl;
    }
    exit(0);
  }
  else if (verbose)
    cout << "15to4 succeeded, sum"<<inputBits<<"="<<sum2<<endl;
}

void testProduct(SecKey& secKey, long bitSize, long bitSize2,
                 long outSize, bool bootstrap)
{
  const Context& context = secKey.getContext();
  const EncryptedArray& ea = *(context.ea);
  long mask = (outSize? ((1L<<outSize)-1) : -1);

  // Choose two random n-bit integers
  long pa = RandomBits_long(bitSize);
  long pb = RandomBits_long(bitSize2);

  // Encrypt the individual bits
  NTL::Vec<Ctxt> eProduct, enca, encb;

  resize(enca, bitSize, Ctxt(secKey));
  for (long i=0; i<bitSize; i++) {
    secKey.Encrypt(enca[i], ZZX((pa>>i)&1));
    if (bootstrap) { // put them at a lower level
      enca[i].bringToSet(context.getCtxtPrimes(5));
    }
  }
  resize(encb, bitSize2, Ctxt(secKey));
  for (long i=0; i<bitSize2; i++) {
    secKey.Encrypt(encb[i], ZZX((pb>>i)&1));
    if (bootstrap) { // put them at a lower level
      encb[i].bringToSet(context.getCtxtPrimes(5));
    }
  }
  if (verbose) {
    cout << "\n  bits-size "<<bitSize<<'+'<<bitSize2;
    if (outSize>0) cout << "->"<<outSize;
    CheckCtxt(encb[0], "b4 multiplication");
  }
  // Test positive multiplication
  vector<long> slots;
  {CtPtrs_VecCt eep(eProduct);  // A wrappers around the output vector
  multTwoNumbers(eep,CtPtrs_VecCt(enca),CtPtrs_VecCt(encb),/*negative=*/false,
                 outSize, &unpackSlotEncoding);
  decryptBinaryNums(slots, eep, secKey, ea);
  } // get rid of the wrapper
  if (verbose)
    CheckCtxt(eProduct[lsize(eProduct)-1], "after multiplication");
  long pProd = pa*pb;
  if (slots[0] != ((pa*pb)&mask)) {
    cout << "BAD\n";
    if (verbose)
      cout << "Positive product error: pa="<<pa<<", pb="<<pb
           << ", but product="<<slots[0]
           << " (should be "<<pProd<<'&'<<mask<<'='<<(pProd&mask)<<")\n";
    exit(0);
  }
  else if (verbose) {
    cout << "positive product succeeded: ";
    if (outSize) cout << "bottom "<<outSize<<" bits of ";
    cout << pa<<"*"<<pb<<"="<<slots[0]<<endl;
  }
  // Test negative multiplication
  secKey.Encrypt(encb[bitSize2-1], ZZX(1));
  decryptBinaryNums(slots, CtPtrs_VecCt(encb), secKey, ea, /*negative=*/true);
  pb = slots[0];
  eProduct.kill();
  {CtPtrs_VecCt eep(eProduct);  // A wrappers around the output vector
  multTwoNumbers(eep,CtPtrs_VecCt(enca),CtPtrs_VecCt(encb),/*negative=*/true,
                 outSize, &unpackSlotEncoding);
  decryptBinaryNums(slots, eep, secKey, ea, /*negative=*/true);
  } // get rid of the wrapper
  if (verbose)
    CheckCtxt(eProduct[lsize(eProduct)-1], "after multiplication");
  pProd = pa*pb;
  if ((slots[0]&mask) != (pProd&mask)) {
    cout << "BAD\n";
    if (verbose)
      cout << "Negative product error: pa="<<pa<<", pb="<<pb
           << ", but product="<<slots[0]
           << " (should be "<<pProd<<'&'<<mask<<'='<<(pProd&mask)<<")\n";
    exit(0);
  }
  else if (verbose) {
    cout << "negative product succeeded: ";
    if (outSize) cout << "bottom "<<outSize<<" bits of ";
    cout << pa<<"*"<<pb<<"="<<slots[0]<<endl;
  }

#ifdef HELIB_DEBUG
  const Ctxt* minCtxt = nullptr;
  long minLvl=1000;
  for (const Ctxt& c: eProduct) {
    long lvl = c.logOfPrimeSet();
    if (lvl < minLvl) {
      minCtxt = &c;
      minLvl = lvl;
    }
  }
  decryptAndPrint((cout<<" after multiplication: "), *minCtxt, secKey, ea,0);
  cout << endl;
#endif
}


void testAdd(SecKey& secKey, long bitSize1, long bitSize2,
             long outSize, bool bootstrap)
{
  const Context& context = secKey.getContext();
  const EncryptedArray& ea = *(context.ea);
  long mask = (outSize? ((1L<<outSize)-1) : -1);

  // Choose two random n-bit integers
  long pa = RandomBits_long(bitSize1);
  long pb = RandomBits_long(bitSize2);

  // Encrypt the individual bits
  NTL::Vec<Ctxt> eSum, enca, encb;

  resize(enca, bitSize1, Ctxt(secKey));
  for (long i=0; i<bitSize1; i++) {
    secKey.Encrypt(enca[i], ZZX((pa>>i)&1));
    if (bootstrap) { // put them at a lower level
      enca[i].bringToSet(context.getCtxtPrimes(5));
    }
  }
  resize(encb, bitSize2, Ctxt(secKey));
  for (long i=0; i<bitSize2; i++) {
    secKey.Encrypt(encb[i], ZZX((pb>>i)&1));
    if (bootstrap) { // put them at a lower level
      encb[i].bringToSet(context.getCtxtPrimes(5));
    }
  }
  if (verbose) {
    cout << "\n  bits-size "<<bitSize1<<'+'<<bitSize2;
    if (outSize>0) cout << "->"<<outSize;
    cout <<endl;
    CheckCtxt(encb[0], "b4 addition");
  }

  // Test addition
  vector<long> slots;
  {CtPtrs_VecCt eep(eSum);  // A wrapper around the output vector
  addTwoNumbers(eep, CtPtrs_VecCt(enca), CtPtrs_VecCt(encb),
                outSize, &unpackSlotEncoding);
  decryptBinaryNums(slots, eep, secKey, ea);
  } // get rid of the wrapper
  if (verbose) CheckCtxt(eSum[lsize(eSum)-1], "after addition");
  long pSum = pa+pb;
  if (slots[0] != ((pa+pb)&mask)) {
    cout << "BAD\n";
    if (verbose)
      cout << "addTwoNums error: pa="<<pa<<", pb="<<pb
           << ", but pSum="<<slots[0]
           << " (should be ="<<(pSum&mask)<<")\n";
    exit(0);
  }
  else if (verbose) {
    cout << "addTwoNums succeeded: ";
    if (outSize) cout << "bottom "<<outSize<<" bits of ";
    cout << pa<<"+"<<pb<<"="<<slots[0]<<endl;
  }

#ifdef HELIB_DEBUG
  const Ctxt* minCtxt = nullptr;
  long minLvl=1000;
  for (const Ctxt& c: eSum) {
    long lvl = c.logOfPrimeSet();
    if (lvl < minLvl) {
      minCtxt = &c;
      minLvl = lvl;
    }
  }
  decryptAndPrint((cout<<" after addition: "), *minCtxt, secKey, ea,0);
  cout << endl;
#endif
}
