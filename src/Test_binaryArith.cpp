/* Copyright (C) 2012-2017 IBM Corp.
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
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <NTL/BasicThreadPool.h>
NTL_CLIENT

#include "EncryptedArray.h"
#include "FHE.h"

#include "intraSlot.h"
#include "binaryArith.h"

#ifdef DEBUG_PRINTOUT
#include "debugging.h"
#endif
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
  {  2, 27000, 32767, 15, 31,  7, 151, 11628, 28087,25824, 30,  6, -10, 28, 4}
};

void test15for4(FHESecKey& secKey);
void testAddTwo(FHESecKey& secKey, long bitSize, long outSize);
void testAddMany(FHESecKey& secKey, long bitSize, long outSize);
void testProduct(FHESecKey& secKey, long bitSize, long outSize);

int main(int argc, char *argv[])
{
  ArgMapping amap;
  long prm=1;
  amap.arg("prm", prm, "parameter size (0-tiny,...,3-huge)");
  long bitSize = 5;
  amap.arg("bitSize", bitSize, "bitSize of input integers (<=10)");
  long outSize = 0;
  amap.arg("outSize", outSize, "bitSize of output integers", "as many as needed");
  long nTests = 3;
  amap.arg("nTests", nTests, "number of tests to run");
  long seed=0;
  amap.arg("seed", seed, "PRG seed");
  amap.arg("verbose", verbose, "print more information");

  amap.parse(argc, argv);
  assert(prm >= 0 && prm < 4);
  if (seed) NTL::SetSeed(ZZ(seed));
  if (bitSize<=0) bitSize=5;
  else if (bitSize>10) bitSize=10;

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

  if (verbose)
    cout << "computing key-independent tables..." << std::flush;
  FHEcontext context(m, p, /*r=*/1, gens, ords);
  context.bitsPerLevel = B;
  buildModChain(context, /*L=*/bitSize+5, c,/*extraBits=*/8);
  context.makeBootstrappable(mvec, /*t=*/0,
                             /*flag=*/false, /*cacheType=DCRT*/2);
  buildUnpackSlotEncoding(unpackSlotEncoding, *context.ea);
  if (verbose) {
    cout << " done.\n";
    context.zMStar.printout();
    cout << "\ncomputing key-dependent tables..." << std::flush;
  }
  FHESecKey secKey(context);
  secKey.GenSecKey(/*Hweight=*/128);
  addSome1DMatrices(secKey); // compute key-switching matrices
  addFrbMatrices(secKey);
  secKey.genRecryptData();
  if (verbose) cout << " done\n";

  activeContext = &context; // make things a little easier sometimes
#ifdef DEBUG_PRINTOUT
  dbgEa = (EncryptedArray*) context.ea;
  dbgKey = &secKey;
#endif

  for (long i=0; i<nTests; i++)
    test15for4(secKey);

  /*
  for (long i=0; i<nTests; i++)
    testAddTwo(secKey, bitSize, outSize);

  for (long i=0; i<nTests; i++)
    testAddMany(secKey, bitSize, outSize);
  */
  for (long i=0; i<nTests; i++)
    testProduct(secKey, bitSize, outSize);

  if (verbose) printAllTimers(cout);
  cout << "=== all tests pass ===\n";
  return 0;
}

void test15for4(FHESecKey& secKey)
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
  long numOutputs
    = fifteenOrLess4Four(CtPtrs_vectorCt(outBuf), CtPtrs_vectorPt(inPtrs));

  // Check the result
  long sum2=0;
  for (int i=0; i<numOutputs; i++) {
    ZZX poly;
    secKey.Decrypt(poly, outBuf[i]);
    sum2 += to_long(ConstTerm(poly)) << i;
  }
  if (sum != sum2) {
    cout << "15to4 error: inputs="<<inputBits<<", sum="<<sum;
    cout << " but sum2="<<sum2<<endl;
    exit(0);
  }
  else if (verbose)
    cout << "15to4 succeeded, sum"<<inputBits<<"="<<sum2<<endl;
}

void testAddTwo(FHESecKey& secKey, long bitSize, long outSize)
{
}

void testAddMany(FHESecKey& secKey, long bitSize, long outSize)
{
}

void testProduct(FHESecKey& secKey, long bitSize, long outSize)
{
  const EncryptedArray& ea = *(secKey.getContext().ea);
  long mask = (outSize? ((1L<<outSize)-1) : -1);

  // Choose two random n-bit integers
  long pa = RandomBits_long(bitSize);
  long pb = RandomBits_long(bitSize);

  // Encrypt the individual bits
  NTL::Vec<Ctxt> eProduct, enca, encb;

  resize(enca, bitSize, Ctxt(secKey));
  resize(encb, bitSize, Ctxt(secKey));
  for (long i=0; i<bitSize; i++) {
    secKey.Encrypt(enca[i], ZZX((pa>>i)&1));
    secKey.Encrypt(encb[i], ZZX((pb>>i)&1));
  }

  // Test positive multiplication
  vector<long> slots;
  {CtPtrs_VecCt eep(eProduct);  // A wrappers around the output vector
  multTwoNumbers(eep,CtPtrs_VecCt(enca),CtPtrs_VecCt(encb),/*negative=*/false,
                 outSize, &unpackSlotEncoding);
  decryptBinaryNums(slots, eep, secKey, ea);
  } // get rid of the wrapper
  long pProd = pa*pb;
  if (slots[0] != ((pa*pb)&mask)) {
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
  secKey.Encrypt(encb[bitSize-1], ZZX(1));
  decryptBinaryNums(slots, CtPtrs_VecCt(encb), secKey, ea, /*negative=*/true);
  pb = slots[0];
  eProduct.kill();
  {CtPtrs_VecCt eep(eProduct);  // A wrappers around the output vector
  multTwoNumbers(eep,CtPtrs_VecCt(enca),CtPtrs_VecCt(encb),/*negative=*/true,
                 outSize, &unpackSlotEncoding);
  decryptBinaryNums(slots, eep, secKey, ea, /*negative=*/true);
  } // get rid of the wrapper
  pProd = pa*pb;
  if ((slots[0]&mask) != (pProd&mask)) {
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

#ifdef DEBUG_PRINTOUT
  const Ctxt* minCtxt = nullptr;
  long minLvl=1000;
  for (const Ctxt& c: eProduct) {
    long lvl = c.findBaseLevel();
    if (lvl < minLvl) {
      minCtxt = &c;
      minLvl = lvl;
    }
  }
  decryptAndPrint(cout, *minCtxt, secKey, ea,0);
#endif
}


#if 0 // OLD CODE
void testEncryptedBits(const FHEcontext& context, const FHESecKey& secKey)
{
  const int bitSize = 7;
  const int sizeLimit = 5;
  const long mask = (1L << sizeLimit) -1;

  long numInputs = 20+NTL::RandomBnd(16); // random # between 20 and 35
  vector<vector<Ctxt> > out;
  vector<vector<Ctxt> > inBuf(numInputs, vector<Ctxt>(bitSize,Ctxt(secKey)) );

  vector<long> inputs(numInputs, 0);
  vector<Ctxt*> inPtrs[35];
  long sum=0;
  for (int i=0; i<numInputs; i++) {
    resize(inPtrs[i], lsize(inBuf[i]), (Ctxt*)nullptr);
    for (int j=0; j<lsize(inPtrs[i]); j++) {
      if (NTL::RandomBnd(bitSize)>0) { // leave empty with small probability
        inPtrs[i][j] = &(inBuf[i][j]);
        long bit = NTL::RandomBnd(2);  // a random bit
        inputs[i] += (bit << j);
        secKey.Encrypt(inBuf[i][j], ZZX(bit));
      }
    }
    sum += inputs[i];
  }

  // Add these numbers

  // long numOutputs = fifteenOrLess4Four(out[0], out[1], out[2], out[3],
  //                                      inPtrs, numInputs, sizeLimit);

  std::vector<Ctxt> eSum;
  addManyNums(eSum, inPtrs, numInputs, sizeLimit);

  // Check the result
  long sum2 = 0;

  // vector<long> outputs(numOutputs, 0);
  // for (int i=0; i<numOutputs; i++) {
  //   vector<long> nums;
  //   decryptBinaryNums(nums, out[i]);
  //   outputs[i] = nums[0];
  //   sum2 += outputs[i];
  // }

  vector<long> nums;
  decryptBinaryNums(nums, eSum);
  sum2 = nums[0];

  if ((sum&mask) != (sum2&mask)) {
    cout << "sum error: inputs="<<inputs<<", sum="<<sum<<endl;
    //    cout << "    but outputs="<<outputs<<", sum2="<<sum2<<endl;
    cout << "    but sum2="<<sum2<<endl;
    printLevelsVV(cout, out);
    exit(0);
  }
  // else
  //   cout << "success, inputs="<<inputs<<", output="<<outputs<<endl;
}
#endif
