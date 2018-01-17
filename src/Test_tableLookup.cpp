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
// Test_tableLookup.cpp - test the table-lookup module
#include <iostream>
#include <NTL/BasicThreadPool.h>
#include "intraSlot.h"
#include "tableLookup.h"

#ifdef DEBUG_PRINTOUT
#include "debugging.h"
#endif

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

// test function declerations
void testLookup(const FHESecKey& sKey, long insize, long outsize);
void testWritein(const FHESecKey& sKey, long insize, long nTests);


int main(int argc, char *argv[])
{
  ArgMapping amap;
  long prm=1;
  amap.arg("prm", prm, "parameter size (0-tiny,...,4-huge)");
  long bitSize = 5;
  amap.arg("bitSize", bitSize, "bitSize of input integers (<=32)");
  long outSize = 0;
  amap.arg("outSize", outSize, "bitSize of output integers", "as many as needed");
  long nTests = 3;
  amap.arg("nTests", nTests, "number of tests to run");
  bool bootstrap = false;
  amap.arg("bootstrap", bootstrap, "test multiplication with bootstrapping");
  long seed=0;
  amap.arg("seed", seed, "PRG seed");
  long nthreads=1;
  amap.arg("nthreads", nthreads, "number of threads");
  amap.arg("verbose", verbose, "print more information");

  amap.parse(argc, argv);
  assert(prm >= 0 && prm < 5);
  assert(bitSize<=7);
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
  if (bootstrap) L=30; // that should be enough
  else           L = 3 +bitSize;
  
  if (verbose) {
    cout <<"input bitSize="<<bitSize<<", output size bound="<<outSize
         <<", running "<<nTests<<" tests for each function\n";
    if (nthreads>1) cout << "  using "<<NTL::AvailableThreads()<<" threads\n";
    cout << "computing key-independent tables..." << std::flush;
  }
  FHEcontext context(m, p, /*r=*/1, gens, ords);
  context.bitsPerLevel = B;
  buildModChain(context, L, c,/*extraBits=*/8);
  if (bootstrap) {
    context.makeBootstrappable(mvec, /*t=*/0,
                               /*flag=*/false, /*cacheType=DCRT*/2);
  }
  buildUnpackSlotEncoding(unpackSlotEncoding, *context.ea);
  if (verbose) {
    cout << " done.\n";
    context.zMStar.printout();
    cout << " L="<<L<<", B="<<B<<endl;
    cout << "\ncomputing key-dependent tables..." << std::flush;
  }
  FHESecKey secKey(context);
  secKey.GenSecKey(/*Hweight=*/128);
  addSome1DMatrices(secKey); // compute key-switching matrices
  addFrbMatrices(secKey);
  if (bootstrap) secKey.genRecryptData();
  if (verbose) cout << " done\n";

  activeContext = &context; // make things a little easier sometimes
#ifdef DEBUG_PRINTOUT
  dbgEa = (EncryptedArray*) context.ea;
  dbgKey = &secKey;
#endif

  testLookup(secKey, bitSize, outSize);
  cout << "  *** testLookup PASS ***\n";

  testWritein(secKey, bitSize, nTests);
  cout << "  *** testWritein PASS ***\n";

  if (verbose) printAllTimers(cout);
  return 0;
}

static void encryptIndex(std::vector<Ctxt>& ei, long index,
                         const FHESecKey& sKey)
{
  for (long i=0; i<lsize(ei); i++)
    sKey.Encrypt(ei[i], to_ZZX((index>>i) &1)); // i'th bit of index
}

static long decryptIndex(std::vector<Ctxt>& ei, const FHESecKey& sKey)
{
  long num=0;
  for (long i=0; i<lsize(ei); i++) {
    ZZX poly;
    sKey.Decrypt(poly, ei[i]);
    num += to_long(NTL::ConstTerm(poly)) << i;
  }
  return num;
}

void testLookup(const FHESecKey& sKey, long insize, long outsize)
{
  // Build a table s.t. T[i] = 2^{outsize -1}/(i+1), i=0,...,2^insize -1
  std::vector<zzX> T;
  buildLookupTable(T, [](double x){ return 1/(x+1.0);},
                   insize,  /*scale_in=*/0, /*sign_in=*/0,
                   outsize, /*scale_out=*/1-outsize, /*sign_out=*/0,
                   *(sKey.getContext().ea));

  assert(lsize(T)==(1L<<insize));
  for (long i=0; i<lsize(T); i++) {
    Ctxt c(sKey);
    std::vector<Ctxt> ei(insize, c);
    encryptIndex(ei, i, sKey);              // encrypt the index
    tableLookup(c, T, CtPtrs_vectorCt(ei)); // get the encrypted entry
    // decrypt and compare
    ZZX poly;  sKey.Decrypt(poly, c); // decrypt
    zzX poly2; convert(poly2, poly);  // convert to zzX
    if (poly2 != T[i]) {
      cout << "testLookup error: decrypted T["<<i<<"]\n";
      exit(0);
    }
  }
}

void testWritein(const FHESecKey& sKey, long size, long nTests)
{
  const EncryptedArray& ea = *(sKey.getContext().ea);
  long tSize = 1L << size; // table size

  // encrypt a random table
  vector<long> pT(tSize, 0);         // plaintext table
  vector<Ctxt> T(tSize, Ctxt(sKey)); // encrypted table
  for (long i=0; i<size; i++) {
    long bit = RandomBits_long(1);   // a random bit
    sKey.Encrypt(T[i], to_ZZX(bit));
    pT[i] = bit;
  }

  // Add 1 to 20 random entries in the table
  for (long count=0; count<nTests; count++) {
    // encrypt a random index into the table
    long index = RandomBnd(tSize); // 0 <= index < tSize
    vector<Ctxt> I(size, Ctxt(sKey));
    encryptIndex(I, index, sKey);

    // do the table write-in
    tableWriteIn(CtPtrs_vectorCt(T), CtPtrs_vectorCt(I), &unpackSlotEncoding);
    pT[index]++; // add 1 to entry 'index' in the plaintext table
  }

  // Check that the ciphertext and plaintext tables still match
  for (int i=0; i<tSize; i++) {
    ZZX poly;
    sKey.Decrypt(poly, T[i]);
    long decrypted = to_long(NTL::ConstTerm(poly));
    long p = T[i].getPtxtSpace();
    if ((pT[i] - decrypted) % p) { // not equal mod p
      cout << "testWritein error: decrypted T["<<i<<"]="<<decrypted
           <<" but should be "<<pT[i]<<" (mod "<<p<<")\n";
      exit(0);
    }
  }
}
