/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

/* Test_Permutations.cpp - Applying plaintext permutation to encrypted vector
 */
#include <NTL/ZZ.h>
NTL_CLIENT

#include "NumbTh.h"
#include "timing.h"
#include "permutations.h"
#include "EncryptedArray.h"

void testCtxt(long m, long p, long widthBound=0, long L=0, long r=1);

void usage(char *prog) 
{
  cout << "Usage: "<<prog<<" [test=? [optional parameters...]]\n";
  cout << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cout << "  e.g, 'test=1 m=108 p=2 r=1\n";
  cout << "  test is either 0 (plaintext) or 1 (ciphertext)[default=1]\n\n";
  cout << "test=0, permuting plaintext hypercubes (dimension upto 4):\n";
  cout << "  ord1,ord2,ord3,ord4 size of dimensions 1..4 [default ord1=30, ord2,3,4=0]\n";
  cout << "  good1,good2,good3,good4 native rotation flags (0/1) [default=1]\n";
  cout << "  depth bounds the depth of permutation network [default=5]\n";
  cout << "\ntest=1, permuting ciphertext slots:\n";
  cout << "  m is the cyclotomic field [default=4369]\n";
  cout << "  p,r define the plaintext space p^r [default p=2,r=1]\n";
  cout << "  depth bounds the depth of permutation network [default=5]\n";
  cout << "  L is number of levels in chain [default=depth]\n\n";
  cout << "dry=1 for dry run [default=0]\n";
  exit(0);
}

void testCube(Vec<GenDescriptor>& vec, long widthBound)
{
  GeneratorTrees trees;
  long cost = trees.buildOptimalTrees(vec, widthBound);
  cout << "@TestCube: trees=" << trees << endl;
  cout << " cost =" << cost << endl;
  Vec<long> dims;
  trees.getCubeDims(dims);
  CubeSignature sig(dims);

  for (long cnt=0; cnt<3; cnt++) {
    Permut pi;
    randomPerm(pi, trees.getSize());
    //    if (pi.length()<100)  cout << "pi="<<pi<<endl;

    PermNetwork net;
    net.buildNetwork(pi, trees);
    //    if (pi.length()<100) {
    //      cout << "permutations network {[gIdx,e,isID,shifts]} = " << endl;
    //      cout << net << endl;
    //    }

    HyperCube<long> cube1(sig), cube2(sig);
    for (long i=0; i<cube1.getSize(); i++) cube1[i] = i;
    HyperCube<long> cube3 = cube1;
    applyPermToVec(cube2.getData(), cube1.getData(), pi); // direct application
    net.applyToCube(cube3); // applying permutation netwrok
    if (cube2==cube3) cout << "GOOD\n";
    else {
      cout << "BAD\n";
      if (cube1.getSize()<100) {
	cout << "in="<<cube1.getData() << endl;
	cout << "out1="<<cube2.getData()<<", out2="
	     << cube3.getData()<<endl<<endl;
      }
    }
  }
}

void testCtxt(long m, long p, long widthBound, long L, long r)
{
  cout << "@testCtxt(m="<<m<<",p="<<p<<",depth="<<widthBound<< ",r="<<r<<")";

  FHEcontext context(m,p,r);
  EncryptedArray ea(context); // Use G(X)=X for this ea object

  // Some arbitrary initial plaintext array
  vector<long> in(ea.size());
  for (long i=0; i<ea.size(); i++) in[i] = i % p;

  // Setup generator-descriptors for the PAlgebra generators
  Vec<GenDescriptor> vec(INIT_SIZE, ea.dimension());
  for (long i=0; i<ea.dimension(); i++)
    vec[i] = GenDescriptor(/*order=*/ea.sizeOfDimension(i),
			   /*good=*/ ea.nativeDimension(i), /*genIdx=*/i);

  // Some default for the width-bound, if not provided
  if (widthBound<=0) widthBound = 1+log2((double)ea.size());

  // Get the generator-tree structures and the corresponding hypercube
  GeneratorTrees trees;
  long cost = trees.buildOptimalTrees(vec, widthBound);
  cout << ": trees=" << trees << endl;
  cout << " cost =" << cost << endl;

  //  Vec<long> dims;
  //  trees.getCubeDims(dims);
  //  CubeSignature sig(dims);

  // 1/2 prime per level should be more or less enough, here we use 1 per layer
  if (L<=0) L = 1+trees.numLayers();
  buildModChain(context, /*nLevels=*/L, /*nDigits=*/3);
  cout << "**Using "<<L<<" primes (of which "
       << context.ctxtPrimes.card() << " are Ctxt-primes)\n";

  // Generate a sk/pk pair
  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(64); // A Hamming-weight-64 secret key
  Ctxt ctxt(publicKey);

  for (long cnt=0; cnt<3; cnt++) {
    resetAllTimers();
    // Choose a random permutation
    Permut pi;
    randomPerm(pi, trees.getSize());
    //    if (pi.length()<100)  cout << "pi="<<pi<<endl;

    // Build a permutation network for pi
    PermNetwork net;
    net.buildNetwork(pi, trees);
    //    if (pi.length()<100) {
    //      cout << "permutations network {[gIdx,e,isID,shifts]} = " << endl;
    //      cout << net << endl;
    //    }

    // make sure we have the key-switching matrices needed for this network
    cout << "  ** generating matrices... " << flush;
    addMatrices4Network(secretKey, net);
    cout << "done" << endl;

    // Apply the permutation pi to the plaintext
    vector<long> out1(ea.size());
    vector<long> out2(ea.size());
    applyPermToVec(out1, in, pi); // direct application

    // Encrypt plaintext array, then apply permutation network to ciphertext
    ea.encrypt(ctxt, publicKey, in);
    cout << "  ** applying permutation network to ciphertext... " << flush;
    double t = GetTime();
    net.applyToCtxt(ctxt, ea); // applying permutation netwrok
    t = GetTime() -t;
    cout << "done in " << t << " seconds" << endl;
    ea.decrypt(ctxt, secretKey, out2);

    if (out1==out2) cout << "GOOD\n";
    else {
      cout << "************ BAD\n";
    }
    // printAllTimers();
  }
}


/* m = 31, p = 2, phi(m) = 30
  ord(p)=5
  generator 6 has order (== Z_m^*) of 6
  T = [1 6 5 30 25 26 ]

  m = 61, p = 3, phi(m) = 60
  ord(p)=10
  generator 13 has order (== Z_m^*) of 3
  generator 2 has order (!= Z_m^*) of 2
  T = [1 2 13 26 47 33 ]

  m = 683, p = 2, phi(m) = 682
  ord(p)=22
  generator 3 has order (== Z_m^*) of 31

  m = 47127, p = 2, phi(m) = 30008
  ord(p)=22
  generator 5 has order (== Z_m^*) of 682
  generator 13661 has order (== Z_m^*) of 2
*/

int main(int argc, char *argv[])
{
  argmap_t argmap;
  argmap["test"] = "1";
  argmap["m"] = "4369";
  argmap["p"] = "2";
  argmap["r"] = "1";
  argmap["depth"] = "5";
  argmap["L"] = "0";
  argmap["ord1"] = "30";
  argmap["ord2"] = "0";
  argmap["ord3"] = "0";
  argmap["ord4"] = "0";
  argmap["good1"] = "1";
  argmap["good2"] = "1";
  argmap["good3"] = "1";
  argmap["good4"] = "1";
  argmap["dry"] = "1";

  // get parameters from the command line

  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long test = atoi(argmap["test"]);
  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long m = atoi(argmap["m"]);
  long depth = atoi(argmap["depth"]);
  long L = atoi(argmap["L"]);

  long ord1 = atoi(argmap["ord1"]);
  long ord2 = atoi(argmap["ord2"]);
  long ord3 = atoi(argmap["ord3"]);
  long ord4 = atoi(argmap["ord4"]);
  long good1 = atoi(argmap["good1"]);
  long good2 = atoi(argmap["good2"]);
  long good3 = atoi(argmap["good3"]);
  long good4 = atoi(argmap["good4"]);

  bool dry = atoi(argmap["dry"]);

  setDryRun(dry);
  if (test==0) {
    Vec<GenDescriptor> vec;
    long nGens;
    if (ord2<=1) nGens=1;
    else if (ord3<=1) nGens=2;
    else if (ord4<=1) nGens=3;
    else nGens=4;
    vec.SetLength(nGens);

    switch (nGens) {
    case 4:  vec[3] = GenDescriptor(ord4, good4, /*genIdx=*/3);
    case 3:  vec[2] = GenDescriptor(ord3, good3, /*genIdx=*/2);
    case 2:  vec[1] = GenDescriptor(ord2, good2, /*genIdx=*/1);
    default: vec[0] = GenDescriptor(ord1, good1, /*genIdx=*/0);
    }
    cout << "***Testing ";
    if (isDryRun()) cout << "(dry run) ";
    for (long i=0; i<vec.length(); i++)
      cout << "("<<vec[i].order<<","<<vec[i].good<<")";
    cout << ", depth="<<depth<<"\n";
    testCube(vec, depth);
  }
  else {
    setTimersOn();
    cout << "***Testing m="<<m<<", p="<<p<<", depth="<<depth<< endl;
    testCtxt(m,p,depth,L,r);
  }
}



#if 0
  cout << "***Testing m=31, p=2, width=3\n"; // (6 good)
  testCtxt(/*m=*/31, /*p=*/2, /*width=*/3);

  cout << "\n***Testing m=61, p=3, width=3\n"; // (3 good), (2, bad)
  testCtxt(/*m=*/61, /*p=*/3, /*width=*/3);

  cout << "\n***Testing m=683, p=2, width=5\n"; // (31, good)
  testCtxt(/*m=*/683, /*p=*/2, /*width=*/5);

  //  cout << "\n***Testing m=47127, p=2, width=11\n"; // (682,good),(2,good)
  //  testCtxt(/*m=*/47127, /*p=*/2, /*width=*/11);

  // Test 1: a single good small prime-order generator (3)
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 1);
  vec[0] = GenDescriptor(/*order=*/3, /*good=*/true, /*genIdx=*/0);
  cout << "***Testing (3,good), width=1\n";
  testCube(vec, /*width=*/1);
  }

  // Test 2: a single bad larger prime-order generator (31)
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 1);
  vec[0] = GenDescriptor(/*order=*/31, /*good=*/false, /*genIdx=*/0);
  cout << "\n***Testing (31,bad), width=5\n";
  testCube(vec, /*width=*/5);
  }

  // Test 3: two generators with small prime orders (2,3), both bad
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 2);
  vec[0] = GenDescriptor(/*order=*/2, /*good=*/false, /*genIdx=*/0);
  vec[1] = GenDescriptor(/*order=*/3, /*good=*/false, /*genIdx=*/1);
  cout << "\n***Testing [(2,bad),(3,bad)], width=3\n";
  testCube(vec, /*width=*/3);
  }

  // Test 4: two generators with small prime orders (2,3), one good
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 2);
  vec[0] = GenDescriptor(/*order=*/3, /*good=*/true, /*genIdx=*/0);
  vec[1] = GenDescriptor(/*order=*/2, /*good=*/false, /*genIdx=*/1);
  cout << "\n***Testing [(3,good),(2,bad)], width=3\n";
  testCube(vec, /*width=*/3);
  }

  // Test 5: a single good composite-order generator (6)
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 1);
  vec[0] = GenDescriptor(/*order=*/6, /*good=*/true, /*genIdx=*/0);
  cout << "\n***Testing (6,good), width=3\n";
  testCube(vec, /*width=*/3);
  }

  // Test 6: (6,good),(2,bad)
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 2);
  vec[0] = GenDescriptor(/*order=*/6,/*good=*/true, /*genIdx=*/0);
  vec[1] = GenDescriptor(/*order=*/ 2, /*good=*/false,/*genIdx=*/1);
  cout << "\n**Testing [(6,good),(2,bad)], width=5\n";
  testCube(vec, /*width=*/5);
  }

  // Test 7: the "general case", (682,good),(2,bad)
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 2);
  vec[0] = GenDescriptor(/*order=*/682,/*good=*/true, /*genIdx=*/0);
  vec[1] = GenDescriptor(/*order=*/ 2, /*good=*/false,/*genIdx=*/1);
  cout << "\n**Testing [(682,good),(2,bad)], width=11\n";
  testCube(vec, /*width=*/11);
  }
#endif
