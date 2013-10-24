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

#include "NumbTh.h"
#include "PAlgebra.h"
#include "permutations.h"

void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'n1= L1=10 good1=0\n\n";
  cerr << "  n is [default=108]\n";
  cerr << "  L is [default=7]\n";
  cerr << "  good is the good-generator flag [default=0]\n";
  exit(0);
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

void testIt(Vec<GenDescriptor>& vec, long widthBound)
{
  GeneratorTrees trees;
  trees.buildOptimalTrees(vec, widthBound);
  Vec<long> dims;
  trees.getCubeDims(dims);
  CubeSignature sig(dims);

  for (long cnt=0; cnt<3; cnt++) {
    Permut pi;
    randomPerm(pi, trees.getSize());
    if (pi.length()<100)  cout << "pi="<<pi<<endl;

    PermNetwork net;
    net.buildNetwork(pi, trees);
    if (pi.length()<100) {
      cout << "permutations network {[gIdx,e,isID,shifts]} = " << endl;
      cout << net << endl;
    }

    HyperCube<long> cube1(sig), cube2(sig);
    for (long i=0; i<cube1.getSize(); i++) cube1[i] = i+1;
    HyperCube<long> cube3 = cube1;
    applyPermToVec(cube2.getData(), cube1.getData(), pi); // direct application
    net.applyToCube(cube3); // applying permutation netwrok
    if (cube2==cube3) cout << "yay\n\n";
    else {
      cout << "blech\n";
      cout << "in="<<cube1.getData() << endl;
      cout << "out1="<<cube2.getData()<<", out2="<<cube3.getData()<<endl<<endl;
    }
  }
}

int main(int argc, char *argv[])
{
  /*  argmap_t argmap;
  argmap["n"] = "108";
  argmap["L"] = "5";
  argmap["good"] = "0";
  if (!parseArgs(argc, argv, argmap)) {
    cerr << "bad args\n";
    exit(0);
  }
  long n = atoi(argmap["n"]);
  long L = atoi(argmap["L"]);
  bool good = !!atoi(argmap["good"]);
  */

#if 0
  // Test 1: a single good small prime-order generator (3)
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 1);
  vec[0] = GenDescriptor(/*order=*/3, /*good=*/true, /*genIdx=*/0);
  cout << "***Testing (3,good), width=1\n";
  testIt(vec, /*width=*/1);
  }

  // Test 2: a single bad larger prime-order generator (31)
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 1);
  vec[0] = GenDescriptor(/*order=*/31, /*good=*/false, /*genIdx=*/0);
  cout << "\n***Testing (31,bad), width=5\n";
  testIt(vec, /*width=*/5);
  }
#endif

  // Test 3: two generators with small prime orders (2,3), both bad
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 2);
  vec[0] = GenDescriptor(/*order=*/3, /*good=*/false, /*genIdx=*/0);
  vec[1] = GenDescriptor(/*order=*/2, /*good=*/false, /*genIdx=*/1);
  cout << "\n***Testing [(3,bad),(2,bad)], width=3\n";
  testIt(vec, /*width=*/3);
  }

#if 0
  // Test 4: two generators with small prime orders (2,3), one good
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 2);
  vec[0] = GenDescriptor(/*order=*/3, /*good=*/true, /*genIdx=*/0);
  vec[1] = GenDescriptor(/*order=*/2, /*good=*/false, /*genIdx=*/1);
  cout << "\n***Testing [(3,good),(2,bad)], width=3\n";
  testIt(vec, /*width=*/3);
  }

  // Test 5: a single bad composite-order generator (6)
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 1);
  vec[0] = GenDescriptor(/*order=*/6, /*good=*/false, /*genIdx=*/0);
  cout << "\n***Testing (6,bad), width=3\n";
  testIt(vec, /*width=*/3);
  }

  // Test 6: a single good composite-order generator (6)
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 1);
  vec[0] = GenDescriptor(/*order=*/6, /*good=*/true, /*genIdx=*/0);
  cout << "\n***Testing (6,good), width=3\n";
  testIt(vec, /*width=*/3);
  }

  // Test 7: the "general case", (682,good),(2,bad)
  {
  Vec<GenDescriptor> vec(INIT_SIZE, 2);
  vec[0] = GenDescriptor(/*order=*/682,/*good=*/true, /*genIdx=*/0);
  vec[1] = GenDescriptor(/*order=*/ 2, /*good=*/false,/*genIdx=*/1);
  cout << "**Testing [(682,good),(2,bad)], width=11\n";
  testIt(vec, /*width=*/11);
  }
#endif
}
