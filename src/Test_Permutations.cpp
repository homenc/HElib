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

void testIt(const Vec<GenDescriptor>& vec, long width);

void test1()
{
  Vec<GenDescriptor> vec(INIT_SIZE, 1);
  vec[0] = GenDescriptor(/*order=*/6, /*good=*/true, /*genIdx=*/0);
  cout << "**Testing (6,good), width=3\n";
  testIt(vec, /*width=*/3);
}

void test2()
{
  Vec<GenDescriptor> vec(INIT_SIZE, 2);
  vec[0] = GenDescriptor(/*order=*/3, /*good=*/true, /*genIdx=*/0);
  vec[1] = GenDescriptor(/*order=*/2, /*good=*/false, /*genIdx=*/1);
  cout << "**Testing [(3,good),(2,bad)], width=3\n";
  testIt(vec, /*width=*/3);
}

void test3()
{
  Vec<GenDescriptor> vec(INIT_SIZE, 1);
  vec[0] = GenDescriptor(/*order=*/31, /*good=*/true, /*genIdx=*/0);
  cout << "**Testing (31,good), width=5\n";
  testIt(vec, /*width=*/5);
}

void test4()
{
  Vec<GenDescriptor> vec(INIT_SIZE, 2);
  vec[0] = GenDescriptor(/*order=*/682,/*good=*/true, /*genIdx=*/0);
  vec[1] = GenDescriptor(/*order=*/ 2, /*good=*/false,/*genIdx=*/1);
  cout << "**Testing [(682,good),(2,bad)], width=11\n";
  testIt(vec, /*width=*/11);
}

void testIt(const Vec<GenDescriptor>& vec, long width)
{
  GeneratorTrees trees;
  trees.buildOptimalTrees(vec,width);
  cout << trees << endl;
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

  test1();
  test2();
  test3();
  test4();
}
