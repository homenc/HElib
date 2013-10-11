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
#include "PAlgebra.h"
#include "permutations.h"

// Type-1 is a 2x3 good tree, type-2 is a 2(good)x3(good)x2(bad)
void BuildTree(OneGeneratorTree& T, long type)
{
  T.CollapseToRoot(); // ensure the tree has nothing except perhaps the root

  long nodeIdx = 0;
  long size = 6;
  SubDimension n6(/*genIdx=*/0, /*size=*/6, /*e=*/1, /*good=*/true);
  SubDimension n3(/*genIdx=*/0, /*size=*/3, /*e=*/0, /*good=*/true);
  SubDimension n2g(/*genIdx=*/0, /*size=*/2, /*e=*/0, /*good=*/true);
  if (type == 2) {
    size = 12;
    SubDimension n12(/*genIdx=*/0, /*size=*/2, /*e=*/1, /*good=*/true);
    SubDimension n2b(/*genIdx=*/0, /*size=*/2, /*e=*/0, /*good=*/false);
    T.PutDataInRoot(n12);
    nodeIdx = T.AddChildren(0, n6, n2b);
  } else
    T.PutDataInRoot(n6);

  T.AddChildren(nodeIdx, n2g, n3);
  computeEvalues(T, 0, size);
}

void TestIt(Vec<OneGeneratorTree>& ts)
{
  GeneratorTrees trees(ts);
  trees.ComputeCubeMapping();

  Permut p;
  RandomPerm(p, trees.getSize());
  PermNetwork ntwrk;
  ntwrk.BuildNetwork(p, trees);
}

int main()
{
  // Build a "good" 2x3 tree for the order-6 generator g=5 in Z_31^*/(2)
  // PAlgebra al(/*m=*/31, /*p=*/2);
  // Build a 2x3x3 tree for the order-18 generator g=24 in Z_127^*/(2)
  // PAlgebra al(/*m=*/127, /*p=*/2);
  // SubDimension rtDim(/*genIdx=*/0, /*size=*/18, /*e=*/1, /*good=*/true);

  // Test #1: a single 2x3 tree
  {
  Vec<OneGeneratorTree> trees(INIT_SIZE, 1);
  BuildTree(trees[0], 1);
  TestIt(trees);
  }
  // Test #2: a single 2x3x2 tree
  {
  Vec<OneGeneratorTree> trees(INIT_SIZE, 1);
  BuildTree(trees[0], 2);
  TestIt(trees);
  }
  // Test #3: both trees
  {
  Vec<OneGeneratorTree> trees(INIT_SIZE, 2);
  BuildTree(trees[0], 2);
  BuildTree(trees[1], 1);
  TestIt(trees);
  }

#if 0
  Vec<long> dims(INIT_SIZE,3);
  {
    long dd[3] = {2,3,2};
    for (long i=0; i<3; i++) dims[i] = dd[i];
  }
  CubeSignature cs(dims);

#define psize 12
  Permut p(INIT_SIZE, psize);
  { long pp[psize] = { 0, 3, 9, 4, 2, 7, 8, 11, 1, 6, 5, 10 };
    for (long i=0; i<psize; i++) p[i] = pp[i];
  }
  std::cout << "input="<<p<<endl<<endl;

  // break it into three permutations over the colums and rows
  vector<ColPerm> out;
  breakPermByDim(out,p,cs);

  for (long i=0; i<(long)out.size(); i++) {
    out[i].printout(std::cout);
    long nPerms = out[i].getSize()/dims[out[i].getPermDim()];
    for (long j=0; j<nPerms; j++) {
      out[i].extractSlice(p,j);
      std::cout <<"  "<<j<<": "<<p<<endl;
    }
    out[i].makeExplicit(p);
    std::cout << " expl: " << p << endl;
    out[i].getShiftAmounts(p);
    std::cout << " shft: " << p << endl << endl;
  }
#endif
}
