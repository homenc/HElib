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
#include <NTL/ZZ.h>
#include "NumbTh.h"
#include "permutations.h"

// Apply a permutation to a function, out[i]=in[p1[i]]
void ApplyPermutation(Vec<long>& out, const Vec<long>& in, const Permut& p1)
{
  assert(&out != &in); // NOT an in-place procedure
  out.SetLength(p1.length());
  for (long i=0; i<p1.length(); i++)
    out[i] = in.at(p1[i]);
}

// Apply two permutations to a function out[i]=in[p2[p1[i]]]
void ApplyPermutations(Vec<long>& out, const Vec<long>& in,
		       const Permut& p2, const Permut& p1)
{
  assert(&out != &in); // NOT an in-place procedure
  assert(p1.length() == p2.length());
  out.SetLength(p1.length());
  for (long i=0; i<p1.length(); i++)
    out[i] = in.at(p2.at(p1[i]));
}


// Write a slice (row) permutation explicitly
void SlicePerm::makeExplicit(Permut& out) const
{
  // pi consists of separate contiguous permutations, each over [0,n-1]
  long n = getProd(dim);
  long nPerms = getSize()/n; // how many separate permutations

  out.SetLength(getSize());  // allocate space and initialize
  long offset = 0;
  while (nPerms--) {
    for (long i=0; i<n; i++)
      out[offset+i] = (*this)[offset+i] + offset;
    offset += n;
  }
}

// Extrtact the permutation over the i'th slice
void SlicePerm::extractSlice(Permut& out, long i) const
{
  // *this consists of separate contiguous permutations, each over [0,n-1]
  long n = getProd(dim);
  long offset = i*n; // The beginning of the i'th slice
  assert (i>=0 && offset < getSize());

  out.SetLength(n); // allocate space and initialize
  for (long j=0; j<n; j++) out[j] = (*this)[offset+j];
}

// When dim=m (the last dimension), ColPerm and SlicePerm are equivalent
SlicePerm& SlicePerm::operator=(const ColPerm &other)
{
  assert(&getSig()==&other.getSig());// Ensure the same cube dimensions
  assert(other.dim==getNumDims()-1); //   and dim = the last dimension,
  getData() = ((HyperCube<long>&) other).getData(); // then copy the raw data.
  dim = other.dim;
  return *this;
}

// Write a column permutation explicitly
void ColPerm::makeExplicit(Permut& out) const
{
  /* the data vector consists of interleaved permutations over [0,n-1],
   * n=getDim(dim), and the "step" between elements of any single permutation
   * is s=getProd(dim+1). Hence we have m=getSize()/(n*s) "chunks" of
   * interleaved permutations. 
   *
   * For example, with a 2x3x2 cube and dim=1, we have size=12, n=4 and s=2,
   * and hence two "chunks". Data vector [1  1  0  2  2  0  2  0  1  1  0  2]
   * in this example is interpreted as four permutations:
   *                                     [1     0     2                     ] 
   *                                     [   1     2     0                  ] 
   *                                     [                  2     1     0   ] 
   *                                     [                     0     1     2] 
   * Written explicitly, we then get:    [2  3  0  5  4  1 10  7  8  9  6 11].
   */
  long n = getDim(dim);     // the permutations are over [0,n-1]
  long s = getProd(dim+1);  // step between elements of a permutation
  long m = getSize()/(n*s); // how many chunks of interleaved permutations

  out.SetLength(getSize()); // allocate space and initialize
  long offset = 0;
  while (m--) { // go over the chunks one at a time
    for (long i=0; i<n; i++) for (long j=0; j<s; j++) {
	long ind = i*s +j;
	out[ind+offset] = (*this)[ind+offset]*s +j +offset;
      }
    offset += n*s;
  }
}

// For each position in the data vector, compute how many slots it should be
// shifted inside its small permutation. Returns zero if all the shift amount
// are zero, nonzero otherwise.
long ColPerm::getShiftAmounts(Vec<long>& out) const
{
  /* the data vector consists of interleaved permutations over [0,n-1],
   * n=getDim(dim), and the "step" between elements of any single permutation
   * is s=getProd(dim+1). Hence we have m=getSize()/(n*s) "chunks" of
   * interleaved permutations. 
   *
   * For example, with a 2x3x2 cube and dim=1, we have size=12, n=4 and s=2,
   * and hence two "chunks". Data vector [1  1  0  2  2  0  2  0  1  1  0  2]
   * in this example is interpreted as four permutations:
   *                                     [1     0     2                     ]
   *                                     [   1     2     0                  ]
   *                                     [                  2     1     0   ]
   *                                     [                     0     1     2]
   * Hence the shift amount are:         [1     -1    0                     ]
   *                                     [   1     1    -2                  ]
   *                                     [                  2     0    -2   ]
   *                                     [                     0     0     0]
   * so the output vector is [1  1 -1  1  0 -2  2  0  0  0 -2  0].
   */
  long n = getDim(dim);     // the permutations are over [0,n-1]
  long s = getProd(dim+1);  // step between elements of a permutation
  long m = getSize()/(n*s); // how many chunks of interleaved permutations
  long nonZero = 0;

  out.SetLength(getSize()); // allocate space and initialize
  long offset = 0;
  while (m--) { // go over the chunks one at a time
    for (long i=0; i<n; i++) for (long j=0; j<s; j++) {
	long ind = i*s +j +offset;
	out[ind] = (*this)[ind] -i;
	nonZero |= out[ind];
      }
    offset += n*s;
  }
  return nonZero;
}

// Extrtact the permutation over the i'th column
void ColPerm::extractSlice(Permut& out, long i) const
{
  long n = getDim(dim);    // the permutations are over [0,n-1]
  assert (i>=0 && i*n<getSize());

  long s = getProd(dim+1); // step between elements of a permutation

  // every chunk consists of s permutations, which chunk contains the i'th
  long chunk = i / s;
  long offset = (i % s) + (chunk*n*s); // offset from beginning of pi

  out.SetLength(n); // allocate space and initialize
  for (long j=0; j<n; j++)
    out[j] = (*this)[offset + j*s];
}

// When dim=m (the last dimension), ColPerm and SlicePerm are equivalent
ColPerm& ColPerm::operator=(const SlicePerm &other)
{
  assert(&getSig()==&other.getSig());// Ensure the same cube dimensions
  assert(other.dim==getNumDims()-1); //   and dim = the last dimension,
  getData() = ((HyperCube<long>&) other).getData(); // then copy the raw data.
  dim = other.dim;
  return *this;
}


/* Takes a permutation pi over an m-dimensional cube C=Z_{n1} x ... x Z_{nm}
 * and expresses pi as a product pi = rho_{2m-1} o ... o rho_2 o rho_1 where
 * each rho_i is a column permutation along one dimension. Specifically for
 * i<m, the permutations rho_i and rho_{2(m-1)-i} permute the i'th dimension
 ************************************************************************/
void breakPermByDim(vector<ColPerm>& out, 
		    const Permut &pi, const CubeSignature& sig)
{
  assert(sig.getSize()==pi.length());

  // Allocate two temporary SlicePerm's 
  SlicePerm tmp1(sig,pi); // initialize to pi
  SlicePerm tmp2(sig);    // empty permutation
  SlicePerm* tp1 = &tmp1;
  SlicePerm* tp2 = &tmp2;

  // Allocate the output permutations
  long m = sig.getNumDims();
  ColPerm dummy(sig);
  out.assign(2*m -1, dummy); // allocate space and initialize
  if (m == 1) { // special case, no need to break
    out[0] = tmp1;
    return;
  }

  for (long i=0; i<m-2; i++) {
    tp1->breakPermTo3(out[i], *tp2, out[2*m-i-2]);
    std::swap(tp1,tp2);
  }
  // use pointer hack to cast between SlicePerm and ColPerm
  tp2 = (SlicePerm*) &out[m-1];
  tp1->breakPermTo3(out[m-2], *tp2, out[m]);
}

// Break a permutation into column-row-column format. The input pi permutes
// each dimension-i subcube, and in the output rho1,rho3 permute only along
// the i'th dimension and rho2 permutes each dimension-i+1 subcube.
// This routine cannot permute in-place, it is assumed that pi and rho2 point
// to disjoint vectors.
void SlicePerm::breakPermTo3(ColPerm& rho1, 
			     SlicePerm& rho2, ColPerm& rho3) const
{
  assert(&rho1.getSig()==&getSig());
  assert(&rho2.getSig()==&getSig());
  assert(&rho3.getSig()==&getSig());

  // *this consists of separate permutations over [0,n-1], and each
  // of these is viewed as a permutation over an n1 x n2 cebe

  long n1 = getDim(dim); // Size of this dimension, inherited from HyperCube
  long n2 = getProd(dim+1); 
  long n = getProd(dim); // = n1*n2;

  // representing I_n as I_n1 x I_n2: i == n2*rep[i].first + rep[i].second
  pair<long,long> rep[n];
  for (long ind=0,i=0; i<n1; i++) for (long j=0; j<n2; j++,ind++) {
      rep[ind].first = i;
      rep[ind].second = j;
    }

  long offset = 0;
  long nPerms = getSize()/n;  // how many separate permutations
  while (nPerms > 0) { // break the permutations one at a time

    // Construct a bipartite n2-by-n2 graph for pi (cf. Lemma 1 in [GHS12a]).
    // For each j=pi(i) with representations i=(i1,i2) and j=(j1,j2), we put
    // in the bipartite graph an edge i2->j2 and label it by i.
    BipartitleGraph bg;
    for (long i=0; i<n; i++) {
      long j = (*this)[offset+i]; // the image of i under the permutation
      // when i = (i1,i2) and j=(j1,j2), add an edge from i2 to j2 labeled i
      bg.addEdge(rep[i].second, rep[j].second, i);
      // cerr <<"  "<<rep[i].second<<"->"<<rep[j].second<<" label "<<i<<endl;
    }
    //  cerr << endl;

    // The bipartite graph is n1-regular, so we can break its edges into
    // n1 perfect matchings, which are numbered 1,2,...,n1.
    bg.partitionToMatchings();

    // The output permutations are defined by the representation i<->(i1,i2),
    // the target permutation pi, and the coloring of the bipartite graph.
    // Denote by sigma(i) the color of the edge labeled by i, sigma(i) \in
    // {1,...,n1}. Also denote by (pi^1(i),pi^2(i)) the representation of
    // pi(i). Then:
    //
    // + rho_1 is defined by           (i1, i2) -> (sigma(i), i2)
    // + rho_2 is defined by     (sigma(i), i2) -> (sigma(i),pi^2(i))
    // + rho_3 is defined by (sigma(i),pi^2(i)) -> (pi^1(i), pi^2(i))
    //
    // rho_1 is a permutation because for every value of i2 (corresponding to
    // a left node in the graph), all the edges leaving that node have different
    // colors.
    //
    // rho_2 is a permutation since the edges of each color form a perfect,
    // so for every left-node j2 and color c there is a single c-colored edge
    // going into j2. The label of that edge determines a unique origin index
    // i=(i1,i2), and therefore also the pre-image of (c,j2) under rho_2, which
    // is (sigma(i),i2)=(c,i2).
    //
    // rho_3 is a permutation because rho_1,rho_2,pi are permutations, and
    // pi = rho_3 o rho_2 o rho_1.
    //
    // Note that the edges are colored 1..n2 while our represenation above
    // has the second digits in the range 0..n2-1, so below we use sigma(i)-1
    // rather than sigma(i).

    for (long i2=0; i2<n2; i2++) // go over all edges in the bipartite graph
      for (LNeighborList::iterator it=bg.left[i2].neighbors.begin(); 
	   it!=bg.left[i2].neighbors.end(); ++it) {
	LabeledEdge& e = it->second; // An edge in the bipartite graph;
	long i = e.label; // labeled by i
	long c = e.color -1; // colored by c (after the -1 adjustment)
	// i2 = e.from = rep[i].second;
	long j = (*this)[offset+i];
	long j1 = rep[j].first;
	long j2 = e.to;   // = it->first = rep[j].second;

	long tmp1 = c*n2 + i2;  // the image of i under rho1 = (c,i2)
	rho1[offset+i] = c;
	long tmp2 = c*n2 + j2;  // the image of tmp1 under rho2 =(c,j2)
	rho2[offset+tmp1] = j2;
	rho3[offset+tmp2] = j1; // image of tmp2 under rho3 =(j1,j2)=pi(i)
	// cerr <<"  "<< i <<" -> "<< tmp1 <<" -> "<< tmp2 <<" ->  "
	//	<< (*this)[offset+i] << endl;
      }
    --nPerms;
    offset += n;
  }
  rho1.dim = rho3.dim = dim;
  rho2.dim = dim+1;
}

void randomPerm(Permut& perm, long n)
{
  perm.SetLength(n);
  for (long j = 0; j < n; j++)
     perm[j] = j;
   
  // random shuffle
  for (long m = n; m > 0; m--) {
     long p = RandomBnd(m);
     // swap positions p and m-1 of perm
     long tmp = perm[p];
     perm[p] = perm[m-1];
     perm[m-1] = tmp;
  }
}

// A recursive procedure for computing the e exponent values
// FIXME: This code should be part of the logic that builds the tree
void computeEvalues(const GeneratorTree &T, long idx, long genOrd)
{
  // if either child is missing then we are at a leaf (this is a full tree)
  long left = T[idx].getLeftChild();
  long right = T[idx].getRightChild();
  if (left < 0 || right < 0) return;

  SubDimension& lData = (SubDimension&) T[left].getData();
  SubDimension& rData = (SubDimension&) T[right].getData();

  // The entire tree size is sz1 * sz2;
  long sz1 = lData.size;  // size of left subtree
  long sz2 = rData.size; // size of right subtree

  long ee = T[idx].getData().e;
  // If right tree is bad, copy e from parent to right and scale left
  if (!rData.good) {
    rData.e = ee;
    lData.e = ee * sz2;
  }
  // If right tree is good and left is bad, copy parent to left and scale right
  else if (!lData.good) {
    lData.e = ee;
    rData.e = ee * sz1;
  }
  // If both are good, scale both using CRT coefficients (assuming GCD==1)
  else {
    long f1 = CRTcoeff(sz1, sz2); // f1 = 0 mod sz1, 1 mod sz2
    long f2 = sz1*sz2 +1 - f1;    // f2 = 1 mod sz1, 0 mod sz2
    lData.e = MulMod(ee, f2, genOrd);
    rData.e = MulMod(ee, f1, genOrd);
  }
  // Recurse on the two subtrees
  computeEvalues(T, left, genOrd);
  computeEvalues(T, right, genOrd);
}

// Get the cube dimensions corresponding to a vector of trees,
// the ordered vector with one dimension per leaf in any of the trees.
void GeneratorTree::GetCubeDims(Vec<long>& dims,
				const vector<GeneratorTree>& trees)
{
  // how many dimensions do we need
  long nDims = 0;
  for (long i=0; i<(long)trees.size(); i++)
    nDims += trees[i].getNleaves();
  dims.SetLength(nDims); // set the size

  // copy dims from the leaves in all the trees
  long idx = 0;
  for (long i=0; i<(long)trees.size(); i++) {
    const GeneratorTree& T = trees[i];
    for (long leaf=T.FirstLeaf(); leaf>=0; leaf=T.NextLeaf(leaf))
      dims[idx++] = T[leaf].getData().size;
  }
}

void PermNetwork::BuildNetwork(const Permut& pi, 
			       const vector<GeneratorTree>& trees,
			       const PAlgebra& zmStar)
{
  // FIXME: For now we only implement the single-generator case
  assert(trees.size()==1);

  Vec<long> dims;
  GeneratorTree::GetCubeDims(dims, trees);
  CubeSignature sig(dims); // make a cube-signature object
  std::cerr << "\ndims="<<dims<<endl;

  const GeneratorTree &T = trees[0];
  long genIdx = T[0].getData().genIdx;
  long genOrd = zmStar.OrderOf(genIdx);

  // Traverse the tree and compute the e exponent in each node
  // FIXME: This should be part of the logic that builds the tree
  ((TreeNode<SubDimension>&)T[0]).getData().e = 1; // The root e value is 1
  computeEvalues(T,0,genOrd); // compute the exponents in the rest of the tree

  // Compute the total number of levels in the target permutation network
  long nLvls=0;
  for (long i=0, leaf=T.FirstLeaf(); leaf>=0; i++, leaf=T.NextLeaf(leaf)) {
    long lvls=1; // how many levels to add for this leaf

    const SubDimension& leafData = T[leaf].getData();
    if (leafData.benes != NULL) // a benes leaf
      lvls = 3; /* FIXME: leafData.benes->getNLvls(); */

    if (T.NextLeaf(leaf)>=0) lvls *= 2; // not last leaf, so not middle level
    nLvls += lvls;
  }

  // Compute the mapping from linear array to cube and back
  Permut map2cube, map2array;
  T.CubeMapping(map2cube, map2array);

  std::cerr << "pi =      "<<pi<<endl;
  std::cerr << "map2cube ="<<map2cube<<endl;
  std::cerr << "map2array="<<map2array<<endl;

  // Compute the permutation on the cube, rho = map2cube o pi o map2array
  Permut rho;
  ApplyPermutations(rho, map2cube, pi, map2array);
  std::cerr << "rho =     "<<rho<<endl;

  // Break rho along the different dimensions
  vector<ColPerm> perms;
  breakPermByDim(perms, rho, sig);

  for (long i=0; i<(long)perms.size(); i++) {
    Permut tmp;
    perms[i].makeExplicit(tmp);
    std::cerr << " prems["<<i<<"]="<<tmp<<endl;
  }

  // Go over the different permutations and build the corresponding levels
  this->SetLength(nLvls); // allocate space
  Vec<long> shifts;
  for (long i=0,j=0,leaf=T.FirstLeaf(); leaf>=0; leaf=T.NextLeaf(leaf),i++) {
    const SubDimension& leafData = T[leaf].getData();
    const ColPerm& p1 = perms[i]; // leaf corresponding to permutations i,n-i
    const ColPerm& p2 = perms[perms.size()-i-1];

    if (leafData.benes != NULL) { // A Benes leaf
      // FIXME: do something here..
      j += 3; /* FIXME: leafData.benes->getNLvls(); */
      continue;
    }

    // A naive permutation leaf, corresponding to 1 or 2 levels
    PermNetLayer& l1 = (*this)[j];
    l1.isID = !p1.getShiftAmounts(shifts);
    if (!l1.isID) {
      std::cerr << "layer "<<i<<": "<<shifts<<endl;
      if (leafData.good) // For good leaves, shift by -x is the same as size-x
	for (long k=0; k<shifts.length(); k++)
	  if (shifts[k]<0) shifts[k] += leafData.size;
      // li.shifts[i] := shifts[map2cube[i]]
      ApplyPermutation(l1.shifts, shifts, map2cube);
      std::cerr << "       : "<<l1.shifts<<endl;
    }
    else std::cerr << "layer "<<i<<"= identity\n";

    if (j == this->length()-j-1) break; // the middle layer (last leaf)

    PermNetLayer& l2 = (*this)[this->length()-j-1];
    l2.isID = !p2.getShiftAmounts(shifts);
    if (!l2.isID) {
      std::cerr << "layer "<<(perms.size()-i-1)<<": "<<shifts<<endl;
      if (leafData.good) // For good leaves, shift by -x is the same as size-x
	for (long k=0; k<shifts.length(); k++)
	  if (shifts[k]<0) shifts[k] += leafData.size;
      // li.shifts[i] := shifts[map2cube[i]]
      ApplyPermutation(l2.shifts, shifts, map2cube);
      std::cerr << "       : "<<l2.shifts<<endl;
    }
    else std::cerr << "layer "<<(perms.size()-i-1)<<"= identity\n";

    j++;
  }
}

// Adds one to the little-endian representation of an integer in base digits,
// returns true if there was an overflow
static bool addOne(Vec<long>& rep, const Vec<long> digits)
{
  for (long i=rep.length()-1; i>=0; --i) {
    rep[i]++;
    if (rep[i] >= digits[i])
      rep[i] -= digits[i];
    else
      return false;
  }
  return true;
}

void GeneratorTree::CubeMapping(Permut& map2cube, Permut& mapBack) const
{
  Vec<long> dims(INIT_SIZE, getNleaves());
  Vec<long> coefs(INIT_SIZE,getNleaves());
  for (long i=getNleaves()-1, leaf=LastLeaf(); i>=0; i--,leaf=PrevLeaf(leaf)) {
    dims[i] = (*this)[leaf].getData().size;
    coefs[i] = (*this)[leaf].getData().e;
  }
  std::cerr << "coefs="<<coefs<<endl;

  // A representation of an integer with digits from dims
  Vec<long> rep(INIT_SIZE, getNleaves());
  for (long i=0; i<rep.length(); i++) rep[i]=0; // initialize to zero

  // initialize to all zero
  long sz = (*this)[0].getData().size;
  map2cube.SetLength(sz);
  mapBack.SetLength(sz);
  for (long i=0; i<map2cube.length(); i++) map2cube[i]=mapBack[i]=0;

  // compute the permutations
  for (long i=1; i<sz; i++) {
    addOne(rep, dims); // representation of i in base dims
    for (long j=0; j<coefs.length(); j++) {
      long tmp = MulMod(rep[j], coefs[j], sz);
      map2cube[i] = AddMod(map2cube[i], tmp, sz);
    }
    mapBack[map2cube[i]] = i;
  }
}


int main()
{
  // Build a "good" 2x3 tree for the order-6 generator g=5 in Z_31^*/(2)
  PAlgebra al(/*m=*/31, /*p=*/2);

  // Build a 2x3x3 tree for the order-18 generator g=24 in Z_127^*/(2)
  //PAlgebra al(/*m=*/127, /*p=*/2);
  //SubDimension rtDim(/*genIdx=*/0, /*size=*/18, /*e=*/1, /*good=*/true);

  for (long iii=0; iii<3; iii++) {
  // Initialize the tree with just the size-6 root
  SubDimension rtDim(/*genIdx=*/0, /*size=*/6, /*e=*/1, /*good=*/true);
  GeneratorTree T(rtDim);

  // A random size-6 permutation
  Permut p;
  randomPerm(p, T[0].getData().size);

  // Add the size-2 and size-3 children
  SubDimension lftDim(/*genIdx=*/0, /*size=*/2, /*e=*/3, /*good=*/false);
  SubDimension rgtDim(/*genIdx=*/0, /*size=*/3, /*e=*/4, /*good=*/true);
  T.AddChildren(0, rgtDim, lftDim);

  // A vector of trees, with only this one tree in it
  vector<GeneratorTree> trees;
  trees.push_back(T);

  PermNetwork ntwrk;
  ntwrk.BuildNetwork(p, trees, al);
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
