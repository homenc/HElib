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
#include "permutations.h"

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
  // the data vector consists of interleaved permutations over [0,n-1],
  // n=getDim(dim), and the "step" between elements of any single permutation
  // is s=getProd(dim+1). Hence we have m=getSize()/(n*s) "chunks" of
  // interleaved permutations.
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

  // Allocate the output permutations
  long m = sig.getNumDims();
  ColPerm dummy(sig);
  out.assign(2*m -1, dummy); // allocate space and initialize

  // Allocate two temporary SlicePerm's 
  SlicePerm tmp1(sig,pi); // initialize to pi
  SlicePerm tmp2(sig);    // empty permutation
  SlicePerm* tp1 = &tmp1;
  SlicePerm* tp2 = &tmp2;

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

#if 0
int main()
{
  // A 2x3x2 cube
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
    std::cout << " expl: " << p << endl << endl;
  }
}
#endif
