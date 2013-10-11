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
/**
 * @file permutations.h
 * @brief Permutations over Hypercubes and their slices.
 **/
#ifndef _PERMUTATIONS_H_
#define _PERMUTATIONS_H_

#include <iostream>
#include "matching.h"
#include "hypercube.h"
#include "PAlgebra.h"

using namespace std;
using namespace NTL;

//! A simple permutation is just a vector with p[i]=\pi_i
typedef Vec<long> Permut;

//! Apply a permutation to a function, out[i]=in[p1[i]] (NOT in-place)
void ApplyPermToFunc(Vec<long>& out, const Vec<long>& in, const Permut& p1);

//! Apply two permutations to a function out[i]=in[p2[p1[i]]] (NOT in-place)
void ApplyPermsToFunc(Vec<long>& out, const Vec<long>& in,
		      const Permut& p2, const Permut& p1);

//! @brief A random size-n permutation
void RandomPerm(Permut& perm, long n);

class SlicePerm; // a "row permutation"
class ColPerm;   // a "column permutation"

/**
 * @class SlicePerm
 * @brief Permuting subcubes (rows/slices) of a hypercube.
 *
 * SlicePerm is derived from a HyperCube<long>, and it uses the cube object
 * to store the actual permutation data. The interpretation of this data,
 * however, depends on the data member dim.
 * 
 * A SlicePerm is breaking the cube into consecutive subcubes and permuting
 * each subcube separately. Specifically, the cube is partitioned into 
 * n=getSize()/getProd(dim) subcubes, so we have n consequtive permutations,
 * all over the domain [0,getProd(dim)-1]. A SlicePerm with dim=0 is just a
 * regular permutation over the entire cube.
 *
 * For example, for a 2x3x2 cube and dim=1 (i.e., the 2nd dimension), the data
 * vector could be something like [1 0 2 5 4 3 3 0 5 4 1 2]. This means the
 * permutations [1 0 2 5 4 3] over the 1st 2x3 subcube and [3 0 5 4 1 2] over
 * the 2nd subcube. The same permutation written explicitly would be
 * [1 0 2 5 4 3 9 6 11 10 7 8].
 **/
class SlicePerm : public HyperCube<long> {
private:
  long dim;
  SlicePerm(); // disabled

public:
  friend class ColPerm;

  explicit SlicePerm(const CubeSignature& _sig): HyperCube<long>(_sig)
  { dim = -1; }
  SlicePerm(const CubeSignature& _sig, const Permut& p): HyperCube<long>(_sig)
  { *this = p; } // use assignment operator to copy data

  long getPermDim() const {return dim;}

  void makeExplicit(Permut& out) const;   //! Write the permutation explicitly
  void extractSlice(Permut& out, long i) const;//! The perm over the ith subcube

  //! When dim=0 (the first dimension), SlicePerm and Permut are equivalent
  SlicePerm& operator=(const Permut &other) {
    dim = 0;           // Set the dimension to 0
    getData() = other; // Copy the raw data
    return *this;
  }

  //! When dim=m (the last dimension), ColPerm and SlicePerm are equivalent
  SlicePerm& operator=(const ColPerm &other);

  //! @brief Break a row permutation into column-row-column format.
  //! The input permutes each dimension-i subcube, and in the output rho1,rho3
  //! permute only along the i'th dimension and rho2 permutes each dimension-i+1
  //! subcube. This method cannot permute in-place, it is assumed that rho2 and
  //! *this point to disjoint objects.
  void breakPermTo3(ColPerm& rho1, SlicePerm& rho2, ColPerm& rho3) const;

 //! A test/debugging method
  void printout(ostream& s) { // a test/debugging method
    s << "Cube signature: " << getSig() << endl;
    s << "  dim="<<dim<<endl;
    s << "  data="<<getData()<<endl;
  }
};

/**
 * @class ColPerm
 * @brief Permuting a single dimension (column) of a hypercube
 * 
 * ColPerm is derived from a HyperCube<long>, and it uses the cube object
 * to store the actual permutation data. The interpretation of this data,
 * however, depends on the data member dim.
 * 
 * The cube is partitioned into n=getSize()/getDim(dim) subcubes, and each of
 * them is permuted separatly. Hence pi consists of n interleaved permutations,
 * each over the domain [0,getDim(dim)-1], and the "step" between elements of
 * any single permutation inside the data vector is getProd(dim+1). Permuting
 * along the last dimension is the same as SlicePerm for that dimension.
 * 
 * For example, permuting a 2x3x2 cube along dim=1 (the 2nd dimention), we
 * could have the data vector as  [ 1  1  0  2  2  0  2  0  1  1  0  2 ]. (In
 * this example we have n=4 and step=2.) This means the four subcubes are
 * permuted by the permutations     [ 1     0     2                      ]
 *                                  [    1     2     0                   ]
 *                                  [                   2     1     0    ]
 *                                  [                      0     1     2 ].
 * Written explicitly, we get:      [ 2  3  0  5  4  1 10  7  8  9  6 11 ].
 * 
 * Another representation that we provide is by "shift amount": how many
 * slots each element needs to move inside its small permutation. For the
 * example above, this will be:     [ 1     -1    0                      ]
 *                                  [    1     1    -2                   ]
 *                                  [                   2     0    -2    ]
 *                                  [                      0     0     0 ]
 * so we write the permutatation as [ 1  1 -1  1  0 -2  2  0  0  0 -2  0 ].
 **/
class ColPerm : public HyperCube<long> {
private:
  long dim;
  ColPerm(); // disabled

public:
  friend class SlicePerm;

  explicit ColPerm(const CubeSignature& _sig): HyperCube<long>(_sig)
  { dim = -1; }

  long getPermDim() const {return dim;}

  void makeExplicit(Permut& out) const;    // Write the permutation explicitly
  void extractSlice(Permut& out, long i) const; // The perm over the ith subcube

  //! For each position in the data vector, compute how many slots it should be
  //! shifted inside its small permutation. Returns zero if all the shift amount
  //! are zero, nonzero values otherwise.
  long getShiftAmounts(Permut& out) const;

  //! When dim=m (the last dimension), ColPerm and SlicePerm are equivalent
  ColPerm& operator=(const SlicePerm &other);

 //! A test/debugging method
  void printout(ostream& s) { // a test/debugging method
    s << "Cube signature: " << getSig() << endl;
    s << "  dim="<<dim<<endl;
    s << "  data="<<getData()<<endl;
  }
};

/**
 * @brief Takes a permutation pi over m-dimensional cube C=Z_{n1} x...x Z_{nm}
 * and expresses pi as a product pi = rho_{2m-1} o ... o rho_2 o rho_1 where
 * each rho_i is a column permutation along one dimension. Specifically for
 * i<m, the permutations rho_i and rho_{2(m-1)-i} permute the i'th dimension
 **/
void breakPermByDim(vector<ColPerm>& out, 
		    const Permut &pi, const CubeSignature& sig);

/****
 * A simple implementation of full binary trees (each non-leaf has 2 children)
 ****/
template<class T> class FullBinaryTree;

// A node in a tree, these nodes will be kept in a vector, so we use
// indexes rather than pointers
template<class T> class TreeNode {
  T data;
  long parent;
  long leftChild, rightChild;
  long prev, next; // useful, e.g., to connect all leaves in a list

  void makeNullIndexes() {parent = leftChild = rightChild = prev = next = -1;}

public:
  TreeNode() { makeNullIndexes(); }
  explicit TreeNode(const T& d): data(d) { makeNullIndexes(); }

  T& getData() { return data; }
  const T& getData() const { return data; }

  long getParent() const { return parent; }
  long getLeftChild() const { return leftChild; }
  long getRightChild() const { return rightChild; }
  long getPrev() const { return prev; }
  long getNext() const { return next; }

  friend class FullBinaryTree<T>;
};


// A binary tree, the root is always the node at index 0
template<class T> class FullBinaryTree {
  vector< TreeNode<T> > nodes;
  long nLeaves;             // how many leaves in this tree
  long firstLeaf, lastLeaf; // index of the first/last leaves

public:
  FullBinaryTree() { nLeaves=0; firstLeaf = lastLeaf = -1; } // empty tree

  explicit FullBinaryTree(const T& d)  // tree with only a root
  {
    nLeaves = 1;
    TreeNode<T> n(d);
    nodes.push_back(n);
    firstLeaf = lastLeaf = 0;
  }

  void PutDataInRoot(const T& d)
  {
    if (nodes.size()==0) { // make new root
      TreeNode<T> n(d);
      nodes.push_back(n);
      firstLeaf = lastLeaf = 0;
      nLeaves = 1;
    }
    else nodes[0].data = d; // Root exists, just update data
  }

  // Provide some of the interfaces of the underlying vector
  long size() { return (long) nodes.size(); }

  TreeNode<T>& operator[](long i) { return nodes[i]; }
  const TreeNode<T>& operator[](long i) const { return nodes[i]; }

  TreeNode<T>& at(long i) { return nodes.at(i); }
  const TreeNode<T>& at(long i) const { return nodes.at(i); }

  T& DataOfNode(long i) const { return nodes.at(i).data; }

  long getNleaves() const { return nLeaves; }
  long FirstLeaf() const { return firstLeaf; }
  long NextLeaf(long i) const { return nodes.at(i).next; }
  long PrevLeaf(long i) const { return nodes.at(i).prev; }
  long LastLeaf() const { return lastLeaf; }

  long RootIdx() const { return 0; }
  long ParentIdx(long i) const { return nodes.at(i).parent; }
  long LeftChildIdx(long i) const { return nodes.at(i).leftChild; }
  long RightChildIdx(long i) const { return nodes.at(i).rightChild; }

  //! If the parent is a leaf, add to it tho children with the given data,
  //! else just update the data of the two children of this parent.
  //! Returns the index of the left child, the right-child index is one
  //! more than the left-child index.
  long AddChildren(long parentIdx, const T& leftData, const T& rightData);

  //! Remove all nodes in the tree except for the root
  void CollapseToRoot()
  {
    if (nodes.size() > 1) {
      nodes.resize(1);
      firstLeaf = lastLeaf = 0;
      nLeaves = 1;
    }
  }
};
template <class T>
long FullBinaryTree<T>::AddChildren(long parentIdx, 
				    const T& leftData, const T& rightData)
{
  assert (parentIdx >= 0 && parentIdx < (long)(nodes.size()));

  // If parent is a leaf, add to it two children
  if (nodes[parentIdx].leftChild==-1 && nodes[parentIdx].rightChild==-1) {
    long childIdx = nodes.size();
    TreeNode<T> n1(leftData);
    nodes.push_back(n1); // add left child to vector
    TreeNode<T> n2(rightData);
    nodes.push_back(n2);// add right child to vector

    TreeNode<T>& parent = nodes[parentIdx];
    TreeNode<T>& left = nodes[childIdx];
    TreeNode<T>& right = nodes[childIdx+1];

    parent.leftChild = childIdx;            // point to children from parent
    parent.rightChild= childIdx+1;
    left.parent = right.parent = parentIdx; // point to parent from children

    // remove parent and insert children to the linked list of leaves
    left.prev = parent.prev;
    left.next = childIdx+1;
    right.prev = childIdx;
    right.next = parent.next;
    if (parent.prev>=0) { // parent was not the 1st leaf
      nodes[parent.prev].next = childIdx;
      parent.prev = -1;
    }
    else // parent was the first leaf, now its left child is 1st
      firstLeaf = childIdx;

    if (parent.next>=0) { // parent was not the last leaf
      nodes[parent.next].prev = childIdx+1;
      parent.next = -1;
    }
    else // parent was the last leaf, now its left child is last
      lastLeaf = childIdx+1;

    nLeaves++; // we replaced a leaf by a parent w/ two leaves
  }
  else { // parent is not a leaf, update the two children
    TreeNode<T>& parent = nodes[parentIdx];
    assert(parent.leftChild>=0 && parent.rightChild>=0);

    TreeNode<T>& left = nodes[parent.leftChild];
    TreeNode<T>& right = nodes[parent.rightChild];
    left.data = leftData;
    right.data = rightData;
  }
  return nodes[parentIdx].leftChild;
}


class BenesData; // information on a generalized Benes network

// A node in a tree relative to some generator
class SubDimension {
 public:
  long genIdx; // sub-dimension of what generator
  long size;   // Size of cube slice
  long e;      // shift-by-1 in this sub-dim is done via X -> X^{g^e}
  bool good;   // good or bad

  // If this is a Benes leaf, a description of its network (else NULL)
  BenesData* benes;

  explicit SubDimension(long idx=0, long sz=0, 
			long ee=0, bool gd=false, BenesData* bns=NULL)
  { genIdx=idx; size=sz; e=ee; good=gd; benes=bns; }

  /*  SubDimension& operator=(const SubDimension& other)
    { genIdx=other.genIdx; size=other.size; 
      e=other.e; good=other.good; benes=other.benes;
      return *this;
      } 
  */
};
typedef FullBinaryTree<SubDimension> OneGeneratorTree;// tree for one generator

//! A recursive procedure for computing the e exponent values
void computeEvalues(const OneGeneratorTree &T, long idx, long genOrd);

//! A vector of generator trees, one per generator in Zm*/(p)
class GeneratorTrees  {
  Vec<OneGeneratorTree> trees;
  Permut map2cube, map2array;

 public:
  GeneratorTrees() {} // default constructor

  GeneratorTrees(const Vec<OneGeneratorTree>& _trees): trees(_trees) {}

  // Initialze trees with only the roots.
  GeneratorTrees(const Vec<SubDimension>& dims);

  long length() const { return trees.length(); }
  OneGeneratorTree& operator[](long i) { return trees[i]; }
  const OneGeneratorTree& operator[](long i) const { return trees[i]; }

  OneGeneratorTree& at(long i) { return trees.at(i); }
  const OneGeneratorTree& at(long i) const { return trees.at(i); }

  long getSize() const { return map2cube.length(); }
  OneGeneratorTree& getGenTree(long i) { return trees.at(i); }
  const OneGeneratorTree& getGenTree(long i) const { return trees.at(i); }

  const Vec<long>& Map2Cube() const { return map2cube; }
  const Vec<long>& Map2Array() const { return map2array; }
  Vec<long>& Map2Cube() { return map2cube; }
  Vec<long>& Map2Array() { return map2array; }

  long Map2Cube(long i) const { return map2cube[i]; }
  long Map2Array(long i) const { return map2array[i]; }

  //! Get the cube dimensions corresponding to the vector of trees,
  //! the ordered vector with one dimension per leaf in all the trees.
  void getCubeDims(Vec<long>& dims) const;

  // Returns coordinates of i relative to leaves of the tree
  //  void getCoordinates(Vec<long>&, long i) const;

  //! Compute the trees corresponding to the "optimal" way of breaking
  //! a permutation into dimensions, subject to some constraints
  void BuildOptimalTrees(long widthBound);

  /**
   * @brief Computes permutations mapping between linear array and the cube.
   *
   * If the cube dimensions (i.e., leaves of tree) are n1,n2,...,nt and
   * N=\prod_j n_j is the size of the cube, then an integer i can be
   * represented in either the mixed base of the n_j's or in "CRT basis"
   * relative to the leaves: Namely either
   *                            i = \sum_{j<=t}  i_j  * \prod_{k>j} n_k,
   *                         or i = \sum_leaf i'_leaf * leaf.e mod N.
   *
   * The breakPermByDim procedure expects its input in the mixed-base
   * representation, and the maps are used to convert back and forth.
   * Specifically, let (i'_1,...,i'_t) be the CRT representation of i in
   * this cube, and j = \sum_{j=1}^t i'_j * \prod_{k>j} n_k, then we have
   * map2cube[i]=j and map2array[j]=i.
   **/
  void ComputeCubeMapping();
};


// Permutation networks

class Ctxt;
class EncryptedArray;
class PermNetwork;

// The information needed to apply one layer of a permutation network
class PermNetLayer {
  long genIdx; // shift-by-1 in this layer is done via X -> X^{g^e}
  long e;
  Vec<long> shifts; // shifts[i] is how much to shift slot i
  bool isID; // a silly optimization, does this layer copmute the identity?

  friend class PermNetwork;
};

class PermNetwork {
  Vec<PermNetLayer> layers;

public:
  PermNetwork() {}; // empty network
  PermNetwork(const Permut& pi,const GeneratorTrees& trees)
    { BuildNetwork(pi, trees); }

  long Width() const { return layers.length(); }

  // Take as input a permutation pi and the trees of all the generators,
  // and prepares the permutation network for this pi
  void BuildNetwork(const Permut& pi, const GeneratorTrees& trees);

  //! Apply network to permute a ciphertext
  void ApplyToCtxt(Ctxt& c);

  //! Apply network to array, used mostly for debugging
  void ApplyToArray(HyperCube<long>& v);

  //! Apply network to plaintext polynomial, used mostly for debugging
  void ApplyToPtxt(ZZX& p, const EncryptedArray& ea);
};

#endif /* ifndef _PERMUTATIONS_H_ */
