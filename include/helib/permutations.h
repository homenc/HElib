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
/**
 * @file permutations.h
 * @brief Permutations over Hypercubes and their slices.
 **/
#ifndef HELIB_PERMUTATIONS_H
#define HELIB_PERMUTATIONS_H

#include <helib/PAlgebra.h>
#include <helib/matching.h>
#include <helib/hypercube.h>
#include <helib/apiAttributes.h>

namespace helib {

//! A simple permutation is just a vector with p[i]=\pi_i
typedef NTL::Vec<long> Permut;

//! Apply a permutation to a std::vector, out[i]=in[p1[i]] (NOT in-place)
template <typename T>
void applyPermToVec(NTL::Vec<T>& out, const NTL::Vec<T>& in, const Permut& p1);
template <typename T>
void applyPermToVec(std::vector<T>& out,
                    const std::vector<T>& in,
                    const Permut& p1);

//! Apply two permutations to a std::vector out[i]=in[p2[p1[i]]] (NOT in-place)
template <typename T>
void applyPermsToVec(NTL::Vec<T>& out,
                     const NTL::Vec<T>& in,
                     const Permut& p2,
                     const Permut& p1);
template <typename T>
void applyPermsToVec(std::vector<T>& out,
                     const std::vector<T>& in,
                     const Permut& p2,
                     const Permut& p1);

//! @brief A random size-n permutation
void randomPerm(Permut& perm, long n);

/**
 * @class ColPerm
 * @brief Permuting a single dimension (column) of a hypercube
 *
 * ColPerm is derived from a HyperCube<long>, and it uses the cube object
 * to store the actual permutation data. The interpretation of this data,
 * however, depends on the data member dim.
 *
 * The cube is partitioned into columns of size n = getDim(dim):
 * a single column consists of the n entries whose indices i have the
 * same coordinates in all dimensions other than dim.  The entries
 * in any such column form a permutation on [0..n).
 *
 * For a given ColPerm perm, one way to access each column is as follows:
 *   for slice_index = [0..perm.getProd(0, dim))
 *     CubeSlice slice(perm, slice_index, dim)
 *     for col_index = [0..perm.getProd(dim+1))
 *        getHyperColumn(column, slice, col_index)
 *
 * Another way is to use the getCoord and addCoord methods.
 *
 * For example, permuting a 2x3x2 cube along dim=1 (the 2nd dimension), we
 * could have the data std::vector as  [ 1  1  0  2  2  0  2  0  1  1  0  2 ].
 * This means the four columns are
 * permuted by the permutations     [ 1     0     2                      ]
 *                                  [    1     2     0                   ]
 *                                  [                   2     1     0    ]
 *                                  [                      0     1     2 ].
 * Written explicitly, we get:      [ 2  3  0  5  4  1 10  7  8  9  6 11 ].
 *
 * Another representation that we provide is by "shift amount": how many
 * slots each element needs to move inside its small permutation. For the
 * example above, this will be:     [ 1    -1     0                      ]
 *                                  [    2    -1    -1                   ]
 *                                  [                   2     0    -2    ]
 *                                  [                      0     0     0 ]
 * so we write the permutation as [ 1  1 -1  1  0 -2  2  0  0  0 -2  0 ].
 **/
class ColPerm : public HyperCube<long>
{
private:
  long dim;
  ColPerm(); // disabled

public:
  explicit ColPerm(const CubeSignature& _sig) : HyperCube<long>(_sig)
  {
    dim = -1;
  }

  long getPermDim() const { return dim; }
  void setPermDim(long _dim)
  {
    assertInRange(_dim,
                  0l,
                  getNumDims(),
                  "Algebra does not have a dimension of index _dim: " +
                      std::to_string(_dim));
    dim = _dim;
  }

  void makeExplicit(Permut& out) const; // Write the permutation explicitly

  //! For each position in the data std::vector, compute how many slots it
  //! should be shifted inside its small permutation. Returns zero if all the
  //! shift amounts are zero, nonzero values otherwise.
  long getShiftAmounts(Permut& out) const;

  //! Get multiple layers of a Benes permutation network. Returns in out[i][j]
  //! the shift amount to move item j in the i'th layer. Also isID[i]=true if
  //! the i'th layer is the identity (i.e., contains only 0 shift amounts).
  void getBenesShiftAmounts(NTL::Vec<Permut>& out,
                            NTL::Vec<bool>& idID,
                            const NTL::Vec<long>& benesLvls) const;

  //! A test/debugging method
  void printout(std::ostream& s)
  { // a test/debugging method
    s << "Cube signature: " << getSig() << std::endl;
    s << "  dim=" << dim << std::endl;
    s << "  data=" << getData() << std::endl;
  }
};
std::ostream& operator<<(std::ostream& s, const ColPerm& p);

/**
 * @brief Takes a permutation pi over m-dimensional cube C=Z_{n1} x...x Z_{nm}
 * and expresses pi as a product pi = rho_{2m-1} o ... o rho_2 o rho_1 where
 * each rho_i is a column permutation along one dimension. Specifically for
 * i<m, the permutations rho_i and rho_{2(m-1)-i} permute the i'th dimension
 **/
void breakPermByDim(std::vector<ColPerm>& out,
                    const Permut& pi,
                    const CubeSignature& sig);

/**
 * @class GeneralBenesNetwork
 * @brief Implementation of generalized Benes Permutation Network
 **/
class GeneralBenesNetwork
{
private:
  long n; // size of perm, n > 1, not necessarily power of 2
  long k; // recursion depth, k = least integer k s/t 2^k >= n

  NTL::Vec<NTL::Vec<short>> level;
  // there are 2*k - 1 levels, each with n nodes.
  // level[i][j] is 0, 1, or -1,
  //   which designates an edge from node j at level i
  //   to node j + level[i][j]*shamt(i) at level i+1

  GeneralBenesNetwork(); // default constructor disabled

public:
  //! computes recursion depth k for generalized Benes network of size n.
  //! the actual number of levels in the network is 2*k-1
  static long depth(long n)
  {
    long k = 1;
    while ((1L << k) < n)
      k++;
    return k;
  }

  //! maps a level number i = 0..2*k-2 to a recursion depth d = 0..k-1
  //! using the formula d = (k-1)-|(k-1)-i|
  static long levelToDepthMap(UNUSED long n, long k, long i)
  {
    assertInRange<InvalidArgument>(i,
                                   0l,
                                   2 * k - 1,
                                   "Level number i not in [0, 2 * k - 1)");
    return (k - 1) - labs((k - 1) - i);
  }

  //! shift amount for level number i=0..2*k-2 using the formula
  //! ceil( floor(n/2^d) / 2), where d = levelToDepthMap(i)
  static long shamt(long n, long k, long i)
  {
    long d = levelToDepthMap(n, k, i);
    return ((n >> d) + 1) >> 1;
  }

  long getDepth() const { return k; }
  long getSize() const { return n; }
  long getNumLevels() const { return 2 * k - 1; }

  const NTL::Vec<short>& getLevel(long i) const
  {
    assertInRange<InvalidArgument>(i,
                                   0l,
                                   2 * k - 1,
                                   "Level number i not in [0, 2 * k - 1)");
    return level[i];
  }

  long levelToDepthMap(long i) const { return levelToDepthMap(n, k, i); }
  long shamt(long i) const { return shamt(n, k, i); }

  // constructor
  GeneralBenesNetwork(const Permut& perm);

  // test correctness

  bool testNetwork(const Permut& perm) const;
};

/**
 * @class FullBinaryTree
 * @brief A simple implementation of full binary trees (each non-leaf has 2
 * children)
 */
template <typename T>
class FullBinaryTree;

//! @class TreeNode
//! @brief A node in a full binary tree.
//!
//! These nodes are in a std::vector, so we use indexes rather than pointers
template <typename T>
class TreeNode
{
  T data;
  long parent;
  long leftChild, rightChild;
  long prev, next; // useful, e.g., to connect all leaves in a list

  void makeNullIndexes() { parent = leftChild = rightChild = prev = next = -1; }

public:
  TreeNode() { makeNullIndexes(); }
  explicit TreeNode(const T& d) : data(d) { makeNullIndexes(); }

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
template <typename T>
class FullBinaryTree
{
  long aux; // a "foreign key" for this tree (holds generator index)
  std::vector<TreeNode<T>> nodes;
  long nLeaves;           // how many leaves in this tree
  long frstLeaf, lstLeaf; // index of the first/last leaves

public:
  FullBinaryTree(long _aux = 0)
  {
    aux = _aux;
    nLeaves = 0;
    frstLeaf = lstLeaf = -1;
  } // empty tree

  explicit FullBinaryTree(const T& d, long _aux = 0) // tree with only a root
  {
    aux = _aux;
    nLeaves = 1;
    TreeNode<T> n(d);
    nodes.push_back(n);
    frstLeaf = lstLeaf = 0;
  }

  void putDataInRoot(const T& d)
  {
    if (nodes.size() == 0) { // make new root
      TreeNode<T> n(d);
      nodes.push_back(n);
      frstLeaf = lstLeaf = 0;
      nLeaves = 1;
    } else
      nodes[0].data = d; // Root exists, just update data
  }

  // Provide some of the interfaces of the underlying std::vector
  long size() { return (long)nodes.size(); }

  TreeNode<T>& operator[](long i) { return nodes[i]; }
  const TreeNode<T>& operator[](long i) const { return nodes[i]; }

  TreeNode<T>& at(long i) { return nodes.at(i); }
  const TreeNode<T>& at(long i) const { return nodes.at(i); }

  T& DataOfNode(long i) { return nodes.at(i).data; }
  const T& DataOfNode(long i) const { return nodes.at(i).data; }

  long getAuxKey() const { return aux; }
  void setAuxKey(long _aux) { aux = _aux; }

  long getNleaves() const { return nLeaves; }
  long firstLeaf() const { return frstLeaf; }
  long nextLeaf(long i) const { return nodes.at(i).next; }
  long prevLeaf(long i) const { return nodes.at(i).prev; }
  long lastLeaf() const { return lstLeaf; }

  long rootIdx() const { return 0; }
  long parentIdx(long i) const { return nodes.at(i).parent; }
  long leftChildIdx(long i) const { return nodes.at(i).leftChild; }
  long rightChildIdx(long i) const { return nodes.at(i).rightChild; }

  void printout(std::ostream& s, long idx = 0) const
  {
    s << "[" << aux << " ";
    s << nodes[idx].getData();
    if (leftChildIdx(idx) >= 0)
      printout(s, leftChildIdx(idx));
    if (rightChildIdx(idx) >= 0)
      printout(s, rightChildIdx(idx));
    s << "] ";
  }
  //  friend std::istream& operator>> (std::istream &s, DoubleCRT &d);

  //! If the parent is a leaf, add to it to children with the given data,
  //! else just update the data of the two children of this parent.
  //! Returns the index of the left child, the right-child index is one
  //! more than the left-child index.
  long addChildren(long prntIdx, const T& leftData, const T& rightData);

  //! Remove all nodes in the tree except for the root
  void collapseToRoot()
  {
    if (nodes.size() > 1) {
      nodes.resize(1);
      frstLeaf = lstLeaf = 0;
      nLeaves = 1;
    }
  }
};
template <typename T>
long FullBinaryTree<T>::addChildren(long prntIdx,
                                    const T& leftData,
                                    const T& rightData)
{
  assertInRange(prntIdx, 0l, (long)nodes.size(), "Parent node does not exist");

  // If parent is a leaf, add to it two children
  if (nodes[prntIdx].leftChild == -1 && nodes[prntIdx].rightChild == -1) {
    long childIdx = nodes.size();
    TreeNode<T> n1(leftData);
    nodes.push_back(n1); // add left child to std::vector
    TreeNode<T> n2(rightData);
    nodes.push_back(n2); // add right child to std::vector

    TreeNode<T>& parent = nodes[prntIdx];
    TreeNode<T>& left = nodes[childIdx];
    TreeNode<T>& right = nodes[childIdx + 1];

    parent.leftChild = childIdx; // point to children from parent
    parent.rightChild = childIdx + 1;
    left.parent = right.parent = prntIdx; // point to parent from children

    // remove parent and insert children to the linked list of leaves
    left.prev = parent.prev;
    left.next = childIdx + 1;
    right.prev = childIdx;
    right.next = parent.next;
    if (parent.prev >= 0) { // parent was not the 1st leaf
      nodes[parent.prev].next = childIdx;
      parent.prev = -1;
    } else // parent was the first leaf, now its left child is 1st
      frstLeaf = childIdx;

    if (parent.next >= 0) { // parent was not the last leaf
      nodes[parent.next].prev = childIdx + 1;
      parent.next = -1;
    } else // parent was the last leaf, now its left child is last
      lstLeaf = childIdx + 1;

    nLeaves++; // we replaced a leaf by a parent w/ two leaves
  } else {     // parent is not a leaf, update the two children
    TreeNode<T>& parent = nodes[prntIdx];
    assertTrue(parent.leftChild >= 0, "Left child does not exist");
    assertTrue(parent.rightChild >= 0, "Right child does not exist");

    TreeNode<T>& left = nodes[parent.leftChild];
    TreeNode<T>& right = nodes[parent.rightChild];
    left.data = leftData;
    right.data = rightData;
  }
  return nodes[prntIdx].leftChild;
}

//! A minimal description of a generator for the purpose of building tree
class GenDescriptor
{
public:
  long genIdx; // the index of the corresponding generator in Zm*/(p)
  long order;  // the order of this generator (or a smaller subcube)
  bool good;   // is this a good generator?

  GenDescriptor(long _order, bool _good, long gen = 0) :
      genIdx(gen), order(_order), good(_good)
  {}

  GenDescriptor() {}
};

//! A node in a tree relative to some generator
class SubDimension
{
  static const NTL::Vec<long> dummyBenes; // Useful for initialization
public:
  long size; // size of cube slice, same as 'order' in GenDescriptor
  bool good; // good or bad
  long e;    // shift-by-1 in this sub-dim is done via X -> X^{g^e}

  // A Benes leaf corresponds to either one or two Benes networks, depending
  // on whether or not it is in the middle. If this object is in the middle
  // then scndBenes.length()==0, else scndBenes.length()>=1. Each of the two
  // Benes network can be "trivial", i.e., collapsed to a single layer, which
  // is denoted by benes.length()==1.
  NTL::Vec<long> frstBenes;
  NTL::Vec<long> scndBenes;

  explicit SubDimension(long sz = 0,
                        bool gd = false,
                        long ee = 0,
                        const NTL::Vec<long>& bns1 = dummyBenes,
                        const NTL::Vec<long>& bns2 = dummyBenes)
  {
    size = sz;
    good = gd;
    e = ee;
    frstBenes = bns1;
    scndBenes = bns2;
  }

  long totalLength() const { return frstBenes.length() + scndBenes.length(); }
  /*  SubDimension& operator=(const SubDimension& other)
    { genIdx=other.genIdx; size=other.size;
      e=other.e; good=other.good; benes=other.benes;
      return *this;
      }
  */
  friend std::ostream& operator<<(std::ostream& s, const SubDimension& tree);
};
typedef FullBinaryTree<SubDimension> OneGeneratorTree; // tree for one generator

//! A std::vector of generator trees, one per generator in Zm*/(p)
class GeneratorTrees
{
  long depth; // How many layers in this permutation network
  NTL::Vec<OneGeneratorTree> trees;
  Permut map2cube, map2array;

public:
  GeneratorTrees() { depth = 0; } // default constructor

  // GeneratorTrees(const Vec<OneGeneratorTree>& _trees): trees(_trees)
  //  {depth=0;}
  //
  // Initialize trees with only the roots.
  //  GeneratorTrees(const Vec<SubDimension>& dims);

  long numLayers() const { return depth; } // depth of permutation network
  long numTrees() const { return trees.length(); }   // how many trees
  long getSize() const { return map2cube.length(); } // hypercube size

  OneGeneratorTree& operator[](long i) { return trees[i]; }
  const OneGeneratorTree& operator[](long i) const { return trees[i]; }
  OneGeneratorTree& at(long i) { return trees.at(i); }
  const OneGeneratorTree& at(long i) const { return trees.at(i); }

  OneGeneratorTree& getGenTree(long i) { return trees.at(i); }
  const OneGeneratorTree& getGenTree(long i) const { return trees.at(i); }

  const Permut& mapToCube() const { return map2cube; }
  const Permut& mapToArray() const { return map2array; }
  Permut& mapToCube() { return map2cube; }
  Permut& mapToArray() { return map2array; }

  long mapToCube(long i) const { return map2cube[i]; }
  long mapToArray(long i) const { return map2array[i]; }

  //! Get the "crude" cube dimensions corresponding to the vector of trees,
  //! the ordered vector with one dimension per tree
  void getCubeDims(NTL::Vec<long>& dims) const;

  //! Get the "fine" cube dimensions corresponding to the vector of trees,
  //! the ordered vector with one dimension per leaf in all the trees.
  void getCubeSubDims(NTL::Vec<long>& dims) const;

  // Returns coordinates of i relative to leaves of the tree
  //  void getCoordinates(Vec<long>&, long i) const;

  //! Compute the trees corresponding to the "optimal" way of breaking
  //! a permutation into dimensions, subject to some constraints. Returns
  //! the cost (# of 1D shifts) of this solution.
  //! Returns NTL_MAX_LONG if no solution
  long buildOptimalTrees(const NTL::Vec<GenDescriptor>& vec, long depthBound);

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

  friend std::ostream& operator<<(std::ostream& s, const GeneratorTrees& t);
};

// Permutation networks

class Ctxt;
class EncryptedArray;
class PermNetwork;

//! @class PermNetLayer
//! @brief The information needed to apply one layer of a permutation network
class PermNetLayer
{
  long genIdx; // shift-by-1 in this layer is done via X -> X^{g^e}
  long e;
  NTL::Vec<long> shifts; // shifts[i] is how much to shift slot i
  bool isID; // a silly optimization, does this layer compute the identity?

  friend class PermNetwork;
  friend std::ostream& operator<<(std::ostream& s, const PermNetwork& net);

public:
  long getGenIdx() const { return genIdx; }
  long getE() const { return e; }
  const NTL::Vec<long>& getShifts() const { return shifts; }
  bool isIdentity() const { return isID; }
};

//! A full permutation network
class PermNetwork
{
  NTL::Vec<PermNetLayer> layers;

  //! Compute one or more layers corresponding to one network of a leaf
  void setLayers4Leaf(long lyrIdx,
                      const ColPerm& p,
                      const NTL::Vec<long>& benesLvls,
                      long gIdx,
                      const SubDimension& leafData,
                      const Permut& map2cube);

public:
  PermNetwork(){}; // empty network
  PermNetwork(const Permut& pi, const GeneratorTrees& trees)
  {
    buildNetwork(pi, trees);
  }

  long depth() const { return layers.length(); }

  //! Take as input a permutation pi and the trees of all the generators,
  //! and prepares the permutation network for this pi
  void buildNetwork(const Permut& pi, const GeneratorTrees& trees);

  //! Apply network to permute a ciphertext
  void applyToCtxt(Ctxt& c, const EncryptedArray& ea) const;

  //! Apply network to array, used mostly for debugging
  void applyToCube(HyperCube<long>& v) const;

  //! Apply network to plaintext polynomial, used mostly for debugging
  void applyToPtxt(NTL::ZZX& p, const EncryptedArray& ea) const;

  const PermNetLayer& getLayer(long i) const { return layers[i]; }

  friend std::ostream& operator<<(std::ostream& s, const PermNetwork& net);
};

} // namespace helib

#endif // ifndef HELIB_PERMUTATIONS_H
