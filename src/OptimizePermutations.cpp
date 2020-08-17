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
/* @file OptimizePermutations.cpp
 * @brief Implementation of optimized permutation networks
 */

#include <cstdlib>
#include <list>
#include <sstream>
#include <memory>

#include <NTL/vector.h>
#include <helib/NumbTh.h>
#include <helib/EncryptedArray.h>
#include <helib/permutations.h>
#include <helib/apiAttributes.h>

namespace helib {

//! \cond FALSE (make doxygen ignore these classes)
template <typename T>
class ClassHash
{ // helper class to make defining hash functions easier
public:
  size_t operator()(const T& t) const { return t.hash(); }
};
//! \endcond

// routines for finding optimal level-collapsing strategies for Benes networks

// Removes duplicates from list x. Uses vector aux as scratch space
//   invariants:
//     * all elements in x lie in a range a..b
//     * aux[a..b] is false before and after removeDups is called
void removeDups(std::list<long>& x, bool* aux)
{
  for (std::list<long>::iterator i = x.begin(); i != x.end();) {
    if (aux[*i])
      i = x.erase(i);
    else {
      aux[*i] = true;
      i++;
    }
  }

  for (std::list<long>::iterator i = x.begin(); i != x.end(); i++)
    aux[*i] = false;
}

// Creates a new list with the old values and old values +/- offset.
// results outside the range -n+1 .. n-1 are discarded
// and all resulting duplicates are removed
void addOffset(std::list<long>& x,
               long offset,
               long n,
               bool* aux,
               UNUSED bool good = false)
{
  for (std::list<long>::iterator i = x.begin(); i != x.end(); i++) {
    long val = *i;
    long val1 = val + offset;
    long val2 = val - offset;
    if (val1 > -n && val1 < n)
      x.push_front(val1);
    if (val2 > -n && val2 < n)
      x.push_front(val2);
    /* FIXME: Alternative Impl:
     *        Replace the above two lines by
     *        if (good) {
     *           if (val1 < 0 && val1 > -n) val1 += n;
     *           if (val2 < 0 && val2 > -n) val2 += n;
     *        }
     *        if (val1 > -n && val1 < n && !aux[val1]) {
     *           x.push_front(val1);
     *           aux[val1] = true;
     *        }
     *        if (val2 > -n && val2 < n && !aux[val2]) {
     *           x.push_front(val2);
     *           aux[val2] = true;
     *        }
     * Then we can eliminate both removeDups and reducedCount
     */
  }

  removeDups(x, aux);
}

// Counts the number of unique elements mod n in x
long reducedCount(const std::list<long>& x, long n, bool* aux)
{
  long res = 0;

  for (std::list<long>::const_iterator i = x.begin(); i != x.end(); i++) {
    long val = *i;
    if (val < 0)
      val += n;

    if (!aux[val]) {
      res++;
      aux[val] = true;
    }
  }

  for (std::list<long>::const_iterator i = x.begin(); i != x.end(); i++) {
    long val = *i;
    if (val < 0)
      val += n;
    aux[val] = false;
  }

  return res;
}

// Compute the cost for all (n choose 2) possible ways to collapse levels.
// For j in [0..nlev-i) tab[i][j] is the cost of collapsing levels i..i+j.
// I.e., how many different shift amounts would we need to implement for
// a permutation-network-layer constructed by collapsing these levels.
void buildBenesCostTable(long n,
                         long k,
                         bool good,
                         NTL::Vec<NTL::Vec<long>>& tab)
{
  long nlev = 2 * k - 1;
  tab.SetLength(nlev);
  for (long i = 0; i < nlev; i++)
    tab[i].SetLength(nlev - i);

  NTL::Vec<bool> aux_vec;
  aux_vec.SetLength(2 * n - 1);
  bool* aux = &aux_vec[n - 1];
  for (long i = 0; i < 2 * n - 1; i++)
    aux_vec[i] = false;

  for (long i = 0; i < nlev; i++) {
    std::list<long> x;

    x.push_front(0L);
    for (long j = 0; j < nlev - i; j++) {
      long shamt = GeneralBenesNetwork::shamt(n, k, i + j);
      // The shift amount for this level

      addOffset(x, shamt, n, aux);
      if (good)
        tab[i][j] = reducedCount(x, n, aux) - 1;
      else
        tab[i][j] = x.size() - 1;
      // FIXME: Alternative Impl:
      //        Replace the 5 lines above by
      //        addOffset(x, shamt, n, aux, good);
      //        tab[i][j] = x.size() - 1;
      //        Also initialize aux to false in every iteration
    }
  }
}

//! \cond FALSE (make doxygen ignore these classes)
class LongNode;
typedef std::shared_ptr<LongNode> LongNodePtr;
// A "shared_ptr" is a pointer with (some) garbage collection

// A LongNode is an implementation of std::list<long>, i.e. a linked list of
// counters, representing a particular way of "collapsing levels" in a Benes
// network. Each LongNode holds a count of collapsed levels, and the sum of all
// counter in the list must equal the number of levels in the Benes network.
class LongNode
{
public:
  long count;       // number of levels collapsed
  LongNodePtr next; // next node in the list

  LongNode(long _count, LongNodePtr _next)
  {
    count = _count;
    next = _next;
  }
};
//! \endcond

// Returns the length of the linked list, starting with this pointer
static long length(LongNodePtr ptr)
{
  long res = 0;
  for (LongNodePtr p = ptr; p != nullptr; p = p->next)
    res++;
  return res;
}

// Converts std::list<long> to NTL::Vec<long>
static long listToVec(NTL::Vec<long>& vec, LongNodePtr ptr)
{
  long len = length(ptr);

  vec.SetLength(len);
  long i = 0;
  for (LongNodePtr p = ptr; p != nullptr; p = p->next) {
    vec[i] = p->count;
    i++;
  }
  return len;
}

// Prints out a list of integers
std::ostream& operator<<(std::ostream& s, LongNodePtr p)
{
  if (p == nullptr)
    return s << "[]";

  s << "[" << p->count;
  for (p = p->next; p != nullptr; p = p->next)
    s << " " << p->count;
  return s << "]";
}

// Data structures to hold the memoization table for the dynamic-programming
// computation of the level collapsing in a Benes network. An entry in the
// table is specified by the first-level index i and the allotted budget
// (depth). Once the computation is over, this entry will contain the optimal
// way of breaking this network (encoded as a LongNode list) and the cost of
// this network.
//! \cond FALSE (make doxygen ignore these classes)
class BenesMemoKey
{
public:
  long i;      // Index of 1st level to consider
  long budget; // max-depth of network after collapsing

  BenesMemoKey(long _i, long _budget)
  {
    i = _i;
    budget = _budget;
  }

  size_t hash() const
  {
    std::stringstream s;
    s << i << " " << budget;
    return std::hash<std::string>()(s.str());
  }

  bool operator==(const BenesMemoKey& other) const
  {
    return i == other.i && budget == other.budget;
  }
};

class BenesMemoEntry
{
public:
  long cost;            // Cost of the optimal solution
  LongNodePtr solution; // The solution itself, as a list of LongNode's

  BenesMemoEntry(long _cost, LongNodePtr _solution)
  {
    cost = _cost;
    solution = _solution;
  }

  BenesMemoEntry() : cost(0) {}
};

typedef std::
    unordered_map<BenesMemoKey, BenesMemoEntry, ClassHash<BenesMemoKey>>
        BenesMemoTable;
//! \endcond

// A dynamic program (implemented as a recursive routine with memoization) for
// computing the optimal collapsing of layers in a generalized Benes network.
//
// The parameter i indicates the sub-problem at hand, namely we are trying
// to find the optimal way of collapsing the sub-network from level i to
// the end (level nlev).
//
// The budget is the largest number of levels that we can afford after
// collapsing. (For example, if budget >= nlev-i, then there is no need to
// collapse anything.)

BenesMemoEntry optimalBenesAux(long i,
                               long budget,
                               long nlev,
                               const NTL::Vec<NTL::Vec<long>>& costTab,
                               BenesMemoTable& memoTab)
{
  assertInRange<InvalidArgument>(i,
                                 0l,
                                 nlev,
                                 "Level to collapse index out of bound",
                                 true);
  assertTrue<InvalidArgument>(budget > 0, "No budget left");

  // Did we already solve this problem? If so just return the solution.
  BenesMemoTable::iterator find = memoTab.find(BenesMemoKey(i, budget));
  if (find != memoTab.end()) {
    return find->second;
  }

  // a new subproblem to process...

  long cost;
  LongNodePtr solution;

  if (i == nlev) { // An empty network, nothing to collapse, trivial solution
    cost = 0;
    solution = LongNodePtr();
  } else if (budget == 1) {
    // Almost no budget is left, the only possible solution is to collapse
    // all remaining levels together.

    cost = costTab[i][nlev - i - 1];
    solution = LongNodePtr(new LongNode(nlev - i, LongNodePtr()));
  } else {
    // We consider collapsing levels i..i+j, for all j in [0..nlev-i),
    // then choose the option that gives the best cost.

    long bestCost = NTL_MAX_LONG;
    long bestJ = 0;
    BenesMemoEntry bestT;
    for (long j = 0; j < nlev - i; j++) {
      // If we collapse levels i..i+j, then we need to compute the optimal
      // way of collapsing the rest, i.e. level i+j+1 to the end.
      BenesMemoEntry t =
          optimalBenesAux(i + j + 1, budget - 1, nlev, costTab, memoTab);

      // The total cost of this solution is the cost of the collapsed level
      // (i..i+j) itself, plus the cost of the optimal solution for the rest.
      if (t.cost + costTab[i][j] < bestCost) {
        bestCost = t.cost + costTab[i][j];
        bestJ = j;
        bestT = t;
      }
    }

    cost = bestCost;
    solution = LongNodePtr(new LongNode(bestJ + 1, bestT.solution));
  }

  return memoTab[BenesMemoKey(i, budget)] = BenesMemoEntry(cost, solution);
}

// Computes an optimal level-collapsing strategy for a Benes network
//   n = the size of the network
//   budget = an upper bound on the number of levels in the collapsed network
//   good = flag indicating whether this is with respect to a "good" generator,
//     for which shifts by i and i-n correspond to the same rotation
//   cost = total number of shifts needed by the collapsed network
//   solution = list indicating how levels in a standard benes network
//      are collapsed: if solution = [s_1 s_2 ... s_k], then k <= budget,
//      and the first s_1 levels are collapsed, the next s_2 levels
//      are collapsed, etc.
void optimalBenes(long n,
                  long budget,
                  bool good,
                  long& cost,
                  LongNodePtr& solution)
{
  long k = GeneralBenesNetwork::depth(n); // k = ceiling(log_2 n)
  long nlev = 2 * k - 1; // before collapsing, we have 2k-1 levels

  NTL::Vec<NTL::Vec<long>> costTab;
  // costTab[i][j] to holds the cost of collapsing levels i..i+j

  buildBenesCostTable(n, k, good, costTab);
  // Compute the cost for all (n choose 2) possible ways to collapse levels.

  BenesMemoTable memoTab;
  BenesMemoEntry t = optimalBenesAux(0, budget, nlev, costTab, memoTab);
  // Compute the optimal collapsing of layers in a width-n Benes network

  cost = t.cost;
  solution = t.solution;
}

/********************************************************************/
/***** Routines for finding the optimal splits among generators *****/
/********************************************************************/

//! \cond FALSE (make doxygen ignore these classes)

// A binary tree data structure for splitting generators
class SplitNode;
typedef std::shared_ptr<SplitNode> SplitNodePtr;
// A "std::shared_ptr" is a pointer with (some) garbage collection

class SplitNode
{
public:
  long order; // order associated with this node
  long mid;   // the "middle" token, 0 or 1
  bool good;  // the "good" flag
  LongNodePtr solution1, solution2;
  // The benes solution(s).
  // For non-middle nodes there may be two different solutions,
  // depending on how the budget is split
  SplitNodePtr left, right; // the children (for internal nodes)

  SplitNode(long _order,
            long _mid,
            bool _good,
            LongNodePtr _solution1,
            LongNodePtr _solution2)
  {
    // constructor for leaves
    order = _order;
    mid = _mid;
    good = _good;
    solution1 = _solution1;
    solution2 = _solution2;
    left = right = SplitNodePtr();
  }

  SplitNode(long _order,
            long _mid,
            bool _good,
            SplitNodePtr _left,
            SplitNodePtr _right)
  {
    // constructor for internal nodes
    order = _order;
    mid = _mid;
    good = _good;
    solution1 = solution2 = LongNodePtr();
    left = _left;
    right = _right;
  }

  bool isLeaf() const { return left == nullptr && right == nullptr; }
};
//! \endcond

// Routines for printing the leaves of a generator-tree
void print(std::ostream& s, SplitNodePtr p, bool first)
{
  if (p->isLeaf()) {
    if (!first)
      s << " ";
    s << "[";
    if (p->mid == 1)
      s << "*";
    if (p->good)
      s << "g ";
    else
      s << "b ";
    s << p->order << " " << p->solution1 << " " << p->solution2 << "]";
  } else {
    print(s, p->left, first);
    print(s, p->right, false);
  }
}
std::ostream& operator<<(std::ostream& s, SplitNodePtr p)
{
  s << "[";
  print(s, p, true);
  s << "]";
  return s;
}

// Data structures to hold the memory table for the dynamic-programming
// computation optimizing a generator tree. An entry in the table is specified
// by (order,good-flag,budget,middle-flag). When the computation is over, this
// entry will contain the optimal tree for a single generator (encoded as a
// SplitNode tree) and the cost of this solution.
//! \cond FALSE (make doxygen ignore these classes)
class LowerMemoKey
{
public:
  long order;  // size of hypercube corresponding to this sub-tree
  bool good;   // is this a good subtree
  long budget; // max-depth of network(s) associated to this node
  long mid;    // whether or not this is the middle layer of the network

  LowerMemoKey(long _order, bool _good, long _budget, long _mid)
  {
    order = _order;
    good = _good;
    budget = _budget;
    mid = _mid;
  }

  size_t hash() const
  {
    std::stringstream s;
    s << order << " " << good << " " << budget << " " << mid;
    return std::hash<std::string>()(s.str());
  }

  bool operator==(const LowerMemoKey& other) const
  {
    return order == other.order && good == other.good &&
           budget == other.budget && mid == other.mid;
  }
};

class LowerMemoEntry
{
public:
  long cost;
  SplitNodePtr solution;

  LowerMemoEntry(long _cost, SplitNodePtr _solution)
  {
    cost = _cost;
    solution = _solution;
  }

  LowerMemoEntry() : cost(0) {}
};

typedef std::
    unordered_map<LowerMemoKey, LowerMemoEntry, ClassHash<LowerMemoKey>>
        LowerMemoTable;

// list structure for managing generators

class GenNode;
typedef std::shared_ptr<GenNode> GenNodePtr;
// A "shared_ptr" is a pointer with (some) garbage collection

class GenNode
{
public:
  SplitNodePtr solution; // the solution tree for this generator
  GenNodePtr next;       // next node in the list

  GenNode(SplitNodePtr _solution, GenNodePtr _next)
  {
    solution = _solution;
    next = _next;
  }
};
//! \endcond

// Compute the length of a list
long length(GenNodePtr ptr)
{
  long res = 0;
  for (GenNodePtr p = ptr; p != nullptr; p = p->next)
    res++;
  return res;
}

std::ostream& operator<<(std::ostream& s, GenNodePtr p)
{
  if (p == nullptr) {
    s << "[]";
    return s;
  }

  s << "[" << p->solution;
  p = p->next;
  while (p != nullptr) {
    s << " " << p->solution;
    p = p->next;
  }
  s << "]";

  return s;
}

// upper level memo table

// Data structures to hold the memory table for optimizing the partition to
// separate generator tree. An entry in the table is specified by the index
// of a generator, the budget allocated to this tree, and the middle flag.
// When the computation is over, this entry will contain the optimal list
// of trees (encoded as a GenNode list) and the cost of this solution.
//! \cond FALSE (make doxygen ignore these classes)
class UpperMemoKey
{
public:
  long i;
  long budget;
  long mid;

  UpperMemoKey(long _i, long _budget, long _mid)
  {
    i = _i;
    budget = _budget;
    mid = _mid;
  }

  size_t hash() const
  {
    std::stringstream s;
    s << i << " " << budget << " " << mid;
    return std::hash<std::string>()(s.str());
  }

  bool operator==(const UpperMemoKey& other) const
  {
    return i == other.i && budget == other.budget && mid == other.mid;
  }
};

class UpperMemoEntry
{
public:
  long cost;
  GenNodePtr solution;

  UpperMemoEntry(long _cost, GenNodePtr _solution)
  {
    cost = _cost;
    solution = _solution;
  }

  UpperMemoEntry() : cost(0) {}
};

typedef std::
    unordered_map<UpperMemoKey, UpperMemoEntry, ClassHash<UpperMemoKey>>
        UpperMemoTable;
//! \endcond

// Optimize a single tree: try all possible ways of splitting the order into
// order1*order2 (and also the solution of not splitting at all). For every
// possible split, try all budget allocations and allocations of good and mid.
LowerMemoEntry optimalLower(long order,
                            bool good,
                            long budget,
                            long mid,
                            LowerMemoTable& lowerMemoTable)
{
  assertTrue<InvalidArgument>(order > 1, "Order must be greater than 1");
  assertTrue<InvalidArgument>(mid == 0 || mid == 1, "mid value is not 1 or 2");
  assertTrue<InvalidArgument>(budget > 0, "No budget left");

  // Did we already solve this problem? If so just return the solution.
  LowerMemoTable::iterator find =
      lowerMemoTable.find(LowerMemoKey(order, good, budget, mid));

  if (find != lowerMemoTable.end()) {
    return find->second;
  }

  long cost;
  SplitNodePtr solution;

  if (mid == 0 && budget == 1) {
    // Insufficient budget, no solution is possible

    cost = NTL_MAX_LONG;
    solution = SplitNodePtr();
  } else {
    // First calculate a solution without splitting, making this a leaf node
    LongNodePtr benesSolution1, benesSolution2;

    if (mid == 1) {
      // this is the middle node, so just one Benes network

      optimalBenes(order, budget, good, cost, benesSolution1);
      benesSolution2 = LongNodePtr();
    } else {
      // not the middle node, so we need two Benes networks.
      // if budget is odd, we split it unevenly

      long cost1, cost2;
      optimalBenes(order, budget / 2, good, cost1, benesSolution1);
      if (budget % 2 == 0) { // both networks have the same budget
        cost2 = cost1;
        benesSolution2 = benesSolution1;
      } else { // one network has budget larger by one than the other
        optimalBenes(order, budget - budget / 2, good, cost2, benesSolution2);
      }

      cost = cost1 + cost2;
    }

    // The initial solution corresponds to this being a leaf node
    solution = SplitNodePtr(
        new SplitNode(order, mid, good, benesSolution1, benesSolution2));

    // Recursively try all the ways of splitting order = order1 * order2

    for (long order1 = 2; order1 < order; order1++) {
      // try all factors of order
      // there is some redundancy here, since we consider
      // both the splits (d, n/d) and (n/d, d); however,
      // we utilize this redundancy (see below) to streamline
      // allocation of the "good" token (if we have it)

      // continue;  // temporary hack to force benes for testing
      if (order % order1 != 0)
        continue;

      bool good1 = good;
      bool good2 = good;

      if (good && NTL::GCD(order1, order / order1) != 1)
        good2 = false;

      // The logic is that if the problem is "good"
      // but the split is not relatively prime,
      // then only one of the subproblems is good.
      // since order1 ranges over all factors of order,
      // it suffices to choose the left subproblem to be
      // the good one.

      // Try all possible ways of splitting the budget between the two nodes
      for (long budget1 = 1; budget1 < budget; budget1++) {

        // A slick way of giving the middle token to one of the
        // nodes if we have it, and to none of the nodes if we don't
        for (long mid1 = 0; mid1 <= mid; mid1++) {
          LowerMemoEntry s1 =
              optimalLower(order1, good1, budget1, mid1, lowerMemoTable);
          // FIXME: If s1.cost==NTL_MAX_LONG we do not need to compute
          //        the cost of s2
          LowerMemoEntry s2 = optimalLower(order / order1,
                                           good2,
                                           budget - budget1,
                                           mid - mid1,
                                           lowerMemoTable);
          if (s1.cost != NTL_MAX_LONG && s2.cost != NTL_MAX_LONG &&
              s1.cost + s2.cost < cost) {
            cost = s1.cost + s2.cost;
            // Record the new best solution
            solution = SplitNodePtr(
                new SplitNode(order, mid, good, s1.solution, s2.solution));
          }
        }
      }
    }
  }

  // Record the best solution in the lowerMemoTable and return it
  return lowerMemoTable[LowerMemoKey(order, good, budget, mid)] =
             LowerMemoEntry(cost, solution);
}

// Optimizing a list of trees, trying all the ways of allocating the budget
// and mid token between the trees. This procedure splits the "current
// remaining budget" between trees i through vec.length()-1.
UpperMemoEntry optimalUpperAux(const NTL::Vec<GenDescriptor>& vec,
                               long i,
                               long budget,
                               long mid,
                               UpperMemoTable& upperMemoTable,
                               LowerMemoTable& lowerMemoTable)
{
  assertInRange<InvalidArgument>(i,
                                 0l,
                                 (long)vec.length(),
                                 "Index i does not point to "
                                 "a tree (index out of range)",
                                 true);
  assertTrue<InvalidArgument>(budget >= 0, "Negative budget");
  assertTrue<InvalidArgument>(mid == 0 || mid == 1, "mid value is not 1 or 2");

  // FIXME: Did we already solve this problem? If so just return the solution.
  UpperMemoTable::iterator find =
      upperMemoTable.find(UpperMemoKey(i, budget, mid));
  if (find != upperMemoTable.end()) {
    return find->second;
  }

  long cost;
  GenNodePtr solution;

  if (i == vec.length()) {
    // An empty list, recursion stops here with trivial solution

    cost = 0;
    solution = GenNodePtr();
  } else if (budget == 0) {
    // recursion stops here with no solution

    cost = NTL_MAX_LONG;
    solution = GenNodePtr();
  } else {
    // allocate resources (budget, mid) between generator i and the rest

    long bestCost = NTL_MAX_LONG;
    LowerMemoEntry bestS;
    UpperMemoEntry bestT;

    for (long budget1 = 1; budget1 <= budget; budget1++) {
      for (long mid1 = 0; mid1 <= mid; mid1++) {
        // Optimize the first tree (index i) with the allotted budget1, mid1
        LowerMemoEntry s = optimalLower(vec[i].order,
                                        vec[i].good,
                                        budget1,
                                        mid1,
                                        lowerMemoTable);
        // FIXME: If s.cost==NTL_MAX_LONG we do not need to compute
        //        the cost of t

        // Optimize the rest of the list with the remaining budget and mid
        UpperMemoEntry t = optimalUpperAux(vec,
                                           i + 1,
                                           budget - budget1,
                                           mid - mid1,
                                           upperMemoTable,
                                           lowerMemoTable);
        if (s.cost != NTL_MAX_LONG && t.cost != NTL_MAX_LONG &&
            s.cost + t.cost < bestCost) {
          bestCost = s.cost + t.cost;
          bestS = s;
          bestT = t;
        }
      }
    }

    cost = bestCost;
    if (cost == NTL_MAX_LONG)
      solution = GenNodePtr();
    else
      solution = GenNodePtr(new GenNode(bestS.solution, bestT.solution));
  }

  // Record the best solution in the upperMemoTable and return it
  return upperMemoTable[UpperMemoKey(i, budget, mid)] =
             UpperMemoEntry(cost, solution);
}

/**
 * optimalUpper is the high-level routine used to find
 * an optimal strategy for evaluating a permutation
 * The inputs are:
 *   * vec: a vector describing the generators, each entry
 *       consist of a long/bool pair (order, good)
 *   * budget: the total budget for levels in the network,
 *       which corresponds to the (multiplicative) circuit depth
 * The outputs are:
 *   * cost: the number of rotations needed (a specific
 *       permutation could use fewer)
 *   * solution: the optimal strategy, described below
 *
 * The solution is a linked list of length equal to the length of vec.
 * The ith entry of the list corresponds to the ith entry in vec.
 * Each entry in the list is a binary tree (SplitNodePtr).
 * Every node in such a tree is either a leaf or an internal node
 *   with two children.
 * Every node contains the order associated with that node, along with
 *   a "mid" flag: 1 means this is a "middle" node, 0 means not.
 * In the top-level list of trees, only one tree root will be a middle node,
 *   and for every internal middle node in any tree, only one child
 *   will be a middle node.
 * A leaf node in a tree will contain two "benes network collapsing solutions",
 *   solution1 and solution2.
 *   If this leaf node is a middle node, solution2 is empty.
 *   Otherwise, there are actually two solutions: this is because
 *   the benes network corresponding to this leaf will appear twice
 *   in the overall network, and if the budget allocated to this leaf
 *   is *odd*, we may want to split the budget unevenly between the two
 *   networks (floor(budget/2) and floor(budget/2)+1.
 *   NOTE: one could consider more general splits of the budget
 *   between these two occurrences of the network, but it seems doubtful
 *   to be of benefit.
 * A "benes network collapsing solution" is a list of longs (longNodePtr).
 *   If the list is [n_1 n_2 ... n_k], this means we are to collapse
 *   the first n_1 levels of the Benes network, then the next n_2 levels,
 *   and so on.  The resulting "collapsed" Benes network will have k
 *   levels.
 *********************************************************************/

// Returns the total number of layers in the corresponding (sub)network
static long recursiveCopy2Tree(OneGeneratorTree& gTree,
                               long nodeIdx,
                               SplitNodePtr& solution)
{
  if (solution->isLeaf()) { // Copy the Benes levels
    return listToVec(gTree[nodeIdx].getData().frstBenes, solution->solution1) +
           listToVec(gTree[nodeIdx].getData().scndBenes, solution->solution2);
  }

  // Copy left and right children into the tree. A child with mid-token
  // is copied last, if none has the mid-token then copy them in order.

  SplitNodePtr leftChld, rightChld;
  if (solution->left->mid) {
    rightChld = solution->left;
    leftChld = solution->right;
  } else {
    leftChld = solution->left;
    rightChld = solution->right;
  }
  SubDimension leftData(leftChld->order, leftChld->good);
  SubDimension rightData(rightChld->order, rightChld->good);
  gTree.addChildren(nodeIdx, leftData, rightData);
  return recursiveCopy2Tree(gTree, gTree.leftChildIdx(nodeIdx), leftChld) +
         recursiveCopy2Tree(gTree, gTree.rightChildIdx(nodeIdx), rightChld);
}

// A recursive procedure for computing the e exponent values
static void computeEvalues(const OneGeneratorTree& T, long idx, long genOrd)
{
  // if either child is missing then we are at a leaf (this is a full tree)
  long left = T[idx].getLeftChild();
  long right = T[idx].getRightChild();
  if (left < 0 || right < 0)
    return;

  SubDimension& lData = (SubDimension&)T[left].getData();
  SubDimension& rData = (SubDimension&)T[right].getData();

  // The entire tree size is sz1 * sz2;
  long sz1 = lData.size; // size of left subtree
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
    long f2 = sz1 * sz2 + 1 - f1; // f2 = 1 mod sz1, 0 mod sz2
    lData.e = NTL::MulMod(ee, f2, genOrd);
    rData.e = NTL::MulMod(ee, f1, genOrd);
  }
  // Recurse on the two subtrees
  computeEvalues(T, left, genOrd);
  computeEvalues(T, right, genOrd);
}

// Copy data from a SplitNode tree to a OneGeneratorTree
static long copyToGenTree(OneGeneratorTree& gTree, SplitNodePtr& solution)
{
  SubDimension rootData(solution->order, solution->good, /*e=*/1);
  gTree.putDataInRoot(rootData);
  long len = recursiveCopy2Tree(gTree, gTree.rootIdx(), solution);
  computeEvalues(gTree, gTree.rootIdx(), solution->order);
  return len;
}

// Compute the trees corresponding to the "optimal" way of breaking
// a permutation into dimensions, subject to some constraints
long GeneratorTrees::buildOptimalTrees(const NTL::Vec<GenDescriptor>& gens,
                                       long depthBound)
{
  // TODO: is this check necessary?
  assertTrue<InvalidArgument>(gens.length() >= 0, "negative gens size");
  trees.SetLength(gens.length()); // allocate space if needed

  if (gens.length() == 0) {
    map2cube.SetLength(1, 0);
    map2array.SetLength(1, 0);
    return 0;
  }
  assertTrue<InvalidArgument>(depthBound > 0, "Zero or negative depthBound");

  // reset the trees, starting from only the roots
  for (long i = 0; i < trees.length(); i++) {
    if (trees[i].getNleaves() > 1) // tree is not empty/trivial
      trees[i].collapseToRoot();
  }

  UpperMemoTable upperMemoTable;
  LowerMemoTable lowerMemoTable;

  // Compute a solution in { t.cost, t.solution }
  UpperMemoEntry t =
      optimalUpperAux(gens, 0, depthBound, 1, upperMemoTable, lowerMemoTable);

  // Copy the solution into the trees
  GenNodePtr midPtr;
  long i = 0, treeIdx = 0, midIdx = 0;
  depth = 0; // Also compute the depth of the permutation network
  for (GenNodePtr genPtr = t.solution; genPtr != nullptr;
       genPtr = genPtr->next, i++) {
    if (genPtr->solution->mid) { // Keep the "middle tree" for last
      midPtr = genPtr;
      midIdx = i;
      continue;
    }
    depth += copyToGenTree(trees[treeIdx], genPtr->solution);
    trees[treeIdx].setAuxKey(gens[i].genIdx);
    treeIdx++;
  }
  if (!midPtr) {
    // no solution, undo initialization and return
    // NTL_MAX_LONG

    depth = 0;
    trees.kill();
    map2cube.kill();
    map2array.kill();
    return NTL_MAX_LONG;
  }

  depth += copyToGenTree(trees[treeIdx], midPtr->solution);
  trees[treeIdx].setAuxKey(gens[midIdx].genIdx);

  // Compute the mapping from array to cube and back
  ComputeCubeMapping();

#ifdef HELIB_DEBUG
  NTL::Vec<long> dims; // The "crude" cube dimensions, one dimension per tree
  getCubeDims(dims);
  std::cerr << " dims=" << dims << std::endl;
  std::cerr << " trees=" << *this << std::endl;
  if (map2cube.length() < 100) {
    std::cerr << " map2cube=" << map2cube << std::endl;
    std::cerr << " map2array=" << map2array << std::endl;
  }
  std::cerr << std::endl;
#endif

  return t.cost;
}

} // namespace helib
