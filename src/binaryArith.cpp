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
 * @file binaryArith.cpp
 * @brief Implementing integer addition, multiplication in binary representation
 */
#include <numeric>
#include <climits>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <atomic>
#include <mutex> // std::mutex, std::unique_lock

#include <NTL/BasicThreadPool.h>
#include <helib/binaryArith.h>

#ifdef HELIB_DEBUG
#include <cstdio>
#include <helib/debugging.h>
#endif

#define BPL_ESTIMATE (30)
// FIXME: this should really be dynamic

namespace helib {

#ifdef HELIB_DEBUG
void decryptAndSum(std::ostream& s,
                   const CtPtrMat& numbers,
                   bool twosComplement = false);
#endif

typedef std::pair<long, long> NodeIdx; // nodes are indexed by a pair (i,j)

/**
 * @class DAGnode
 * @brief A node in an addition-DAG structure.
 *
 **/
class DAGnode
{
public:
  NodeIdx idx; // the indexes of this node
  bool isQ;    // if not then isP
  long level;  // The level at the time of computation

  std::atomic_long childrenLeft; // how many children were not computed yet
  DAGnode *parent1, *parent2;

  std::mutex ct_mtx; // controls access to ctxt pointer (and the ctxt itself)
  Ctxt* ct;          // points to the actual ciphertext (or nullptr)

  DAGnode(NodeIdx ii,
          bool qq,
          long lvl,
          long chl = 0,
          DAGnode* pt1 = nullptr,
          DAGnode* pt2 = nullptr) :
      idx(ii),
      isQ(qq),
      level(lvl),
      childrenLeft(chl),
      parent1(pt1),
      parent2(pt2),
      ct(nullptr)
  {}

  DAGnode(DAGnode&& other) :
      // move constructor
      idx(other.idx),
      isQ(other.isQ),
      level(other.level),
      childrenLeft(long(other.childrenLeft)), // copy value of atomic_long
      parent1(other.parent1),
      parent2(other.parent2),
      ct(other.ct)
  {}

  std::string nodeName() const
  {
    return (std::string(isQ ? "Q(" : "P(") + std::to_string(idx.first) + ',' +
            std::to_string(idx.second) + ')');
  }
};

//! A class to help manage the allocation of temporary Ctxt objects
class ScratchCell
{
public:
  std::atomic_bool used;
  std::unique_ptr<Ctxt> ct; // scratch space owns this pointer
  ScratchCell(const Ctxt& c) : used(true), ct(new Ctxt(ZeroCtxtLike, c)) {}
  ScratchCell(ScratchCell&& other) :
      // move constructor
      used(bool(other.used)),
      ct(std::move(other.ct))
  {}
};

/**
 * @class AddDAG
 * @brief A class representing the logic of the order of bit products when
 *        adding two integers.
 *
 * Given two input arrays a[], b[], we build a DAG with each node representing
 * a term either of the form p_{i,j} = \prod_{t=j}^i (a[t]+b[t]), or of the
 * form q_{i,j} = (a[j]*b[j]) * \prod_{t=j+1}^i (a[t]+b[t]). The source nodes
 * are of the forms (a[i]*b[i]) and (a[i]+b[i]), and each non-source node has
 * exactly two parents, whose product yields that node.
 *
 * When building the DAG, we keep the level of each node as high as possible.
 * For example we can set q_{i,j}=p_{i,k}*q_{k-1,j} or q_{i,j}=p_{i,k+1}*q_{k,j}
 * (among other options), and we choose the option that results in the highest
 * level. In addition, we try to minimize the number of nodes in the DAG that
 * actually need to be computed while adding the two numbers (subject to still
 * consuming as few levels as possible).
 **/
class AddDAG
{
  std::mutex scratch_mtx;           // controls access to scratch vector
  std::vector<ScratchCell> scratch; // scratch space for ciphertexts
  std::map<NodeIdx, DAGnode> p;     // p[i,j]= prod_{t=j}^i (a[t]+b[t])
  std::map<NodeIdx, DAGnode> q; // q[i,j]= a[j]b[j]*prod_{t=j+1}^i (a[t]+b[t])
  long aSize, bSize;

  Ctxt* allocateCtxtLike(const Ctxt& c); // Allocate a new ciphertext if needed
  void markAsAvailable(DAGnode* node);   // Mark temporary Ctxt object as unused
  const Ctxt& getCtxt(DAGnode* node,     // Compute a new Ctxt if need
                      const CtPtrs& a,
                      const CtPtrs& b);

  // Add to c the Ctxt from the given node
  void addCtxtFromNode(Ctxt& c, DAGnode* node, const CtPtrs& a, const CtPtrs& b)
  {
    std::unique_lock<std::mutex> lck(node->ct_mtx);
    c += getCtxt(node, a, b);
    if (--(node->childrenLeft) == 0)
      markAsAvailable(node);
  }

public:
  //! Build a plan to add a and b
  void init(const CtPtrs& a, const CtPtrs& b);

  // Build the addition DAG
  AddDAG(const CtPtrs& a, const CtPtrs& b) { init(a, b); }

  //! Perform the actual addition
  void apply(CtPtrs& sum, const CtPtrs& a, const CtPtrs& b, long sizeLimit = 0);

  //! Returns the lowest level in this DAG
  long lowLvl() const
  {
    if (aSize < 1)
      return 0;
    return findQ(bSize - 1, 0)->level;
  }
  //! Returns a pointer to the a 'p' node of index (i,j)
  DAGnode* findP(long i, long j) const
  { // returns nullptr if not exists
    auto it = p.find(NodeIdx(i, j));
    if (it == p.end()) {
      std::cerr << "  findP(" << i << ',' << j << ") not found" << std::endl;
      return nullptr; // not found
    }
    return (DAGnode*)&(it->second);
  }
  //! Returns a pointer to the a 'q' node of index (i,j)
  DAGnode* findQ(long i, long j) const
  { // returns nullptr if not exists
    auto it = q.find(NodeIdx(i, j));
    if (it == q.end()) {
      std::cerr << "  findQ(" << i << ',' << j << ") not found" << std::endl;
      return nullptr; // not found
    }
    return (DAGnode*)&(it->second);
  }
#ifdef HELIB_DEBUG
  void printAddDAG(bool printCT = false);
#endif
};

// When searching for a good middle point, we use a "good default", so that
// in the cases covered by that default solution we get the smallest number
// of nodes with childrenLeft>0.
// When initializing p[i,j] the default should be p[i,i+1-2^e]*p[i-2^e,j]
// where e is the largest exponent with 2^e <= i-j.
inline long defaultPmiddle(long delta)
{
  return 1 << (NTL::NumBits(delta) - 1);
}
// When initializing q[i,j] the default should be p[i,i+1-2^e]*q[i-2^e,j]
// where e is the largest exponent with 2^e+2^{e-1} <= i-j.
inline long defaultQmiddle(long delta)
{
  delta = (2 * delta + 2) / 3; // ceil(2 delta / 3)
  return 1 << (NTL::NumBits(delta) - 1);
}

//! Build a plan to add a and b
void AddDAG::init(const CtPtrs& aa, const CtPtrs& bb)
{
  // make sure that lsize(b) >= lsize(a)
  const CtPtrs& a = (lsize(bb) >= lsize(aa)) ? aa : bb;
  const CtPtrs& b = (lsize(bb) >= lsize(aa)) ? bb : aa;

  aSize = lsize(a);
  bSize = lsize(b);
  assertTrue<InvalidArgument>(aSize >= 1, "a must not be empty");

  // Initialize the p[i,i]'s and q[i,i]'s
  p.clear();
  q.clear();
  for (long i = 0; i < bSize; i++) {
    NodeIdx idx(i, i);
    // The level of b[i]
    long lvl =
        (b.isSet(i) && !(b[i]->isEmpty())) ? b[i]->bitCapacity() : LONG_MAX;
    if (i < aSize) {
      // The level of a[i]
      long aLvl =
          (a.isSet(i) && !(a[i]->isEmpty())) ? a[i]->bitCapacity() : LONG_MAX;
      lvl = std::min(lvl, aLvl);
      if (lvl == LONG_MAX ||
          aLvl == LONG_MAX) // is either a[i] or b[i] is empty
        q.emplace(idx, DAGnode(idx, true, LONG_MAX, 1));
      else
        q.emplace(idx, DAGnode(idx, true, lvl - 1, 1));
    }
    p.emplace(idx, DAGnode(idx, false, lvl, 1));
  }

  // Initialize p[i,j] for bSize>=i>j>0
  for (long delta = 1; delta < bSize; delta++)
    for (long i = bSize - 1; i >= delta; --i) {
      long j = i - delta;
      long mid = i - defaultPmiddle(delta); // initialize to a "good default"
      DAGnode* prnt2 = findP(i, mid + 1);
      DAGnode* prnt1 = findP(mid, j);
      long maxLvl = std::min(prnt1->level, prnt2->level) - BPL_ESTIMATE;
      if (prnt1->level == LONG_MAX ||
          prnt2->level == LONG_MAX) // parent is empty
        maxLvl = LONG_MAX;
      long maxN =
          std::min(long(prnt1->childrenLeft), long(prnt2->childrenLeft));
      for (long m = j; m < i; m++) { // find middle point maximizing lvl(p[i,j])
        if (m == mid)
          continue;
        DAGnode* p2 = findP(i, m + 1);
        DAGnode* p1 = findP(m, j);
        long lvl = std::min(p1->level, p2->level) - BPL_ESTIMATE;
        if (p1->level == LONG_MAX || p2->level == LONG_MAX) // parent is empty
          lvl = LONG_MAX;
        long n = std::min(long(p1->childrenLeft), long(p2->childrenLeft));
        if (lvl > maxLvl || (lvl == maxLvl && n > maxN)) {
          maxLvl = lvl;
          maxN = n;
          prnt1 = p1;
          prnt2 = p2;
        }
      }
      NodeIdx idx(i, j);
      p.emplace(idx, DAGnode(idx, false, maxLvl, 0, prnt1, prnt2));
      prnt1->childrenLeft++;
      prnt2->childrenLeft++;
    }

  // Initialize q[i,j] for bSize>=i>j>=0
  for (long delta = 1; delta < bSize; delta++)
    for (long i = bSize - 1; i >= delta; --i) {
      long j = i - delta;
      if (j >= aSize)
        continue;
      long maxLvl = 0, maxN = 0;
      long mid = i - defaultQmiddle(delta); // initialize to a "good default"
      DAGnode* prnt2 = findP(i, mid + 1);
      DAGnode* prnt1 = findQ(mid, j);
      if (prnt1 != nullptr) {
        maxLvl = std::min(prnt1->level, prnt2->level) - BPL_ESTIMATE;
        if (prnt1->level == LONG_MAX ||
            prnt2->level == LONG_MAX) // parent is empty
          maxLvl = LONG_MAX;
        maxN = long(prnt2->childrenLeft);
      }
      for (long m = j; m < i; m++) { // find middle point maximizing lvl(p[i,j])
        if (m == mid)
          continue;
        DAGnode* p2 = findP(i, m + 1);
        DAGnode* p1 = findQ(m, j);
        if (p1 == nullptr)
          continue;
        long lvl = std::min(p1->level, p2->level) - BPL_ESTIMATE;
        if (p1->level == LONG_MAX || p2->level == LONG_MAX) // parent is empty
          lvl = LONG_MAX;
        long n = long(p2->childrenLeft);
        if (lvl > maxLvl || (lvl == maxLvl && n > maxN)) {
          maxLvl = lvl;
          maxN = n;
          prnt1 = p1;
          prnt2 = p2;
        }
      }
      if (prnt1 == nullptr)
        continue; // cannot create node
      NodeIdx idx(i, j);
      q.emplace(idx, DAGnode(idx, true, maxLvl, 1, prnt1, prnt2));
      prnt1->childrenLeft++;
      prnt2->childrenLeft++;
    }
}

//! Apply the DAG to actually compute the sum
void AddDAG::apply(CtPtrs& sum,
                   const CtPtrs& aa,
                   const CtPtrs& bb,
                   long sizeLimit)
{
  // make sure that lsize(b) >= lsize(a)
  const CtPtrs& a = (lsize(bb) >= lsize(aa)) ? aa : bb;
  const CtPtrs& b = (lsize(bb) >= lsize(aa)) ? bb : aa;
  if (aSize != lsize(a) || bSize != lsize(b))
    throw LogicError("DAG applied to wrong vectors");

  if (sizeLimit == 0)
    sizeLimit = bSize + 1;
  if (lsize(sum) != sizeLimit)
    sum.resize(sizeLimit, &b); // allocate space for the output
  for (long i = 0; i < lsize(sum); i++)
    sum[i]->clear();

  // Allow multi-threading in this loop
  NTL_EXEC_RANGE(sizeLimit, first, last)
  for (long i = first; i < last; i++) { //  for (long i=0; i<sizeLimit; i++) {
    if (i < bSize)
      addCtxtFromNode(*(sum[i]), this->findP(i, i), a, b);
    for (long j = std::min(i - 1, aSize - 1); j >= 0; --j) {
      DAGnode* node = this->findQ(i - 1, j);
      if (node != nullptr)
        addCtxtFromNode(*(sum[i]), node, a, b);
    }
  }
  NTL_EXEC_RANGE_END
}

//! Get the ciphertext for a node, computing it as needed
const Ctxt& AddDAG::getCtxt(DAGnode* node, const CtPtrs& a, const CtPtrs& b)
{
  // NOTE: node->ct_mtx should be locked before calling this function

  if (node->ct == nullptr) { // ciphertext not computed yet, do it now
    if (node->parent1 != nullptr && node->parent2 != nullptr) { // internal node
      // Obtain locks and ciphertexts for both parents. Also reduce the
      // number of children of that parents that still need to be computed
      std::unique_lock<std::mutex> pt1_lck(node->parent1->ct_mtx);
      const Ctxt& c1 = getCtxt(node->parent1, a, b);
      long n1 = --(node->parent1->childrenLeft);

      std::unique_lock<std::mutex> pt2_lck(node->parent2->ct_mtx);
      const Ctxt& c2 = getCtxt(node->parent2, a, b);
      long n2 = --(node->parent2->childrenLeft);

      if (n1 == 0) { // reuse space from parent1
        node->parent1->ct = nullptr;
        node->ct = (Ctxt*)&c1;
        if (c1.isEmpty() || c2.isEmpty())
          node->ct->clear(); // ct is zero if any of the parents is
        else
          node->ct->multiplyBy(c2);
        if (n2 == 0)
          markAsAvailable(node->parent2);
      } else if (n2 == 0) { // reuse space from parent2
        node->parent2->ct = nullptr;
        node->ct = (Ctxt*)&c2;
        if (c1.isEmpty() || c2.isEmpty())
          node->ct->clear(); // ct is zero if any of the parents is
        else
          node->ct->multiplyBy(c1);
      } else { // allocate new space
        node->ct = allocateCtxtLike(c2);
        if (c1.isEmpty() || c2.isEmpty())
          node->ct->clear(); // ct is zero if any of the parents is
        else {
          *(node->ct) = c2;
          node->ct->multiplyBy(c1);
        }
      }
    } else { // no parents, either a[i]+b[i] or a[i]*b[i]
      long i = node->idx.first;
      long j = node->idx.second; // we expect i==j
      const Ctxt* ct_ptr = b.ptr2nonNull();
      assertNotNull(ct_ptr, "ct_ptr must not be null");
      node->ct = allocateCtxtLike(*ct_ptr);

      if (node->isQ) { // This is b[i]*a[j]
        if (b.isSet(i) && !(b[i]->isEmpty()) && a.isSet(j) &&
            !(a[j]->isEmpty())) {
          *(node->ct) = *(b[i]);
          node->ct->multiplyBy(*(a[j]));
        } // if a[j] or b[i] is empty then node->ct is a zero ciphertext
        else
          node->ct->clear();
      } else { // This is b[i]+a[i]
        if (!b.isSet(i) || b[i]->isEmpty())
          node->ct->clear();
        else
          *(node->ct) = *(b[i]);
        if (a.isSet(j) && !(a[j]->isEmpty()))
          *(node->ct) += *(a[j]);
      }
    } // end of no-parents case
  }
  return *(node->ct);
}

//! Adds another cell to scratch space, or use an existing one that's free
Ctxt* AddDAG::allocateCtxtLike(const Ctxt& c)
{
  // look for an unused cell in the scratch array
  for (long i = 0; i < lsize(scratch); i++)
    if (scratch[i].used == false) { // found a free one, try to use it
      bool used = scratch[i].used.exchange(true); // mark it as used
      if (used == false) // make sure no other thread got there first
        return scratch[i].ct.get();
    }

  // If not found, allocate a new cell
  ScratchCell sc(c);      // cell points to new ctxt, with used=true
  Ctxt* pt = sc.ct.get(); // remember the raw pointer
  std::unique_lock<std::mutex> lck(scratch_mtx); // protect scratch vector
  scratch.emplace_back(std::move(sc));           // scratch now owns the pointer
  return pt;                                     // return the raw pointer
}

// Mark a scratch ciphertext as unused. We assume that no two nodes
// ever share a ciphertext object, so if this node is done with the
// object then the object is unused.
void AddDAG::markAsAvailable(DAGnode* node)
{
  // NOTE: node->ct_mtx should be locked before calling this function
  // NOTE: somewhat inefficient, use linear search for the raw pointer
  for (long i = 0; i < (long)scratch.size(); i++)
    if (scratch[i].ct.get() == node->ct)
      scratch[i].used = false;
  node->ct = nullptr;
}

/********************************************************************/
/********************************************************************/

// Use packed bootstrapping, so we can bootstrap all in just one go.
void packedRecrypt(const CtPtrs& a,
                   const CtPtrs& b,
                   std::vector<zzX>* unpackSlotEncoding)
{
  const Ctxt* ct = b.ptr2nonNull(); // find some non-null Ctxt
  if (ct == nullptr)
    ct = a.ptr2nonNull();
  if (ct == nullptr)
    return; // nothing to do

  assertNotNull<InvalidArgument>(unpackSlotEncoding,
                                 "unpackSlotEncoding must not be null");
  assertTrue(ct->getPubKey().isBootstrappable(),
             "public key must be bootstrappable for recryption");

  struct CtPtrs_pair : CtPtrs
  {
    const CtPtrs& a;
    const CtPtrs& b;
    CtPtrs_pair(const CtPtrs& _a, const CtPtrs& _b) : a(_a), b(_b) {}
    Ctxt* operator[](long i) const override
    {
      return (i < lsize(a)) ? a[i] : b[i - lsize(a)];
    }
    long size() const override { return lsize(a) + lsize(b); }
  };
  const CtPtrs_pair ab(a, b);

  packedRecrypt(ab, *unpackSlotEncoding, *(ct->getContext().ea));
}

// Return a number as a vector of bits with little endian-ness
std::vector<long> longToBitVector(long num, long bitSize)
{
  assertTrue<InvalidArgument>(bitSize >= 0, "bitSize must be non-negative.");
  std::vector<long> result;
  for (long i = 0; i < bitSize; num >>= 1, ++i)
    result.push_back(num & 1);
  return result;
}

//! Apply mask across the vector of bits slot-wise.
void binaryMask(CtPtrs& bits, const Ctxt& mask)
{
  for (long i = 0; i < bits.size(); ++i)
    bits[i]->multiplyBy(mask);
}

//! Implementation of output = cond ? trueValue : falseValue
void binaryCond(CtPtrs& output,
                const Ctxt& cond,
                const CtPtrs& trueValue,
                const CtPtrs& falseValue)
{
  assertEq(trueValue.size(),
           falseValue.size(),
           "trueValue and falseValue must have the same size.");
  assertEq(output.size(),
           falseValue.size(),
           "output and input vectors must have the same size.");
  vecCopy(output, trueValue);
  std::vector<Ctxt> falseCopy;
  vecCopy(falseCopy, falseValue);
  binaryMask(output, cond);
  Ctxt negated_cond(cond);
  negated_cond.addConstant(NTL::ZZX(1L));
  CtPtrs_vectorCt falseCopyWrapper(falseCopy);
  binaryMask(falseCopyWrapper, negated_cond);
  // TODO: Change this to use bit-wise XOR
  for (long i = 0; i < output.size(); ++i)
    *output[i] += falseCopy[i];
}

//! Concatenate two binary numbers into a single `CtPtrs` object.
void concatBinaryNums(CtPtrs& output, const CtPtrs& a, const CtPtrs& b)
{
  assertEq(output.size(),
           a.size() + b.size(),
           "output must be of size a.size() + b.size()");
  for (long i = 0; i < a.size(); ++i)
    *output[i] = *a[i];
  for (long i = 0; i < b.size(); ++i)
    *output[i + a.size()] = *b[i];
}

//! Split a binary number into two separate binary numbers.
void splitBinaryNums(CtPtrs& leftSplit, CtPtrs& rightSplit, const CtPtrs& input)
{
  assertEq(leftSplit.size() + rightSplit.size(),
           input.size(),
           "Output sizes must sum to input.size()");
  for (long i = 0; i < leftSplit.size(); ++i)
    *leftSplit[i] = *input[i];
  for (long i = 0; i < rightSplit.size(); ++i)
    *rightSplit[i] = *input[i + leftSplit.size()];
}

//! Shift binary numbers to the left by `shamt`
void leftBitwiseShift(CtPtrs& output, const CtPtrs& input, const long shamt)
{
  assertTrue(shamt >= 0, "Shift amount must be positive.");
  assertEq(output.size(),
           input.size(),
           "output and input must have the same size.");
  for (long i = 0; i < output.size() - shamt; ++i)
    *output[i + shamt] = *input[i];
  for (long i = 0; i < shamt; ++i)
    output[i]->clear();
}

//! Rotate binary numbers by `rotamt`.
void bitwiseRotate(CtPtrs& output, const CtPtrs& input, long rotamt)
{
  assertEq(output.size(),
           input.size(),
           "output and input must be the same size.");
  long bitSize = input.size();
  rotamt = mcMod(rotamt, bitSize);
  for (long i = 0; i < output.size(); ++i)
    *output[i] = *input[mcMod(i - rotamt, bitSize)];
}

//! Compute a bitwise XOR between `lhs` and `rhs`.
void bitwiseXOR(CtPtrs& output, const CtPtrs& lhs, const CtPtrs& rhs)
{
  assertEq(output.size(), lhs.size(), "output and lhs must be the same size.");
  assertEq(lhs.size(), rhs.size(), "lhs and rhs must be the same size.");
  vecCopy(output, lhs);
  for (long i = 0; i < rhs.size(); ++i)
    *output[i] += *rhs[i];
}

//! Compute a bitwise OR between `lhs` and `rhs`.
void bitwiseOr(CtPtrs& output, const CtPtrs& lhs, const CtPtrs& rhs)
{
  assertEq(output.size(), lhs.size(), "output and lhs must be the same size.");
  assertEq(lhs.size(), rhs.size(), "lhs and rhs must be the same size.");
  // a OR b is equivalent to (ab + a + b)
  // Start output off as lhs & rhs
  bitwiseAnd(output, lhs, rhs);
  // Now add on lhs and rhs
  for (long i = 0; i < rhs.size(); ++i) {
    *output[i] += *lhs[i];
    *output[i] += *rhs[i];
  }
}

//! Compute a bitwise AND between `lhs` and `rhs`.
void bitwiseAnd(CtPtrs& output, const CtPtrs& lhs, const CtPtrs& rhs)
{
  assertEq(output.size(), lhs.size(), "output and lhs must be the same size.");
  assertEq(lhs.size(), rhs.size(), "lhs and rhs must be the same size.");
  vecCopy(output, lhs);
  for (long i = 0; i < rhs.size(); ++i)
    output[i]->multiplyBy(*rhs[i]);
}

//! Compute a bitwise AND between `input` and `mask`.
void bitwiseAnd(CtPtrs& output,
                const CtPtrs& input,
                const std::vector<long> mask)
{
  assertEq(output.size(),
           input.size(),
           "output and input must be the same size.");
  vecCopy(output, input);
  for (long i = 0; i < output.size(); ++i)
    if (!mask[i])
      output[i]->clear();
}

//! Compute a bitwise NOT of `input`.
void bitwiseNot(CtPtrs& output, const CtPtrs& input)
{
  assertEq(output.size(),
           input.size(),
           "input and output must have the same size");
  vecCopy(output, input);
  for (long i = 0; i < output.size(); ++i)
    output[i]->addConstant(NTL::ZZ(1L));
}

//! Add two integers in binary representation
void addTwoNumbers(CtPtrs& sum,
                   const CtPtrs& lhs,
                   const CtPtrs& rhs,
                   long sizeLimit,
                   std::vector<zzX>* unpackSlotEncoding)
{
  HELIB_TIMER_START;
  if (lsize(lhs) < 1) {
    vecCopy(sum, rhs, sizeLimit);
    return;
  } else if (lsize(rhs) < 1) {
    vecCopy(sum, lhs, sizeLimit);
    return;
  }

  // Work out the order of multiplications to compute all the carry bits
  AddDAG addPlan(lhs, rhs);

#ifdef HELIB_DEBUG // print plan
  addPlan.printAddDAG();
#endif

  // Ensure that we have enough levels to compute everything,
  // bootstrap otherwise
  if (addPlan.lowLvl() < BPL_ESTIMATE) {
    packedRecrypt(lhs, rhs, unpackSlotEncoding);
    addPlan.init(lhs, rhs);                // Re-compute the DAG
    if (addPlan.lowLvl() < BPL_ESTIMATE) { // still not enough levels
      throw LogicError("not enough levels for addition DAG");
    }
  }
  addPlan.apply(sum, lhs, rhs, sizeLimit); // perform the actual addition
}

// Negate a binary number that is already in 2's complement. Note: input must
// not alias negation.
void negateBinary(CtPtrs& negation, const CtPtrs& input)
{
  assertEq(negation.size(), input.size(), "Arguments must have matching size.");
  std::vector<Ctxt> bitFlippedInput;
  vecCopy(bitFlippedInput, input);
  // First flip all bits of the input.
  for (auto& bit : bitFlippedInput)
    bit.addConstant(NTL::ZZX(1L));
  // Deep copy of input into negation.
  vecCopy(negation, bitFlippedInput);
  // Now add one.
  negation[0]->addConstant(NTL::ZZX(1L));
  // Calculate the resultant carry bits.
  std::vector<Ctxt>& carryBits = bitFlippedInput;
  incrementalProduct(carryBits);
  for (std::size_t i = 1; i < bitFlippedInput.size(); ++i)
    *(negation[i]) += carryBits[i - 1];
}

// Subtract rhs from lhs and put the result in difference.
void subtractBinary(CtPtrs& difference,
                    const CtPtrs& lhs,
                    const CtPtrs& rhs,
                    std::vector<zzX>* unpackSlotEncoding)
{
  assertEq(lhs.size(), rhs.size(), "Size of lhs and rhs must be the same.");
  assertEq(difference.size(),
           rhs.size(),
           "Size of output vector must equal the size of the input vectors.");
  // Negate the rhs and then use the existing add function.
  std::vector<Ctxt> negated_rhs(rhs.size(), *rhs[0]);
  CtPtrs_vectorCt negated_wrapper(negated_rhs);
  negateBinary(negated_wrapper, rhs);
  addTwoNumbers(difference,
                lhs,
                negated_wrapper,
                lhs.size(),
                unpackSlotEncoding);
}

// Return pointers to the three inputs, ordered by size
static std::tuple<const CtPtrs*, const CtPtrs*, const CtPtrs*> orderBySize(
    const CtPtrs& a,
    const CtPtrs& b,
    const CtPtrs& c)
{
  if (lsize(a) <= lsize(b)) {
    if (lsize(b) <= lsize(c))
      return std::make_tuple(&a, &b, &c); // a <= b <= c
    else if (lsize(a) <= lsize(c))
      return std::make_tuple(&a, &c, &b); // a <= c < b
    else
      return std::make_tuple(&c, &a, &b); // c < a <= b
  } else {                                // lsize(b) < lsize(a)
    if (lsize(a) <= lsize(c))
      return std::make_tuple(&b, &a, &c); // b < a <= c
    else if (lsize(b) <= lsize(c))
      return std::make_tuple(&b, &c, &a); // b <= c < a
    else
      return std::make_tuple(&c, &b, &a); // c < b < a
  }
}

// Implementing the basic 3-for-2 trick: u,v,w encrypt bits, return two bits
// x,y such that x+2y = u+v+w over the integers. Outputs can alias the inputs.
static void three4Two(Ctxt& lsb,
                      Ctxt& msb,
                      const Ctxt& u,
                      const Ctxt& v,
                      const Ctxt& w)
{
  Ctxt tmp_v = v;
  Ctxt tmp_w = w;
  lsb = u;
  msb = u;

  lsb += tmp_v;          // u+v
  msb.multiplyBy(tmp_v); // u*v

  tmp_v = lsb; // u+v

  tmp_v.multiplyBy(tmp_w); // (u+v)*w

  lsb += tmp_w; // u+v+w
  msb += tmp_v; // u*v + (u+v)*w = u*v + u*w + v*w
}

// Same as three4Two above, but some of the inputs could be null.
// Returns the number of output bits that are not identically zero.
static long three4Two(Ctxt* lsb, Ctxt* msb, Ctxt* u, Ctxt* v, Ctxt* w)
{
  if (u != nullptr && !u->isEmpty() && v != nullptr && !v->isEmpty() &&
      w != nullptr && !w->isEmpty()) { // if none are empty
    three4Two(*lsb, *msb, *u, *v, *w); // call the function above
    return 2;
  }
  if ((u == nullptr || u->isEmpty()) && (v == nullptr || v->isEmpty()) &&
      (w == nullptr || w->isEmpty())) { // if all are empty
    lsb->clear();                       // result is empty too
    msb->clear();
    return 0;
  }

  // Some are empty, others are not, arrange so that empty are at the end
  if (u == nullptr || u->isEmpty()) {
    if (v == nullptr || v->isEmpty())
      u = w; // only w was non-empty
    else {
      u = v;
      v = w;
    } // v,w were non-empty
  } else if (v == nullptr || v->isEmpty())
    v = w;     // u is non-empty, v was empty
  w = nullptr; // we don't use w anymore

  if (v == nullptr || v->isEmpty()) { // only u is non-empty
    *lsb = *u;
    msb->clear();
    return 1;
  }

  // both u,v are non-empty
  Ctxt tmp = *v;
  *lsb = *u;
  *msb = *u;
  *lsb += tmp;
  msb->multiplyBy(tmp);
  return 2;
}

// Apply the 3-for-2 routine to integers (i.e., an array of bits). The
// inputs need not be of the same size, and size of the output x is
// equal to the largest of them, and the size of the output y is one
// larger. This is safe even when the outputs alias some of the inputs
static void three4Two(CtPtrs& lsb,
                      CtPtrs& msb,
                      const CtPtrs& u,
                      const CtPtrs& v,
                      const CtPtrs& w,
                      long sizeLimit)
{
  HELIB_TIMER_START;
  // Arrange u,v,w by size from smallest to largest
  const CtPtrs *p1, *p2, *p3;
  std::tie(p1, p2, p3) = orderBySize(u, v, w); // size(p3)>=size(p2)>=size(p1)

  if (p3->size() <= 0) { // empty input
    setLengthZero(lsb);
    setLengthZero(msb);
    return;
  }
  if (p1->size() <= 0) { // two or less inputs
    std::vector<Ctxt> tmp;
    vecCopy(tmp, *p2, sizeLimit); // just in case p2, msb share pointers
    vecCopy(msb, *p3, sizeLimit);
    vecCopy(lsb, tmp);
    return;
  }
  if (sizeLimit == 0)
    sizeLimit = p3->size() + 1;

  // Allocate space in the output vectors

  const Ctxt* ctptr = p3->ptr2nonNull();
  std::vector<Ctxt> tmpMsb, tmpLsb;

  long lsbSize = std::min(sizeLimit, lsize(*p3));
  long msbSize = lsbSize;
  if (lsize(*p2) == lsize(*p3) && lsbSize < sizeLimit)
    msbSize++; // possible carry out of last position

  resize(tmpLsb, lsbSize, Ctxt(ZeroCtxtLike, *ctptr));
  resize(tmpMsb, msbSize, Ctxt(ZeroCtxtLike, *ctptr));

  NTL_EXEC_RANGE(msbSize - 1, first, last)
  for (long i = first; i < last; i++) {
    if (i < lsize(*p1))
      three4Two(&tmpLsb[i], &tmpMsb[i + 1], (*p1)[i], (*p2)[i], (*p3)[i]);
    else if (i < lsize(*p2)) {
      three4Two(&tmpLsb[i], &tmpMsb[i + 1], (*p2)[i], (*p3)[i], nullptr);
    } else if (p3->isSet(i))
      tmpLsb[i] = *((*p3)[i]);
  }
  NTL_EXEC_RANGE_END

  if (msbSize == lsbSize) { // we only computed upto lsbSize-1, do the last LSB
    if (p1->isSet(lsbSize - 1))
      tmpLsb[lsbSize - 1] = *((*p1)[lsbSize - 1]);
    if (p2->isSet(lsbSize - 1))
      tmpLsb[lsbSize - 1] += *((*p2)[lsbSize - 1]);
    if (p3->isSet(lsbSize - 1))
      tmpLsb[lsbSize - 1] += *((*p3)[lsbSize - 1]);
  }
  vecCopy(lsb, tmpLsb);
  vecCopy(msb, tmpMsb);
}

//! @brief An implementation of PtrMatrix using vector< PtrVector<T>* >
template <typename T>
struct PtrMatrix_PtPtrVector : PtrMatrix<T>
{
  std::vector<PtrVector<T>*>& rows;
  PtrMatrix_PtPtrVector(std::vector<PtrVector<T>*>& mat) : rows(mat) {}
  PtrVector<T>& operator[](long i) override // returns a row
  {
    return *rows[i];
  }
  const PtrVector<T>& operator[](long i) const override // returns a row
  {
    return *rows[i];
  }
  long size() const override { return lsize(rows); } // How many rows
};

// Calculates the sum of many numbers using the 3-for-2 method
void addManyNumbers(CtPtrs& sum,
                    CtPtrMat& numbers,
                    long sizeLimit,
                    std::vector<zzX>* unpackSlotEncoding)
{
#ifdef HELIB_DEBUG
  std::cout << " addManyNumbers: " << numbers.size()
            << " numbers with size-limit=" << sizeLimit << std::endl;
#endif
  HELIB_TIMER_START;
  const Ctxt* ct_ptr = numbers.ptr2nonNull();
  if (lsize(numbers) < 1 || ct_ptr == nullptr) { // nothing to add
    setLengthZero(sum);
    return;
  }
  if (lsize(numbers) == 1) {
    vecCopy(sum, numbers[0]);
    return;
  }

  bool bootstrappable = ct_ptr->getPubKey().isBootstrappable();
  const EncryptedArray& ea = *(ct_ptr->getContext().ea);

  long leftInQ = lsize(numbers);
  std::vector<CtPtrs*> numPtrs(leftInQ);
  for (long i = 0; i < leftInQ; i++)
    numPtrs[i] = &(numbers[i]);

  // use 3-for-2 repeatedly until only two numbers are leff to add
  while (leftInQ > 2) {
    // If any number is too low level, then bootstrap everything
    PtrMatrix_PtPtrVector<Ctxt> wrapper(numPtrs);
    if (findMinBitCapacity(wrapper) < 3 * ct_ptr->getContext().BPL()) {
      assertNotNull<InvalidArgument>(unpackSlotEncoding,
                                     "unpackSlotEncoding must not be null");
      assertTrue(bootstrappable,
                 "public key must be bootstrappable for recryption");

      packedRecrypt(wrapper, *unpackSlotEncoding, ea, /*belowLvl=*/10);
    }
    // Prepare a vector for pointers to the output of this iteration
    long nTriples = leftInQ / 3;
    long leftOver = leftInQ - (3 * nTriples);
    std::vector<CtPtrs*> numPtrs2(2 * nTriples + leftOver);

    if (leftOver > 0) { // copy the leftover pointers
      numPtrs2[0] = numPtrs[3 * nTriples];
      if (leftOver > 1)
        numPtrs2[1] = numPtrs[3 * nTriples + 1];
    }
    // Allow multi-threading in this loop
    //    NTL_EXEC_RANGE(nTriples, first, last)
    //    for (long i=first; i<last; i++) {   // call the three-for-two
    //    procedure
    for (long i = 0; i < nTriples; i++) { // call the three-for-two procedure
      three4Two(*numPtrs[3 * i],
                *numPtrs[3 * i + 1], // three4Two works in-place
                *numPtrs[3 * i],
                *numPtrs[3 * i + 1],
                *numPtrs[3 * i + 2],
                sizeLimit);

      numPtrs2[leftOver + 2 * i] = numPtrs[3 * i]; // copy the output pointers
      numPtrs2[leftOver + 2 * i + 1] = numPtrs[3 * i + 1];
    }
    //    NTL_EXEC_RANGE_END
    numPtrs.swap(numPtrs2);   // swap input/output vectors
    leftInQ = lsize(numPtrs); // update the size
  }
  // final addition
  addTwoNumbers(sum, *numPtrs[0], *numPtrs[1], sizeLimit, unpackSlotEncoding);
}

// Multiply a positive a by a potentially negative b, we need to sign-extend b
static void multByNegative(CtPtrs& product,
                           const CtPtrs& a,
                           const CtPtrs& b,
                           long sizeLimit,
                           std::vector<zzX>* unpackSlotEncoding)
{
  HELIB_TIMER_START;
  long resSize = lsize(a) + lsize(b);
  if (sizeLimit > 0 && sizeLimit < resSize)
    resSize = sizeLimit;

  NTL::Vec<NTL::Vec<Ctxt>> numbers(NTL::INIT_SIZE, std::min(lsize(a), resSize));
  long nNums = lsize(numbers);
  for (long i = 0; i < nNums; i++)
    numbers[i].SetLength(resSize, Ctxt(ZeroCtxtLike, *(a[0])));

  std::vector<std::pair<long, long>> pairs;
  for (long i = 0; i < nNums; i++)
    for (long j = i; j < resSize; j++)
      if (j < i + lsize(b) && a.isSet(i) && !a[i]->isEmpty() &&
          b.isSet(j - i) && !b[j - i]->isEmpty()) {
        pairs.push_back(std::pair<long, long>(i, j));
      }
  long nPairs = lsize(pairs);

  NTL_EXEC_RANGE(nPairs, first, last)
  for (long idx = first; idx < last; idx++) {
    long i, j;
    std::tie(i, j) = pairs[idx];
    numbers[i][j] = *(b[j - i]);
    numbers[i][j].multiplyBy(*(a[i])); // multiply by the bit of a
  }
  NTL_EXEC_RANGE_END

  // sign extension
  for (long i = 0; i < nNums; i++)
    for (long j = i + lsize(b); j < resSize; j++) {
      numbers[i][j] = numbers[i][i + lsize(b) - 1]; // sign extension
    }

  CtPtrMat_VecCt nums(numbers); // Wrapper around numbers
#ifdef HELIB_DEBUG
  long pa, pb;
  std::vector<long> slots;
  decryptBinaryNums(slots, a, *dbgKey, *dbgEa, false);
  pa = slots[0];
  decryptBinaryNums(slots, b, *dbgKey, *dbgEa, true);
  pb = slots[0];
  decryptAndSum((std::cout << " multByNegative: " << pa << '*' << pb << " = "),
                nums,
                true);
#endif
  addManyNumbers(product, nums, resSize, unpackSlotEncoding);
}

// Multiply two integers (i.e. an array of bits) lhs, rhs.
// Computes the pairwise products x_{i,j} = lhs_i * rhs_j
// then sums the prodcuts using the 3-for-2 method.
void multTwoNumbers(CtPtrs& product,
                    const CtPtrs& lhs,
                    const CtPtrs& rhs,
                    bool rhsTwosComplement,
                    long sizeLimit,
                    std::vector<zzX>* unpackSlotEncoding)
{
  HELIB_TIMER_START;
  long lhsSize = lsize(lhs);
  long rhsSize = lsize(rhs);
  long resSize = lhsSize + rhsSize;
  if (sizeLimit > 0 && sizeLimit < resSize)
    resSize = sizeLimit;

  if (lhs.numNonNull() < 1 || rhs.numNonNull() < 1) {
    setLengthZero(product);
    return; // return 0
  }

#ifdef HELIB_DEBUG
  std::cout << " before multiplication, capacity="
            << findMinBitCapacity({&lhs, &rhs}) << std::endl;
#endif
  // Edge case, if lhs or rhs is 1 bit
  if (lhsSize == 1) {
    if (lhs[0]->isEmpty()) {
      setLengthZero(product);
      return;
    }
    vecCopy(product, rhs, resSize);
    for (long i = 0; i < resSize; i++)
      product[i]->multiplyBy(*(lhs[0]));
    return;
  }
  if (rhsTwosComplement) { // somewhat different implementation for 2s
                           // complement
    multByNegative(product, lhs, rhs, sizeLimit, unpackSlotEncoding);
    return;
  }
  if (rhsSize == 1) {
    if (rhs[0]->isEmpty()) {
      setLengthZero(product);
      return;
    }
    vecCopy(product, lhs, resSize);
    for (long i = 0; i < resSize; i++)
      lhs[i]->multiplyBy(*(rhs[0]));
    return;
  }

  // We make sure temp_lhs is the larger of the two integers
  // to keep the number of additions to a minimum
  const CtPtrs& temp_lhs = (lhsSize >= rhsSize) ? lhs : rhs;
  const CtPtrs& temp_rhs = (lhsSize >= rhsSize) ? rhs : lhs;
  lhsSize = lsize(temp_lhs);
  rhsSize = lsize(temp_rhs);

  NTL::Vec<NTL::Vec<Ctxt>> numbers(NTL::INIT_SIZE,
                                   std::min(lsize(rhs), resSize));
  const Ctxt* ct_ptr = lhs.ptr2nonNull();
  long nNums = lsize(numbers);
  for (long i = 0; i < nNums; i++)
    numbers[i].SetLength(std::min((i + lhsSize), resSize),
                         Ctxt(ZeroCtxtLike, *ct_ptr));
  std::vector<std::pair<long, long>> pairs;
  for (long i = 0; i < nNums; i++)
    for (long j = i; j < lsize(numbers[i]); j++) {
      if (lhs.isSet(j - i) && !(lhs[j - i]->isEmpty()) && rhs.isSet(i) &&
          !(rhs[i]->isEmpty()))
        pairs.push_back(std::pair<long, long>(i, j));
    }
  long nPairs = lsize(pairs);
  NTL_EXEC_RANGE(nPairs, first, last)
  for (long idx = first; idx < last; idx++) {
    long i, j;
    std::tie(i, j) = pairs[idx];
    numbers[i][j] = *(lhs[j - i]);
    numbers[i][j].multiplyBy(*(rhs[i])); // multiply by the bit of rhs
  }
  NTL_EXEC_RANGE_END

  CtPtrMat_VecCt nums(numbers); // A wrapper around numbers
#ifdef HELIB_DEBUG
  long plaintext_lhs, plaintext_rhs;
  std::vector<long> slots;
  decryptBinaryNums(slots, lhs, *dbgKey, *dbgEa, false);
  plaintext_lhs = slots[0];
  decryptBinaryNums(slots, rhs, *dbgKey, *dbgEa, false);
  plaintext_rhs = slots[0];
  decryptAndSum((std::cout << " multTwoNumbers: " << plaintext_lhs << '*'
                           << plaintext_rhs << " = "),
                nums,
                false);
#endif
  addManyNumbers(product, nums, resSize, unpackSlotEncoding);
}

/* seven4Three: adding seven input bits, getting a 3-bit counter
 *
 * input: in[6..0]
 * ----------------
 *         in[6]
 *       b2 b1 = sum of in[2..0] (b2=msb, b1=lsb)
 *       b4 b3 = sum of in[5..3]
 * ------------
 *       c2 c1 = sum of in[6],b1,b3
 *    c4 c3    = sum of b2,b4
 * ------------
 *    d2 d1 c1 = out[2..0]
 */
// The output Ctxts[0..2] must be initialized, can alias the inputs
static void seven4Three(const CtPtrs& out, const CtPtrs& in, long sizeLimit)
{
  // we need 4 scratch ciphertexts
  std::vector<Ctxt> tmp(4, *out[0]);

  // Aliases, referring to the scheme above. Aliases for temporary
  // vars chosen so that inputs, outputs of three4two are distinct

  Ctxt& c1 = *out[0];
  Ctxt& d1 = *out[1];
  Ctxt& d2 = *out[2];

  Ctxt& b1 = tmp[0];
  Ctxt& b2 = tmp[1];
  Ctxt& b3 = tmp[2];
  Ctxt& b4 = tmp[3];
  Ctxt& c2 = d1;
  Ctxt& c3 = b1;
  Ctxt& c4 = b3;

  three4Two(&b1, &b2, in[0], in[1], in[2]); // b2 b1 = 3for2(in[0..2])
  three4Two(&b3, &b4, in[3], in[4], in[5]); // b4 b3 = 3for2(in[3..5])

  three4Two(&c1, &c2, in[6], &b1, &b3); // c2 c1 = 3for2(in[6],b1,b3)
  if (sizeLimit < 2)
    return;
  c3 = b2;

  c3 += b4; // c3 = b2 ^ b4
  c4 = b2;

  c4.multiplyBy(b4); // c4 = b2 * b4
  d2 = c2;

  d1 += c3; // d1 = c2 ^ c3 (d1 alias c2)
  if (sizeLimit < 3)
    return;
  d2.multiplyBy(c3);

  d2 += c4; // d2 = c4 ^ (c2*c3)
}

/* fifteen4Four: adding fifteen input bits, getting a 4-bit counter
 *
 * input: in[14..0]
 * ----------------
 *       b2 b1 = sum of in[2..0] (b2=msb, b1=lsb)
 *       b4 b3 = sum of in[5..3]
 *       b6 b5 = sum of in[8..6]
 *       b8 b7 = sum of in[11..9]
 *      b10 b9 = sum of in[14..12]
 * ------------
 *       b8 b7
 *      b10 b9
 *       c2 c1 = sum of b1,b3,b5
 *    c4 c3    = sum of b2,b4,b6
 * ------------
 *    c4 c3
 *       d2 d1 = sum of b7,b9,c1
 *    d4 d3    = sum of b8,b10,c2
 * ------------
 *          d1
 *    e2 e1    = sum of c3,d2,d3
 * e4 e3       = sum of c4,d4
 * ------------
 * f2 f1 e1 d1 = out[3..0]
 */
// The output Ctxts[0..3] must be initialized, can alias the inputs
static void fifteen4Four(const CtPtrs& out, const CtPtrs& in, long sizeLimit)
{
  // we need 6 scratch ciphertexts
  std::vector<Ctxt> tmp(8, *out[0]);

  // Aliases, referring to the scheme above.

  Ctxt& d1 = *out[0];
  Ctxt& e1 = *out[1];
  Ctxt& f1 = *out[2];
  Ctxt& f2 = *out[3];

  Ctxt& b1 = tmp[0];
  Ctxt& b2 = tmp[1];
  Ctxt& b3 = tmp[2];
  Ctxt& b4 = tmp[3];
  Ctxt& b5 = tmp[4];
  Ctxt& b6 = tmp[5];
  Ctxt& c1 = tmp[6];
  Ctxt& c2 = tmp[7];
  Ctxt& c3 = b1;
  Ctxt& c4 = b3;
  Ctxt& b7 = b5;
  Ctxt& b8 = b2;
  Ctxt& b9 = b4;
  Ctxt& b10 = b6;
  Ctxt& d2 = b7;
  Ctxt& d3 = b9;
  Ctxt& d4 = f2;
  Ctxt& e2 = c1;
  Ctxt& e3 = c2;
  Ctxt& e4 = f2;

  long nThreads = std::min(NTL::AvailableThreads(), 3L);
  NTL_EXEC_INDEX(nThreads, index) // run these three lines in parallel
  switch (index) {
  case 0:
    three4Two(&b1, &b2, in[0], in[1], in[2]); // b2 b1 = 3for2(in[0..2])
    if (nThreads > 1)
      break;
    // FALLTHROUGH
  case 1:
    three4Two(&b3, &b4, in[3], in[4], in[5]); // b4 b3 = 3for2(in[3..5])
    if (nThreads > 2)
      break;
    // FALLTHROUGH
  default:
    three4Two(&b5, &b6, in[6], in[7], in[8]); // b6 b5 = 3for2(in[6..8])
  }
  NTL_EXEC_INDEX_END

  three4Two(c1, c2, b1, b3, b5); // c2 c1 = 3for2(b1,b3,b5)

  three4Two(c3, c4, b2, b4, b6); // c4 c3 = 3for2(b2,b4,b6)

  nThreads = std::min(NTL::AvailableThreads(), 2L);
  NTL_EXEC_INDEX(nThreads, index) // run these two lines in parallel
  switch (index) {
  case 0:
    three4Two(&b7, &b8, in[9], in[10], in[11]); // b8 b7 = 3for2(in[9..11])
    if (nThreads > 1)
      break;
    // FALLTHROUGH
  default:
    three4Two(&b9, &b10, in[12], in[13], in[14]); // b10 b9 = 3for2(in[12..14])
  }
  NTL_EXEC_INDEX_END

  NTL_EXEC_INDEX(nThreads, index) // run these two lines in parallel
  switch (index) {
  case 0:
    three4Two(d1, d2, b7, b9, c1); // d2 d1 = 3for2(b7,b9,c1)
    if (nThreads > 1)
      break;
    // FALLTHROUGH
  default:
    if (sizeLimit >= 2)
      three4Two(d3, d4, b8, b10, c2); // d4 d3 = 3for2(b8,b10,c2)
  }
  NTL_EXEC_INDEX_END
  if (sizeLimit < 2)
    return;

  NTL_EXEC_INDEX(nThreads, index) // run these two blocks in parallel
  switch (index) {
  case 0:
    three4Two(e1, e2, c3, d2, d3); // e2 e1 = 3for2(c3,d2,d3)
    if (nThreads > 1)
      break;
    // FALLTHROUGH
  default:
    if (sizeLimit >= 3) {
      e3 = c4;
      e3 += d4;          // e3 = c4 ^ d4
      e4.multiplyBy(c4); // e4 = c4 * d4 (e4 alias d4)
    }
  }
  NTL_EXEC_INDEX_END
  if (sizeLimit < 3)
    return;

  f1 = e2;
  f1 += e3; // f1 = e2 ^ e3
  if (sizeLimit < 4)
    return;
  e2.multiplyBy(e3);
  f2 += e2; // f2 = e4^(e2*e3)  (f2 alias e4)
}

// Same as above, but some of the pointers may be null.
// Returns number of output bits that are not identically zero.
long fifteenOrLess4Four(const CtPtrs& out, const CtPtrs& in, long sizeLimit)
{
  HELIB_TIMER_START;
  long numNonNull = in.numNonNull();
  if (numNonNull > 7) {
    fifteen4Four(out, in, sizeLimit);
    return 4;
  }

  // At most 7 non-null pointers, collect them in the first entries of a vector
  long lastNonNull = -1;
  std::vector<Ctxt*> inPtrs(7, nullptr);
  for (long i = 0; i < 15; i++)
    if (in.isSet(i))
      inPtrs[++lastNonNull] = in[i];

  if (numNonNull > 3) {
    seven4Three(out, CtPtrs_vectorPt(inPtrs), sizeLimit);
    out[3]->clear(); // msb is zero
    return 3;
  }
  numNonNull = three4Two(out[0], out[1], inPtrs[0], inPtrs[1], inPtrs[2]);
  out[3]->clear(); // msb is zero
  out[2]->clear(); // 2nd msb is zero
  return numNonNull;
}

/********************************************************************/
/***************** test/debugging functions *************************/

// Decrypt the binary numbers that are encrypted in eNums. The bits
// are encrypted in a bit-sliced manner. Namely, encNums[0] contains
// the LSB of all the numbers, encNums[1] the next bits from all, etc.
// If allSlots==false then we only return the subcube with index i=0
// in the last dimension within each ciphertext. Namely, the bit for
// the j'th counter is found in slot of index j*sizeOf(lastDim).
void decryptBinaryNums(std::vector<long>& pNums,
                       const CtPtrs& eNums,
                       const SecKey& sKey,
                       const EncryptedArray& ea,
                       bool twosComplement,
                       bool allSlots)
{
  int offset = 1, size = ea.size();
  if (!allSlots) { // only slots of index i=0 in the last dimension
    offset = ea.sizeOfDimension(ea.dimension() - 1);
    size /= offset;
  }
  pNums.assign(size, 0); // initialize to zero

  for (int i = 0; i < lsize(eNums); i++)
    if (eNums.isSet(i)) {
      std::vector<long> slots;
      ea.decrypt(*eNums[i], sKey, slots);
      for (int j = 0; j < lsize(pNums); j++)
        if (twosComplement && i == lsize(eNums) - 1)
          pNums[j] -= (slots[j * offset] << i);
        else
          pNums[j] += (slots[j * offset] << i);
    }
}

/********************************************************************/
#ifdef HELIB_DEBUG

void AddDAG::printAddDAG(bool printCT)
{
  std::cout << "aSize=" << aSize << ", bSize=" << bSize << std::endl;
  std::cout << "The p[i,j]'s\n============\n";
  for (long delta = 0; delta < bSize; delta++) {
    std::cout << "delta=" << delta << std::endl;
    for (long j = 0; j < bSize - delta; j++) {
      long i = j + delta;
      DAGnode* node = findP(i, j);
      if (node == nullptr)
        continue;
      std::cout << node->nodeName() << ":{ lvl=";
      if (node->level == LONG_MAX)
        std::cout << "XX";
      else
        std::cout << node->level;
      std::cout << ", chLeft=" << int(node->childrenLeft)
                << ", ct=" << node->ct;
      if (node->parent2)
        std::cout << ", prnt2=" << node->parent2->nodeName();
      if (node->parent1)
        std::cout << ", prnt1=" << node->parent1->nodeName();
      std::cout << " }\n";
      if (printCT && node->ct != nullptr)
        decryptAndPrint(std::cout,
                        *(node->ct),
                        *dbgKey,
                        *dbgEa,
                        FLAG_PRINT_VEC);
    }
  }
  std::cout << "\nThe q[i,j]'s\n============\n";
  for (long delta = 0; delta < bSize; delta++) {
    std::cout << "delta=" << delta << std::endl;
    for (long j = 0; j < std::min(aSize, bSize - delta); j++) {
      long i = j + delta;
      DAGnode* node = findQ(i, j);
      if (node == nullptr)
        continue;
      std::cout << node->nodeName() << ":{ lvl=";
      if (node->level == LONG_MAX)
        std::cout << "XX";
      else
        std::cout << node->level;
      std::cout << ", chLeft=" << long(node->childrenLeft)
                << ", ct=" << node->ct;
      if (node->parent2)
        std::cout << ", prnt2=" << node->parent2->nodeName();
      if (node->parent1)
        std::cout << ", prnt1=" << node->parent1->nodeName();
      std::cout << " }\n";
      if (printCT && node->ct != nullptr)
        decryptAndPrint(std::cout,
                        *(node->ct),
                        *dbgKey,
                        *dbgEa,
                        FLAG_PRINT_VEC);
    }
  }
  std::cout << std::endl;
}

void decryptAndSum(std::ostream& s,
                   const CtPtrMat& numbers,
                   bool twosComplement)
{
  s << "sum(";
  long sum = 0;
  for (long i = 0; i < numbers.size(); i++) {
    std::vector<long> slots;
    const CtPtrs& num = numbers[i];
    decryptBinaryNums(slots, num, *dbgKey, *dbgEa, twosComplement);
    s << slots[0] << ' ';
    sum += slots[0];
  }
  s << ")=" << sum << std::endl;
}

#endif // ifdef HELIB_DEBUG

} // namespace helib
