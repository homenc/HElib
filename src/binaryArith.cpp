/* Copyright (C) 2012-2017 IBM Corp.
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
/* binaryArith.cpp
 * Implementing integer addition and multiplication in binary representation.
 */
#include <numeric>
#include <climits>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <atomic>
#include <mutex>          // std::mutex, std::unique_lock

#include <NTL/BasicThreadPool.h>
#include "binaryArith.h"

#ifdef DEBUG_PRINTOUT
#include "debugging.h"
#endif

typedef std::pair<long,long> NodeIdx; // nodes are indexed by a pair (i,j)

/**
 * @class DAGnode
 * @brief A node in an addition-DAG structure.
 **/
class DAGnode {
public:
  NodeIdx idx; // the indexes of this node
  bool isQ;    // if not then isP
  long level;   // The level at the time of computation

  std::atomic_long childrenLeft; // how many children were not computed yet
  DAGnode *parent1, *parent2;

  std::mutex ct_mtx; // controls access to ctxt pointer (and the ctxt itself)
  Ctxt* ct;          // points to the actual ciphertext (or NULL)

  DAGnode(NodeIdx ii, bool qq, long lvl, long chl=0,
           DAGnode* pt1=nullptr, DAGnode* pt2=nullptr):
    idx(ii), isQ(qq), level(lvl),
    childrenLeft(chl), parent1(pt1), parent2(pt2), ct(nullptr) {}

  DAGnode(DAGnode &&other): // move constructor
    idx(other.idx), isQ(other.isQ), level(other.level),
    childrenLeft(long(other.childrenLeft)),// copy value of atomic_long
    parent1(other.parent1), parent2(other.parent2), ct(other.ct) {}

  std::string nodeName() const
  { return (std::string(isQ? "Q(":"P(")
            +std::to_string(idx.first)+','+ std::to_string(idx.second) +')');}
};

class ScratchCell {
public:
  std::atomic_bool used;
  std::unique_ptr<Ctxt> ct; // scratch space owns this pointer
  ScratchCell(const Ctxt& c): used(true), ct(new Ctxt(ZeroCtxtLike,c)) {}
  ScratchCell(ScratchCell&& other): // move constructor
    used(bool(other.used)), ct(std::move(other.ct)) {}
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
 * exactly two parents, whose product yeilds that node.
 *
 * When building the DAG, we attempt to keep the level of each node as high as
 * possible. So for example q_{i,j} can be computed either as p_{i,k}*q_{k-1,j}
 * or p_{i,k+1}*q_{k,j} (among other options), and we choose the option that
 * results in the largest possible level.
 **/
class AddDAG {
  std::vector<ScratchCell> scratch; // scratch space for ciphertexts
  std::map<NodeIdx,DAGnode> p; // p[i,j]= prod_{t=j}^i (a[t]+b[t])
  std::map<NodeIdx,DAGnode> q; // q[i,j]= a[j]b[j]*prod_{t=j+1}^i (a[t]+b[t])
  long aSize, bSize;

  Ctxt* allocateCtxtLike(const Ctxt& c);
  void markAsAvailable(DAGnode* node);
  const Ctxt& getCtxt(DAGnode* node,
                      const CtPtrs& a, const CtPtrs& b);

  void addCtxtFromNode(Ctxt& c, DAGnode* node,
                       const CtPtrs& a, const CtPtrs& b) {
    std::unique_lock<std::mutex> lck(node->ct_mtx);
    c += getCtxt(node, a, b);
    if (--(node->childrenLeft) == 0) markAsAvailable(node);
  }

public:
  long getAsize() const { return aSize; }
  long getBsize() const { return bSize; }

  void init(const CtPtrs& a, const CtPtrs& b);
  void apply(CtPtrs& sum, const CtPtrs& a, const CtPtrs& b, long sizeLimit=0);

  // build the addition DAG
  AddDAG(const CtPtrs& a, const CtPtrs& b) { init(a,b); }

  long lowLvl() const
  {
    if (aSize<1) return 0;
    return findQ(bSize-1,0)->level;
  }
  DAGnode* findP(long i, long j) const { // returns NULL if not exists
    auto it = p.find(NodeIdx(i,j));
    if (it == p.end()) {
      cerr << "  findP("<<i<<','<<j<<") not found\n";
      return nullptr;  // not found
    }
    return (DAGnode*)&(it->second);
  }
  DAGnode* findQ(long i, long j) const { // returns NULL if not exists
    auto it = q.find(NodeIdx(i,j));
    if (it == q.end()) {
      cerr << "  findQ("<<i<<','<<j<<") not found\n";
      return nullptr;  // not found
    }
    return (DAGnode*)&(it->second);
  }
};


// When searching for a good middle point, we use a "good default", so that
// in the cases covered by that default solution we get the smallest number
// of nodes with childrenLeft>0.
// When initializing p[i,j] the default should be p[i,i+1-2^e]*p[i-2^e,j]
// where e is the largest exponent with 2^e <= i-j.
inline long defaultPmiddle(long delta)
{
  return 1 << (NTL::NumBits(delta)-1);
}
// When initializing q[i,j] the default should be p[i,i+1-2^e]*q[i-2^e,j]
// where e is the largest exponent with 2^e+2^{e-1} <= i-j.
inline long defaultQmiddle(long delta)
{
  delta = (2*delta+2)/3; // ceil(2 delta / 3)
  return 1 << (NTL::NumBits(delta)-1);
}

// Initialize the DAG
// FIXME: We often have at the bottom of a/b bits that are known to be zero,
//    the procedures below do not take any advantage of this.
void AddDAG::init(const CtPtrs& aa, const CtPtrs& bb)
{
  // make sure that lsize(b) >= lsize(a)
  const CtPtrs& a = (bb.size()>=aa.size())? aa : bb;
  const CtPtrs& b = (bb.size()>=aa.size())? bb : aa;

  aSize = lsize(a);
  bSize = lsize(b);
  assert (aSize>=1);

  // Initialize the p[i,i]'s and q[i,i]'s
  p.clear();
  q.clear();
  for (long i=0; i<bSize; i++) {
    NodeIdx idx(i,i);
    long lvl = (b.isSet(i) && !(b[i]->isEmpty()))?
      b[i]->findBaseLevel() : LONG_MAX;
    if (i<aSize) {
      long aLvl = (a.isSet(i) && !(a[i]->isEmpty()))?
        a[i]->findBaseLevel() : LONG_MAX;
      lvl = std::min(lvl, aLvl);
      q.emplace(idx,DAGnode(idx, true, lvl-1, 1));
    }
    p.emplace(idx,DAGnode(idx, false, lvl, 1));
  }

  // Initialize p[i,j] for bSize>=i>j>0
  for (long delta=1; delta<bSize; delta++) for (long i=bSize-1; i>=delta; --i) {
      long j = i-delta;
      long mid= i-defaultPmiddle(delta); // initialize to a "good default"
      DAGnode* prnt2 = findP(i,mid+1);
      DAGnode* prnt1 = findP(mid,j);
      long maxLvl = std::min(prnt1->level, prnt2->level) -1;
      if (prnt1->level==LONG_MAX || prnt2->level==LONG_MAX) // parent is empty
        maxLvl = LONG_MAX;
      long maxN =std::min(long(prnt1->childrenLeft), long(prnt2->childrenLeft));
      for (long m=j; m<i; m++) { // find middle point maximizing lvl(p[i,j])
        if (m==mid) continue;
        DAGnode* p2 = findP(i,m+1);
        DAGnode* p1 = findP(m,j);
        long lvl = std::min(p1->level, p2->level) -1;
        if (p1->level==LONG_MAX || p2->level==LONG_MAX) // parent is empty
          lvl = LONG_MAX;
        long n = std::min(long(p1->childrenLeft), long(p2->childrenLeft));
        if (lvl > maxLvl || (lvl == maxLvl && n > maxN)) {
          maxLvl = lvl;
          maxN = n;
          prnt1 = p1;
          prnt2 = p2;
        }
      }
      NodeIdx idx(i,j);
      p.emplace(idx, DAGnode(idx,false,maxLvl,0,prnt1,prnt2));
      prnt1->childrenLeft++;
      prnt2->childrenLeft++;
  }

  // Initialize q[i,j] for bSize>=i>j>=0
  for (long delta=1; delta<bSize; delta++) for (long i=bSize-1; i>=delta; --i) {
      long j = i-delta;
      if (j>=aSize) continue;
      long maxLvl=0, maxN=0;
      long mid= i-defaultQmiddle(delta); // initialize to a "good default"
      DAGnode* prnt2 = findP(i,mid+1);
      DAGnode* prnt1 = findQ(mid,j);
      if (prnt1!=nullptr) {
        maxLvl = std::min(prnt1->level, prnt2->level) -1;
        if (prnt1->level==LONG_MAX || prnt2->level==LONG_MAX)// parent is empty
          maxLvl = LONG_MAX;
        maxN = long(prnt2->childrenLeft);
      }
      for (long m=j; m<i; m++) { // find middle point maximizing lvl(p[i,j])
        if (m==mid) continue;
        DAGnode* p2 = findP(i,m+1);
        DAGnode* p1 = findQ(m,j);
        if (p1==nullptr) continue;
        long lvl = std::min(p1->level, p2->level) -1;
        if (p1->level==LONG_MAX || p2->level==LONG_MAX) // parent is empty
          lvl = LONG_MAX;
        long n = long(p2->childrenLeft);
        if (lvl > maxLvl || (lvl == maxLvl && n > maxN)) {
          maxLvl = lvl;
          maxN = n;
          prnt1 = p1;
          prnt2 = p2;
        }
      }
      if  (prnt1==nullptr) continue; // cannot create node
      NodeIdx idx(i,j);
      q.emplace(idx, DAGnode(idx,true,maxLvl,1,prnt1,prnt2));
      prnt1->childrenLeft++;
      prnt2->childrenLeft++;
    }
}

// Apply the DAG
void AddDAG::apply(CtPtrs& sum,
                   const CtPtrs& aa, const CtPtrs& bb, long sizeLimit)
{
  // make sure that lsize(b) >= lsize(a)
  const CtPtrs& a = (bb.size()>=aa.size())? aa : bb;
  const CtPtrs& b = (bb.size()>=aa.size())? bb : aa;
  if (aSize != a.size() || bSize != b.size())
    throw std::logic_error("DAG applied to wrong vectors");

  if (sizeLimit==0) sizeLimit = bSize+1;
  sum.resize(sizeLimit, &b); // allocate space for the output

  // Allow multi-threading in this loop
  NTL_EXEC_RANGE(sizeLimit, first, last)
  for (long i=first; i<last; i++) {
    if (i<bSize)
      addCtxtFromNode(*(sum[i]), this->findP(i,i), a, b);
    for (long j=std::min(i-1, aSize-1); j>=0; --j) {
      DAGnode *node = this->findQ(i-1,j);
      if (node!=nullptr) addCtxtFromNode(*(sum[i]), node, a, b);
    }
  }
  NTL_EXEC_RANGE_END
}

const Ctxt& AddDAG::getCtxt(DAGnode* node,
                            const CtPtrs& a, const CtPtrs& b)
{
  // NOTE: node->ct_mtx should be locked before calling this function

  if (node->ct == nullptr) { // ciphertext not computed yet, do it now
    if (node->parent1!=nullptr && node->parent2!=nullptr) { // internal node
      std::unique_lock<std::mutex> pt1_lck(node->parent1->ct_mtx);
      const Ctxt& c1 = getCtxt(node->parent1, a, b); // get the parents
      std::unique_lock<std::mutex> pt2_lck(node->parent2->ct_mtx);
      const Ctxt& c2 = getCtxt(node->parent2, a, b);
      long n1 = --(node->parent1->childrenLeft);
      long n2 = --(node->parent2->childrenLeft);
      if (n1==0) {                 // reuse space from parent1
        node->parent1->ct = nullptr;
        node->ct = (Ctxt*) &c1;
        if (c1.isEmpty() || c2.isEmpty())
             node->ct->clear(); // ct is zero if any of the parents is
	else node->ct->multiplyBy(c2);
        if (n2==0)
          markAsAvailable(node->parent2);
      }
      else if (n2==0) {            // reuse space from parent2
        node->parent2->ct = nullptr;
        node->ct = (Ctxt*) &c2;
        if (c1.isEmpty() || c2.isEmpty())
             node->ct->clear(); // ct is zero if any of the parents is
	else node->ct->multiplyBy(c1);
      }
      else {                       // allocate new space
        node->ct = allocateCtxtLike(c2);
        *(node->ct) = c2;
        if (c1.isEmpty() || c2.isEmpty())
             node->ct->clear(); // ct is zero if any of the parents is
	else node->ct->multiplyBy(c1);
      }
    }
    else { // no parents, either a[i]+b[i] or a[i]*b[i]
      long i = node->idx.first;
      long j = node->idx.second; // we expect i==j
      node->ct = allocateCtxtLike(*(b[0]));

      if (node->isQ) { // This is b[i]*a[j]
        if (b.isSet(i) && !(b[i]->isEmpty())
            && a.isSet(j) && !(a[j]->isEmpty())) {
          *(node->ct) = *(b[i]);
          node->ct->multiplyBy(*(a[j]));
        } // if a[j] or b[i] is empty then node->ct is a zero ciphertext
	else node->ct->clear();
      }
      else {           // This is b[i]+a[i]
        if (!b.isSet(i) || b[i]->isEmpty())
	     node->ct->clear();
        else *(node->ct) = *(b[i]);
        if (a.isSet(j) && !(a[j]->isEmpty()))
          *(node->ct) += *(a[j]);
      }
    } // end of no-parents case
  }
  return *(node->ct);
}


// Adds another cell to scratch space, or use an existing one that's free
Ctxt* AddDAG::allocateCtxtLike(const Ctxt& c)
{
  // look for an unused cell in the scratch array
  for (long i=0; i<(long)scratch.size(); i++)
    if (scratch[i].used == false) { // found a free one, try to use it
      bool used = scratch[i].used.exchange(true); // mark it as used
      if (used==false)     // make sure no other thread got there first
        return scratch[i].ct.get();
    }

  // If not found, allocate a new cell
  ScratchCell sc(c);      // cell points to new ctxt, with used=true
  Ctxt* pt = sc.ct.get(); // remember the raw pointer
  scratch.push_back(std::move(sc));  // scratch now owns the pointer
  return pt;         // return the raw pointer
}

// Mark a scratch ciphertext as unused. We assume that no two nodes
// ever share a ciphertext object, so if this node is done with the
// object then the object is unused.
void AddDAG::markAsAvailable(DAGnode* node)
{
  // NOTE: node->ct_mtx should be locked before calling this function
  // NOTE: somewhat inefficient, use linear search for the raw pointer
  for (long i=0; i<(long)scratch.size(); i++)
    if (scratch[i].ct.get()==node->ct)
      scratch[i].used = false;
  node->ct = nullptr;
}

/********************************************************************/
/********************************************************************/

// Use packed bootstrapping, so we can bootstrap all in just one go.
static void packedRecrypt(const CtPtrs& a, const CtPtrs& b,
                          std::vector<zzX>* unpackSlotEncoding)
{
  const Ctxt* ct = b.ptr2nonNull(); // find some non-null Ctxt
  if (ct==nullptr) ct = a.ptr2nonNull();
  if (ct==nullptr) return;    // nothing to do

  assert(unpackSlotEncoding!=nullptr && ct->getPubKey().isBootstrappable());

  struct CtPtrs_pair : CtPtrs {
    const CtPtrs& a;
    const CtPtrs& b;
    CtPtrs_pair(const CtPtrs &_a, const CtPtrs &_b): a(_a), b(_b) {}
    Ctxt* operator[](long i) const override
    { return (i<a.size())? a[i] : b[i-lsize(a)]; }
    long size() const override { return a.size() + b.size(); }
  };
  const CtPtrs_pair ab(a,b);

  packedRecrypt(ab, *unpackSlotEncoding, *(ct->getContext().ea));
}

void addTwoNumbers(CtPtrs& sum, const CtPtrs& a, const CtPtrs& b,
                   long sizeLimit, std::vector<zzX>* unpackSlotEncoding)
{
  if (a.size()<1)      { vecCopy(sum,b); return; }
  else if (b.size()<1) { vecCopy(sum,a); return; }

  // Work out the order of multiplications to compute all the carry bits
  AddDAG mults(a,b);

  // Ensure that we have enough levels to compute everything,
  // bootstrap otherwise
  if (mults.lowLvl()<1) {
    packedRecrypt(a,b,unpackSlotEncoding);
    mults.init(a,b); // Re-compute the DAG
    if (mults.lowLvl()<1) { // still not enough levels
      throw std::logic_error("not enough levels for multiplication DAG");
    }
  }

  //  printAddDAG(mults);
  mults.apply(sum, a, b, sizeLimit);    // perform the actual addition
}

// An object of type T must have an lsize() function
static std::tuple<const CtPtrs*,const CtPtrs*,const CtPtrs*>
orderBySize(const CtPtrs& a, const CtPtrs& b, const CtPtrs& c)
{
  if (a.size() <= b.size()) {
    if (b.size() <= c.size())
      return std::make_tuple(&a,&b,&c); // a <= b <= c
    else if (a.size() <= c.size())
      return std::make_tuple(&a,&c,&b); // a <= c < b
    else
      return std::make_tuple(&c,&a,&b); // c < a <= b
  } else { // lsize(b) < lsize(a)
    if (a.size() <= c.size())
      return std::make_tuple(&b,&a,&c); // b < a <= c
    else if (b.size() <= c.size())
      return std::make_tuple(&b,&c,&a); // b <= c < a
    else
      return std::make_tuple(&c,&b,&a); // c < b < a
  }
}


// Implementing the basic 3-for-2 trick: u,v,w encrypt bits, the
// procedure computes x,y such that x+2y = u+v+w over the integers.
// outputs can alias the inputs
static void three4Two(Ctxt& lsb, Ctxt& msb,
                      const Ctxt& u, const Ctxt& v, const Ctxt& w)
{
  Ctxt tmp_v = v;
  Ctxt tmp_w = w;
  lsb = u;
  msb = u;

  lsb += tmp_v;          // u+v
  msb.multiplyBy(tmp_v); // u*v

  tmp_v = lsb;             // u+v

  tmp_v.multiplyBy(tmp_w); // (u+v)*w

  lsb += tmp_w;            // u+v+w
  msb += tmp_v;   // u*v + (u+v)*w = u*v + u*w + v*w
}

// Same as three4Two above, but some of the inputs could be null.
// Returns the number of output bits that are not identically zero.
static long three4Two(Ctxt* lsb, Ctxt* msb, Ctxt* u, Ctxt* v, Ctxt* w)
{
  if (u!=nullptr && !u->isEmpty()
      && v!=nullptr && !v->isEmpty()
      && w!=nullptr && !w->isEmpty()) {
    three4Two(*lsb,*msb, *u, *v, *w);
    return 2;
  }
  if ((u==nullptr||u->isEmpty())
      && (v==nullptr||v->isEmpty())
      && (w==nullptr||w->isEmpty())) {
    lsb->clear();
    msb->clear();
    return 0;
  }

  // Some are null, others are not, arrange so that nulls are at the end
  if (u==nullptr||u->isEmpty()) {
    if (v==nullptr||v->isEmpty()) u=w;
     else { u=v; v=w; }
  }
  else if (v==nullptr||v->isEmpty()) v=w;
  w = nullptr;

  if (v==nullptr||v->isEmpty()) { // only u is non-null
    *lsb = *u;
    msb->clear();
    return 1;
  }

  // both u,v are non-null
  Ctxt tmp = *v;
  *lsb = *u;
  *msb = *u;
  *lsb += tmp;
  msb->multiplyBy(tmp);
  return 2;
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

  // Aliasas, referring to the scheme above. Aliases for temporary
  // vars chosen so that inputs, outputs of three4two are distinct

  Ctxt& c1=*out[0]; Ctxt& d1=*out[1]; Ctxt& d2=*out[2];

  Ctxt& b1=tmp[0];  Ctxt& b2=tmp[1];  Ctxt& b3=tmp[2];  Ctxt& b4=tmp[3];
  Ctxt& c2=d1;  Ctxt& c3=b1;  Ctxt& c4=b3;

  three4Two(&b1,&b2,in[0],in[1],in[2]); // b2 b1 = 3for2(in[0..2])
  three4Two(&b3,&b4,in[3],in[4],in[5]); // b4 b3 = 3for2(in[3..5])

  three4Two(&c1,&c2, in[6], &b1, &b3);  // c2 c1 = 3for2(in[6],b1,b3)
  if (sizeLimit < 2) return;
  c3 = b2;

  c3 += b4;                           // c3 = b2 ^ b4
  c4 = b2;

  c4.multiplyBy(b4);                  // c4 = b2 * b4
  d2 = c2;

  d1 += c3;                           // d1 = c2 ^ c3 (d1 alias c2)
  if (sizeLimit < 3) return;
  d2.multiplyBy(c3);

  d2 += c4;                           // d2 = c4 ^ (c2*c3)
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

  // Aliasas, referring to the scheme above.

  Ctxt& d1=*out[0]; Ctxt& e1=*out[1]; Ctxt& f1=*out[2]; Ctxt& f2=*out[3];

  Ctxt& b1=tmp[0];  Ctxt& b2=tmp[1];  Ctxt& b3=tmp[2];  Ctxt& b4=tmp[3];
  Ctxt& b5=tmp[4];  Ctxt& b6=tmp[5];  Ctxt& c1=tmp[6];  Ctxt& c2=tmp[7];
  Ctxt& c3=b1;  Ctxt& c4=b3;
  Ctxt& b7=b5;  Ctxt& b8=b2;  Ctxt& b9=b4;  Ctxt& b10=b6;
  Ctxt& d2=b7;  Ctxt& d3=b9;  Ctxt& d4=f2;
  Ctxt& e2=c1;  Ctxt& e3=c2;  Ctxt& e4=f2;

  three4Two(&b1,&b2,in[0],in[1],in[2]); // b2 b1 = 3for2(in[0..2])
  three4Two(&b3,&b4,in[3],in[4],in[5]); // b4 b3 = 3for2(in[3..5])
  three4Two(&b5,&b6,in[6],in[7],in[8]); // b6 b5 = 3for2(in[6..8])

  three4Two(c1,c2, b1, b3, b5);        // c2 c1 = 3for2(b1,b3,b5)

  three4Two(c3,c4, b2, b4, b6);          // c4 c3 = 3for2(b2,b4,b6)

  three4Two(&b7,&b8,in[9],in[10],in[11]);  // b8 b7 = 3for2(in[9..11])
  three4Two(&b9,&b10,in[12],in[13],in[14]);// b10 b9 = 3for2(in[12..14])

  three4Two(d1,d2, b7, b9, c1);        // d2 d1 = 3for2(b7,b9,c1)
  if (sizeLimit < 2) return;
  three4Two(d3,d4, b8, b10, c2);       // d4 d3 = 3for2(b8,b10,c2)

  three4Two(e1,e2, c3, d2, d3);        // e2 e1 = 3for2(c3,d2,d3)
  if (sizeLimit < 3) return;
  e3 = c4;
  e3 += d4;                            // e3 = c4 ^ d4
  e4.multiplyBy(c4);                   // e4 = c4 * d4 (e4 alias d4)

  f1 = e2;
  f1 += e3;                            // f1 = e2 ^ e3
  if (sizeLimit < 4) return;
  e2.multiplyBy(e3);
  f2 += e2;                            // f2 = e4^(e2*e3)  (f2 alias e4)
}
// Same as above, but some of the pointers may be null.
// Returns number of output bits that are not identically zero.
long fifteenOrLess4Four(const CtPtrs& out, const CtPtrs& in, long sizeLimit)
{
  long numNonNull = in.numNonNull();
  if (numNonNull>7) {
    fifteen4Four(out, in, sizeLimit);
    return 4;
  }

  // At most 7 non-null pointers, collect them in the first entires of a vector
  long lastNonNull = -1;
  std::vector<Ctxt*> inPtrs(7,nullptr);
  for (long i=0; i<15; i++)
    if (in.isSet(i)) inPtrs[++lastNonNull] = in[i];

  if (numNonNull>3) {
    seven4Three(out, CtPtrs_vectorPt(inPtrs), sizeLimit);
    out[3]->clear(); // msb is zero
    return 3;
  }
  numNonNull = three4Two(out[0],out[1], inPtrs[0],inPtrs[1],inPtrs[2]);
  out[3]->clear(); // msb is zero
  out[2]->clear(); // 2nd msb is zero
  return numNonNull;
}

// Apply the 3-for-2 routine to integers (i.e., an array of bits). The
// inputs need not be of the same size, and size of the output x is
// equal to the largest of them, and the size of the output y is one
// larger. This is safe even when the outputs alias some of the inputs
static void three4Two(CtPtrs& lsb, CtPtrs& msb,
                      const CtPtrs& u, const CtPtrs& v, const CtPtrs& w,
                      long sizeLimit)
{
  // Arrange u,v,w by size from smallest to largest
  const CtPtrs *p1, *p2, *p3;
  std::tie(p1,p2,p3) = orderBySize(u,v,w); // size(p3)>=size(p2)>=size(p1)

  if (p3->size() <= 0) { // empty input
    lsb.resize(0);
    msb.resize(0);
    return;
  }
  if (p1->size()<=0) { // two or less inputs
    std::vector<Ctxt> tmp;
    vecCopy(tmp,*p2); // just in case p2, msb share pointers
    vecCopy(msb,*p3);
    vecCopy(lsb,tmp);
    return;
  }
  if (sizeLimit==0) sizeLimit = p3->size()+1;

  // Allocate space in the output vectors

  const Ctxt* ctptr = p3->ptr2nonNull();
  std::vector<Ctxt> tmpMsb, tmpLsb;

  long lsbSize = std::min(sizeLimit, p3->size());
  long msbSize = lsbSize;
  if (p2->size()==p3->size() && lsbSize<sizeLimit)
    msbSize++;                  // possible carry out of last position

  resize(tmpLsb, lsbSize, Ctxt(ZeroCtxtLike,*ctptr));
  resize(tmpMsb, msbSize, Ctxt(ZeroCtxtLike,*ctptr));

  for (long i=0; i<std::min(lsize(*p1),lsbSize); i++)
    if (i<lsize(tmpMsb)-1)
      three4Two(&tmpLsb[i], &tmpMsb[i+1], (*p1)[i], (*p2)[i], (*p3)[i]);
    else {
      if (p1->isSet(i)) tmpLsb[i] =  *((*p1)[i]);
      if (p2->isSet(i)) tmpLsb[i] += *((*p2)[i]);
      if (p3->isSet(i)) tmpLsb[i] += *((*p3)[i]);
    }

  for (long i=p1->size(); i<std::min(lsize(*p2),lsbSize); i++) {
    if (p2->isSet(i)) tmpLsb[i] =  *((*p2)[i]);
    if (p3->isSet(i)) tmpLsb[i] += *((*p3)[i]);
    if (i<lsize(tmpMsb)-1 && p2->isSet(i) && p3->isSet(i)) {
      tmpMsb[i+1] = *((*p2)[i]);
      tmpMsb[i+1].multiplyBy(*((*p3)[i]));
    }
  }
  for (long i=p2->size(); i<std::min(lsize(*p3),lsbSize); i++)
    if (p3->isSet(i)) tmpLsb[i] = *((*p3)[i]);

  vecCopy(lsb, tmpLsb);
  vecCopy(msb, tmpMsb);
}


// Calculates the sum of many numbers using the 3-for-2 method
void addManyNumbers(CtPtrs& sum, CtPtrMat& numbers, long sizeLimit,
                    std::vector<zzX>* unpackSlotEncoding)
{
  FHE_TIMER_START;
  if (lsize(numbers)<1) return;
  if (lsize(numbers)==1) { vecCopy(sum, numbers[0]); return; }

  // if just 2 numbers to add then use normal binary addition
  // else enter loop below. We view numbers as a FIFO queue, each
  // time removing the first three entries at the head and adding
  // two new ones at the tail.
  long head=0, tail=0;
  for (long leftInQ=lsize(numbers); leftInQ>2; leftInQ--) {
    long h2 = (head+1) % lsize(numbers);
    long h3 = (head+2) % lsize(numbers);
    long t2 = (tail+1) % lsize(numbers);

    // If any of head,h1,h2 are too low level, then bootstrap everything
    if (findMinLevel({&(numbers[head]),&(numbers[h2]),&(numbers[h3])}) < 2) {
      assert(unpackSlotEncoding!=nullptr);

      const Ctxt* ct = numbers[head].ptr2nonNull(); // find some non-null Ctxt
      if (ct==nullptr) ct = numbers[h2].ptr2nonNull();
      if (ct==nullptr) ct = numbers[h3].ptr2nonNull();
      assert(ct!=nullptr && ct->getPubKey().isBootstrappable());

      const EncryptedArray& ea = *(ct->getContext().ea);

      packedRecrypt(numbers, *unpackSlotEncoding, ea, /*belowLvl=*/10);
    }

    three4Two(numbers[tail], numbers[t2], // three4Two can work in-place
              numbers[head], numbers[h2], numbers[h3], sizeLimit);

    head = (head+3) % lsize(numbers);    
    tail = (tail+2) % lsize(numbers);
  }
  // final addition
  long h2 = (head+1) % lsize(numbers);
  addTwoNumbers(sum, numbers[head], numbers[h2], sizeLimit, unpackSlotEncoding);
} // NOTE: It'd be a little challenging to parallelize this


// Multiply a positive a by a potentially negative b, we need to sign-extend b
static void multByNegative(CtPtrs& product, const CtPtrs& a, const CtPtrs& b,
                           long sizeLimit,std::vector<zzX>* unpackSlotEncoding)
{
  long resSize = lsize(a)+lsize(b);
  if (sizeLimit>0 && sizeLimit<resSize) resSize=sizeLimit;

  NTL::Vec< NTL::Vec<Ctxt> > numbers(INIT_SIZE, std::min(lsize(a),resSize));
  for (long i=0; i<lsize(numbers); i++) {
    numbers[i].SetLength(resSize, Ctxt(ZeroCtxtLike,*(a[0])));
    for (long j=i; j<resSize; j++)
      if (j<i+lsize(b)) {
        if (b.isSet(j-i) && !b[j-i]->isEmpty()
            && a.isSet(i) && !a[i]->isEmpty()) {
          numbers[i][j] = *(b[j-i]);
          numbers[i][j].multiplyBy(*(a[i]));   // multiply by the bit of a
        }
      }
      else numbers[i][j] = numbers[i][i+lsize(b)-1]; // sign extension
  }
#ifdef DEBUG_PRINTOUT
  vector<long> slots;
  cout << " multByNegative: ";
  decryptBinaryNums(slots, a, *dbgKey, *dbgEa, false);
  cout << slots[0] << " * ";
  decryptBinaryNums(slots, b, *dbgKey, *dbgEa, true);
  cout << slots[0] << "=sum(";
  long sum=0;
  for (NTL::Vec<Ctxt>& num : numbers) {
    decryptBinaryNums(slots, CtPtrs_VecCt(num), *dbgKey, *dbgEa, true);
    cout << slots[0] << ' ';
    sum += slots[0];
  }
  cout << ")="<<sum<<endl;
#endif
  CtPtrMat_VecCt nums(numbers);
  addManyNumbers(product, nums, resSize);
}

// Multiply two integers (i.e. an array of bits) a, b.
// Computes the pairwise products x_{i,j} = a_i * b_j
// then sums the prodcuts using the 3-for-2 method.
void multTwoNumbers(CtPtrs& product, const CtPtrs& a, const CtPtrs& b,
                    bool bNegative, long sizeLimit,
                    std::vector<zzX>* unpackSlotEncoding)
{
  long aSize = lsize(a);
  long bSize = lsize(b);
  long resSize = aSize+bSize-1;
  if (sizeLimit>0 && sizeLimit<resSize) resSize=sizeLimit;

  if (a.numNonNull()<1 || b.numNonNull()<1) {
    product.resize(0);
    return; // return 0
  }
  // Edge case, if a or b is 1 bit
  if (aSize==1) {
    if (a[0]->isEmpty()) {
      product.resize(0);
      return;
    }
    vecCopy(product,b,resSize);
    for (long i=0; i<resSize; i++)
      product[i]->multiplyBy(*(a[0]));
    return;
  }
  if (bNegative) { // somewhat different implementation for 2s complement
    multByNegative(product, a, b, sizeLimit, unpackSlotEncoding);
    return;
  }
  if (bSize==1) {
    if (b[0]->isEmpty()) {
      product.resize(0);
      return;
    }
    vecCopy(product,a,resSize);
    for (long i=0; i<resSize; i++)
      a[i]->multiplyBy(*(b[0]));
    return;
  }

  // We make sure aa is the larger of the two integers
  // to keep the number of additions to a minimum
  const CtPtrs& aa = (aSize>=bSize)? a : b;
  const CtPtrs& bb = (aSize>=bSize)? b : a;
  aSize = lsize(aa);
  bSize = lsize(bb);

  NTL::Vec< NTL::Vec<Ctxt> > numbers(INIT_SIZE, std::min(lsize(b),resSize));
  const Ctxt* ct_ptr = a.ptr2nonNull();
  for (long i=0; i<lsize(numbers); i++) {
    numbers[i].SetLength(std::min((i+aSize),resSize),
                         Ctxt( ZeroCtxtLike,*ct_ptr ) );
    for (long j=i; j<lsize(numbers[i]); j++) {
      if (a.isSet(j-i) && !(a[j-i]->isEmpty())
          &&  b.isSet(i) && !(b[i]->isEmpty()) ) {
        numbers[i][j] = *(a[j-i]);
        numbers[i][j].multiplyBy(*(b[i])); // multiply by the bit of b
      }
    }
  }
#ifdef DEBUG_PRINTOUT
  vector<long> slots;
  cout << " multTwoNumbers: ";
  decryptBinaryNums(slots, a, *dbgKey, *dbgEa, false);
  cout << slots[0] << " * ";
  decryptBinaryNums(slots, b, *dbgKey, *dbgEa, false);
  cout << slots[0] << "=sum(";
  long sum=0;
  for (NTL::Vec<Ctxt>& num : numbers) {
    decryptBinaryNums(slots, CtPtrs_VecCt(num), *dbgKey, *dbgEa, false);
    cout << slots[0] << ' ';
    sum += slots[0];
  }
  cout << ")="<<sum<<endl;
#endif
  CtPtrMat_VecCt nums(numbers);
  addManyNumbers(product, nums, sizeLimit, unpackSlotEncoding);
}

/********************************************************************/
/***************** test/debugging functions *************************/

// Decrypt the binary numbers that are encrypted in eNums. The bits
// are encrypted in a bit-sliced manner. Namely, encNums[0] contains
// the LSB of all the numbers, encNums[1] the next bits from all, etc.
// If allSlots==false then we only return the subcube with index i=0
// in the last dimension within each ciphertext. Namely, the bit for
// the j'th counter is found in slot of index j*sizeOf(lastDim).
void decryptBinaryNums(vector<long>& pNums, const CtPtrs& eNums,
                       const FHESecKey& sKey, const EncryptedArray& ea,
                       bool negative, bool allSlots)
{
  int offset=1, size=ea.size();
  if (!allSlots) { // only slots of index i=0 in the last dimension
    offset = ea.sizeOfDimension(ea.dimension()-1);
    size /= offset;
  }
  pNums.assign(size, 0); // initialize to zero

  for (int i=0; i<lsize(eNums); i++) if (eNums.isSet(i)) {
    vector<long> slots;
    ea.decrypt(*eNums[i], sKey, slots);
    for (int j=0; j<lsize(pNums); j++)
      if (negative && i==lsize(eNums)-1)
        pNums[j] -= (slots[j*offset] << i);
      else
        pNums[j] += (slots[j*offset] << i);
  }
}

/********************************************************************/
#ifdef DEBUG_PRINTOUT
#include <cstdio>

void printAddDAG(const AddDAG& md, bool printCT)
{
  long aSize = md.getAsize();
  long bSize = md.getBsize();
  cout << "aSize="<<aSize<<", bSize="<<bSize<<endl;
  cout << "The p[i,j]'s\n============\n";  
  for (long delta=0; delta<bSize; delta++) {
    cout << "delta="<<delta<<endl;
    for (long j=0; j<bSize-delta; j++) {
      long i = j+delta;
      DAGnode* node = md.findP(i,j);
      if (node==nullptr) continue;
      cout << node->nodeName()<<":{ lvl=";
      if (node->level==LONG_MAX) cout << "XX";
      else                       cout << node->level;
      cout <<", chLeft="<<int(node->childrenLeft)
           <<", ct="<<node->ct;
      if (node->parent2)
        cout << ", prnt2="<< node->parent2->nodeName();
      if (node->parent1)
        cout << ", prnt1="<< node->parent1->nodeName();
      cout << " }\n";
      if (printCT && node->ct!=nullptr)
        decryptAndPrint(cout, *(node->ct), *dbgKey, *dbgEa, FLAG_PRINT_VEC);
    }
  }
  cout << "\nThe q[i,j]'s\n============\n";  
  for (long delta=0; delta<bSize; delta++) {
    cout << "delta="<<delta<<endl;
    for (long j=0; j<std::min(aSize, bSize-delta); j++) {
      long i = j+delta;
      DAGnode* node = md.findQ(i,j);
      if (node==nullptr) continue;
      cout << node->nodeName() <<":{ lvl=";
      if (node->level==LONG_MAX) cout << "XX";
      else                       cout << node->level;
      cout <<", chLeft="<<long(node->childrenLeft)
           <<", ct="<<node->ct;
      if (node->parent2)
        cout << ", prnt2="<< node->parent2->nodeName();
      if (node->parent1)
        cout << ", prnt1="<< node->parent1->nodeName();
      cout << " }\n";
      if (printCT && node->ct!=nullptr)
        decryptAndPrint(cout, *(node->ct), *dbgKey, *dbgEa, FLAG_PRINT_VEC);
    }
  }
  cout << endl;
}
#endif // ifdef DEBUG_PRINTOUT
