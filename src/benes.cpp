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
/* benes.cpp - Implementation of Benes permutation networks
 */
#include "FHE.h"
#include "EncryptedArray.h"
#include <NTL/lzz_pXFactoring.h>

#include <cstdlib>
#include <cassert>
#include <list>
#include <tr1/memory>
#include <sstream>


#include <NTL/vector.h>
#include "NumbTh.h"
#include "permutations.h"


/**
 * @class BenesNetwork
 * @brief Implementation of Benes Permutation Network
 **/
class BenesNetwork {
private:

  long k; // Benes network of size 2^k
  Vec< Vec<short> > level;
    // A dimension-k Benes network has 2*k - 1 levels, each with 2^k nodes.
    // level[i][j] \in {0,1}; 0 means a straight edge, 1 means a cross edge
    //   (at node j of level i); 
    // a straight points to the same node at the next level, cross edge from
    //   node j points to node j + delta at the next level, where delta = 2^p
    //   if bit p of j is 0, and delta = -2^p if bit p of j is 1;
    //   here, p is the bit position associated with level i, which is |k-1-i|

  BenesNetwork(); // default constructor disabled

public:

  long getDepth() const { return k; }
  long getSize() const { return 1L << k; }
  long getNumLevels() const { return 2*k-1; }

  const Vec<short>& getLevel(long i) const { return level[i]; }

  //! bit position associated with level i
  long bitPos(long i) const { return labs(k - 1 - i); }

  //! Constructs a Benes network that implements the given permutation
  BenesNetwork(long _k, const Permut& perm); // constructor

  //! Test if *this implements teh permutations perm
  bool testNetwork(const Permut& perm) const;
};

// The recursive function for doing the actual initialization
static void 
recursiveBenesInit(long k, long delta_i, long delta_j,
		   const Permut& perm, const Permut& iperm,
		   Vec< Vec<short> >& level, Vec< Vec<short> >& ilevel);

// Constructs a Benes network that implements the given permutation
BenesNetwork::BenesNetwork(long _k, const Permut& perm)
{
  k = _k;
  assert(k > 0);
  assert(k < NTL_BITS_PER_LONG-1);

  // Check that perm is of size 2^k
  long sz = 1L << k;
  assert(perm.length() == sz);

  // construct the inverse perm, which is convenient in the recursive function.
  // As a byproduct, verify that perm is indeed a permutation on {0,...,2^k-1}.

  Permut iperm;
  iperm.SetLength(sz);

  for (long j = 0; j < sz; j++)
    iperm[j] = -1;

  for (long j = 0; j < sz; j++) {
    long j1 = perm[j];
    assert(j1 >= 0 && j1 < sz);
    iperm[j1] = j;
  }


  for (long j = 0; j < sz; j++)
    assert(iperm[j] != -1);
  

  // allocate space for the levels graph

  level.SetLength(2*k-1);
  for (long i = 0; i < 2*k-1; i++)
    level[i].SetLength(sz);

  // allocate space for the reverse levels graph...
  // makes the recursive construction more convenient

  Vec< Vec<short> > ilevel;
  ilevel.SetLength(2*k-1);
  for (long i = 0; i < 2*k-1; i++)
    ilevel[i].SetLength(sz);

  // recursively construct the levels graph

  recursiveBenesInit(k, 0, 0, perm, iperm, level, ilevel);

}

// The recursive function for doing the actual initialization
void recursiveBenesInit(long k, long delta_i, long delta_j,
			const Permut& perm, const Permut& iperm,
			Vec< Vec<short> >& level, Vec< Vec<short> >& ilevel)
{
  if (k == 1) {
    // recursion stops here...
    // only two possibilities for perm: the identity perm or the swap perm

    if (perm[0] == 0) {
      // the identity perm

      level[delta_i][delta_j] = 0;
      level[delta_i][delta_j+1] = 0;
      ilevel[delta_i][delta_j] = 0;
      ilevel[delta_i][delta_j+1] = 0;
    }
    else {
      // the swap perm

      level[delta_i][delta_j] = 1;
      level[delta_i][delta_j+1] = 1;
      ilevel[delta_i][delta_j] = 1;
      ilevel[delta_i][delta_j+1] = 1;
    }
    return;
  }

  long sz = 1L << k;
  long hsz = sz/2;
  long nlev = level.length();

  // id_perm: the identity permutation on {0,...,sz-1}
  Permut id_perm;
  id_perm.SetLength(sz);
  for (long j = 0; j < sz; j++) id_perm[j] = j;
  
  // *xperm[0] is the perm on LHS of network
  // *xperm[1] is the perm on RHS of network
  // *xiperm[0] is the inv perm on LHS of network
  // *xiperm[1] is the inv perm on RHS
  const Permut *xperm[] = { &id_perm, &perm };
  const Permut *xiperm[] = { &id_perm, &iperm };

  // *first_level[0] is the first level when traversing left to right
  // *first_level[1] is the first level when traversing right to left
  // *ifirst_level[0] is the reversed first level when traversing left to right
  // *ifirst_level[1] is the reversed first level when traversing right to left
  Vec<short> *first_level[] = { &level[delta_i], &ilevel[nlev-1-delta_i] };
  Vec<short> *ifirst_level[] = { &ilevel[delta_i], &level[nlev-1-delta_i] };

  // *last_level[0] is the last level when traversing left to right
  // *last_level[1] is the last level when traversing right to left
  // *ilast_level[0] is the reversed last level when traversing left to right
  // *ilast_level[1] is the reversed last level when traversing right to left
  Vec<short> *last_level[] = { &level[nlev-1-delta_i], &ilevel[delta_i] };
  Vec<short> *ilast_level[] = { &ilevel[nlev-1-delta_i], &level[delta_i] };

  // inner_perm[0][0] upper internal perm
  // inner_perm[0][1] upper internal inv perm
  // inner_perm[1][0] lower internal perm
  // inner_perm[1][1] lower internal inv perm
  Permut inner_perm[2][2];

  inner_perm[0][0].SetLength(hsz);
  inner_perm[0][1].SetLength(hsz);
  inner_perm[1][0].SetLength(hsz);
  inner_perm[1][1].SetLength(hsz);

  // marked[0] indicates which nodes on left have been marked
  // marked[1] indicates which nodes on right have been marked
  Vec<bool> marked[2];
  marked[0].SetLength(sz);
  for (long j = 0; j < sz; j++) marked[0][j] = false;
  marked[1].SetLength(sz);
  for (long j = 0; j < sz; j++) marked[1][j] = false;

  // counter[0] is used to count through nodes on left side
  // counter[1] is used to count through nodes on right side
  long counter[2];
  counter[0] = 0;
  counter[1] = 0;


  int d = 0; // initial direction left to right
  long n = -1; // current node, initially undefined
  long e = 0;  // current edge

  for (;;) {
    if (n == -1) { // scan for unmarked node

      while (counter[d] < sz && marked[d][counter[d]]) counter[d]++;
      if (counter[d] >= sz) break; // we're done!!

      n = counter[d];
      e = 0;
    }

    marked[d][n] = true; // mark node n

    // traverse path through graph

    long v = (*xperm[d])[n]; // value associated with this node
    long net = ((n + e*hsz) % sz) / hsz; 
       // net = 0 => upper internal network
       // net = 1 => lower internal network

    long n1 = n % hsz; 
      // node adjacent to n in the internal network

    long n3 = (*xiperm[1-d])[v]; 
      // corresponding node of the other side of the external network

    long n2 = n3 % hsz;
      // node adjacent to n3 in the internal network

    long e2 = ((n3 + net*hsz) % sz) / hsz;
      // edge type for edge connecting n2 to n3

    /* here is the picture:
     *
     *       e               e2
     *  n ------ n1 ---- n2 ----- n3
     *          \           /
     *           \_________/
     *                |
     *         internal network
     */
#if 0  // print debug info
    cout << "d=" << d << " " << "net=" << net << " " << "n=" << n << " ";
    cout << "e=" << e  << " " << "n1=" << n1 << " " << "n2=" << n2 << " ";
    cout << "e2=" << e2 << " " << "n3=" << n3 << endl;
#endif

    // update external edges in level graph

    (*first_level[d])[delta_j+n] = e;
    (*ifirst_level[d])[delta_j+net*hsz+n1] = e;

    (*last_level[d])[delta_j+net*hsz+n2] = e2;
    (*ilast_level[d])[delta_j+n3] = e2;

    // update internal permutations

    inner_perm[net][d][n2] = n1;
    inner_perm[net][1-d][n1] = n2;

    // mark node n3

    marked[1-d][n3] = true;

    // calculate new node and edge

    long n3_sib = (n3 + hsz) % sz; // n3's sibling

    if (!marked[1-d][n3_sib]) {
      // n3's sib is ummarked, so traverse back using
      // starting from that node, using the same edge type as n2 <--> n3

      n = n3_sib;
      e = e2;
    }
    else {
       // start afresh

       n = -1;
    }

    d = 1-d; // reverse direction
  }

  // done constructing the external edges and internal permutations...
  // time to recurse

  // recurse on upper internal network:
  recursiveBenesInit(k-1, delta_i+1, delta_j, inner_perm[0][0], inner_perm[0][1],
                     level, ilevel);

  // recurse on lower intrernal network
  recursiveBenesInit(k-1, delta_i+1, delta_j + hsz, inner_perm[1][0], inner_perm[1][1],
                     level, ilevel);  
}

bool BenesNetwork::testNetwork(const Permut& perm) const
{
  long sz = getSize();
  long nlev = getNumLevels();

  for (long j = 0; j < sz; j++) {
    // find correct position for j

    long j1 = j;
    for (long i = 0; i < nlev; i++) {
      const Vec<short>& lev = getLevel(i);
      if (lev[j1] == 0) continue;

      long p = bitPos(i);
      long delta = (j1 & (1L << p)) ? -(1L << p) : (1L << p);
      j1 += delta;
    }
    if (perm[j1] != j) return false;
  }
  return true;
}


// aux routines for GeneralBenesNetwork

static 
long GB_depth(long n)
// computes recursion depth k for generalized Benes network of size n.
// the actual number of levels in the network is 2*k-1
{
  long k = 1;
  while ((1L << k) < n) k++;
  return k;
}


static inline
long GB_levelToDepthMap(long n, long k, long i) 
// maps a level number i = 0..2*k-2 to a recursion depth d = 0..k-1
// using the formula d = (k-1)-|(k-1)-i|
{
  assert(i >= 0 && i < 2*k-1);
  return (k-1) - labs((k-1)-i);
}


static inline
long GB_shamt(long n, long k, long i) 
// shift amount for level number i=0..2*k-2
// using the formula ceil( floor(n/2^d) / 2), 
//   where d = levelToDepthMap(i)
{
  long d = GB_levelToDepthMap(n, k, i);
  return ((n >> d) + 1) >> 1;
}

class GeneralBenesNetwork {
private:

  long n; // size of perm, n > 1, not necessarily power of 2

  long k; // recursion depth, k = least integer k s/t 2^k >= n

  Vec< Vec<short> > level;
    // there are 2*k - 1 levels, each wity n nodes.
    // level[i][j] is 0, 1, or -1, 
    //   which designates an edge from node j at level i 
    //   to node j + level[i][j]*shamt(i) at level i+1

  GeneralBenesNetwork(); // default constructor disabled

public:

  long getDepth() const { return k; }
  long getSize() const { return n; }
  long getNumLevels() const { return 2*k-1; }

  const Vec<short>& getLevel(long i) const 
  { 
    assert(i >= 0 && i < 2*k-1);
    return level[i];
  }

  long levelToDepthMap(long i) const { return GB_levelToDepthMap(n, k, i); }

  long shamt(long i) const { return GB_shamt(n, k, i); }
   
  // constructor
  GeneralBenesNetwork(const Permut& perm);

  // test correctness

  bool testNetwork(const Permut& perm) const;

};


static void 
recursiveGeneralBenesInit(long n, long k, long d, long delta_j,
                          const Permut& perm, const Permut& iperm,
                          Vec< Vec<short> >& level, Vec< Vec<short> >& ilevel)
{
  long sz = perm.length();

  if (d == k-1) {
    // recursion stops

    if (sz == 1) {
      level[d][delta_j] = 0;
      ilevel[d][delta_j] = 0;
    }
  
    if (sz == 2) {
      // only two possibilities for perm: the identity perm or the swap perm
  
      if (perm[0] == 0) {
        // the identity perm
  
        level[d][delta_j] = 0;
        level[d][delta_j+1] = 0;
        ilevel[d][delta_j] = 0;
        ilevel[d][delta_j+1] = 0;
      }
      else {
        // the swap perm
  
        level[d][delta_j] = 1;
        level[d][delta_j+1] = -1;
        ilevel[d][delta_j] = 1;
        ilevel[d][delta_j+1] = -1;
      }
    }
    return;
  }

  long nlev = 2*k-1;

  long sz0 = GB_shamt(n, k, d); // size of upper internal network
  long sz1 = sz - sz0;

#if 0
  cout << "**** " << d << " " << delta_j << " " << sz0 << " " << sz1 << "\n";
  cout << perm << "\n";
#endif


  // id_perm: the identity permutation on {0,...,sz-1}
  Permut id_perm;
  id_perm.SetLength(sz);
  for (long j = 0; j < sz; j++) id_perm[j] = j;
  
  // *xperm[0] is the perm on LHS of network
  // *xperm[1] is the perm on RHS of network
  // *xiperm[0] is the inv perm on LHS of network
  // *xiperm[1] is the inv perm on RHS
  const Permut *xperm[] = { &id_perm, &perm };
  const Permut *xiperm[] = { &id_perm, &iperm };

  // *first_level[0] is the first level when traversing left to right
  // *first_level[1] is the first level when traversing right to left
  // *ifirst_level[0] is the reversed first level when traversing left to right
  // *ifirst_level[1] is the reversed first level when traversing right to left
  Vec<short> *first_level[] = { &level[d], &ilevel[nlev-1-d] };
  Vec<short> *ifirst_level[] = { &ilevel[d], &level[nlev-1-d] };

  // *last_level[0] is the last level when traversing left to right
  // *last_level[1] is the last level when traversing right to left
  // *ilast_level[0] is the reversed last level when traversing left to right
  // *ilast_level[1] is the reversed last level when traversing right to left
  Vec<short> *last_level[] = { &level[nlev-1-d], &ilevel[d] };
  Vec<short> *ilast_level[] = { &ilevel[nlev-1-d], &level[d] };

  // inner_perm[0][0] upper internal perm
  // inner_perm[0][1] upper internal inv perm
  // inner_perm[1][0] lower internal perm
  // inner_perm[1][1] lower internal inv perm
  Permut inner_perm[2][2];

  inner_perm[0][0].SetLength(sz0);
  inner_perm[0][1].SetLength(sz0);
  inner_perm[1][0].SetLength(sz1);
  inner_perm[1][1].SetLength(sz1);

  // marked[0] indicates which nodes on left have been marked
  // marked[1] indicates which nodes on right have been marked
  Vec<bool> marked[2];
  marked[0].SetLength(sz);
  for (long j = 0; j < sz; j++) marked[0][j] = false;
  marked[1].SetLength(sz);
  for (long j = 0; j < sz; j++) marked[1][j] = false;

  // counter[0] is used to count through nodes on left side
  // counter[1] is used to count through nodes on right side
  long counter[2];
  counter[0] = 0;
  counter[1] = 0;



  long dir = 0; // direction, initially left to right
  long e = 0;  // current edge, initially straight

  long node; // current node
  long stop_node; // stopping node, initially undefined
  long stop_side;  // stopping side
  

  if (sz0 == sz1) {
    // an even split
    node = 0;
    stop_node = sz0;
    stop_side = 0;
  }
  else if (sz0 > sz1) {
    // an odd split, top larger
    node = sz0-1;
    stop_node = sz0-1;
    stop_side = 1;
  }
  else { // sz0 < sz1
    // an odd split, bottom larger
    node = sz-1;
    stop_node = sz-1;
    stop_side = 1;
  }

 

  for (;;) {
    marked[dir][node] = true; // mark current node

    // traverse path through graph

    long v = (*xperm[dir])[node]; // value associated with this node

    long net = long(node >= sz0)^labs(e);
       // net = 0 => upper internal network
       // net = 1 => lower internal network

    long node1 = node - sz0*long(node >= sz0); 
      // node adjacent to node in the internal network

    long node3 = (*xiperm[1-dir])[v]; 
      // corresponding node of the other side of the external network

    long node2 = node3 - sz0*long(node3 >= sz0);
      // node adjacent to node3 in the internal network

    long e2 = (node3 - (node2 + net*sz0))/sz0;
      // edge type for edge connecting node2 to node3

    /* here is the picture:
     *
     *         e                       e2
     *  node ------ node1 ---- node2 ----- node3
     *               \           /
     *                \_________/
     *                     |
     *              internal network
     */
#if 0  // print debug info
    cout << "dir=" << dir << " " << "net=" << net << " " << "node=" << node << " ";
    cout << "e=" << e  << " " << "node1=" << node1 << " " << "node2=" << node2 << " ";
    cout << "e2=" << e2 << " " << "node3=" << node3 << endl;
#endif

    // update external edges in level graph

    (*first_level[dir])[delta_j+node] = e;
    (*ifirst_level[dir])[delta_j+net*sz0+node1] = -e;

    (*last_level[dir])[delta_j+net*sz0+node2] = e2;
    (*ilast_level[dir])[delta_j+node3] = -e2;

    // update internal permutations

    inner_perm[net][dir][node2] = node1;
    inner_perm[net][1-dir][node1] = node2;

    // change direction
    dir = 1-dir;

    // mark node3 on new side
    marked[dir][node3] = true;

    // check for stop node
    if (node3 == stop_node && dir == stop_side) {
       // search for unmarked node 

      while (counter[dir] < sz && marked[dir][counter[dir]]) counter[dir]++;
      if (counter[dir] >= sz) break; // we're done!!


      // update node, e, stop_node, stop_side
      node = counter[dir];
      e = 0;

      if (node < sz0)
        stop_node = node + sz0;
      else
        stop_node = node - sz0; 

      stop_side = dir;
    }
    else {
       // calculate sibling node and continue
      if (node3 < sz0) 
        node = node3 + sz0; 
      else 
        node = node3 - sz0;

      e = e2;
    }
  }

  // done constructing the external edges and internal permutations...
  // time to recurse

  // recurse on upper internal network:
  recursiveGeneralBenesInit(n, k, d+1, delta_j, 
                            inner_perm[0][0], inner_perm[0][1], level, ilevel);

  // recurse on lower intrernal network
  recursiveGeneralBenesInit(n, k, d+1, delta_j + sz0, 
                            inner_perm[1][0], inner_perm[1][1], level, ilevel);  
}


GeneralBenesNetwork::GeneralBenesNetwork(const Permut& perm)
{
  n = perm.length();

  // check that n > 1
  assert(n > 1);

  // compute recursion depth k = least integer k s/t 2^k >= n
  k = GB_depth(n);

  // construct the inverse perm, which is convenient in the recursive function.
  // As a byproduct, verify that perm is indeed a permutation on {0,...,n-1}.

  Permut iperm;
  iperm.SetLength(n);

  for (long j = 0; j < n; j++)
    iperm[j] = -1;

  for (long j = 0; j < n; j++) {
    long j1 = perm[j];
    assert(j1 >= 0 && j1 < n);
    iperm[j1] = j;
  }


  for (long j = 0; j < n; j++)
    assert(iperm[j] != -1);

  // allocate space for the levels graph

  level.SetLength(2*k-1);
  for (long i = 0; i < 2*k-1; i++)
    level[i].SetLength(n);

  // allocate space for the reverse levels graph...
  // makes the recursive construction more convenient

  Vec< Vec<short> > ilevel;
  ilevel.SetLength(2*k-1);
  for (long i = 0; i < 2*k-1; i++)
    ilevel[i].SetLength(n);

  // recursively construct the levels graph

  recursiveGeneralBenesInit(n, k, 0, 0, perm, iperm, level, ilevel);
}


bool GeneralBenesNetwork::testNetwork(const Permut& perm) const
{
  long sz = getSize();
  long nlev = getNumLevels();

  for (long j = 0; j < sz; j++) {
    // find correct position for j

    long j1 = j;
    for (long i = 0; i < nlev; i++) {
      const Vec<short>& lev = getLevel(i);
      j1 += shamt(i)*lev[j1];
    }
    if (perm[j1] != j) return false;
  }
  return true;
}

// ******* optimization code **********

// helper class to make defining hash functions easier

template<class T> 
class ClassHash {
public:
  size_t operator()(const T& t) const { return t.hash(); }
};

// routines for finding optimal level collapsing strategies

void removeDups(list<long>& x, bool *aux)
// removes duplicates from list x
// uses vector aux as scratch space
//   invariants: 
//     * all elements in x lie in a range a..b
//     * aux[a..b] is false before and after removeDups is called
{
  for (list<long>::iterator i = x.begin(); i != x.end(); ) { 
    if (aux[*i])
      i = x.erase(i);
    else {
      aux[*i] = true;
      i++;
    }
  }

  for (list<long>::iterator i = x.begin(); i != x.end(); i++)
    aux[*i] = false;
}

void addOffset(list<long>& x, long offset, long n, bool *aux)
// creates a new list with all the old values, plus
//   all the old values with +/- offset added.
// results outside the range -n+1 .. n-1 are discarded
//   and all resulting duplicates are removed
{
  for (list<long>::iterator i = x.begin(); i != x.end(); i++) {
    long val = *i;
    long val1 = val + offset;
    long val2 = val - offset;
    if (val1 > -n && val1 < n) x.push_front(val1);
    if (val2 > -n && val2 < n) x.push_front(val2);
  }

  removeDups(x, aux);
}

long reducedCount(const list<long>& x, long n, bool *aux)
// counts the number of unique elements mod n in x 
{
  long res = 0;

  for (list<long>::const_iterator i = x.begin(); i != x.end(); i++) {
    long val = *i;
    if (val < 0) val += n;

    if (!aux[val]) {
      res++;
      aux[val] = true;
    }
  }

  for (list<long>::const_iterator i = x.begin(); i != x.end(); i++) {
    long val = *i;
    if (val < 0) val += n;
    aux[val] = false;
  }

  return res;
}


void buildBenesCostTable(long n, long k, bool good, Vec< Vec<long> >& tab)
{
  long nlev = 2*k-1;
  tab.SetLength(nlev);
  for (long i = 0; i < nlev; i++) tab[i].SetLength(nlev-i);

  // for j in [0..nlev-i) tab[i][j] is the cost of collapsing 
  // levels i..i+j


  Vec<bool> aux_vec;
  aux_vec.SetLength(2*n-1);
  bool *aux = &aux_vec[n-1];
  for (long i = 0; i < 2*n-1; i++) aux_vec[i] = false;

  for (long i = 0; i < nlev; i++) {
    list<long> x;

    x.push_front(0L);
    for (long j = 0; j < nlev-i; j++) {
      long shamt = GB_shamt(n, k, i+j);
      addOffset(x, shamt, n, aux);
      if (good)
        tab[i][j] = reducedCount(x, n, aux) - 1;
      else
        tab[i][j] = x.size() - 1;
    }
  }

}



class LongNode;
typedef tr1::shared_ptr<LongNode> LongNodePtr;


class LongNode {
public:
  long count;           // number of levels collapsed
  LongNodePtr next; // next node in the list

  LongNode(long _count, LongNodePtr _next) {
    count = _count; next = _next;
  }
};

long length(LongNodePtr ptr)
{
  long res = 0;
  for (LongNodePtr p = ptr; p != NULL; p = p->next) res++;
  return res;
}


void listToVec(LongNodePtr ptr, Vec<long>& vec)
{
  long len = length(ptr);

  vec.SetLength(len);
  long i = 0;
  for (LongNodePtr p = ptr; p != NULL; p = p->next) {
    vec[i] = p->count; 
    i++;
  }
}


ostream& operator<<(ostream& s, LongNodePtr p)
{
  if (p == NULL) {
    s << "[]";
    return s;
  }

  s << "[" << p->count;
  p = p->next;
  while (p != NULL) {
    s << " " << p->count;
    p = p->next;
  }
  cout << "]";

  return s;
}



class BenesMemoKey {
public:
  long i;
  long budget;

  BenesMemoKey(long _i, long _budget) {
    i = _i; budget = _budget;
  }

  size_t hash() const {
    stringstream s;
    s << i << " " << budget;
    return tr1::hash<string>()(s.str());
  }

  bool operator==(const BenesMemoKey& other) const {
    return i == other.i && budget == other.budget;
  }
};

class BenesMemoEntry {
public:
  long cost;
  LongNodePtr solution;

  BenesMemoEntry(long _cost, LongNodePtr _solution) {
    cost = _cost; solution = _solution;
  }

  BenesMemoEntry() { }
};



typedef tr1::unordered_map< BenesMemoKey, BenesMemoEntry, 
                            ClassHash<BenesMemoKey> > BenesMemoTable;

BenesMemoEntry optimalBenesAux(long i, long budget, long nlev, 
                                     const Vec< Vec<long> >& costTab, 
                                     BenesMemoTable& memoTab)
{
  assert(i >= 0 && i <= nlev);
  assert(budget > 0);

#if 0
  cout << "enter(" << i << "," << budget << ")\n";
#endif

  BenesMemoTable::iterator find = memoTab.find(BenesMemoKey(i, budget));

  if (find != memoTab.end()) {
#if 0
    cout << "found(" << i << "," << budget << ")\n";
#endif
    return find->second;
  }

  // a new subproblem to process...

  long cost;
  LongNodePtr solution;

  if (i == nlev) {
    // nothing to collapse, trivial solution

    cost = 0;
    solution = LongNodePtr();
    
  }
  else if (budget == 1) {
    // only one possible solution, which is to collapse
    // all remaining levels into a single level

    cost = costTab[i][nlev-i-1];
    solution = LongNodePtr(new LongNode(nlev-i, LongNodePtr()));
  }
  else {
    // budget > 1, so we consider collapsing levels i..i+j, for j in [0..nlev-i)

    long bestCost = NTL_MAX_LONG;
    long bestJ = 0;
    BenesMemoEntry bestT;
    for (long j = 0; j < nlev-i; j++) {
      BenesMemoEntry t = optimalBenesAux(i+j+1, budget-1, nlev, costTab, memoTab);
      if (t.cost + costTab[i][j] < bestCost) {
        bestCost = t.cost + costTab[i][j];
        bestJ = j;
        bestT = t;
      }
    }

    cost = bestCost;
    solution = LongNodePtr(new LongNode(bestJ+1, bestT.solution));
  }

#if 0
  cout << "computed(" << i << "," << budget << ") = (" << cost << ",";
  print(cout, solution);
  cout << ")\n";
#endif
  return memoTab[BenesMemoKey(i, budget)] = BenesMemoEntry(cost, solution);
}

void optimalBenes(long n, long budget, bool good, 
                     long& cost, LongNodePtr& solution)
// computes an optimal level-collapsing strategy for a Benes network
//   n = the size of the network
//   budget = an upper bound on the number of levels in the collapsed network
//   good = flag indicating whether this is with respect to a "good" generator,
//     for which shifts by i and i-n correspond to the same rotation
//   cost = total number of shifts needed by the collapsed network
//   solution = list indicating how levels in a standard benes network
//      are collapsed: if solution = [s_1 s_2 ... s_k], then k <= budget,
//      and the first s_1 levels are collapsed, the next s_2 levels are collapsed, etc.
{
  long k = GB_depth(n);
  long nlev = 2*k - 1;

  Vec< Vec<long> > costTab;

  buildBenesCostTable(n, k, good, costTab);

  // cout << good << "\n";
  // cout << costTab << "\n";

  BenesMemoTable memoTab;

  BenesMemoEntry t = optimalBenesAux(0, budget, nlev, costTab, memoTab);

  cost = t.cost;
  solution = t.solution;
  
}


// routines for finding the optimal splits among generators

// binary tree data structure for spliting generators

class SplitNode;
typedef tr1::shared_ptr<SplitNode> SplitNodePtr;


class SplitNode {
public:
  long order; // order associated with this node
  long mid;   // the "middle" token, 0 or 1
  LongNodePtr solution1, solution2; 
    // the benes solution(s)...for non-middle nodes,
    // there may be two different solutions, depending
    // on how the budget is split
  SplitNodePtr left, right;  // the children (for internal nodes)
  
  SplitNode(long _order, long _mid, LongNodePtr _solution1, LongNodePtr _solution2) {
  // constructor for leaves
    order = _order; mid = _mid; solution1 = _solution1; solution2 = _solution2;
    left = right = SplitNodePtr();
  }

  SplitNode(long _order, long _mid, SplitNodePtr _left, SplitNodePtr _right) {
  // constructor for internal nodes
    order = _order; mid = _mid; solution1 = solution2 = LongNodePtr();
    left = _left; right = _right;
  }

  bool isLeaf() const { return left == NULL && right == NULL; }
    
};

void print(ostream& s, SplitNodePtr p, bool first)
{
  if (p->isLeaf()) {
    if (!first) s << " ";
    s << "[";
    if (p->mid == 1) s << "* ";
    s << p->order << " " << p->solution1 << " " << p->solution2 << "]";
  }
  else {
    print(s, p->left, first);
    print(s, p->right, false);
  }
}


ostream& operator<<(ostream& s, SplitNodePtr p)
{
  s << "[";
  print(s, p, true);
  s << "]";
  return s;
}


// lower level memo table


class LowerMemoKey {
public:
  long order;
  bool good;
  long budget;
  long mid;

  LowerMemoKey(long _order, bool _good, long _budget, long _mid) {
    order = _order; good = _good; budget = _budget; mid = _mid;
  }

  size_t hash() const {
    stringstream s;
    s << order << " " << good << " " << budget << " " << mid;
    return tr1::hash<string>()(s.str());
  }

  bool operator==(const LowerMemoKey& other) const {
    return order == other.order && good == other.good &&
           budget == other.budget && mid == other.mid;
  }
};

class LowerMemoEntry {
public:
  long cost;
  SplitNodePtr solution;

  LowerMemoEntry(long _cost, SplitNodePtr _solution) {
    cost = _cost; solution = _solution;
  }

  LowerMemoEntry() { }
};



typedef tr1::unordered_map< LowerMemoKey, LowerMemoEntry, 
                            ClassHash<LowerMemoKey> > LowerMemoTable;


// list structure for managing generators

class GenNode;
typedef tr1::shared_ptr<GenNode> GenNodePtr;


class GenNode {
public:
  SplitNodePtr solution;  // the solution tree for this generator
  GenNodePtr next; // next node in the list

  GenNode(SplitNodePtr _solution, GenNodePtr _next) {
    solution = _solution; next = _next;
  }
};

long length(GenNodePtr ptr)
{
  long res = 0;
  for (GenNodePtr p = ptr; p != NULL; p = p->next) res++;
  return res;
}


ostream& operator<<(ostream& s, GenNodePtr p)
{
  if (p == NULL) {
    s << "[]";
    return s;
  }

  s << "[" << p->solution;
  p = p->next;
  while (p != NULL) {
    s << " " << p->solution;
    p = p->next;
  }
  cout << "]";

  return s;
}


// upper level memo table

class UpperMemoKey {
public:
  long i;
  long budget;
  long mid;

  UpperMemoKey(long _i, long _budget, long _mid) {
    i = _i; budget = _budget; mid = _mid;
  }

  size_t hash() const {
    stringstream s;
    s << i << " " << budget << " " << mid;
    return tr1::hash<string>()(s.str());
  }

  bool operator==(const UpperMemoKey& other) const {
    return i == other.i && budget == other.budget && mid == other.mid;
  }
};

class UpperMemoEntry {
public:
  long cost;
  GenNodePtr solution;

  UpperMemoEntry(long _cost, GenNodePtr _solution) {
    cost = _cost; solution = _solution;
  }

  UpperMemoEntry() { }
};



typedef tr1::unordered_map< UpperMemoKey, UpperMemoEntry, 
                            ClassHash<UpperMemoKey> > UpperMemoTable;


class GenDescriptor {
public:
  long order;
  bool good;

  GenDescriptor(long _order, bool _good) {
    order = _order; good = _good;
  }

  GenDescriptor() { }
};

LowerMemoEntry optimalLower(long order, bool good, long budget, long mid, 
                            LowerMemoTable& lowerMemoTable)
{
  assert(order > 1);
  assert(mid == 0 || mid == 1);
  assert(budget > 0);

  LowerMemoTable::iterator find = 
    lowerMemoTable.find(LowerMemoKey(order, good, budget, mid));

  if (find != lowerMemoTable.end()) {
    return find->second;
  }

  long cost;
  SplitNodePtr solution;

  if (mid == 0 && budget == 1) {
    // no solution possible

    cost = NTL_MAX_LONG;
    solution = SplitNodePtr();
  }
  else {
    // we first calculate a solution without splitting

    LongNodePtr benesSolution1, benesSolution2;

    if (mid == 1) { 
      // this is the middle node, so just one Benes network

      optimalBenes(order, budget, good, cost, benesSolution1);
      benesSolution2 = LongNodePtr();
    }
    else {
      // not the middle node, so we need two Benes networks.
      // if budget is odd, we split it unevenly

      long cost1, cost2;
      optimalBenes(order, budget/2, good, cost1, benesSolution1);
      if (budget % 2 == 0) {
        cost2 = cost1;
        benesSolution2 = benesSolution1;
      }
      else {
        optimalBenes(order, budget - budget/2, good, cost2, benesSolution2);
      }

      cost = cost1 + cost2;
    }

    solution = SplitNodePtr(new SplitNode(order, mid, 
                                              benesSolution1, benesSolution2));


    // now that we have a candidate solution that has no splits,
    // we recursively try splitting the order

    for (long order1 = 2; order1 < order; order1++) {
      if (order % order1 != 0) continue;

      bool good1 = good;
      bool good2 = good;

      if (good && GCD(order1, order/order1) != 1)
        good2 = false;

      // The logic is that if the problem is "good"
      // but the split is not relatively prime,
      // then only one of the subproblems is good.
      // since order1 ranges over all factors of order,
      // it suffices to choose the left subproblem to be
      // the good one.

      for (long budget1 = 1; budget1 < budget; budget1++) {
        for (long mid1 = 0; mid1 <= mid; mid1++) {
          LowerMemoEntry s1 = optimalLower(order1, good1, budget1, mid1, 
                                           lowerMemoTable);
          LowerMemoEntry s2 = optimalLower(order/order1, good2, budget-budget1, 
                                           mid-mid1, lowerMemoTable);
          if (s1.cost != NTL_MAX_LONG && s2.cost != NTL_MAX_LONG &&
              s1.cost + s2.cost < cost) {
            cost = s1.cost + s2.cost;
            solution = SplitNodePtr(new SplitNode(order, mid, s1.solution,
                                                      s2.solution));
          }
        }
      }
    }
  }

  return lowerMemoTable[LowerMemoKey(order, good, budget, mid)] =
           LowerMemoEntry(cost, solution);

  
}

UpperMemoEntry optimalUpperAux(const Vec<GenDescriptor>& vec, long i, long budget,
                               long mid, UpperMemoTable& upperMemoTable,
                               LowerMemoTable& lowerMemoTable) 
{
  assert(i >= 0 && i <= vec.length());
  assert(budget >= 0);
  assert(mid == 0 || mid == 1);

  UpperMemoTable::iterator find = upperMemoTable.find(UpperMemoKey(i, budget, mid));
  if (find != upperMemoTable.end()) {
    return find->second;
  }

  long cost;
  GenNodePtr solution;

  if (i == vec.length()) {
    // recursion stops here with trivial solution

    cost = 0;
    solution = GenNodePtr();
  }
  else if (budget == 0) {
    // recursion stops here with no solution

    cost = NTL_MAX_LONG;
    solution = GenNodePtr();
  }
  else {
    // allocate resources (budget, mid) between generator i and the rest

    long bestCost = NTL_MAX_LONG;
    LowerMemoEntry bestS;
    UpperMemoEntry bestT;

    for (long budget1 = 1; budget1 <= budget; budget1++) {
      for (long mid1 = 0; mid1 <= mid; mid1++) {
        LowerMemoEntry s = optimalLower(vec[i].order, vec[i].good, budget1, mid1,
                                        lowerMemoTable);
        UpperMemoEntry t = optimalUpperAux(vec, i+1, budget-budget1, mid-mid1,
                                           upperMemoTable, lowerMemoTable);
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

  return upperMemoTable[UpperMemoKey(i, budget, mid)] 
           = UpperMemoEntry(cost, solution);
}
                     

// optimalUpper is the high-level routine used to find
// an optimal strategy for evaluating a permutation
// The inputs are:
//   * vec: a vector describing the generators, each entry
//       consist of a long/bool pair (order, good)
//   * budget: the total budget for levels in the network,
//       which corresponds to the (multiplicative) circuit depth
// The outputs are:
//   * cost: the number of rotations needed (a specific
//       permutation could use fewer)
//   * solution: the optimal strategy, described below
//
// The solution is a linked list of length equal to the length of vec.
// The ith entry of the list corresponds to the ith entry in vec.
// Each entry in the list is a binary tree (SplitNodePtr).
// Every node in such a tree is either a leaf or an internal node
//   with two children.
// Every node contains the order associated with that node, along with
//   a "mid" flag: 1 means this is a "middle" node, 0 means not. 
// In the top-level list of trees, only one tree root will be a middle node,
//   and for every internal middle node in any tree, only one child
//   will be a middle node.
// A leaf node in a tree will contain two "benes network collapsing solutions",
//   solution1 and solution2.
//   If this leaf node is a middle node, solution2 is empty.
//   Otherwise, there are actually two solutions: this is because
//   the benes network corresponding to this leaf will appear twice
//   in the overall network, and if the budget allocated to this leaf
//   is *odd*, we may want to split the budget unevenly between the two
//   networks (floor(budget/2) and floor(budget/2)+1.
//   NOTE: one could consider more general splits of the budget
//   between these two occurrences of the network, but it seems doubtful
//   to be of benefit.
// A "benes network collapsing solution" is a list of longs (longNodePtr).
//   If the list is [n_1 n_2 ... n_k], this means we are to collapse
//   the first n_1 levels of the Benes network, then the next n_2 levels,
//   and so on.  The resulting "collapsed" Benes network will have k
//   levels.


  

void optimalUpper(const Vec<GenDescriptor>& vec, long budget,
                  long& cost, GenNodePtr& solution)
{
  assert(budget > 0);

  UpperMemoTable upperMemoTable;
  LowerMemoTable lowerMemoTable;

  UpperMemoEntry t = 
    optimalUpperAux(vec, 0, budget, 1, upperMemoTable, lowerMemoTable);

  cost = t.cost;
  solution = t.solution;
}




static ZZX makeIrredPoly(long p, long d)
{
	assert(d >= 1);
  assert(ProbPrime(p));

  if (d == 1) return ZZX(1, 1); // the monomial X

  zz_pBak bak; bak.save();
  zz_p::init(p);
  return to_ZZX(BuildIrred_zz_pX(d));
}

void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'R=4 L=9 k=80'\n\n";
  cerr << "  R is log_2 of the logical array size [default=1]\n";
  cerr << "  p is the plaintext base [default=2]" << endl;
  cerr << "  r is the lifting [default=1]" << endl;
  cerr << "  d is the degree of the field extension [default==1]\n";
  cerr << "    (d == 0 => factors[0] defined the extension)\n";
  cerr << "  c is number of columns in the key-switching matrices [default=2]\n";
  cerr << "  k is the security parameter [default=80]\n";
  cerr << "  L is the # of primes in the modulus chai [default=10]\n";
  cerr << "  s is the minimum number of slots [default=4]\n";
  cerr << "  m is a specific modulus\n";
  exit(0);
}

// divc(), rotate() defined in NumbTh.h

void rotateSlots(const EncryptedArray& ea, Vec< copied_ptr<Ctxt> >& v, long amt)
// copied_ptr is "smart pointer" with shallow-copy semantics, see cloned_ptr.h
{
  long nblocks = v.length();
  long nslots = ea.size();
  long N = nblocks * nslots;

  if (N == 0) return;

  amt %= N;
  if (amt < 0) amt += N;

  if (amt == 0) return;

  long q = amt / nslots;
  long r = amt % nslots;

 
  rotate(v, q);

  if (r != 0) {
    // construct appropriate mask: first r slots 0,
    // remaining slots are 1

    vector<long> mask;
    mask.resize(nslots);
    for (long j = 0; j < r; j++) mask[j] = 0;
    for (long j = r; j < nslots; j++) mask[j] = 1;

    ZZX pmask;
    ea.encode(pmask, mask); 

    const FHEPubKey& publicKey = v[0]->getPubKey();
    Ctxt c0(publicKey), c1(publicKey), carry(publicKey);

    for (long i = 0; i < nblocks; i++) {
      c0 = *v[i];
      ea.rotate(c0, r);
      c1 = c0;
      c1.multByConstant(pmask);
      c0 -= c1;
   
      c1 += carry;

      *v[i] = c1;
      carry = c0;

    }
    *v[0] += carry;
  }
}


void applyNetwork(const BenesNetwork& net, 
                  const EncryptedArray& ea, 
                  Vec< copied_ptr<Ctxt> >& v)
{
  long sz = net.getSize();
  long nlev = net.getNumLevels();

  long nblocks = v.length();
  long nslots = ea.size();
  long N = nblocks * nslots;

  assert(sz <= N);
  assert(sz > N - nslots);

  for (long i = 0; i < nlev; i++) {
    const Vec<short>& lev = net.getLevel(i);
    long p = net.bitPos(i);
    long shamt1 = -(1L << p);
    long shamt2 = +(1L << p);

    Vec< copied_ptr<Ctxt> > v1, v2;
    v1 = v; 
    v2 = v;

    for (long b = 0; b < nblocks; b++) {
      vector<long> mask1, mask2;
      mask1.resize(nslots);
      mask2.resize(nslots);

      for (long s = 0; s < nslots; s++) {
        long j = b*nslots + s;
        mask1[s] = mask2[s] = 0;
        if (j < sz && lev[j] != 0) {
          if (j & (1L << p)) 
            mask1[s] = 1;
          else
            mask2[s] = 1;
        }
      }

      ZZX pmask1, pmask2;
      ea.encode(pmask1, mask1);
      ea.encode(pmask2, mask2);

      v1[b]->multByConstant(pmask1);
      v2[b]->multByConstant(pmask2);
      *v[b] -= *v1[b];
      *v[b] -= *v2[b];
    }

    rotateSlots(ea, v1, shamt1);
    rotateSlots(ea, v2, shamt2);

    for (long b = 0; b < nblocks; b++) {
      *v[b] += *v1[b];
      *v[b] += *v2[b];
    }
  }
}





int main(int argc, char *argv[])
{

#if 1
  argmap_t argmap;
  argmap["n"] = "2";
  argmap["L"] = "10";
  argmap["good"] = "0";
  if (!parseArgs(argc, argv, argmap)) {
    cerr << "bad args\n";
    exit(0);
  }
  long n = atoi(argmap["n"]);
  long L = atoi(argmap["L"]);
  bool good = !!atoi(argmap["good"]);

  for (long iter=0; iter < 20; iter++) {
    Permut perm;
    RandomPerm(perm,n);      
    // cout << perm << "\n";
    GeneralBenesNetwork net(perm);   
    if (net.testNetwork(perm))  cout << ".";
    else                        cout << "X";
  }
  cout << "\n";

  long cost;
  LongNodePtr solution;

  optimalBenes(n, L, good, cost, solution);

  cout << cost << "\n";
  cout << solution << "\n";
  cout << length(solution) << "\n";

  cout << "*************\n";

  long cost1;
  GenNodePtr solution1;

  Vec<GenDescriptor> vec;

  vec.SetLength(2);
  vec[0] = GenDescriptor(n, good);
  vec[1] = GenDescriptor(101, good);
  optimalUpper(vec, L, cost1, solution1);

  cout << cost1 << "\n";
  cout << solution1 << "\n";
  
  

  
#endif

#if 0
  argmap_t argmap;
  argmap["R"] = "1";

  argmap["p"] = "2";  // plaintext space characteristic
  argmap["r"] = "1";  // if r>1 plaintext space is Z_{p^r}
  argmap["d"] = "1";  // if d>1 plaintext space is GF(p^d)
  argmap["c"] = "2";  // number of columns in key-switching matrices
  argmap["k"] = "80"; // security parameter
  argmap["L"] = "10"; // number of primes in modulus-chain
  argmap["s"] = "0";  // if s>0, need to have at least s slots
  argmap["m"] = "0";  // if m>0, use m'th cyclotomic polynomial

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long R = atoi(argmap["R"]);
  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long d = atoi(argmap["d"]);
  long c = atoi(argmap["c"]);
  long k = atoi(argmap["k"]);
  long s = atoi(argmap["s"]);
  long L = atoi(argmap["L"]);
  long chosen_m = atoi(argmap["m"]);

  long w = 64; // Hamming weight of secret key

  long m = FindM(k, L, c, p, d, s, chosen_m, true);

  long sz = 1L << R; // size of the permutation

  // Build the FHEContext for this instance of the cryptosystem
  FHEcontext context(m, p, r);
  buildModChain(context, L, c);
  context.zMStar.printout();

  // Generate a key-pair
  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(w); // A Hamming-weight-w secret key

  // Compute tables for the plaintext space
  ZZX G;
  if (d == 0)
    G = context.alMod.getFactorsOverZZ()[0];
  else
    G = makeIrredPoly(p, d); 

  cerr << "G = " << G << "\n";
  cerr << "generating key-switching matrices... ";
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  cerr << "done\n";

  cerr << "computing masks and tables for rotation...";
  EncryptedArray ea(context, G);
  cerr << "done\n";

  long nslots = ea.size(); // how many plaintext slots in each ciphertext
  long nblocks = divc(sz, nslots); // nblocks = ceiling(sz/nslots)
    // nblocks is # of ciphertexts needed to hold a permutation of zise sz
  long remslots = nslots*nblocks - sz; // # of leftover slots

  vector<long> mask; // a mask to zero-out the extra slots in the last ctxt
  mask.resize(nslots);
  for (long j = 0; j < nslots-remslots; j++) mask[j] = 1;
  for (long j = nslots-remslots; j < nslots; j++) mask[j] = 0;
  PlaintextArray pmask(ea);
  pmask.encode(mask); // encode the mask as a plaintext polynomial

  Vec< copied_ptr<Ctxt> > cvec; // a vector of pointers to ciphertexts
  cvec.SetLength(nblocks);      // allocate space and initialize to empty
  for (long i = 0; i < nblocks; i++) cvec[i].set_ptr(new Ctxt(publicKey));

  Vec< copied_ptr<PlaintextArray> > pvec; // vector of pointers to ptxt arrays
  pvec.SetLength(nblocks);
  for (long i = 0; i < nblocks; i++) pvec[i].set_ptr(new PlaintextArray(ea));

  for (long i = 0; i < nblocks; i++)
    pvec[i]->random();

  pvec[nblocks-1]->mul(pmask); // zero out leftover slots in last ptxt array
  // Shai: why do we care about these slots?

  // Print the plaintext before the permutation
  for (long i = 0; i < nblocks; i++) {
    pvec[i]->print(cout);
    cout << "\n";
  }
  cout << "\n";

  // Encrypt the plaintext arrays, then permute them
  for (long i = 0; i < nblocks; i++)
    ea.encrypt(*cvec[i], publicKey, *pvec[i]); // encrypt all the blocks

  // Choose a random permutation of size sz
  Permut perm;
  RandomPerm(perm, sz);
  cout << "perm = " << perm << "\n";

  // Setup a Benes network for this permutation
  BenesNetwork net(R, perm);

  // Apply the Benes network to the encrypted arrays
  applyNetwork(net, ea, cvec);

  // Decrypt the permuted arrays
  for (long i = 0; i < nblocks; i++)
    ea.decrypt(*cvec[i], secretKey, *pvec[i]);

  // Print the permuted plaintext
  for (long i = 0; i < nblocks; i++) {
    pvec[i]->print(cout);
    cout << "\n";
  }
  cout << "\n";

#endif
}

// benes_x R=4 p=47 m=46
