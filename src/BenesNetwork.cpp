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
#include <NTL/lzz_pXFactoring.h>
#include <helib/EncryptedArray.h>

#include <cstdlib>
#include <list>
#include <sstream>
#include <memory>

#include <NTL/vector.h>
#include <helib/NumbTh.h>
#include <helib/permutations.h>

namespace helib {

static void recursiveGeneralBenesInit(long n,
                                      long k,
                                      long d,
                                      long delta_j,
                                      const Permut& perm,
                                      const Permut& iperm,
                                      NTL::Vec<NTL::Vec<short>>& level,
                                      NTL::Vec<NTL::Vec<short>>& ilevel)
{
  long sz = perm.length();

  if (d == k - 1) {
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
        level[d][delta_j + 1] = 0;
        ilevel[d][delta_j] = 0;
        ilevel[d][delta_j + 1] = 0;
      } else {
        // the swap perm

        level[d][delta_j] = 1;
        level[d][delta_j + 1] = -1;
        ilevel[d][delta_j] = 1;
        ilevel[d][delta_j + 1] = -1;
      }
    }
    return;
  }

  long nlev = 2 * k - 1;

  // size of upper internal network
  long sz0 = GeneralBenesNetwork::shamt(n, k, d);
  long sz1 = sz - sz0;

  assertTrue(labs(sz0 - sz1) <= 1l, "sz1 must be within 1 of sz0");

  // id_perm: the identity permutation on {0,...,sz-1}
  Permut id_perm;
  id_perm.SetLength(sz);
  for (long j = 0; j < sz; j++)
    id_perm[j] = j;

  // *xperm[0] is the perm on LHS of network
  // *xperm[1] is the perm on RHS of network
  // *xiperm[0] is the inv perm on LHS of network
  // *xiperm[1] is the inv perm on RHS
  const Permut* xperm[] = {&id_perm, &perm};
  const Permut* xiperm[] = {&id_perm, &iperm};

  // *first_level[0] is the first level when traversing left to right
  // *first_level[1] is the first level when traversing right to left
  // *ifirst_level[0] is the reversed first level when traversing left to right
  // *ifirst_level[1] is the reversed first level when traversing right to left
  NTL::Vec<short>* first_level[] = {&level[d], &ilevel[nlev - 1 - d]};
  NTL::Vec<short>* ifirst_level[] = {&ilevel[d], &level[nlev - 1 - d]};

  // *last_level[0] is the last level when traversing left to right
  // *last_level[1] is the last level when traversing right to left
  // *ilast_level[0] is the reversed last level when traversing left to right
  // *ilast_level[1] is the reversed last level when traversing right to left
  NTL::Vec<short>* last_level[] = {&level[nlev - 1 - d], &ilevel[d]};
  NTL::Vec<short>* ilast_level[] = {&ilevel[nlev - 1 - d], &level[d]};

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
  NTL::Vec<bool> marked[2];
  marked[0].SetLength(sz);
  for (long j = 0; j < sz; j++)
    marked[0][j] = false;
  marked[1].SetLength(sz);
  for (long j = 0; j < sz; j++)
    marked[1][j] = false;

  // counter[0] is used to count through nodes on left side
  // counter[1] is used to count through nodes on right side
  long counter[2];
  counter[0] = 0;
  counter[1] = 0;

  long dir = 0; // direction, initially left to right
  long e = 0;   // current edge, initially straight

  long node;      // current node
  long stop_node; // stopping node, initially undefined
  long stop_side; // stopping side

  if (sz0 == sz1) {
    // an even split
    node = 0;
    stop_node = sz0;
    stop_side = 0;
  } else if (sz0 > sz1) {
    // an odd split, top larger
    node = sz0 - 1;
    stop_node = sz0 - 1;
    stop_side = 1;
  } else { // sz0 < sz1
    // an odd split, bottom larger
    node = sz - 1;
    stop_node = sz - 1;
    stop_side = 1;
  }

  for (;;) {
    marked[dir][node] = true; // mark current node

    // traverse path through graph

    long v = (*xperm[dir])[node]; // value associated with this node

    long net = long(node >= sz0) ^ labs(e);
    // net = 0 => upper internal network
    // net = 1 => lower internal network

    long node1 = node - sz0 * long(node >= sz0);
    // node adjacent to node in the internal network

    long node3 = (*xiperm[1 - dir])[v];
    // corresponding node of the other side of the external network

    long node2 = node3 - sz0 * long(node3 >= sz0);
    // node adjacent to node3 in the internal network

    long e2 = (node3 - (node2 + net * sz0)) / sz0;
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

    // update external edges in level graph

    (*first_level[dir])[delta_j + node] = e;
    (*ifirst_level[dir])[delta_j + net * sz0 + node1] = -e;

    (*last_level[dir])[delta_j + net * sz0 + node2] = e2;
    (*ilast_level[dir])[delta_j + node3] = -e2;

    // update internal permutations

    inner_perm[net][dir][node2] = node1;
    inner_perm[net][1 - dir][node1] = node2;

    // change direction
    dir = 1 - dir;

    // mark node3 on new side
    marked[dir][node3] = true;

    // check for stop node
    if (node3 == stop_node && dir == stop_side) {
      // search for unmarked node

      while (counter[dir] < sz && marked[dir][counter[dir]])
        counter[dir]++;
      if (counter[dir] >= sz)
        break; // we're done!!

      // update node, e, stop_node, stop_side
      node = counter[dir];
      e = 0;

      if (node < sz0)
        stop_node = node + sz0;
      else
        stop_node = node - sz0;

      stop_side = dir;
    } else {
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
  recursiveGeneralBenesInit(n,
                            k,
                            d + 1,
                            delta_j,
                            inner_perm[0][0],
                            inner_perm[0][1],
                            level,
                            ilevel);

  // recurse on lower internal network
  recursiveGeneralBenesInit(n,
                            k,
                            d + 1,
                            delta_j + sz0,
                            inner_perm[1][0],
                            inner_perm[1][1],
                            level,
                            ilevel);
}

GeneralBenesNetwork::GeneralBenesNetwork(const Permut& perm)
{
  n = perm.length();

  // check that n > 1
  assertTrue<InvalidArgument>(n > 1l,
                              "permutation length must be greater than one");

  // compute recursion depth k = least integer k s/t 2^k >= n
  k = GeneralBenesNetwork::depth(n);

  // construct the inverse perm, which is convenient in the recursive function.
  // As a byproduct, verify that perm is indeed a permutation on {0,...,n-1}.

  Permut iperm;
  iperm.SetLength(n);

  for (long j = 0; j < n; j++)
    iperm[j] = -1;

  for (long j = 0; j < n; j++) {
    long j1 = perm[j];
    assertInRange(j1, 0l, n, "permutation element out of range");
    iperm[j1] = j;
  }

  for (long j = 0; j < n; j++) {
    assertTrue(iperm[j] != -1l, "permutation element not processed");
  }

  // allocate space for the levels graph

  level.SetLength(2 * k - 1);
  for (long i = 0; i < 2 * k - 1; i++)
    level[i].SetLength(n);

  // allocate space for the reverse levels graph...
  // makes the recursive construction more convenient

  NTL::Vec<NTL::Vec<short>> ilevel;
  ilevel.SetLength(2 * k - 1);
  for (long i = 0; i < 2 * k - 1; i++)
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
      const NTL::Vec<short>& lev = getLevel(i);
      j1 += shamt(i) * lev[j1];
    }
    if (perm[j1] != j)
      return false;
  }
  return true;
}

} // namespace helib
