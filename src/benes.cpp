
#include <cstdlib>
#include <cassert>
#include <NTL/vector.h>
#include "NumbTh.h"

using namespace std;
using namespace NTL;



class BenesNetwork {
private:

  long k; // Benes network of size 2^k
  Vec< Vec<short> > level;
    // 2*k - 1 levels
    // each level has 2^k nodes
    // level[i][j] = 0, 1; 0 means a straight edge,
    //   1 means a slanted edge (at node j of level i); 
    // a slanted edge from node j points to node j + delta
    //   at the next level, where delta = 2^p if bit p
    //   of j is 0, and delta = -2^p if bit p of j is 1;
    //   here, p is the bit position associated with level i,
    //   which is |k-1-i|

  BenesNetwork(); // disabled
  

public:

  long getK() const { return k; }
  long getSize() const { return 1L << k; }
  long getNumLevels() const { return 2*k-1; }

  const Vec<short>& getLevel(long i) const { return level[i]; }

  long bitPos(long i) const { return labs(k - 1 - i); }
    // bit position associated with level i

  BenesNetwork(long _k, const Vec<long>& perm); // constructor



};




void recursiveBenesInit(long k, long delta_i, long delta_j,
                        const Vec<long>& perm, 
                        const Vec<long>& iperm,
                        Vec< Vec<short> >& level,
                        Vec< Vec<short> >& ilevel)
{

  if (k == 1) {
    // recursion stops here...only two possibolities for
    // perm: the identity perm or the swap perm

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
  Vec<long> id_perm;
  id_perm.SetLength(sz);
  for (long j = 0; j < sz; j++) id_perm[j] = j;
  
  // *xperm[0] is the perm on LHS of network
  // *xperm[1] is the perm on RHS of network
  // *xiperm[0] is the inv perm on LHS of network
  // *xiperm[1] is the inv perm on RHS
  const Vec<long> *xperm[] = { &id_perm, &perm };
  const Vec<long> *xiperm[] = { &id_perm, &iperm };

  // *first_level[0] is the first level when traversing left to right
  // *ifirst_level[1] is the reversed first level when traversing right to left
  Vec<short> *first_level[] = { &level[delta_i], &ilevel[nlev-1-delta_i] };
  Vec<short> *ifirst_level[] = { &ilevel[delta_i], &level[nlev-1-delta_i] };

  // *last_level[0] is the last level when traversing left to right
  // *ilast_level[1] is the reversed last level when traversing right to left
  Vec<short> *last_level[] = { &level[nlev-1-delta_i], &ilevel[delta_i] };
  Vec<short> *ilast_level[] = { &ilevel[nlev-1-delta_i], &level[delta_i] };

  // inner_perm[0][0] upper internal perm
  // inner_perm[0][1] upper internal inv perm
  // inner_perm[1][0] lower internal perm
  // inner_perm[1][1] lower internal inv perm
  Vec<long> inner_perm[2][2];

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
    if (n == -1) { 
      // scan for unmarked node

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

             e               e2
        n ------ n1 ---- n2 ----- n3
                \           /
                 \_________/
                      |
               internal network

    */

    // print debug info

#if 0

    cout << "d=" << d << " ";
    cout << "net=" << net << " ";
    cout << "n=" << n << " ";
    cout << "e=" << e  << " ";
    cout << "n1=" << n1 << " ";
    cout << "n2=" << n2 << " ";
    cout << "e2=" << e2 << " ";
    cout << "n3=" << n3 << " ";
    cout << "\n";
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
       

BenesNetwork::BenesNetwork(long _k, const Vec<long>& perm)
{
  k = _k;
  assert(k > 0);
  assert(k < NTL_BITS_PER_LONG-1);

  long sz = 1L << k;
  assert(perm.length() == sz);

  // construct the inverse perm, which is convenient in
  // the recursive construction...as a byproduct, this
  // also verifies that perm is indeed a permutation on
  // {0, .., 2^k-1 }

  Vec<long> iperm;
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


bool testNetwork(const BenesNetwork& net, const Vec<long>& perm)
{
  long sz = net.getSize();
  long nlev = net.getNumLevels();

  for (long j = 0; j < sz; j++) {
    // find correct position for j

    long j1 = j;

    for (long i = 0; i < nlev; i++) {
      const Vec<short>& lev = net.getLevel(i);
      if (lev[j1] == 0) continue;

      long p = net.bitPos(i);
      long delta = (j1 & (1L << p)) ? -(1L << p) : (1L << p);
      j1 += delta;
    }

    if (perm[j1] != j) return false;
  }

  return true;

}

int main(int argc, char *argv[])
{
   argmap_t argmap;
   argmap["k"] = "2";

   if (!parseArgs(argc, argv, argmap)) {
      cerr << "bad args\n";
      exit(0);
   }


   long k = atoi(argmap["k"]);
   long n = 1L << k;

   for (long iter=0; iter < 10; iter++) {

      Vec<long> perm;
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
   
      BenesNetwork net(k, perm);
   
      if (testNetwork(net, perm))
         cout << ".";
      else
         cout << "X";
   }

   cout << "\n";

}
