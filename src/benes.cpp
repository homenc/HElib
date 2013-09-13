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

  long getK() const { return k; }
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

  // Pad the permutation to size 2^k. If it is smaller then the extra
  // entries are set to the identity permutation. 
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
    // only two possibolities for perm: the identity perm or the swap perm

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
  argmap_t argmap;
  argmap["R"] = "1";

#if 0
  if (!parseArgs(argc, argv, argmap)) {
    cerr << "bad args\n";
    exit(0);
  }
  long R = atoi(argmap["R"]);
  long n = 1L << R;

  for (long iter=0; iter < 10; iter++) {
    Permut perm;
    randomPerm(perm,n);      
    BenesNetwork net(R, perm);   
    if (testNetwork(net, perm)) cout << ".";
    else                        cout << "X";
  }
  cout << "\n";
#endif

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
  randomPerm(perm, sz);
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
}

// benes_x R=4 p=47 m=46
