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
#include <cassert>
#include <cstdio>

#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"

#define pSize (NTL_SP_NBITS/2) /* The size of levels in the chain */

long rotationAmount(const EncryptedArray& ea, const FHEPubKey& publicKey,
	       bool withMatrix=false);
void timeOps(const EncryptedArray& ea, const FHEPubKey& publicKey, Ctxt& ret,
	     const vector<Ctxt>& c, ZZX& p, long nTests, long nPrimes=0);

// Returns either a random automorphism amount or an amount
// for which we have a key-switching matrix s^k -> s.
long rotationAmount(const EncryptedArray& ea, const FHEPubKey& publicKey,
	       bool onlyWithMatrix)
{
  const PAlgebra& pa = ea.getContext().zMStar;
  long nSlots = pa.getNSlots();
  long r = RandomBnd(nSlots);
  long k = pa.ith_rep(r);
  if (onlyWithMatrix) { // return the 1st step in the path to k
    const KeySwitch& matrix = publicKey.getNextKSWmatrix(k,0);
    k = matrix.fromKey.getPowerOfX();
  }
  return k;
}

void timeOps(const EncryptedArray& ea, const FHEPubKey& publicKey, Ctxt& ret,
	     const vector<Ctxt>& c, ZZX& p, long nTests, long nPrimes)
{
  assert(c.size()>=3);
  vector<Ctxt> cc = c;
  // perform operations at a lower level
  if (nPrimes>0) for (long i=0; i<(long)cc.size(); i++)
    cc[i].modDownToLevel(nPrimes);
  else 
    nPrimes = cc[0].findBaseLevel();

  // inner-product of size-5 vectors
  if (nPrimes > 2) {
    cerr << "." << std::flush;
    startFHEtimer("dimension-5-innerProduct");
    innerProduct(ret,cc,cc);
    ret.modDownToLevel(ret.findBaseLevel());
    stopFHEtimer("dimension-5-innerProduct");
  }

  // Multiplication with 2,3 arguments
  cerr << "." << std::flush;
  for (long i=0; i<nTests; i++) {
    Ctxt c0 = cc[0];
    startFHEtimer("multiplyBy");
    c0.multiplyBy(cc[1]);
    c0.modDownToLevel(c0.findBaseLevel());      // mod-down if needed
    stopFHEtimer("multiplyBy");
    ret += c0; // Just so the compiler doesn't optimize it away
  }

  if (nPrimes > 2) {
    cerr << "." << std::flush;
    for (long i=0; i<nTests; i++) {
      Ctxt c0 = cc[0];
      startFHEtimer("multiplyBy2");
      c0.multiplyBy2(cc[1],cc[2]);
      c0.modDownToLevel(c0.findBaseLevel());      // mod-down if needed
      stopFHEtimer("multiplyBy2");
      ret += c0; // Just so the compiler doesn't optimize it away
    }
  }

  // Multiply by constant
  cerr << "." << std::flush;
  for (long i=0; i<3*nTests; i++) {
    Ctxt c0 = cc[0];
    startFHEtimer("multByConstant");
    c0.multByConstant(p);
    stopFHEtimer("multByConstant");
    ret -= c0; // Just so the compiler doesn't optimize it away
  }

  // Add constant
  cerr << "." << std::flush;
  for (long i=0; i<10*nTests; i++) {
    Ctxt c0 = cc[0];
    startFHEtimer("addConstant");
    c0.addConstant(p);
    stopFHEtimer("addConstant");
    ret += c0; // Just so the compiler doesn't optimize it away
  }

  // Rotation by an amount k for which we have a key-switching matrix
  cerr << "." << std::flush;
  for (long i=0; i<nTests; i++) {
    Ctxt c0 = cc[0];
    long k = rotationAmount(ea,publicKey,/*withMatrix=*/true);
    startFHEtimer("automorph-with-matrix");
    c0.smartAutomorph(k);
    c0.modDownToLevel(c0.findBaseLevel());      // mod-down if needed
    stopFHEtimer("automorph-with-matrix");    
    ret += c0; // Just so the compiler doesn't optimize it away
  }
  // Rotation by a random amount k
  cerr << "." << std::flush;
  for (long i=0; i<nTests; i++) {
    Ctxt c0 = cc[0];
    long k = rotationAmount(ea,publicKey,/*withMatrix=*/false);
    startFHEtimer("automorph");
    c0.smartAutomorph(k);
    c0.modDownToLevel(c0.findBaseLevel());      // mod-down if needed
    stopFHEtimer("automorph");
    ret += c0; // Just so the compiler doesn't optimize it away
  }
}

void  TimeIt(long m, long p, long r, bool d_eq_1)
{
  cerr << "\n\n******** TimeIt: m=" << m
       << ", p=" << p
       << ", r=" << r
       << ", d=" << d_eq_1
       << endl;
  resetAllTimers();

  // Build the context, with as many primes as possible wrt to secparam=80
  startFHEtimer("Initialize context");
  FHEcontext context(m, p, r);
  stopFHEtimer("Initialize context");

  long phim = context.zMStar.getPhiM();
  long L = floor((7.2*phim)/(pSize* /*cc*/1.33* (110+/*k*/80))) -1;
  assert(L>1); // Make sure we have at least a few primes

  startFHEtimer("buildModChain");
  buildModChain(context, L, /*c=*/3);
  stopFHEtimer("buildModChain");

  startFHEtimer("keyGen");
  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(64); // A Hamming-weight-64 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  stopFHEtimer("keyGen");

  ZZX G;
  if (d_eq_1)
    SetX(G); // set G(X)=X
  else 
    G = context.alMod.getFactorsOverZZ()[0];

  startFHEtimer("Initialize EncryptedArray");
  EncryptedArray ea(context, G);
  stopFHEtimer("Initialize EncryptedArray");
  
  ZZX poly;
  PlaintextArray pp(ea);
  vector<PlaintextArray> vp(5,pp);
  for (long i=0; i<(long)vp.size(); i++)
    vp[i].random();

  Ctxt cc(publicKey);
  vector<Ctxt> vc(5,cc);

  for (long i=0; i<(long)vc.size(); i++) {
    startFHEtimer("encode");
    ea.encode(poly, vp[i]);
    stopFHEtimer("encode");
    startFHEtimer("SK-encrypt");
    secretKey.Encrypt(vc[i], poly);
    stopFHEtimer("SK-encrypt");
    startFHEtimer("PK-encrypt");
    publicKey.Encrypt(vc[i], poly);
    stopFHEtimer("PK-encrypt");
  }

  cerr << "Initialization time:\n";
  printAllTimers();
  resetAllTimers();
  cerr << endl;

  long nTests = 5;

  for (long i=2; i<L; i*=2) {
    cerr << "Operations with "<<i<<" primes in the chain: ";
    timeOps(ea, publicKey, cc,vc, poly, nTests, i);
    cerr << endl;
    printAllTimers();
    resetAllTimers();
    cerr << endl;
  }
  cerr << "Operations with "<<L<<" primes in the chain: ";
  timeOps(ea, publicKey, cc,vc, poly, nTests);
  cerr << endl;
  printAllTimers();
  resetAllTimers();
  cerr << endl;

  for (long i=0; i<(long)vc.size(); i++) {
    startFHEtimer("decrypt");
    secretKey.Decrypt(poly, vc[i]);
    stopFHEtimer("decrypt");
    startFHEtimer("decode");
    ea.decode(vp[i], poly);
    stopFHEtimer("decode");
  }

  cerr << "Decoding/decryption time:\n";
  printAllTimers();
}


void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'm=11441 p=2 r=1'\n\n";
  cerr << "  m determines the cyclotomic ring, defaults to all the set\n";
  cerr << "    m in { 4051, 4369, 4859, 10261,11023,11441,\n";
  cerr << "          18631,20485,21845, 49981,53261       }\n";
  cerr << "  p is the plaintext base [default=2]" << endl;
  cerr << "  r is the lifting [default=1]" << endl;
  cerr << "  d is the degree of the field extension [default==1]\n";
  cerr << "    (d == 0 => factors[0] defined the extension)\n";
  exit(0);
}

int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["p"] = "2";
  argmap["r"] = "1";
  argmap["m"] = "0";
  argmap["d"] = "1";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long m = atoi(argmap["m"]);
  long d = atoi(argmap["d"]);

  long ms[11] = { 4051, 4369, 4859, 10261,11023,11441,
		 18631,20485,21845, 49981,53261};
  if (m>0)
    TimeIt(m, p, r, (d==1));
  else for (long i=0; i<11; i++)
      TimeIt(ms[i], p, r, (d==1));
}
