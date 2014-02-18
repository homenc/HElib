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
void timeOps(const EncryptedArray& ea, const FHEPubKey& publicKey,
	     Ctxt& c0, Ctxt& c1, Ctxt& c2, Ctxt& c3, ZZX& p, long nTests,
	     long nPrimes=0);

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

void timeOps(const EncryptedArray& ea, const FHEPubKey& publicKey, 
	     Ctxt& c0, Ctxt& c1, Ctxt& c2, Ctxt& c3, ZZX& p, long nTests,
	     long nPrimes)
{
  if (nPrimes>0) {  // perform operations at a lower level
    c0.modDownToLevel(nPrimes);
    c1.modDownToLevel(nPrimes);
    c2.modDownToLevel(nPrimes);
  }

  // Multiplication with 2,3 arguments
  cerr << "." << std::flush;
  for (long i=0; i<nTests; i++) {
    Ctxt cc1 = c1;
    startFHEtimer("multiplyBy");
    cc1.multiplyBy(c0);
    stopFHEtimer("multiplyBy");
    c3 += cc1; // Just so the compiler doesn't optimize it away
  }

  cerr << "." << std::flush;
  for (long i=0; i<nTests; i++) {
    Ctxt cc2 = c2;
    startFHEtimer("multiplyBy2");
    cc2.multiplyBy2(c0,c1);
    stopFHEtimer("multiplyBy2");
    c3 += cc2; // Just so the compiler doesn't optimize it away
  }

  // Multiply by constant
  cerr << "." << std::flush;
  for (long i=0; i<3*nTests; i++) {
    Ctxt cc0 = c0;
    startFHEtimer("multByConstant");
    cc0.multByConstant(p);
    stopFHEtimer("multByConstant");
    c3 -= cc0; // Just so the compiler doesn't optimize it away
  }

  // Add constant
  cerr << "." << std::flush;
  for (long i=0; i<10*nTests; i++) {
    Ctxt cc0 = c0;
    startFHEtimer("addConstant");
    cc0.addConstant(p);
    stopFHEtimer("addConstant");
    c3 += cc0; // Just so the compiler doesn't optimize it away
  }

  // Rotation by an amount k for which we have a key-switching matrix
  cerr << "." << std::flush;
  for (long i=0; i<nTests; i++) {
    Ctxt cc0 = c0;
    long k = rotationAmount(ea,publicKey,/*withMatrix=*/true);
    startFHEtimer("automorph-with-matrix");
    cc0.smartAutomorph(k);
    stopFHEtimer("automorph-with-matrix");    
    c3 += cc0; // Just so the compiler doesn't optimize it away
  }
  // Rotation by a random amount k
  cerr << "." << std::flush;
  for (long i=0; i<nTests; i++) {
    Ctxt cc0 = c0;
    long k = rotationAmount(ea,publicKey,/*withMatrix=*/false);
    startFHEtimer("automorph");
    cc0.smartAutomorph(k);
    stopFHEtimer("automorph");
    c3 += cc0; // Just so the compiler doesn't optimize it away
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
  
  PlaintextArray p0(ea);
  PlaintextArray p1(ea);
  PlaintextArray p2(ea);
  Ctxt c0(publicKey), c1(publicKey), c2(publicKey), c3(publicKey);

  p0.random();
  p1.random();
  p2.random();

  ZZX pp;
  startFHEtimer("encode");
  ea.encode(pp, p0);
  stopFHEtimer("encode");
  startFHEtimer("SK-encrypt");
  secretKey.Encrypt(c0, pp);
  startFHEtimer("SK-encrypt");
  startFHEtimer("PK-encrypt");
  publicKey.Encrypt(c0, pp);
  startFHEtimer("PK-encrypt");

  startFHEtimer("encode");
  ea.encode(pp, p1);
  stopFHEtimer("encode");
  startFHEtimer("SK-encrypt");
  secretKey.Encrypt(c1, pp);
  startFHEtimer("SK-encrypt");
  startFHEtimer("PK-encrypt");
  publicKey.Encrypt(c1, pp);
  startFHEtimer("PK-encrypt");

  startFHEtimer("encode");
  ea.encode(pp, p2);
  stopFHEtimer("encode");
  startFHEtimer("SK-encrypt");
  secretKey.Encrypt(c2, pp);
  startFHEtimer("SK-encrypt");
  startFHEtimer("PK-encrypt");
  publicKey.Encrypt(c2, pp);
  startFHEtimer("PK-encrypt");

  cerr << "Initialization time:\n";
  printAllTimers();
  cerr << endl;

  long nTests = 3*(53261/m); // more tests for smaller values of m
  resetAllTimers();
  cerr << "Operations with "<<L<<" primes in the chain: ";
  timeOps(ea, publicKey, c0,c1,c2,c3, pp, nTests);// timing with all primes in the chain
  cerr << endl;
  printAllTimers();
  cerr << endl;

  if (L >= 14) {
    L -= (L-2)/3;
    resetAllTimers();
    cerr << "Operations with "<<L<<" primes in the chain: ";
    timeOps(ea, publicKey, c0,c1,c2,c3, pp, nTests, L); // timing with 2/3 of the primes
    cerr << endl;
    printAllTimers();
    cerr << endl;
  }
  if (L >= 8) {
    L -= (L-2)/2;
    resetAllTimers();
    cerr << "Operations with "<<L<<" primes in the chain: ";
    timeOps(ea, publicKey, c0,c1,c2,c3, pp, nTests, L); // timing with 1/3 of the primes
    cerr << endl;
    printAllTimers();
    cerr << endl;
  }
  if (L >= 4) {
    resetAllTimers();
    cerr << "Operations with 2 primes in the chain: ";
    timeOps(ea, publicKey, c0,c1,c2,c3, pp, nTests, 2); // timing with only two primes
    cerr << endl;
    printAllTimers();
    cerr << endl;
  }
  resetAllTimers();

  startFHEtimer("decrypt");
  secretKey.Decrypt(pp, c0);
  startFHEtimer("decrypt");
  startFHEtimer("decode");
  ea.decode(p0, pp);
  stopFHEtimer("decode");

  startFHEtimer("decrypt");
  secretKey.Decrypt(pp, c1);
  startFHEtimer("decrypt");
  startFHEtimer("decode");
  ea.decode(p1, pp);
  stopFHEtimer("decode");

  startFHEtimer("decrypt");
  secretKey.Decrypt(pp, c2);
  startFHEtimer("decrypt");
  startFHEtimer("decode");
  ea.decode(p2, pp);
  stopFHEtimer("decode");

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
