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
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include <NTL/lzz_pXFactoring.h>

#include <cassert>
#include <cstdio>

#include <tr1/memory> 
// for tr1::shared_ptr<T>

/********

Experimental program for equality testing...

********/

typdef tr1::shared_ptr<Ctxt> CtxtPtr;

void printVector(const vector<ZZX>& v) 
{
}

void computeIt(vector<CtxtPtr>& res, const Ctxt& ctxt, unsigned long w, long n)
{

}


void  TestIt(long L, long m, long c)
{
  FHEcontext context(m, 2, 1); // p = 2, r = 1
  long d = context.zMStar.getOrdP(); 

  buildModChain(context, L, c);

  context.zMStar.printout();
  cerr << endl;
#ifdef DEBUG
  cerr << context << endl;
#endif

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(w); // A Hamming-weight-w secret key


  ZZX G;

  G = makeIrredPoly(p, d); 

  cerr << "generating key-switching matrices... ";
  addFrbMatrices(secretKey);
  cerr << "done\n";


  cerr << "computing masks and tables for rotation...";
  EncryptedArray ea(context, G);
  cerr << "done\n";



  long nslots = ea.size();

  double t = GetTime();
  resetAllTimers();

  vector<ZZX> v;
  v.resize(nslots);
  for (long i = 0; i < nslots; i++) {
    GF2X f;
    random(f, d);
    conv(v[i], f);
  }

  printVector(v);

  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, v);
  // ctxt encrypts a vector where each slots is a random
  // polynomial of degree < d

  unsigned long w = RandomBits_ulong(d);
  // w is a random d-bit number

  long n = d;


  vector<CtxtPtr> res;
  res.resize(n);
  for (long j = 0; j < n; j++) {
    res[j] = CtxtPtr(new Ctxt(publicKey));
  }

  computeIt(res, ctxt, w, n);

  for (long j = 0; j < n; j++) {
    vector<ZZX> v1;
    ea.decrypt(*res[j], secretKey, v1); 
    printVector(v1);
  }

}


void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'R=4 L=9 k=80'\n\n";
  cerr << "  R is the number of rounds\n";
  cerr << "  p is the plaintext base [default=2]" << endl;
  cerr << "  r is the lifting [default=1]" << endl;
  cerr << "  d is the degree of the field extension [default==1]\n";
  cerr << "    (d == 0 => factors[0] defined the extension)\n";
  cerr << "  c is number of columns in the key-switching matrices [default=2]\n";
  cerr << "  k is the security parameter [default=80]\n";
  cerr << "  L is the # of primes in the modulus chai [default=4*R]\n";
  cerr << "  s is the minimum number of slots [default=4]\n";
  cerr << "  m is a specific modulus\n";
  exit(0);
}


int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["R"] = "1";
  argmap["p"] = "2";
  argmap["r"] = "1";
  argmap["d"] = "1";
  argmap["c"] = "2";
  argmap["k"] = "80";
  argmap["L"] = "0";
  argmap["s"] = "0";
  argmap["m"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long R = atoi(argmap["R"]);
  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long d = atoi(argmap["d"]);
  long c = atoi(argmap["c"]);
  long k = atoi(argmap["k"]);
  //  long z = atoi(argmap["z"]);
  long L = atoi(argmap["L"]);
  if (L==0) { // determine L based on R,r
    L = 3*R+3;
    if (p>2 || r>1) { // add some more primes for each round
      long addPerRound = 2*ceil(log((double)p)*r*3)/(log(2.0)*NTL_SP_NBITS) +1;
      L += R * addPerRound;
    }
  }
  long s = atoi(argmap["s"]);
  long chosen_m = atoi(argmap["m"]);

  long w = 64; // Hamming weight of secret key
  //  long L = z*R; // number of levels

  long m = FindM(k, L, c, p, d, s, chosen_m, true);

  setTimersOn();
  TestIt(R, p, r, d, c, k, w, L, m);

  cerr << endl;
  printAllTimers();
  cerr << endl;

}

// call to get our running test case:
// Test_General_x p=23 m=20485 L=20 R=5
