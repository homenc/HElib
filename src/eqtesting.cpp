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
/*
 * @file eqtesting.cpp
 * @brief Experimental program for equality testing...
 */
#include <NTL/lzz_pXFactoring.h>
NTL_CLIENT
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"

#include <cassert>
#include <cstdio>


// computes ctxt^{2^d-1} using a method that takes
// O(log d) automorphisms and multiplications
void fastPower(Ctxt& ctxt, long d) 
{
  if (d == 1) return;

  Ctxt orig = ctxt;

  long k = NumBits(d);
  long e = 1;

  for (long i = k-2; i >= 0; i--) {
    Ctxt tmp1 = ctxt;
    tmp1.smartAutomorph(1L << e);
    ctxt.multiplyBy(tmp1);
    e = 2*e;

    if (bit(d, i)) {
      ctxt.smartAutomorph(2);
      ctxt.multiplyBy(orig);
      e += 1;
    }
  }
}

// incrementalZeroTest sets each res[i], for i=0..n-1, to
// a ciphertext in which each slot is 0 or 1 according
// to whether or not bits 0..i of corresponding slot in ctxt
// is zero (1 if not zero, 0 if zero).
// It is assumed that res and each res[i] is already initialized
// by the caller.
// Complexity: O(d + n log d) smart automorphisms
//             O(n d) 
void incrementalZeroTest(Ctxt* res[], const EncryptedArray& ea,
			 const Ctxt& ctxt, long n)
{
  FHE_TIMER_START;
  long nslots = ea.size();
  long d = ea.getDegree();

  // compute linearized polynomial coefficients

  vector< vector<ZZX> > Coeff;
  Coeff.resize(n);

  for (long i = 0; i < n; i++) {
    // coeffients for mask on bits 0..i
    // L[j] = X^j for j = 0..i, L[j] = 0 for j = i+1..d-1

    vector<ZZX> L;
    L.resize(d);

    for (long j = 0; j <= i; j++) 
      SetCoeff(L[j], j);

    vector<ZZX> C;

    ea.buildLinPolyCoeffs(C, L);

    Coeff[i].resize(d);
    for (long j = 0; j < d; j++) {
      // Coeff[i][j] = to the encoding that has C[j] in all slots
      // FIXME: maybe encrtpted array should have this functionality
      //        built in
      vector<ZZX> T;
      T.resize(nslots);
      for (long s = 0; s < nslots; s++) T[s] = C[j];
      ea.encode(Coeff[i][j], T);
    }
  }

  vector<Ctxt> Conj(d, ctxt);
  // initialize Cong[j] to ctxt^{2^j}
  for (long j = 0; j < d; j++) {
    Conj[j].smartAutomorph(1L << j);
  }

  for (long i = 0; i < n; i++) {
    res[i]->clear();
    for (long j = 0; j < d; j++) {
      Ctxt tmp = Conj[j];
      tmp.multByConstant(Coeff[i][j]);
      *res[i] += tmp;
    }

    // *res[i] now has 0..i in each slot
    // next, we raise to the power 2^d-1

    fastPower(*res[i], d);
  }
  FHE_TIMER_STOP;
}

#ifdef DEBUG_TEST
void printBits(const vector<ZZX>& v, long n) 
{
  long len = v.size();
  if (n>50 || len>32) return;
  for (long i = 0; i < len; i++) {
    for (long j = n-1; j >= 0; j--)
      cout << coeff(v[i], j);
    cout << " ";
  }
  cout << "\n";
}


void  TestIt(long c, long k, long w, long L, long m, long n)
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

  G = makeIrredPoly(2, d); 
  // G = context.alMod.getFactorsOverZZ()[0];

  cerr << "generating key-switching matrices... ";
  addFrbMatrices(secretKey);
  addSome1DMatrices(secretKey);
  cerr << "done\n";


  cerr << "computing masks and tables for rotation...";
  EncryptedArray ea(context, G);
  cerr << "done\n";



  long nslots = ea.size();

  if (n <= 0 || n > d) n = d;
  


  vector<ZZX> v;
  v.resize(nslots);
  for (long i = 0; i < nslots; i++) {
    GF2X f;
    random(f, n);
    conv(v[i], f);
  }

  printBits(v, n);

  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, v);
  // ctxt encrypts a vector where each slots is a random
  // polynomial of degree < n


  Ctxt* res[n];
  for (long j = 0; j < n; j++) res[j] = new Ctxt(publicKey); // allocate

  resetAllTimers();

  incrementalZeroTest(res, ea, ctxt, n);

  for (long j = 0; j < n; j++) {
    vector<ZZX> v1;
    ea.decrypt(*res[j], secretKey, v1); 
    printBits(v1, n);
  }

  for (long j = 0; j < n; j++) delete res[j]; // cleanup
}


void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'R=4 L=9 k=80'\n\n";
  cerr << "  c is number of columns in the key-switching matrices [default=2]\n";
  cerr << "  k is the security parameter [default=80]\n";
  cerr << "  L is the # of primes in the modulus chain [default=20]\n";
  cerr << "  m is a specific modulus\n";
  cerr << "  n is the number of masks (defaults to all)\n";
  exit(0);
}


int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["c"] = "2";
  argmap["k"] = "80";
  argmap["L"] = "20";
  argmap["m"] = "0";
  argmap["n"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long c = atoi(argmap["c"]);
  long k = atoi(argmap["k"]);
  long L = atoi(argmap["L"]);
  long chosen_m = atoi(argmap["m"]);
  long n = atoi(argmap["n"]);

  long w = 64; // Hamming weight of secret key
  //  long L = z*R; // number of levels

  long m = FindM(k, L, c, 2, 1, 0, chosen_m, true);

  setTimersOn();
  TestIt(c, k, w, L, m, n);

  cerr << endl;
  printAllTimers();
  cerr << endl;

}

// call to get our running test case:
// eqtesting_x m=20485 
// eqtesting_x m=105 for quick testing

#endif
