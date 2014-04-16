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
/**
 * @file Test_PolyEval.cpp
 * @brief Homomorphic Polynomial Evaluation
 */
#include "polyEval.h"
#include "EncryptedArray.h"

void testIt(long d, long k, long p, long r, long m, long L,
	    bool isMonic=false)
{
  FHEcontext context(m, p, r);
  buildModChain(context, L, /*c=*/3);
  EncryptedArray ea(context);

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(/*w=*/64);// A Hamming-weight-64 secret key
  //  addSome1DMatrices(secretKey); // compute key-switching matrices
  Ctxt inCtxt(publicKey), outCtxt(publicKey);

#ifdef DEBUG_PRINTOUT
  ea_pt = &ea;        // for debugging purposes
  sk_pt = &secretKey;
#endif

  // evaluate at random points (at least one co-prime with p)
  long p2r = context.alMod.getPPowR();
  vector<long> x;
  ea.random(x);
  while (GCD(x[0],p)!=1) { x[0] = RandomBnd(p2r); }
  ea.encrypt(inCtxt, publicKey, x);

  ZZX poly;
  for (long i=d; i>=0; i--)
    SetCoeff(poly, i, RandomBnd(p)); // coefficients are random
  if (isMonic) SetCoeff(poly, d);    // set top coefficient to 1

  // Evaluate poly on the ciphertext
  cout << "* degree-"<<d<<", m="<<m<<", L="<<L<<", p^r="<<p2r<<endl;
  polyEval(outCtxt, poly, inCtxt, k);

  // Check the result
  vector<long> y;
  ea.decrypt(outCtxt, secretKey, y);
  for (long i=0; i<ea.size(); i++) {
    long ret = polyEvalMod(poly, x[i], p2r);
    if (ret != y[i]) {
      cerr << "ouch, p="<<p2r<<", poly="<<poly<<", x["<<i<<"]="<<x[i]
	   << ", y="<<y[i] << "!=" << ret << endl;
      exit(0);
    }
  }
  cout << " okay\n" << std::flush;
}

void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  p is the plaintext base [default=3]" << endl;
  cerr << "  r is the lifting [default=2]" << endl;
  cerr << "  m is a specific cyclotomic ring\n";
  cerr << "  d is the polynomial degree [default=undefined]" << endl;
  cerr << "    d=undefined means trying a few powers d=1,...,4,25,...,34"<<endl;
  cerr << "  k is the baby-step parameter [default=undefined]" << endl;
  cerr << "    if k is undefined it is computed from d" << endl;
  exit(0);
}

int main(int argc, char *argv[])
{
  argmap_t argmap;
  argmap["p"] = "3";
  argmap["r"] = "2";
  argmap["m"] = "0";
  argmap["d"] = "-1";
  argmap["k"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long m = atoi(argmap["m"]);
  long d = atoi(argmap["d"]);
  long k = atoi(argmap["k"]);

  long max_d = (d<=0)? 35 : d;
  long L = 1+NextPowerOfTwo(max_d);
  if (m<2)
    m = FindM(/*secprm=*/80, L, /*c=*/3, p, 1, 0, m, true);

  // Test both monic and non-monic polynomials of this degree
  if (d>=0) {
    testIt(d, k, p, r, m, L, false);
    testIt(d, k, p, r, m, L, false);
    testIt(d, k, p, r, m, L, false);
    testIt(d, k, p, r, m, L, true);
    return 0;
  }

  // Test degrees 1 to 4 and 25 through 35
  testIt(1, k, p, r, m, L, true);
  testIt(2, k, p, r, m, L, true);
  testIt(3, k, p, r, m, L, true);
  testIt(4, k, p, r, m, L, true);
  for (d=25; d<=35; d++)
    testIt(d, k, p, r, m, L);

  return 0;
}
