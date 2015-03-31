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
 * Test_PolyEval.cpp - Homomorphic Polynomial Evaluation
 */
#include <NTL/ZZ.h>
NTL_CLIENT
#include "polyEval.h"
#include "EncryptedArray.h"

#ifdef DEBUG_PRINTOUT
extern FHESecKey* dbgKey;
extern EncryptedArray* dbgEa;
#endif

bool testEncrypted(long d, const EncryptedArray& ea,
		   const FHESecKey& secretKey)
{
  const FHEcontext& context = ea.getContext();
  const FHEPubKey& publicKey = secretKey;
  long p = publicKey.getPtxtSpace();
  zz_pBak bak; bak.save(); zz_p::init(p);
  zz_pXModulus phimX = conv<zz_pX>(context.zMStar.getPhimX());

  // Choose random plaintext polynomials
  zz_pX pX = random_zz_pX(deg(phimX)-1);
  Vec<zz_pX> ppoly(INIT_SIZE, d);
  for (long i=0; i<ppoly.length(); i++) random(ppoly[i], deg(phimX)-1);

  // Evaluate the non-encrypted polynomial
  zz_pX pres = (ppoly.length()>0)? ppoly[ppoly.length()-1] : zz_pX::zero();
  for (long i=ppoly.length()-2; i>=0; i--) {
    MulMod(pres, pres, pX, phimX);
    pres += ppoly[i];
  }

  // Encrypt the random polynomials
  Ctxt cX(publicKey);
  Vec<Ctxt> cpoly(INIT_SIZE, d, cX);
  secretKey.Encrypt(cX, conv<ZZX>(pX));
  for (long i=0; i<ppoly.length(); i++)
    secretKey.Encrypt(cpoly[i], conv<ZZX>(ppoly[i]));

  // Evaluate the encrypted polynomial
  polyEval(cX, cpoly, cX);

  // Compare the results
  ZZX ret;
  secretKey.Decrypt(ret, cX);
  zz_pX cres = conv<zz_pX>(ret);
  bool success = (cres == pres);
  if (success) cout << " encrypted poly match, ";
  else         cout << " encrypted poly MISMATCH\n";
  return success;
}

void testIt(long d, long k, long p, long r, long m, long L,
	    bool isMonic=false)
{
  FHEcontext context(m, p, r);
  long p2r = context.alMod.getPPowR();
  buildModChain(context, L, /*c=*/3);
  EncryptedArray ea(context);

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(/*w=*/64);// A Hamming-weight-64 secret key
  //  addSome1DMatrices(secretKey); // compute key-switching matrices

#ifdef DEBUG_PRINTOUT
  dbgEa = &ea;        // for debugging purposes
  dbgKey = &secretKey;
#endif

  cout << (isDryRun()? "* dry run, " : "* ")
       << "degree-"<<d<<", m="<<m<<", L="<<L<<", p^r="<<p2r<<endl;

  // evaluate encrypted poly at encrypted point
  if (!testEncrypted(d, ea, secretKey)) exit(0);

  // evaluate at random points (at least one co-prime with p)
  vector<long> x;
  ea.random(x);
  while (GCD(x[0],p)!=1) { x[0] = RandomBnd(p2r); }
  Ctxt inCtxt(publicKey), outCtxt(publicKey);
  ea.encrypt(inCtxt, publicKey, x);

  ZZX poly;
  for (long i=d; i>=0; i--)
    SetCoeff(poly, i, RandomBnd(p2r)); // coefficients are random
  if (isMonic) SetCoeff(poly, d);    // set top coefficient to 1

  // Evaluate poly on the ciphertext
  polyEval(outCtxt, poly, inCtxt, k);

  // Check the result
  vector<long> y;
  ea.decrypt(outCtxt, secretKey, y);
  for (long i=0; i<ea.size(); i++) {
    long ret = polyEvalMod(poly, x[i], p2r);
    if (ret != y[i]) {
      cout << "plaintext poly MISMATCH\n";
      exit(0);
    }
  }
  cout << "plaintext poly match\n" << std::flush;
}

void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cout << "  dry=1 for dry run [default=0]\n";
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
  argmap["dry"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long m = atoi(argmap["m"]);
  long d = atoi(argmap["d"]);
  long k = atoi(argmap["k"]);
  bool dry = atoi(argmap["dry"]);

  long max_d = (d<=0)? 35 : d;
  long L = 5+NextPowerOfTwo(max_d);
  if (m<2)
    m = FindM(/*secprm=*/80, L, /*c=*/3, p, 1, 0, m, true);
  setDryRun(dry);

  // Test both monic and non-monic polynomials of this degree
  if (d>=0) {
    testIt(d, k, p, r, m, L, false);
    testIt(d, k, p, r, m, L, true);
    return 0;
  }

  // Test degrees 1 to 4 and 25 through 35
  testIt(1, k, p, r, m, L, true);
  testIt(3, k, p, r, m, L, true);
  for (d=25; d<=33; d+=2)
    testIt(d, k, p, r, m, L);

  return 0;
}
