/* Copyright (C) 2012-2017 IBM Corp.
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

static bool noPrint = false;

bool testEncrypted(long d, const EncryptedArray& ea,
		   const FHESecKey& secretKey)
{
  const FHEcontext& context = ea.getContext();
  const FHEPubKey& publicKey = secretKey;
  long p = publicKey.getPtxtSpace();
  zz_pBak bak; bak.save(); zz_p::init(p);
  zz_pXModulus phimX = conv<zz_pX>(ea.getPAlgebra().getPhimX());

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
  if (success) std::cout << " encrypted poly match, ";
  else         std::cout << " encrypted poly MISMATCH\n";
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

  if (!noPrint) std::cout << (isDryRun()? "* dry run, " : "* ")
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
      std::cout << "plaintext poly MISMATCH\n";
      exit(0);
    }
  }
  std::cout << "plaintext poly match\n" << std::flush;
}

void usage(char *prog) 
{
  std::cout << "Usage: "<<prog<<" [ optional parameters ]...\n";
  std::cout << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  std::cout << "  dry=1 for dry run [default=0]\n";
  std::cout << "  p is the plaintext base [default=3]" << endl;
  std::cout << "  r is the lifting [default=2]" << endl;
  std::cout << "  m is a specific cyclotomic ring\n";
  std::cout << "  d is the polynomial degree [default=undefined]" << endl;
  std::cout << "    d=undefined means trying a few powers d=1,...,4,25,...,34"<<endl;
  std::cout << "  k is the baby-step parameter [default=undefined]" << endl;
  std::cout << "    if k is undefined it is computed from d" << endl;
  std::cout << "  noPrint suppresses printouts [default=0]" << endl;
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
  argmap["noPrint"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long m = atoi(argmap["m"]);
  long d = atoi(argmap["d"]);
  long k = atoi(argmap["k"]);
  bool dry = atoi(argmap["dry"]);
  noPrint = atoi(argmap["noPrint"]);

  long max_d = (d<=0)? 35 : d;
  long L = 5+NextPowerOfTwo(max_d);
  if (m<2)
    m = FindM(/*secprm=*/80, L, /*c=*/3, p, 1, 0, m, !noPrint);
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
