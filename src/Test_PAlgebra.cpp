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
#include <string>
#include <sstream>
#include <NTL/ZZ.h>
NTL_CLIENT
#include "NumbTh.h"
#include "FHEContext.h"

void usage() 
{
  cout << "Usage: Test_PAlgebra_x m=<int> [ p=<int> ] [ r=<int> ]" << endl;
  cout << "  m is an integer that determines Z_m^*" << endl;
  cout << "  p is an integer that determines the plaintext base [default=2]" << endl;
  cout << "  r is an integer that determines the lifting [default=1]" << endl;
  exit(0);
}

int main(int argc, char *argv[]) 
{

  argmap_t argmap;
  argmap["m"] = "0"; 
  argmap["p"] = "2";
  argmap["r"] = "1";

  if (!parseArgs(argc, argv, argmap)) usage();

  unsigned long m = atoi(argmap["m"]);

  if (!m) usage();

  unsigned long p = atoi(argmap["p"]);

  unsigned long r = atoi(argmap["r"]);

  cout << "m = " << m << ", p = " << p <<  ", r = " << r << endl;

  vector<long> f;
  factorize(f,m);
  cout << "factoring "<<m<<" gives [";
  for (unsigned long i=0; i<f.size(); i++)
    cout << f[i] << " ";
  cout << "]\n";

  PAlgebra al(m, p);
  al.printout();
  cout << "\n";

  PAlgebraMod almod(al, r);

  FHEcontext context(m, p, r);
  buildModChain(context, 5, 2);

  stringstream s1;
  writeContextBase(s1, context);
  s1 << context;

  string s2 = s1.str();

  cout << s2 << endl;

  stringstream s3(s2);

  unsigned long m1, p1, r1;
  vector<long> gens, ords;
  readContextBase(s3, m1, p1, r1, gens, ords);

  FHEcontext c1(m1, p1, r1, gens, ords);
  s3 >> c1;

  if (context == c1)
    cout << "equal\n";
  else
    cout << "not equal\n";

  return 0;

#if 0


  PAlgebraModTwo al2(al);
  PAlgebraMod2r al2r(r,al2);

  if (false) {
    long nslots = al.NSlots();
    long ngens = al.numOfGens();

    for (long i = 0; i < nslots; i++) {
      const int* v = al.dLog(al.ith_rep(i));
      cout << i << " ";
      for (long j = 0; j < ngens; j++)
         cout << v[j] << " ";
      cout << "\n";
      
    }

    return 0;
  }

  //  GF2X::HexOutput = 1; // compact hexadecimal printout
  //  cout << "Phi_m(X) = " << al.PhimX()  << endl;
  //  vec_GF2X Fs = al2.Factors();
  //  cout << "Factors = " << Fs << endl;

  //  GF2X Q = Fs[0];
  //  for (long i=1; i<Fs.length(); i++) Q *= Fs[i];
  //  cout << "mult of factors = " << Q << endl;

  GF2X G;          // G is the AES polynomial, G(X)= X^8 +X^4 +X^3 +X +1
  SetCoeff(G,8); SetCoeff(G,4); SetCoeff(G,3); SetCoeff(G,1); SetCoeff(G,0);
  //  GF2X G; SetCoeff(G,4); SetCoeff(G,1); SetCoeff(G,0); // X^4 +X +1
  //  cout << "G = " << G << ", deg(G)=" << deg(G) << endl;

  vector<GF2X> maps;
  al2.mapToSlots(maps, G);
  //  cout << "maps = [";
  //  for (unsigned long i=0; i<maps.size(); i++)
  //    cout << maps[i] << " ";
  //  cout << "]\n";


  GF2X X; SetX(X); // The polynomial X
  GF2X ptxt;       // plaintext has X in all the slots
  cout << "embedding X in plaintext slots... ";
  al2.embedInAllSlots(ptxt, X, maps);
  cout << "done\n";
  //  cout << "ptxt = " << ptxt << endl;

  // Debugging printout: p modulo all the factors
  //  vector<GF2X> crt;
  //  al2.CRT_decompose(crt,ptxt);
  //  cout << "ptxt mod factors = [";
  //  for (unsigned long i=0; i<crt.size(); i++) cout << crt[i] << " ";
  //  cout << "]\n";

  // Decode the plaintext back to a vector of elements,
  // and check that they are all equal to X
  vector<GF2X> alphas;
  cout << "decoding plaintext slots... ";
  al2.decodePlaintext(alphas, ptxt, G, maps);
  cout << "done\n";
  //  cout << "alphas = [";
  //  for (unsigned long i=0; i<alphas.size(); i++)
  //    cout << alphas[i] << " ";
  //  cout << "]\n";

  cout << "comparing " << alphas.size() << " plaintext slots to X... ";
  for (unsigned long i=0; i<alphas.size(); i++) if (alphas[i] != X) {
      cout << "\n  alphas["<<i<<"] = "<<alphas[i]<<" != X\n\n";
      exit(0);
  }
  cout << "all tests completed successfully\n\n";

  // encode and decode random polynomials
  for (unsigned long i=0; i<alphas.size(); i++) 
    random(alphas[i], 8); // random degree-7 polynomial mod 2

  cout << "embedding random GF(2^8) elements in plaintext slots... ";
  al2.embedInSlots(ptxt, alphas, maps);
  cout << "done\n";

  // Compute p^2 mod Phi_m(X) and also p(X^2) mod Phi_m(X) and
  // verify that they both decode to a vector of X^2 in all the slots

  cout << "squaring and decoding plaintext slots... ";

  X *= X;                            // X^2
  //  GF2X ptxt2;
  //  SqrMod(ptxt2,ptxt,al2.PhimXMod());      // ptxt2 = ptxt^2 mod Phi_m(X)
  CompMod(ptxt, ptxt, X, al2.PhimXMod()); // ptxt = ptxt(X^2) mod Phi_m(X)

  //  // sanity chack: these should be the same mod 2 (but not mod 2^r)
  //  if (ptxt != ptxt2) cout << "ptxt^2 != ptxt(X^2) mod Phi_m(X)\n";

  vector<GF2X> betas;
  al2.decodePlaintext(betas, ptxt, G, maps);
  cout << "done\n";

  if (alphas.size() != betas.size()) Error("wrong number of slots decoded");

  cout << "comparing decoded plaintext slots... ";
  for (unsigned long i=0; i<alphas.size(); i++) {
    SqrMod(alphas[i],alphas[i],G); // should get alpha[i]^2 mod G
    if (alphas[i] != betas[i]) {
      cout << "\n  alphas["<<i<<"] = "<<alphas[i]
           <<" != " << "betas["<<i<<"] = " << betas[i] << "\n\n";
      exit(0);
    }
  }
  cout << "all tests completed successfully\n\n";
  // return 0;

  al2r.restoreContext();

  vector<zz_pX> maps1;
  zz_pX X1; SetX(X1);

  cerr << "HERE1\n";
  al2r.mapToSlots(maps1, X1);
  cerr << "HERE1a\n";

  vector<zz_pX> alphas1;
  alphas1.resize(maps.size());

  for (long i = 0; i < lsize(alphas1); i++)
    random(alphas1[i], 1);

  zz_pX ptxt1;

  cerr << "HERE2\n";
  al2r.embedInSlots(ptxt1, alphas1, maps1);

  cerr << "HERE3\n";
  vector<zz_pX> betas1;
  al2r.decodePlaintext(betas1, ptxt1, X1, maps1);

  assert(alphas1 == betas1);
  

  return 0;
#endif
}
