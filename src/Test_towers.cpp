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

namespace std {} using namespace std;
namespace NTL {} using namespace NTL;

#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include "powerful.h"
#include <NTL/lzz_pXFactoring.h>

#include <cassert>


void buildLinPolyMatrix(mat_zz_pE& M, long p, long d1)
{
   long d = zz_pE::degree();
   assert(d % d1 == 0);
   long d2 = d/d1;
   long q = power_long(p, d1);

   M.SetDims(d2, d2);

   for (long j = 0; j < d2; j++) 
      conv(M[0][j], zz_pX(j, 1));

   for (long i = 1; i < d2; i++)
      for (long j = 0; j < d2; j++)
         M[i][j] = power(M[i-1][j], q);
}



void buildLinPolyCoeffs(vec_zz_pE& C_out, const vec_zz_pE& L, long p, long r, long d1)
{
   mat_zz_pE M;
   buildLinPolyMatrix(M, p, d1);

   vec_zz_pE C;
   ppsolve(C, M, L, p, r);

   C_out = C;
}

void applyLinPoly(zz_pE& beta, const vec_zz_pE& C, const zz_pE& alpha, long p, long d1)
{
   long d = zz_pE::degree();
   assert(d % d1 == 0);
   long d2 = d/d1;
   assert(d2 == C.length());
   long q = power_long(p, d1);

   zz_pE gamma, res;

   gamma = to_zz_pE(zz_pX(1, 1));
   res = C[0]*alpha;
   for (long i = 1; i < d2; i++) {
      gamma = power(gamma, q);
      res += C[i]*to_zz_pE(CompMod(rep(alpha), rep(gamma), zz_pE::modulus()));
   }

   beta = res;
}

zz_pE convert2to1(const Mat<zz_p>& M2, const Mat<zz_p>& M2i, long d1,
                  const Vec<zz_pX>& v)
{
  long d = M2.NumRows();
  assert(d % d1 == 0);
  long d2 = d/d1;

  assert(d == zz_pE::degree());
  assert(v.length() <= d2);

  Vec<zz_p> w;
  w.SetLength(d);

  for (long i = 0; i < v.length(); i++) {
    assert(deg(v[i]) < d1);
    for (long j = 0; j <= deg(v[i]); j++) 
      w[i*d1 + j] = v[i][j];
  }

  Vec<zz_p> z = w*M2;
  return conv<zz_pE>( conv<zz_pX>( z ) );
}

Vec<zz_pX> convert1to2(const Mat<zz_p>& M2, const Mat<zz_p>& M2i, long d1,
                       const zz_pE& beta)
{
  long d = M2.NumRows();
  assert(d % d1 == 0);
  long d2 = d/d1;
  assert(d == zz_pE::degree());

  Vec<zz_p> z = VectorCopy(rep(beta), d);
  Vec<zz_p> w = z*M2i;

  Vec<zz_pX> res;
  res.SetLength(d2);

  for (long i = 0; i < d2; i++) {
    Vec<zz_p> tmp;
    tmp.SetLength(d1);
    for (long j = 0; j < d1; j++)
      tmp[j] = w[i*d1+j];

    res[i] = conv<zz_pX>(tmp);
  }

  return res;
}


void  TestIt(long R, long p, long r, long d, long c, long k, long w, 
               long L, long m1, long m2)
{
  cerr << "*** TestIt: R=" << R 
       << ", p=" << p
       << ", r=" << r
       << ", d=" << d
       << ", c=" << c
       << ", k=" << k
       << ", w=" << w
       << ", L=" << L
       << ", m1=" << m1
       << ", m2=" << m2
       << endl;

  long m = m1*m2;

  assert(GCD(p, m) == 1);


  FHEcontext context(m, p, r);

  long phim = context.zMStar.getPhiM();
  long phim1 = phi_N(m1);
  long phim2 = phi_N(m2);

  long nslots = context.zMStar.getNSlots();

  buildModChain(context, L, c);

  context.zMStar.printout();
  cerr << endl;

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(w); // A Hamming-weight-w secret key

  assert(d == 0); // that's the assumption here

  ZZX GG;

  if (d == 0)
    GG = context.alMod.getFactorsOverZZ()[0];
  else
    GG = makeIrredPoly(p, d); 

  d = deg(GG);
  long d1 = multOrd(p, m1);
  long d2 = d/d1;

  cout << "d1=" << d1 << ", d2=" << d2 << ", d=" << d << "\n"; 


  zz_p::init(power_long(p, r));
  zz_pX G = conv<zz_pX>(GG);

  zz_pX alpha = zz_pX(m2, 1) % G; // X^m2 has order m1, and should
                                   // have min poly of degree d1

  // compute min poly of alpha by linear algebra 

  Mat<zz_p> M1;
  M1.SetDims(d1, d);

  for (long i = 0; i < d1; i++) {
    VectorCopy(M1[i], PowerMod(alpha, i, G), d);
  }

  Vec<zz_p> V1;
  VectorCopy(V1, PowerMod(alpha, d1, G), d);

  Mat<zz_p> M1sq;

  Mat<zz_p> R1;
  R1.SetDims(d, d1);

  for (;;) {
    cout << ".";

    for (long i = 0; i < d; i++)
      for (long j = 0; j < d1; j++)
        random(R1[i][j]);

    M1sq = M1*R1;
  
    Mat<long> M1sqInt = conv< Mat<long> >(M1sq);
    {
       zz_pBak bak; bak.save();
       zz_p::init(p);
       Mat<zz_p> M1sq_modp = conv< Mat<zz_p> >(M1sqInt);
       if (determinant(M1sq_modp) != 0) break;
    }
 
  }

  cout << "\n";

  Vec<zz_p> V1sq = V1*R1;

  Mat<zz_p> M1sqi;
  ppInvert(M1sqi, M1sq, p, r);

  Vec<zz_p> W1 = V1sq * M1sqi;

  assert(W1*M1 == V1);

  zz_pX H = zz_pX(d1, 1) - conv<zz_pX>(W1);
  // H is the min poly of alpha

  assert(CompMod(H, alpha, G) == 0);

  // now construct the min poly of X over ZZ_{p^r}[alpha]

  Mat<zz_p> M2;
  M2.SetDims(d, d);

  for (long i = 0; i < d2; i++) {
    // construct rows [i..i+d1)
    for (long j = 0; j < d1; j++) {
       VectorCopy(M2[i*d1+j], (PowerMod(alpha, j, G) << i) % G, d);
    }
  }

  Mat<zz_p> M2i;
  ppInvert(M2i, M2, p, r);

  // M2 is the matrix that takes us from the two-step tower
  // to the one-step tower, and M2i is its inverse.

  Vec<zz_p> V2, W2;
  VectorCopy(V2, zz_pX(d, 1) % G, d);

  W2 = V2 * M2i;

  assert(W2 * M2 == V2);



  zz_pE::init(G);

  Vec<zz_pE> C_out, map1;

  Vec< Vec<zz_pX> > map2;
  map2.SetLength(d2);
  long idx_in = 1;
  long idx_out = d2-1;
  map2[idx_in].SetLength(idx_out+1);
  map2[idx_in][idx_out] = 1;
  // so map2 projects idx_in onto idx_out

  map1.SetLength(d2);
  for (long i = 0; i < d2; i++) 
    map1[i] = convert2to1(M2, M2i, d1, map2[i]); 

  buildLinPolyCoeffs(C_out, map1, p, r, d1);

  zz_pE testval1;
  random(testval1);

  Vec<zz_pX> testval2 = convert1to2(M2, M2i, d1, testval1);
  zz_pE testval1a = convert2to1(M2, M2i, d1, testval2);

  cout << testval1 << "\n";
  cout << testval1a << "\n";

  cout << testval2 << "\n";

  zz_pE resval1;
  applyLinPoly(resval1, C_out, testval1, p, d1);
  Vec<zz_pX> resval2 = convert1to2(M2, M2i, d1, resval1);

  cout << resval2;

}




void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'R=1 p=2 k=80'\n\n";
  cerr << "  R is the number of rounds\n";
  cerr << "  p is the plaintext base [default=2]" << endl;
  cerr << "  r is the lifting [default=1]" << endl;
  cerr << "  d is the degree of the field extension [default==0]\n";
  cerr << "    (d == 0 => factors[0] defined the extension)\n";
  cerr << "  c is number of columns in the key-switching matrices [default=2]\n";
  cerr << "  k is the security parameter [default=80]\n";
  cerr << "  L is the # of primes in the modulus chai [default=4*R]\n";
  cerr << "  s is the minimum number of slots [default=4]\n";
  cerr << "  m defined the cyclotomic polynomial Phi_m(X)\n";
  cerr << "  seed is the PRG seed\n";
  exit(0);
}


int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["R"] = "1";
  argmap["p"] = "2";
  argmap["r"] = "1";
  argmap["d"] = "0";
  argmap["c"] = "2";
  argmap["k"] = "80";
  argmap["L"] = "0";
  argmap["s"] = "0";
  argmap["m1"] = "0";
  argmap["m2"] = "0";
  argmap["seed"] = "0";

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
    if (r==1) L = 2*R+2;
    else      L = 4*R;
  }
  long s = atoi(argmap["s"]);
  long m1 = atoi(argmap["m1"]);
  long m2 = atoi(argmap["m2"]);
  long seed = atoi(argmap["seed"]);

  long w = 64; // Hamming weight of secret key
  //  long L = z*R; // number of levels
  //

  assert(m1 != 0 && m2 != 0);

  if (seed) SetSeed(conv<ZZ>(seed));

  TestIt(R, p, r, d, c, k, w, L, m1, m2);

  cerr << endl;
  printAllTimers();
  cerr << endl;

}

// m1=5 m2=17   -- d1=4, d2=2, d=8
// m1=7 m2=13   -- d1=3, d2=4, d=12
