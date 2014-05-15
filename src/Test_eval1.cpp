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
#include "powerful.h"
#include <NTL/lzz_pXFactoring.h>

#include <cassert>

void init_representatives(Vec<long>& representatives, long m, long p)
{
  Vec<bool> available;
  available.SetLength(m);

  long num_available = 0;

  for (long i = 0; i < m; i++) {
    if (GCD(i, m) == 1) {
      available[i] = true;
      num_available++;
    }
    else
      available[i] = false;
  }

  representatives.SetLength(0);

  while (num_available > 0) {

    // choose next available at random
    long i;
    do {
      i = RandomBnd(m);
    } while (!available[i]);

    append(representatives, i);

    // mark all conjugates as unavailable
    long j = i;
    do {
      available[j] = false;
      num_available--;
      j = MulMod(j, p, m);
    } while (j != i);
  }
}

void init_slot_mappings(Vec<long>& slot_index, 
                        Vec<long>& slot_rotate, 
                        const Vec<long>& representatives, 
                        long m,
                        long p,
                        const FHEcontext& context)
{
   long nslots = representatives.length();

   assert(nslots == long(context.zMStar.getNSlots()));

   slot_index.SetLength(nslots);
   slot_rotate.SetLength(nslots);

   Vec<bool> used; // for debugging
   used.SetLength(nslots);
   for (long i = 0; i < nslots; i++) used[i] = false;
   
   for (long i = 0; i < nslots; i++) {
     long t = representatives[i];
     long h = 0;
     long idx;
     while ((idx = context.zMStar.indexOfRep(InvMod(t, m))) == -1) {
       t = MulMod(t, p, m);
       h++;
     }

     assert(!used[idx]);
     used[idx] = true;
     slot_index[i] = idx;
     slot_rotate[i] = h;
   }
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

  assert(GCD(p, m) == 1 && multOrd(p, m1) == multOrd(p, m));

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

  cerr << "GG = " << GG << "\n";
  cerr << "generating key-switching matrices... ";
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need
  cerr << "done\n";

  printAllTimers();
  resetAllTimers();



  cerr << "computing masks and tables for rotation...";
  EncryptedArray ea(context, GG);
  cerr << "done\n";

  Vec<long> rep1;
  init_representatives(rep1, m1, p);

  Vec<long> rep2;
  init_representatives(rep2, m2, 1);
  

  Vec<long> representatives;
  for (long i = 0; i < rep1.length(); i++)
    for (long j = 0; j < rep2.length(); j++) {
      // chinese remaindering
      long x1 = rep1[i];
      long x2 = rep2[j];

      long x = mcMod(x1*m2*InvMod(m2, m1) + x2*m1*InvMod(m1, m2), m);

      append(representatives, x);
    }


  cout << representatives << "\n";

  Vec<long> slot_index, slot_rotate;

  init_slot_mappings(slot_index, slot_rotate, representatives, m, p, context);

  cout << slot_index << "\n";
  cout << slot_rotate << "\n";

  zz_p::init(p);

  zz_pX G = conv<zz_pX>(GG);

  Vec<zz_pX> eval_points;
  eval_points.SetLength(nslots);
  for (long i = 0; i < nslots; i++) 
    eval_points[i] = zz_pX(representatives[i], 1) % G; 


  zz_pX F;
  random(F, phim);

  Vec<zz_pX> eval_values;
  eval_values.SetLength(nslots);

  for (long i = 0; i < nslots; i++)
    eval_values[i] = CompMod(F, eval_points[i], G);



  // Crazy code from powerful module

  Vec<long> phiVec;
  phiVec.SetLength(2);
  phiVec[0] = phim1;
  phiVec[1] = phim2;

  Vec<long> powVec;
  powVec.SetLength(2);
  powVec[0] = m1;
  powVec[1] = m2;

  Vec<long> divVec;
  computeDivVec(divVec, m, powVec);

  Vec<long> invVec;
  computeInvVec(invVec, divVec, powVec);

  CubeSignature shortSig(phiVec);
  CubeSignature longSig(powVec);

  Vec<long> polyToCubeMap;
  Vec<long> cubeToPolyMap;
  computePowerToCubeMap(polyToCubeMap, cubeToPolyMap, m, powVec, invVec, longSig);

  Vec<long> shortToLongMap;
  computeShortToLongMap(shortToLongMap, shortSig, longSig);

  Vec<long> longToShortMap;
  computeLongToShortMap(longToShortMap, m, shortToLongMap);
  

  Vec<zz_pX> cycVec;
  computeCycVec(cycVec, powVec);


  ZZX PhimX = Cyclotomic(m);
  zz_pX phimX = conv<zz_pX>(PhimX);

  HyperCube<zz_p> cube(shortSig);
  HyperCube<zz_p> tmpCube(longSig);

  convertPolyToPowerful(cube, tmpCube, F, cycVec, 
                        polyToCubeMap, shortToLongMap);

  zz_pX poly1;

  convertPowerfulToPoly(poly1, cube, m, shortToLongMap, cubeToPolyMap, phimX);

  if (F == poly1)
    cout << "*********** :-)\n";
  else {
    cout << "*********** :-(\n";
    cout << F << "\n";
    cout << poly1 << "\n";
  }

  // end crazy code
  //

  // evaluate F via the cube structure
  
  Vec<zz_pX> alt_values;
  
  for (long i = 0; i < rep1.length(); i++) {
    for (long j = 0; j < rep2.length(); j++) {
      cerr << ".";
      long x1 = rep1[i];
      long x2 = rep2[j];

      zz_pX res;
      for (long l = 0; l < phim; l++) {
        long s1 = cube.getCoord(l, 0);
        long s2 = cube.getCoord(l, 1);
        long s = mcMod(m2*s1*x1 + m1*s2*x2, m);

        res += zz_pX(s, cube[l]) % G; 
      }

      append(alt_values, res);
    }
  }

  if (eval_values == alt_values) 
    cout << "right on!!\n";
  else
    cout << "no way!!!\n";


  zz_pX H = zz_pX(p, 1) % G;

  vector<ZZX> adjusted_values;
  adjusted_values.resize(nslots);

  for (long i = 0; i < nslots; i++) {
    zz_pX V = eval_values[i];
    long h = slot_rotate[i];
    for (long j = 0; j < h; j++) 
      V = CompMod(V, H, G);
    
    adjusted_values[ slot_index[i] ] = conv<ZZX>(V);
  }

  ZZX FF1;
  ea.encode(FF1, adjusted_values);
  
  zz_pX F1 = conv<zz_pX>(FF1);

  if (F1 == F) 
    cout << "yes!!\n";
  else {
    cout << "NO!!!\n";

    ZZX FF = conv<ZZX>(F);
    vector<ZZX> decode_vec;
    ea.decode(decode_vec, FF);

    for (long i = 0; i < nslots; i++) {
      if (adjusted_values[i] != decode_vec[i]) {
        cout << "*** " << i << " ";
        zz_pX U, V;
        U = conv<zz_pX>(adjusted_values[i]);

        for (long j = 0; j < d; j++) {
          for (long l = 0; l < nslots; l++) 
            if (U == conv<zz_pX>(decode_vec[l])) 
              cout << "fixed " << j << " " << (l-i);

          U = CompMod(U, H, G);
        }
        cout << "\n";
      }
    }
  }

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

// working test case: m1=89 m2=23
