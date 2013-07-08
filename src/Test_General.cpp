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

/**************

1. c1.multiplyBy(c0)
2. c0 += random constant
3. c2 *= random constant
4. tmp = c1
5. ea.shift(tmp, random amount in [-nSlots/2, nSlots/2])
6. c2 += tmp
7. ea.rotate(c2, random amount in [1-nSlots, nSlots-1])
8. c1.negate()
9. c3.multiplyBy(c2) 
10. c0 -= c3

**************/


ZZX makeIrredPoly(long p, long d)
{
  assert(d >= 1);
  assert(ProbPrime(p));

  if (d == 1) return ZZX(1, 1); // the monomial X

  zz_pBak bak; bak.save();
  zz_p::init(p);
  return to_ZZX(BuildIrred_zz_pX(d));
}


void  TestIt(long R, long p, long r, long d, long c, long k, long w, 
               long L, long m)
{
  cerr << "*** TestIt: R=" << R 
       << ", p=" << p
       << ", r=" << r
       << ", d=" << d
       << ", c=" << c
       << ", k=" << k
       << ", w=" << w
       << ", L=" << L
       << ", m=" << m
       << endl;

  FHEcontext context(m, p, r);
  buildModChain(context, L, c);

  context.zMStar.printout();
  cerr << endl;

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(w); // A Hamming-weight-w secret key


  ZZX G;

  if (d == 0)
    G = context.alMod.getFactorsOverZZ()[0];
  else
    G = makeIrredPoly(p, d); 

  cerr << "G = " << G << "\n";
  cerr << "generating key-switching matrices... ";
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  cerr << "done\n";

  printAllTimers();
  resetAllTimers();

  cerr << "computing masks and tables for rotation...";
  EncryptedArray ea(context, G);
  cerr << "done\n";


  double t = GetTime();

  long nslots = ea.size();

  PlaintextArray p0(ea);
  PlaintextArray p1(ea);
  PlaintextArray p2(ea);
  PlaintextArray p3(ea);

  p0.random();
  p1.random();
  p2.random();
  p3.random();

  Ctxt c0(publicKey), c1(publicKey), c2(publicKey), c3(publicKey);
  ea.encrypt(c0, publicKey, p0);
  ea.encrypt(c1, publicKey, p1);
  ea.encrypt(c2, publicKey, p2);
  ea.encrypt(c3, publicKey, p3);

  for (long i = 0; i < R; i++) {

    cerr << "*** round " << i << "..."<<endl;

     long shamt = RandomBnd(2*(nslots/2) + 1) - (nslots/2);
                  // random number in [-nslots/2..nslots/2]
     long rotamt = RandomBnd(2*nslots - 1) - (nslots - 1);
                  // random number in [-(nslots-1)..nslots-1]

     // two random constants
     PlaintextArray const1(ea);
     PlaintextArray const2(ea);
     const1.random();
     const2.random();

     p1.mul(p0); // c1.multiplyBy(c0)
     p0.add(const1); // c0 += random constant
     p2.mul(const2); // c2 *= random constant
     PlaintextArray tmp_p(p1); // tmp = c1
     tmp_p.shift(shamt); // ea.shift(tmp, random amount in [-nSlots/2, nSlots/2])
     p2.add(tmp_p); // c2 += tmp
     p2.rotate(rotamt); // ea.rotate(c2, random amount in [1-nSlots, nSlots-1])
     p1.negate(); // c1.negate()
     p3.mul(p2); // c3.multiplyBy(c2) 
     p0.sub(p3); // c0 -= c3

     ZZX const1_poly, const2_poly;
     ea.encode(const1_poly, const1);
     ea.encode(const2_poly, const2);

     c1.multiplyBy(c0);              CheckCtxt(c1, "c1*=c0");
     c0.addConstant(const1_poly);    CheckCtxt(c0, "c0+=k0");
     c2.multByConstant(const2_poly); CheckCtxt(c2, "c2*=k0");
     Ctxt tmp(c1);
     ea.shift(tmp, shamt);           CheckCtxt(tmp, "tmp=c1>>$");
     c2 += tmp;                      CheckCtxt(c2, "c2+=tmp");
     ea.rotate(c2, rotamt);          CheckCtxt(c2, "c2>>>=$");
     c1.negate();                    CheckCtxt(c1, "c1=-c1");
     c3.multiplyBy(c2);              CheckCtxt(c3, "c3*=c2");
     c0 -= c3;                       CheckCtxt(c0, "c0=-c3");

     PlaintextArray pp0(ea);
     PlaintextArray pp1(ea);
     PlaintextArray pp2(ea);
     PlaintextArray pp3(ea);

     ea.decrypt(c0, secretKey, pp0);
     ea.decrypt(c1, secretKey, pp1);
     ea.decrypt(c2, secretKey, pp2);
     ea.decrypt(c3, secretKey, pp3);

     cerr << "\n";

     if (!pp0.equals(p0)) cerr << "oops 0\n";
     if (!pp1.equals(p1)) cerr << "oops 1\n";
     if (!pp2.equals(p2)) cerr << "oops 2\n";
     if (!pp3.equals(p3)) cerr << "oops 3\n";
  }

  t = GetTime() - t;

  cerr << "time for circuit: " << t << "\n";
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
    if (r==1) L = 2*R+2;
    else      L = 4*R;
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

