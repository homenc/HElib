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
#include <iostream>
#include <NTL/ZZX.h>
NTL_CLIENT

#include "NumbTh.h"
#include "DoubleCRT.h"
#include "SingleCRT.h"

ZZX mod(const ZZX& f, const ZZ& n) 
{
  ZZX t;
  t = f;
  PolyRed(t, n);
  return t;
}

int main(int argc, char *argv[]) 
{
  if (argc<2) {
    cout << "\nUsage: " << argv[0] << " m [n]\n\n";
    cout << "  m is an integer that determines Z_m^*\n";
    cout << "  optional n is the number of prime in the chain (in [3,100], default=5)\n\n";
    exit(0);
  }
  long m = atoi(argv[1]);
  if (m > (1L << 20)) Error("Parameter m cannot exceed 2^20");

  long n = 5;
  if (argc>2) {
    n = atoi(argv[2]);
    if (n<3 || n>100) n=5;
  }
  FHEcontext context(m);
  long phim = context.zMstar.phiM();
  double totalsize = AddPrimesByNumber(context,5); // five primes in the chain
  ZZ product = context.productOfPrimes(context.ctxtPrimes);

  cout << "Generated moduli-chain of " << context.numPrimes()
       << " primes, with log(product)=" << totalsize<< ".\n";
  for (long i=0; i<context.numPrimes(); i++) {
    cout << "  prime #" << i << " " << context.ithPrime(i) << endl;
  }
  cout << "Product of primes=" << product << endl;

  // test conversions between ZZX, SingleCRT, and DoubleCRT
  ZZX p1,p2,p3;
  SingleCRT s1(context),s2(context),s3(context);
  DoubleCRT d1(context),d2(context),d3(context);

  // A simple sanity-check, using tyhe constant polynomial f(X)=12
  d1 = 12;
  d1.toPoly(p1);
  cout << "poly(DoubleCRT(12))=" << p1 << endl;

  // choose p1 at random with degree phi(m)-1 modulo the product
  ZZ_p::init(product);
  ZZ_pContext ctxt; ctxt.save();
  {ZZ_pX tmp; random(tmp, phim); // a random polynomial mod p1
  conv(p1, tmp);}                // convert to ZZX format
  PolyRed(p1,product);// make the coefficients of p1 in [-prod/2,+prod/2)

  //  cout << " p1=" << p1 << endl;

  s1 = p1; // convert ZZX to SingleCRT
  d1 = p1; // convert ZZX to DoubleCRT
  d2 = s1; // conver SingleCRT to DoubleCRT
  cout << "\nTest  1 " << ((d1==d2)? "PASS" : "FAIL")
       << ": poly=>singleCRT=>DoubleCRT =?= poly=>DoubleCRT\n";

  d1.toPoly(p2); // Convert DoubleCRT to ZZX
  //  cout << " p2=" << p2 << endl;
  cout << "Test  2 " << ((p1==p2)? "PASS" : "FAIL")
       << ": poly=>DoubleCRT=>poly =?= poly\n";

  s1.toPoly(p2); // Convert SingleCRT to ZZX
  cout << "Test  3 " << ((p1==p2)? "PASS" : "FAIL")
       << ": poly=>SingleCRT=>poly =?= poly\n";

  s2 = d2; // Convert DoubleCRT to SingleCRT
  cout << "Test  4 " << ((s1==s2)? "PASS" : "FAIL")
       << ": SingleCRT=>DoubleCRT=>SingleCRT =?= SingleCRT\n";

 // check that level manipulation does not change a small polynomial

#if 0
  d1.sampleSmall(); // make it a polynomial with coefficients -1,0,1
  d1.toPoly(p1);
  d1.setNlevels(2);
  d1.toPoly(p2);
  d1.setNlevels(n);
  d1.toPoly(p3);
  cout << "\nTest  5 " << ((p1==p2)? "PASS" : "FAIL")
       << ": SmallDoubleCRT.reduceLevels=>poly =?= SmallDoubleCRT=>poly\n";
  cout << "Test  6 " << ((p2==p3)? "PASS" : "FAIL")
       << ": SmallDoubleCRT.increaseLevels=>poly =?= SmallDoubleCRT=>poly\n";

  d1.toPoly(p1, /*fromLvl=*/1, /*numLvls=*/n-1);
  cout << "Test  7 " << ((p1==p3)? "PASS" : "FAIL")
       << ": SmallDoubleCRT=>poly(top-levels) =?= SmallDoubleCRT=>poly\n";

  d1.setNlevels(2);
  d1.randomize();
  d1.toPoly(p1);
  d1.setNlevels(n);
  d1.toPoly(p2);
  cout << "Test  8 " << ((p1==p2)? "PASS" : "FAIL")
       << ": DoubleCRT.increaseLevels=>poly =?= DoubleCRT=>poly\n";
#endif

  // test DoubleCRT arithmetic

  // choose p2 at random with degree phi(m)-1 and coefficients in {-1,0,1}
  sampleSmall(p1, phim);
  sampleSmall(p2, phim);

  d1 = p1;
  d2 = p2;

  // Addition/subtraction

  p2 = mod(p2 + p1, product);
  d2 += d1;


  cout << "\nTest  9 " << ((to_ZZX(d2)==p2)? "PASS" : "FAIL")
       << ": poly(d2+d1) =?= poly(d2)+poly(d1)\n";

  p2 = mod(p2 - p1, product);
  d2 -= d1;

  cout << "Test 10 " << ((to_ZZX(d2)==p2)? "PASS" : "FAIL")
       << ": poly(d2-d1) =?= poly(d2)-poly(d1)\n";


  p2 = mod(p2 - 3, product);
  d2 -= 3;

  cout << "Test 11 " << ((to_ZZX(d2)==p2)? "PASS" : "FAIL")
       << ": poly(d2-const) =?= poly(d2)-const\n";



  // Multiplication by constant

  p2 = mod(p2 * 70, product);
  d2 *= 70;

  cout << "Test 12 " << ((to_ZZX(d2)==p2)? "PASS" : "FAIL")
       << ": poly(d2*const) =?= poly(d2)*const\n";


  // Multiplication by polynomial
  ZZX PhimX = context.zMstar.PhimX();

  p2 = mod((p2 * p1) % PhimX, product);
  d2 *= d1;

  cout << "Test 13 " << ((to_ZZX(d2)==p2)? "PASS" : "FAIL")
       << ": poly(d2*d1) =?= poly(d2)*poly(d1) mod Phi_m\n";

  // same as above, but use d2 *= p1 instead of d2 *= d1
  p2 = mod((p2 * p1) % PhimX, product);
  d2 *= p1;

  cout << "Test 14 " << ((to_ZZX(d2)==p2)? "PASS" : "FAIL")
       << ": poly(d2*p1) =?= poly(d2)*p1 mod Phi_m\n";

  ModComp(p2, p2, ZZX(11, 1), PhimX);
  p2 = mod(p2, product);
  d2.automorph(11);

  cout << "Test 15 " << ((to_ZZX(d2)==p2)? "PASS" : "FAIL")
       << ": automorphism test\n";
  




#if 0

  // Division by constant
  d2.toPoly(p2);
  d2 /= 705;
  d2.toPoly(p3);   // p3 should be p2/705 mod prod
  p3 *= 705;
  PolyRed(p3,product); // map the coefficients to [-prod/2,prod/2)
  cout << "Test 15 " << ((p2==p3)? "PASS" : "FAIL")
       << ": poly(d2/const) =?= poly(d2)/const\n";

  // Automorphism test 
  d1.toPoly(p1);
  ctxt.restore();
  ZZ_pX pp = to_ZZ_pX(p1);
  ZZ_pX x11; SetCoeff(x11,11);  // x11 = X^11
  ZZ_pXModulus F(to_ZZ_pX(PhimX));
  CompMod(pp,pp,x11,F); // pp := p1(X^11) mod Phi_m(X)
  conv(p2,pp);          // convert back to ZZX format
  PolyRed(p2,product);  // map the coefficients to [-prod/2,prod/2)
  d1.automorph(11);            // X -> X^11
  d1.toPoly(p1);       // convert to polynomial representation
  cout << "Test 16 " << ((p2==p1)? "PASS" : "FAIL")
       << ": poly(d1 >> const) =?= poly(d1) >> const\n";

  // test SignleCRT arithmetic

  // Addition/subtraction
  d1.toSingleCRT(s1);
  d2.toSingleCRT(s2);
  d1.toPoly(p1);
  d2.toPoly(p2);
  s2 += s1;
  s2.toPoly(p3);   // p3 should be p1+p2 mod prod
  p2 += p1;
  PolyRed(p2,product); // map the coefficients to [-prod/2,prod/2)
  cout << "\nTest 17 " << ((p2==p3)? "PASS" : "FAIL")
       << ": poly(s2+s1) =?= poly(s2)+poly(s1)\n";

  s2.toPoly(p2);
  s2 -= s1;
  s2.toPoly(p3);   // p3 should be p2-p1 mod prod
  p2 -= p1;
  PolyRed(p2,product); // map the coefficients to [-prod/2,prod/2)

  cout << "Test 18 " << ((p2==p3)? "PASS" : "FAIL")
       << ": poly(s2-s1) =?= poly(s2)-poly(s1)\n";

  //  cout << p2 << "\n";
  //  cout << p3 << "\n";

  s2.toPoly(p2);
  s2 -= 3;
  s2.toPoly(p3);   // p3 should be p2-3 mod prod
  p2 -= 3;
  PolyRed(p2,product); // map the coefficients to [-prod/2,prod/2)
  cout << "Test 19 " << ((p2==p3)? "PASS" : "FAIL")
       << ": poly(s2-const) =?= poly(s2)-const\n";

  // in-place operations
  s1.toPoly(p1,1,1); // p1 = s1.polys[1]
  PolyRed(p1,context.ithPrime(1),true);  // reduce to [0,q1)
  s2 = s1; 
  if (s2 != s1) Error("s2 != s1 after assignment s2=s1");
  s1 -= s1.polys[1];
  s2 -= p1;
  cout << "Test 20 " << ((s2==s1)? "PASS" : "FAIL")
       << ": s -= poly(lvl-1) =?= s - lvl1\n";

  // Multiplication by constant
  d2.toPoly(p2);
  d2 *= 70;
  d2.toPoly(p3);   // p3 should be p2*70 mod prod
  p2 *= 70;
  PolyRed(p2,product); // map the coefficients to [-prod/2,prod/2)
  cout << "Test 21 " << ((p2==p3)? "PASS" : "FAIL")
       << ": poly(d2*const) =?= poly(d2)*const\n";

  // Division by constant
  s2.toPoly(p2);
  s2 /= 705;
  s2.toPoly(p3);   // p3 should be p2/705 mod prod
  p3 *= 705;
  PolyRed(p3,product); // map the coefficients to [-prod/2,prod/2)
  cout << "Test 22 " << ((p2==p3)? "PASS" : "FAIL")
       << ": poly(d2/const) =?= poly(d2)/const\n";

  cout << endl;

#endif
}
