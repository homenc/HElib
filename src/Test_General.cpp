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

/* Test_General.cpp - A general test program that uses a mix of operations over four ciphertexts.
 */
#include <NTL/ZZ.h>
#include <NTL/BasicThreadPool.h>
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include <NTL/lzz_pXFactoring.h>

#include <cassert>
#include <cstdio>

#ifdef DEBUG_PRINTOUT
#define debugCompare(ea,sk,p,c) {\
  NewPlaintextArray pp(ea);\
  ea.decrypt(c, sk, pp);\
  if (!equals(ea, pp, p)) { \
    std::cout << "oops:\n"; std::cout << p << "\n"; \
    std::cout << pp << "\n"; \
    exit(0); \
  }}
#else
#define debugCompare(ea,sk,p,c)
#endif

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

static bool noPrint = false;

void  TestIt(long R, long p, long r, long d, long c, long k, long w, 
               long L, long m, const Vec<long>& gens, const Vec<long>& ords)
{
  char buffer[32];
  if (!noPrint) {
    std::cout << "\n\n******** TestIt" << (isDryRun()? "(dry run):" : ":");
    std::cout << " R=" << R 
	      << ", p=" << p
	      << ", r=" << r
	      << ", d=" << d
	      << ", c=" << c
	      << ", k=" << k
	      << ", w=" << w
	      << ", L=" << L
	      << ", m=" << m
	      << ", gens=" << gens
	      << ", ords=" << ords
	      << endl;
  }
  vector<long> gens1, ords1;
  convert(gens1, gens);
  convert(ords1, ords);

  FHEcontext context(m, p, r, gens1, ords1);
  buildModChain(context, L, c);

  ZZX G;
  if (d == 0)
    G = context.alMod.getFactorsOverZZ()[0];
  else
    G = makeIrredPoly(p, d); 

  if (!noPrint) {
    context.zMStar.printout();
    std::cout << endl;

    std::cout << "security=" << context.securityLevel()<<endl;
    std::cout << "# ctxt primes = " << context.ctxtPrimes.card() << "\n";
    std::cout << "# bits in ctxt primes = " 
	 << long(context.logOfProduct(context.ctxtPrimes)/log(2.0) + 0.5) << "\n";
    std::cout << "# special primes = " << context.specialPrimes.card() << "\n";
    std::cout << "# bits in special primes = " 
	 << long(context.logOfProduct(context.specialPrimes)/log(2.0) + 0.5) << "\n";
    std::cout << "G = " << G << "\n";
  }

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(w); // A Hamming-weight-w secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need

  EncryptedArray ea(context, G);
  long nslots = ea.size();

  NewPlaintextArray p0(ea);
  NewPlaintextArray p1(ea);
  NewPlaintextArray p2(ea);
  NewPlaintextArray p3(ea);

  random(ea, p0);
  random(ea, p1);
  random(ea, p2);
  random(ea, p3);

  Ctxt c0(publicKey), c1(publicKey), c2(publicKey), c3(publicKey);
  ea.encrypt(c0, publicKey, p0);
  // {ZZX ppp0; ea.encode(ppp0, p0); c0.DummyEncrypt(ppp0);} // dummy encryption
  ea.encrypt(c1, publicKey, p1); // real encryption
  ea.encrypt(c2, publicKey, p2); // real encryption
  ea.encrypt(c3, publicKey, p3); // real encryption

  resetAllTimers();

  FHE_NTIMER_START(Circuit);

  for (long i = 0; i < R; i++) {

    if (!noPrint) std::cout << "*** round " << i << "..."<<endl;

     long shamt = RandomBnd(2*(nslots/2) + 1) - (nslots/2);
                  // random number in [-nslots/2..nslots/2]
     long rotamt = RandomBnd(2*nslots - 1) - (nslots - 1);
                  // random number in [-(nslots-1)..nslots-1]

     // two random constants
     NewPlaintextArray const1(ea);
     NewPlaintextArray const2(ea);
     random(ea, const1);
     random(ea, const2);

     ZZX const1_poly, const2_poly;
     ea.encode(const1_poly, const1);
     ea.encode(const2_poly, const2);

     mul(ea, p1, p0);     // c1.multiplyBy(c0)
     c1.multiplyBy(c0);
     if (!noPrint) CheckCtxt(c1, "c1*=c0");
     debugCompare(ea,secretKey,p1,c1);

     add(ea, p0, const1); // c0 += random constant
     c0.addConstant(const1_poly);
     if (!noPrint) CheckCtxt(c0, "c0+=k1");
     debugCompare(ea,secretKey,p0,c0);

     mul(ea, p2, const2); // c2 *= random constant
     c2.multByConstant(const2_poly);
     if (!noPrint) CheckCtxt(c2, "c2*=k2");
     debugCompare(ea,secretKey,p2,c2);

     NewPlaintextArray tmp_p(p1); // tmp = c1
     Ctxt tmp(c1);
     sprintf(buffer, "c2>>=%d", (int)shamt);
     shift(ea, tmp_p, shamt); // ea.shift(tmp, random amount in [-nSlots/2,nSlots/2])
     ea.shift(tmp, shamt);
     if (!noPrint) CheckCtxt(tmp, buffer);
     debugCompare(ea,secretKey,tmp_p,tmp);

     add(ea, p2, tmp_p);  // c2 += tmp
     c2 += tmp;
     if (!noPrint) CheckCtxt(c2, "c2+=tmp");
     debugCompare(ea,secretKey,p2,c2);

     sprintf(buffer, "c2>>>=%d", (int)rotamt);
     rotate(ea, p2, rotamt); // ea.rotate(c2, random amount in [1-nSlots, nSlots-1])
     ea.rotate(c2, rotamt);
     if (!noPrint) CheckCtxt(c2, buffer);
     debugCompare(ea,secretKey,p2,c2);

     ::negate(ea, p1); // c1.negate()
     c1.negate();
     if (!noPrint) CheckCtxt(c1, "c1=-c1");
     debugCompare(ea,secretKey,p1,c1);

     mul(ea, p3, p2); // c3.multiplyBy(c2) 
     c3.multiplyBy(c2);
     if (!noPrint) CheckCtxt(c3, "c3*=c2");
     debugCompare(ea,secretKey,p3,c3);

     sub(ea, p0, p3); // c0 -= c3
     c0 -= c3;
     if (!noPrint) CheckCtxt(c0, "c0=-c3");
     debugCompare(ea,secretKey,p0,c0);

  }

  c0.cleanUp();
  c1.cleanUp();
  c2.cleanUp();
  c3.cleanUp();

  FHE_NTIMER_STOP(Circuit);

  if (!noPrint) {
    std::cout << endl;
    printAllTimers();
    std::cout << endl;
  }
  resetAllTimers();
  FHE_NTIMER_START(Check);
   
  NewPlaintextArray pp0(ea);
  NewPlaintextArray pp1(ea);
  NewPlaintextArray pp2(ea);
  NewPlaintextArray pp3(ea);
   
  ea.decrypt(c0, secretKey, pp0);
  ea.decrypt(c1, secretKey, pp1);
  ea.decrypt(c2, secretKey, pp2);
  ea.decrypt(c3, secretKey, pp3);
   
  if (equals(ea, pp0, p0) && equals(ea, pp1, p1)
      && equals(ea, pp2, p2) && equals(ea, pp3, p3))
       std::cout << "GOOD\n";
  else std::cout << "BAD\n";

  FHE_NTIMER_STOP(Check);
   
  std::cout << endl;
  if (!noPrint) {
    printAllTimers();
    std::cout << endl;
  }

#if 0
  vector<Ctxt> vc(L,c0);            // A vector of L ciphertexts
  vector<NewPlaintextArray> vp(L, p0); // A vector of L plaintexts
  for (long i=0; i<L; i++) {
    vp[i].random();                     // choose a random plaintext 
    ea.encrypt(vc[i], publicKey, vp[i]); // encrypt it
    if (i>0) vp[i].mul(vp[i-1]); // keep a running product of plaintexts
  }
  incrementalProduct(vc); // Compute the same running product homomorphically

  // Check that the products match
  bool fail = false;
  for (long i=0; i<L; i++) {
    ea.decrypt(vc[i], secretKey, p0); // decrypt it
    if (!p0.equals(vp[i])) {
      fail = true;
      std::cout << "incrementalProduct oops "<<i<< endl;
    }
  }
  if (!fail) std::cout << "incrementalProduct works\n";
#endif
}


/* A general test program that uses a mix of operations over four ciphertexts.
 * Usage: Test_General_x [ name=value ]...
 *   R       number of rounds  [ default=1 ]
 *   p       plaintext base  [ default=2 ]
 *   r       lifting  [ default=1 ]
 *   d       degree of the field extension  [ default=1 ]
 *              d == 0 => factors[0] defines extension
 *   c       number of columns in the key-switching matrices  [ default=2 ]
 *   k       security parameter  [ default=80 ]
 *   L       # of levels in the modulus chain  [ default=heuristic ]
 *   s       minimum number of slots  [ default=0 ]
 *   repeat  number of times to repeat the test  [ default=1 ]
 *   m       use specified value as modulus
 *   mvec    use product of the integers as  modulus
 *              e.g., mvec='[5 3 187]' (this overwrite the m argument)
 *   gens    use specified vector of generators
 *              e.g., gens='[562 1871 751]'
 *   ords    use specified vector of orders
 *              e.g., ords='[4 2 -4]', negative means 'bad'
 */
int main(int argc, char **argv) 
{
  setTimersOn();

  ArgMapping amap;

  bool dry=false;
  amap.arg("dry", dry, "dry=1 for a dry-run");

  long R=1;
  amap.arg("R", R, "number of rounds");

  long p=2;
  amap.arg("p", p, "plaintext base");

  long r=1;
  amap.arg("r", r,  "lifting");

  long d=1;
  amap.arg("d", d, "degree of the field extension");
  amap.note("d == 0 => factors[0] defines extension");

  long c=2;
  amap.arg("c", c, "number of columns in the key-switching matrices");

  
  long k=80;
  amap.arg("k", k, "security parameter");

  long L=0;
  amap.arg("L", L, "# of levels in the modulus chain",  "heuristic");

  long s=0;
  amap.arg("s", s, "minimum number of slots");

  long repeat=1;
  amap.arg("repeat", repeat,  "number of times to repeat the test");

  long chosen_m=0;
  amap.arg("m", chosen_m, "use specified value as modulus", NULL);

  Vec<long> mvec;
  amap.arg("mvec", mvec, "use product of the integers as  modulus", NULL);
  amap.note("e.g., mvec='[5 3 187]' (this overwrite the m argument)");

  Vec<long> gens;
  amap.arg("gens", gens, "use specified vector of generators", NULL);
  amap.note("e.g., gens='[562 1871 751]'");

  Vec<long> ords;
  amap.arg("ords", ords, "use specified vector of orders", NULL);
  amap.note("e.g., ords='[4 2 -4]', negative means 'bad'");

  long seed=0;
  amap.arg("seed", seed, "PRG seed");

  long nt=1;
  amap.arg("nt", nt, "num threads");

  amap.arg("noPrint", noPrint, "suppress printouts");

  amap.parse(argc, argv);

  SetSeed(ZZ(seed));
  SetNumThreads(nt);
  
  if (L==0) { // determine L based on R,r
    L = 3*R+3;
    if (p>2 || r>1) { // add some more primes for each round
      long addPerRound = 2*ceil(log((double)p)*r*3)/(log(2.0)*FHE_p2Size) +1;
      L += R * addPerRound;
    }
  }

  long w = 64; // Hamming weight of secret key
  //  long L = z*R; // number of levels

  if (mvec.length()>0)
    chosen_m = computeProd(mvec);
  std::cout << argv[0] << ": ";
  long m = FindM(k, L, c, p, d, s, chosen_m, !noPrint);

  setDryRun(dry);
  for (long repeat_cnt = 0; repeat_cnt < repeat; repeat_cnt++) {
    TestIt(R, p, r, d, c, k, w, L, m, gens, ords);
  }
}

// call to get our running test case:
// Test_General_x p=23 m=20485 L=10 R=5
//
// another call to get an example where phi(m) is very
// close to m:
// Test_General_x m=18631 L=10 R=5
