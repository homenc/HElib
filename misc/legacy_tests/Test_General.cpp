/* Copyright (C) 2012-2020 IBM Corp.
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
#include <helib/helib.h>
#include <NTL/lzz_pXFactoring.h>

#include <cassert>
#include <cstdio>
#include <helib/ArgMap.h>
#include <helib/fhe_stats.h>
#include <helib/log.h>

NTL_CLIENT
using namespace helib;

//#define HELIB_DEBUG

#if 1
#include <helib/debugging.h>
#define debugCompare(ea,sk,p,c) {\
  PtxtArray pp(ea);\
  pp.decrypt(c, sk);\
  if (pp != p) { \
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

static bool noPrint = true;
static long special_bits;

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

  if (!noPrint) fhe_stats = true;

  
  Context context { ContextBuilder<BGV>()

      .m(m).p(p).r(r)
      .gens(gens1).ords(ords1)
      .bits(L).c(c).bitsInSpecialPrimes(special_bits)

  };

#if 0
  // This works in C++17, but not earlier, as it then
  // requires a copy or move constructor
  Context context =
    ContextBuilder<BGV>()
      .m(m).p(p).r(r)
      .gens(gens1).ords(ords1)
      .bits(L).c(c).bitsInSpecialPrimes(special_bits);
#endif

#if 0
  Context context(m, p, r, gens1, ords1);
  buildModChain(context, L, c,
    /*willBeBootstrappable=*/false,
    /*t=*/0,
    /*resolution=*/3,
    /*bitsInSpecialPrimes=*/special_bits);
#endif

  ZZX G;
  if (d == 0)
    G = context.alMod.getFactorsOverZZ()[0];
  else
    G = makeIrredPoly(p, d);

  if (!noPrint) {
    context.zMStar.printout();
    std::cout << endl;

    std::cout << "security=" << context.securityLevel()<<endl;
    std::cout << "# small primes = " << context.smallPrimes.card() << "\n";
    std::cout << "# ctxt primes = " << context.ctxtPrimes.card() << "\n";
    std::cout << "# bits in ctxt primes = "
	 << long(context.logOfProduct(context.ctxtPrimes)/log(2.0) + 0.5) << "\n";
    std::cout << "# special primes = " << context.specialPrimes.card() << "\n";
    std::cout << "# bits in special primes = "
	 << long(context.logOfProduct(context.specialPrimes)/log(2.0) + 0.5) << "\n";
    std::cout << "G = " << G << "\n";

  }

  SecKey secretKey(context);
  const PubKey& publicKey = secretKey;
  secretKey.GenSecKey(); // A +-1/0 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need

  std::shared_ptr<EncryptedArray> ea_ptr(std::make_shared<EncryptedArray>(context, G));
  // Alias to avoid issues with previous code
  EncryptedArray& ea = *ea_ptr;
  long nslots = ea.size();
#ifdef HELIB_DEBUG
  dbgKey = &secretKey;
  dbgEa  = ea_ptr;
#endif

  PtxtArray p0(ea), p1(ea), p2(ea), p3(ea);

  p0.random();
  p1.random();
  p2.random();
  p3.random();

  Ctxt c0(publicKey), c1(publicKey), c2(publicKey), c3(publicKey);
  p0.encrypt(c0);
  p1.encrypt(c1); 
  p2.encrypt(c2); 
  p3.encrypt(c3); 

  resetAllTimers();

  HELIB_NTIMER_START(Circuit);

  for (long i = 0; i < R; i++) {

    if (!noPrint) std::cout << "*** round " << i << "..."<<endl;

     long shamt = RandomBnd(2*(nslots/2) + 1) - (nslots/2);
                  // random number in [-nslots/2..nslots/2]
     long rotamt = RandomBnd(2*nslots - 1) - (nslots - 1);
                  // random number in [-(nslots-1)..nslots-1]

     // two random constants
     PtxtArray const1(ea), const2(ea);
     const1.random();
     const2.random();


     p1 *= p0;     // c1.multiplyBy(c0)
     c1.multiplyBy(c0);
     if (!noPrint) CheckCtxt(c1, "c1*=c0");
     debugCompare(ea,secretKey,p1,c1);

     p0 += const1; // c0 += random constant
     c0.addConstant(const1);
     if (!noPrint) CheckCtxt(c0, "c0+=k1");
     debugCompare(ea,secretKey,p0,c0);

     p2 *= const2; // c2 *= random constant
     c2.multByConstant(const2);
     if (!noPrint) CheckCtxt(c2, "c2*=k2");
     debugCompare(ea,secretKey,p2,c2);

#if 0
     p2 *= 3L;
     c2 *= 3L;
     if (!noPrint) CheckCtxt(c2, "c2*=3");
     debugCompare(ea,secretKey,p2,c2);

     p2 += 3L;
     c2 += 3L;
     if (!noPrint) CheckCtxt(c2, "c2+=3");
     debugCompare(ea,secretKey,p2,c2);


#endif

     PtxtArray tmp_p(p1); // tmp = c1
     Ctxt tmp(c1);
     sprintf(buffer, "tmp=c1>>=%d", (int)shamt);
     shift(tmp_p, shamt); // ea.shift(tmp, random amount in [-nSlots/2,nSlots/2])
     ea.shift(tmp, shamt);
     if (!noPrint) CheckCtxt(tmp, buffer);
     debugCompare(ea,secretKey,tmp_p,tmp);

     p2 += tmp_p;  // c2 += tmp
     c2 += tmp;
     if (!noPrint) CheckCtxt(c2, "c2+=tmp");
     debugCompare(ea,secretKey,p2,c2);

     sprintf(buffer, "c2>>>=%d", (int)rotamt);
     rotate(p2, rotamt); // ea.rotate(c2, random amount in [1-nSlots, nSlots-1])
     ea.rotate(c2, rotamt);
     if (!noPrint) CheckCtxt(c2, buffer);
     debugCompare(ea,secretKey,p2,c2);

     p1.negate(); // c1.negate()
     c1.negate();
     if (!noPrint) CheckCtxt(c1, "c1=-c1");
     debugCompare(ea,secretKey,p1,c1);

     p3 *= p2; // c3.multiplyBy(c2)
     c3.multiplyBy(c2);
     if (!noPrint) CheckCtxt(c3, "c3*=c2");
     debugCompare(ea,secretKey,p3,c3);

     p0 -= p3; // c0 -= c3
     c0 -= c3;
     if (!noPrint) CheckCtxt(c0, "c0=-c3");
     debugCompare(ea,secretKey,p0,c0);

  }

  c0.cleanUp();
  c1.cleanUp();
  c2.cleanUp();
  c3.cleanUp();

  HELIB_NTIMER_STOP(Circuit);

  if (!noPrint) {
    std::cout << endl;
    printAllTimers();
    std::cout << endl;
  }
  resetAllTimers();
  HELIB_NTIMER_START(Check);

  PtxtArray pp0(ea);
  PtxtArray pp1(ea);
  PtxtArray pp2(ea);
  PtxtArray pp3(ea);

  pp0.decrypt(c0, secretKey);
  pp1.decrypt(c1, secretKey);
  pp2.decrypt(c2, secretKey);
  pp3.decrypt(c3, secretKey);

  if (pp0 == p0 && pp1 == p1
      && pp2 == p2 && pp3 == p3)
       std::cout << "GOOD\n";
  else std::cout << "BAD\n";

  HELIB_NTIMER_STOP(Check);

  if (!noPrint) {
    printAllTimers();
    print_stats(cout);
    std::cout << endl;
  }
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
 *   L       # of bits in the modulus chain  [ default=heuristic ]
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
  helog.setLogToStderr();

  setTimersOn();

  ArgMap amap;

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

  long L=500;
  amap.arg("L", L, "# of bits in the modulus chain");

  long s=0;
  amap.arg("s", s, "minimum number of slots");

  long repeat=1;
  amap.arg("repeat", repeat,  "number of times to repeat the test");

  long chosen_m=0;
  amap.arg("m", chosen_m, "use specified value as modulus", nullptr);

  Vec<long> mvec;
  amap.arg("mvec", mvec, "use product of the integers as  modulus", nullptr);
  amap.note("e.g., mvec='[5 3 187]' (this overwrite the m argument)");

  Vec<long> gens;
  amap.arg("gens", gens, "use specified vector of generators", nullptr);
  amap.note("e.g., gens='[562 1871 751]'");

  Vec<long> ords;
  amap.arg("ords", ords, "use specified vector of orders", nullptr);
  amap.note("e.g., ords='[4 2 -4]', negative means 'bad'");

  long seed=0;
  amap.arg("seed", seed, "PRG seed");

  long nt=1;
  amap.arg("nt", nt, "num threads");

  amap.arg("noPrint", noPrint, "suppress printouts");

  amap.arg("special_bits", special_bits, "# of bits in special primes");

  amap.parse(argc, argv);

  SetSeed(ZZ(seed));
  SetNumThreads(nt);


  long w = 64; // Hamming weight of secret key

  if (mvec.length()>0)
    chosen_m = computeProd(mvec);
  long m = FindM(k, L, c, p, d, s, chosen_m, !noPrint);

  setDryRun(dry);
  for (long repeat_cnt = 0; repeat_cnt < repeat; repeat_cnt++) {
    TestIt(R, p, r, d, c, k, w, L, m, gens, ords);
  }
}

