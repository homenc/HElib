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
#include <NTL/vec_long.h>
#include "NumbTh.h"
#include "FHEContext.h"

#define pSize 16   /* empirical average size of small primes */
#define p0Size 29  /* size of first 1-2 small primes */

NTL_CLIENT




long FindM(long k, long L, long c, long p, long d, long s, long chosen_m, bool verbose)
{
  // get a lower-bound on the parameter N=phi(m):
  // 1. Empirically, we use ~20-bit small primes in the modulus chain (the main
  //    constraints is that 2m must divide p-1 for every prime p). The first
  //    prime is larger, a 40-bit prime. (If this is a 32-bit machine then we
  //    use two 20-bit primes instead.)
  // 2. With L levels, the largest modulus for "fresh ciphertexts" has size
  //          q0 ~ p0 * p^{L} ~ 2^{40+20L}
  // 3. We break each ciphertext into upto c digits, do each digit is as large
  //    as    D=2^{(40+20L)/c}
  // 4. The added noise variance term from the key-switching operation is
  //    c*N*sigma^2*D^2, and this must be mod-switched down to w*N (so it is
  //    on part with the added noise from modulus-switching). Hence the ratio
  //    P that we use for mod-switching must satisfy c*N*sigma^2*D^2/P^2<w*N,
  //    or    P > sqrt(c/w) * sigma * 2^{(40+20L)/c}
  // 5. With this extra P factor, the key-switching matrices are defined
  //    relative to a modulus of size
  //          Q0 = q0*P ~ sqrt{c/w} sigma 2^{(40+20L)(1+1/c)}
  // 6. To get k-bit security we need N>log(Q0/sigma)(k+110)/7.2, i.e. roughly
  //          N > (40+20L)(1+1/c)(k+110) / 7.2

  double cc = 1.0+(1.0/(double)c);
  long N = (long) ceil((pSize*L+p0Size)*cc*(k+110)/7.2);

  // pre-computed values of [phi(m),m,d]
  long ms[][4] = {
    //phi(m)  m  ord(2) c_m*1000
    { 1176,  1247, 28,  3736}, // gens=5(42)
    { 2880,  3133, 24,  3254}, // gens=6(60), 7(!2)
    { 4096,  4369, 16,  3422}, // gens=129(16),3(!16)
    { 5292,  5461, 14,  4160}, // gens=3(126),509(3)
    { 5760,  8435, 24,  8935}, // gens=58(60),1686(2),11(!2)
    { 8190,  8191, 13,  1273}, // gens=39(630)
    {10752, 11441, 48,  3607}, // gens=7(112),5(!2)
    {12000, 13981, 20,  2467}, // gens=10(30),23(10),3(!2)
    {11520, 15665, 24, 14916}, // gens=6(60),177(4),7(!2)
    {15004, 15709, 22,  3867}, // gens=5(682)
    {18000, 18631, 25,  4208}, // gens=17(120),1177(6)
    {15360, 20485, 24, 12767}, // gens=6(80),242(4),7(2)
    {16384, 21845, 16, 12798}, // gens=129(16),273(4),3(!16)
    {17280 ,21931, 24, 18387}, // gens=6(60),467(6),11(!2)
    {19200, 21607, 40, 35633}, // gens=13(120),2789(2),3(!2)
    {21168, 27305, 28, 15407}, // gens=6(126),781(6)
    {23040, 23377, 48,  5292}, // gens=35(240),5(!2)
    {24576, 24929, 48,  5612}, // gens=12(256),5(!2)
    {27000, 32767, 15, 20021}, // gens=3(150),873(6),6945(2)
    {31104, 31609, 72,  5149}, // gens=11(216),5(!2)
    {43690, 43691, 34, 0},     // gens=69(1285)
    {49500, 49981, 30, 0},     // gens=3(1650)
    {46080, 53261, 24, 33409}, // gens=3(240),242(4),7(!2)
    {54000, 55831, 25, 0},     // gens=22(360),3529(6)
    {49140, 57337, 39,  2608}, // gens=39(630),40956(2)
    {51840, 59527, 72, 21128}, // gens=58(60),1912(6),7(!2)
    {61680, 61681, 40,  1273}, // gens=33(771),17(!2)
    {65536, 65537, 32,  1273}, // gens=2(32),3(!2048)
    {75264, 82603, 56, 36484}, // gens=3(336),24294(2),7(!2)
    {84672, 92837, 56, 38520}  // gens=18(126),1886(6),3(!2)
  };

  long m = 0;
  size_t i=0;

  // find the first m satisfying phi(m)>=N and d | ord(p) in Z_m^*
  // and phi(m)/ord(p) >= s
  if (chosen_m) {
    if (GCD(p, chosen_m) == 1) {
      long ordP = multOrd(p, chosen_m);
      if (d == 0 || ordP % d == 0) {
        // chosen_m is OK
        m = chosen_m;
      }
    }
  }
  else {
    for (i=0; i<sizeof(ms)/sizeof(long[4]); i++) { 
      if (ms[i][0] < N || GCD(p, ms[i][1]) != 1) continue;
      long ordP = multOrd(p, ms[i][1]);
      long nSlots = ms[i][0]/ordP;
      if (d != 0 && ordP % d != 0) continue;
      if (nSlots < s) continue;

      m = ms[i][1];
      break;
    }
  }

  if (m==0) Error("Cannot support these parameters");
  if (verbose) {
    cerr << "*** Bound N="<<N<<", choosing m="<<m <<", phi(m)="<< ms[i][0]
         << endl;
  }

  return m;
}

// A global variable, pointing to the "current" context
FHEcontext* activeContext = NULL;

void FHEcontext::productOfPrimes(ZZ& p, const IndexSet& s) const
{
  p = 1;
  for (long i = s.first(); i <= s.last(); i = s.next(i))
    p *= ithPrime(i);
}

void FHEcontext::AddPrime(long p, bool special)
{
  long twoM = 2 * zMStar.getM();  // ensure p-1 is divisible by 2m
  assert( ProbPrime(p) && p % twoM == 1 && !inChain(p) );

  long i = moduli.size();
  moduli.push_back( Cmodulus(zMStar, p, 0) );

  if (special)
    specialPrimes.insert(i);
  else
    ctxtPrimes.insert(i);
}

// Adds to the chain primes whose product is at least totalSize bits
double AddPrimesBySize(FHEcontext& context, double totalSize, bool special)
{
  if (!context.zMStar.getM() || context.zMStar.getM() > (1<<20)) // sanity checks
    Error("AddModuli1: m undefined or larger than 2^20");

  long p = (1UL << NTL_SP_NBITS)-1;   // Start from as large prime as possible
  long twoM = 2 * context.zMStar.getM(); // make p-1 divisible by 2m
  p -= (p%twoM); // 0 mod 2m
  p += twoM +1;  // 1 mod 2m, a 2m quantity is subtracted below

  bool lastPrime = false;
  double sizeLeft = totalSize; // how much is left to do
  while (sizeLeft > 0.0) {
    // A bit of convoluted logic attempting not to overshoot totalSize by much
    if (sizeLeft < log((double)p) && !lastPrime) { // decrease p
      lastPrime = true;
      p = ceil(exp(sizeLeft));
      p -= (p%twoM)-1; // make p-1 divisible by 2m
      twoM = -twoM;    // increase p below, rather than decreasing it
    }
    do { p -= twoM; } while (!ProbPrime(p)); // next prime
    if (!context.inChain(p)) {
      context.AddPrime(p, special);
      sizeLeft -= log((double)p);
    }
  }
  return totalSize-sizeLeft;
}

// Adds nPrimes primes to the chain, returns the bitsize of the product of
// all primes in the chain.
double AddPrimesByNumber(FHEcontext& context, long nPrimes, 
			 long p/*=starting point*/, bool special)
{
  if (!context.zMStar.getM() || context.zMStar.getM() > (1<<20))  // sanity checks
    Error("FHEcontext::AddModuli2: m undefined or larger than 2^20");

  long twoM = 2 * context.zMStar.getM();

  // make sure that p>0 and that p-1 is divisible by m
  if (p<1) p = 1;
  p -= (p % twoM) -1;

  double sizeSoFar = 0.0;
  while (nPrimes>0) {
    do { p += twoM; } while (!ProbPrime(p)); // next prime
    if (!context.inChain(p)) {
      context.AddPrime(p, special);
      nPrimes -= 1;
      sizeSoFar += log((double)p);
    }
  }
  return sizeSoFar;
}

void buildModChain(FHEcontext &context, long nLvls, long nDgts)
{
  // The first 1-2 primes of total p0size bits
  #if (NTL_SP_NBITS > p0Size)
    AddPrimesByNumber(context, 1, 1UL<<p0Size); // add a single prime
  #else
    AddPrimesByNumber(context, 2, 1UL<<(p0Size/2)); // add two primes
  #endif

  // The next nLvls primes, as small as possible
  AddPrimesByNumber(context, nLvls);

  // calculate the size of the digits

  context.digits.resize(nDgts); // allocate space

  IndexSet s1;
  double sizeSoFar = 0.0;
  double maxDigitSize = 0.0;
  if (nDgts>1) { // we break ciphetext into a few digits when key-switching
    double dsize = context.logOfProduct(context.ctxtPrimes)/nDgts; // estimate
    double target = dsize-(pSize/3.0);
    long idx = context.ctxtPrimes.first();
    for (long i=0; i<nDgts-1; i++) { // compute next digit
      IndexSet s;
      while (idx <= context.ctxtPrimes.last() && (empty(s)||sizeSoFar<target)) {
        s.insert(idx);
	sizeSoFar += log((double)context.ithPrime(idx));
	idx = context.ctxtPrimes.next(idx);
      }
      assert (!empty(s));
      context.digits[i] = s;
      s1.insert(s);
      double thisDigitSize = context.logOfProduct(s);
      if (maxDigitSize < thisDigitSize) maxDigitSize = thisDigitSize;
      target += dsize;
    }
    IndexSet s = context.ctxtPrimes / s1; // all the remaining primes
    context.digits[nDgts-1] = s;
    double thisDigitSize = context.logOfProduct(s);
    if (maxDigitSize < thisDigitSize) maxDigitSize = thisDigitSize;
  }
  else { 
    maxDigitSize = context.logOfProduct(context.ctxtPrimes);
    context.digits[0] = context.ctxtPrimes;
  }

  // Add primes to the chain for the P factor of key-switching
  double sizeOfSpecialPrimes 
    = maxDigitSize + log(nDgts/32.0)/2 + log(context.stdev *2);

  AddPrimesBySize(context, sizeOfSpecialPrimes, true);
}

bool FHEcontext::operator==(const FHEcontext& other) const
{
  if (zMStar != other.zMStar) return false;
  if (alMod != other.alMod) return false;

  if (moduli.size() != other.moduli.size()) return false;
  for (size_t i=0; i<moduli.size(); i++) {
    const Cmodulus& m1 = moduli[i];
    const Cmodulus& m2 = other.moduli[i];
    if (m1.getQ() != m2.getQ()) return false;
  }

  if (ctxtPrimes != other.ctxtPrimes) return false;
  if (specialPrimes != other.specialPrimes) return false;

  if (digits.size() != other.digits.size()) return false;
  for (size_t i=0; i<digits.size(); i++)
    if (digits[i] != other.digits[i]) return false;

  if (stdev != other.stdev) return false;
  return true;
}

void writeContextBase(ostream& str, const FHEcontext& context)
{
  str << "[" << context.zMStar.getM() 
      << " " << context.zMStar.getP()
      << " " << context.alMod.getR() 
      << "]\n";
}

ostream& operator<< (ostream &str, const FHEcontext& context)
{
  str << "[\n";

  // standard-deviation
  str << context.stdev << "\n";

  // the "special" index 
  str << context.specialPrimes << "\n ";

  // output the primes in the chain
  str << context.moduli.size() << "\n";
  for (long i=0; i<(long)context.moduli.size(); i++)
    str << context.moduli[i].getQ() << " ";
  str << "\n ";

  // output the digits
  str << context.digits.size() << "\n";
  for (long i=0; i<(long)context.digits.size(); i++)
    str << context.digits[i] << " ";

  str <<"\n";

  str << "]\n";

  return str;
}

void readContextBase(istream& str, unsigned long& m, unsigned long& p, unsigned long& r)
{
  // Advance str beyond first '[' 
  seekPastChar(str, '[');  // this function is defined in NumbTh.cpp

  str >> m >> p >> r;

  seekPastChar(str, ']');
}

istream& operator>> (istream &str, FHEcontext& context)
{
  seekPastChar(str, '[');  // this function is defined in NumbTh.cpp

  // Get the standard deviation
  str >> context.stdev;

  IndexSet s;
  str >> s; // read the special set

  context.moduli.clear();
  context.specialPrimes.clear();
  context.ctxtPrimes.clear();

  long nPrimes;
  str >> nPrimes;
  for (long i=0; i<nPrimes; i++) {
    long p;
    str >> p; 
    context.AddPrime(p, s.contains(i));
  }

  // read in the partition to digits
  long nDigits;
  str >> nDigits;
  context.digits.resize(nDigits);
  for (long i=0; i<(long)context.digits.size(); i++)
    str >> context.digits[i];

  seekPastChar(str, ']');
  return str;
}

