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

#include "FHEContext.h"
#include "EvalMap.h"
#include "powerful.h"


long FindM(long k, long L, long c, long p, long d, long s, long chosen_m, bool verbose)
{
  // get a lower-bound on the parameter N=phi(m):
  // 1. Each level in the modulus chain corresponds to pSize=p2Size/2
  //    bits (where we have one prime of this size, and all the others are of
  //    size p2Size).
  //    When using DoubleCRT, we need 2m to divide q-1 for every prime q.
  // 2. With L levels, the largest modulus for "fresh ciphertexts" has size
  //          Q0 ~ p^{L+1} ~ 2^{(L+1)*pSize}
  // 3. We break each ciphertext into upto c digits, do each digit is as large
  //    as    D=2^{(L+1)*pSize/c}
  // 4. The added noise variance term from the key-switching operation is
  //    c*N*sigma^2*D^2, and this must be mod-switched down to w*N (so it is
  //    on par with the added noise from modulus-switching). Hence the ratio
  //    P that we use for mod-switching must satisfy c*N*sigma^2*D^2/P^2<w*N,
  //    or    P > sqrt(c/w) * sigma * 2^{(L+1)*pSize/c}
  // 5. With this extra P factor, the key-switching matrices are defined
  //    relative to a modulus of size
  //          Q0 = q0*P ~ sqrt{c/w} sigma 2^{(L+1)*pSize*(1+1/c)}
  // 6. To get k-bit security we need N>log(Q0/sigma)(k+110)/7.2, i.e. roughly
  //          N > (L+1)*pSize*(1+1/c)(k+110) / 7.2

  // Compute a bound on m, and make sure that it is not too large
  double cc = 1.0+(1.0/(double)c);
  double dN = ceil((L+1)*FHE_pSize*cc*(k+110)/7.2);
  long N = NTL_SP_BOUND;
  if (N > dN) N = dN;
  else {
    cerr << "Cannot support a bound of " << dN;
    Error(", aborting.\n");
  }

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
  else if (p==2) { // use pre-computed table, divisors of 2^n-1 for some n's

    static long ms[][4] = {  // pre-computed values of [phi(m),m,d]
      //phi(m), m, ord(2),c_m*1000 (not used anymore)
      { 1176,  1247, 28,  3736}, // gens=5(42)
      { 2880,  3133, 24,  3254}, // gens=6(60), 7(!2)
      { 4050,  4051, 50, 0},     // gens=130(81)
      { 4096,  4369, 16,  3422}, // gens=129(16),3(!16)
      { 4704,  4859, 28, 0},     // gens=7(42),3(!4)
      { 5292,  5461, 14,  4160}, // gens=3(126),509(3)
      { 5760,  8435, 24,  8935}, // gens=58(60),1686(2),11(!2)
      { 7500,  7781, 50, 0},     // gens=353(30),3(!5)
      { 8190,  8191, 13,  1273}, // gens=39(630)
      { 9900, 10261, 30, 0},     // gens=3(330)
      {10752, 11441, 48,  3607}, // gens=7(112),5(!2)
      {10800, 11023, 45, 0},     // gens=270(24),2264(2),3(!5)
      {12000, 13981, 20,  2467}, // gens=10(30),23(10),3(!2)
      {11520, 15665, 24, 14916}, // gens=6(60),177(4),7(!2)
      {14112, 14351, 18, 0},     // gens=7(126),3(!4)
      {15004, 15709, 22,  3867}, // gens=5(682)
      {18000, 18631, 25,  4208}, // gens=17(120),1177(6)
      {15360, 20485, 24, 12767}, // gens=6(80),242(4),7(2)
      {16384, 21845, 16, 12798}, // gens=129(16),273(4),3(!16)
      {17280 ,21931, 24, 18387}, // gens=6(60),467(6),11(!2)
      {19200, 21607, 40, 35633}, // gens=13(120),2789(2),3(!2)
      {21168, 27305, 28, 15407}, // gens=6(126),781(6)
      {23040, 23377, 48,  5292}, // gens=35(240),5(!2)
      {23310, 23311, 45, 0},     // gens=489(518)
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

  // If m is not set yet, just set it close to N. This may be a lousy
  // choice of m for this p, since you will get a small number of slots.

  if (m==0) {
    // search only for odd values of m, to make phi(m) a little closer to m
    for (long candidate=N|1; candidate<10*N; candidate+=2) {
      if (GCD(p,candidate)!=1) continue;

      long ordP = multOrd(p,candidate); // the multiplicative order of p mod m
      if (d>1 && ordP%d!=0 ) continue;
      if (ordP > 100) continue;  // order too big, we will get very few slots

      long n = phi_N(candidate); // compute phi(m)
      if (n < N) continue;       // phi(m) too small

      m = candidate;  // all tests passed, return this value of m
      break;
    }
  }

  if (verbose) {
    cerr << "*** Bound N="<<N<<", choosing m="<<m <<", phi(m)="<< phi_N(m)
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

// Find the next prime and add it to the chain
long FHEcontext::AddPrime(long initialP, long delta, bool special, 
                          bool findRoot)
{
  // long twoM = 2 * zMStar.getM();
  // assert((initialP % twoM == 1) && (delta % twoM == 0));
  // NOTE: this assertion will fail for the half-prime in ALT_CRT

  long p = initialP;
  do { p += delta; } // delta could be positive or negative
  while (p>initialP/16 && p<NTL_SP_BOUND && !(ProbPrime(p) && !inChain(p)));

  if (p<=initialP/16 || p>=NTL_SP_BOUND) return 0; // no prime found

  long i = moduli.size(); // The index of the new prime in the list
  moduli.push_back( Cmodulus(zMStar, p, findRoot ? 0 : 1) );

  if (special)
    specialPrimes.insert(i);
  else
    ctxtPrimes.insert(i);

  return p;
}

long FHEcontext::AddFFTPrime(bool special)
{
  zz_pBak bak; bak.save(); // Backup the NTL context

  do {
    zz_p::FFTInit(fftPrimeCount);
    fftPrimeCount++;
  } while (inChain(zz_p::modulus()));

  long i = moduli.size(); // The index of the new prime in the list
  long p = zz_p::modulus();

  moduli.push_back( Cmodulus(zMStar, 0, 1) ); // a dummy Cmodulus object

  if (special)
    specialPrimes.insert(i);
  else
    ctxtPrimes.insert(i);

  return p;
}

// Adds several primes to the chain. If byNumber=true then totalSize specifies
// the number of primes to add. If byNumber=false then totalSize specifies the
// target total bitsize of all the added primes.
// The function returns the total bitsize of all the added primes.
double AddManyPrimes(FHEcontext& context, double totalSize, 
		     bool byNumber, bool special)
{
  double nBits = 0.0;     // How many bits added so far
  double sizeSoFar = 0.0;
  if (!context.zMStar.getM() || context.zMStar.getM()>(1<<20))// sanity checks
    Error("AddManyPrimes: m undefined or larger than 2^20");

  if (ALT_CRT) {
    while (sizeSoFar < totalSize) {
      long p = context.AddFFTPrime(special);
      nBits += log((double)p);
      sizeSoFar = byNumber? (sizeSoFar+1.0) : nBits;
    }
  }
  else {
#ifdef NO_HALF_SIZE_PRIME
    long sizeBits = context.bitsPerLevel;
#else
    long sizeBits = 2*context.bitsPerLevel;
#endif
    if (special) {
      long numPrimes = ceil(totalSize/NTL_SP_NBITS);// how many special primes
      sizeBits = ceil(totalSize/numPrimes);         // what's the size of each
    }
    long twoM = 2 * context.zMStar.getM();

    if (sizeBits>NTL_SP_NBITS) sizeBits = NTL_SP_NBITS;
    long sizeBound = 1L << sizeBits;
    if (sizeBound < twoM*log2(twoM)*8) {
      sizeBits = ceil(log2(twoM*log2(twoM)))+3;
      sizeBound = 1L << sizeBits;
    }

    // make p-1 divisible by m*2^k for as large k as possible
    while (twoM < sizeBound/(sizeBits*2)) twoM *= 2;

    long bigP = sizeBound - (sizeBound%twoM) +1; // 1 mod 2m
    long p = bigP+twoM; // The twoM is subtracted in the AddPrime function

    // FIXME: The last prime could be smaller
    while (sizeSoFar < totalSize) {
      if ((p = context.AddPrime(p,-twoM,special))) { // found a prime
        nBits += log((double)p);
        sizeSoFar = byNumber? (sizeSoFar+1.0) : nBits;
      }
      else { // we ran out of primes, try a lower power of two
        twoM /= 2;
        assert(twoM > (long)context.zMStar.getM()); // can we go lower?
        p = bigP;
      }
    }
  }
  return nBits;
}

void buildModChain(FHEcontext &context, long nLevels, long nDgts)
{
#ifdef NO_HALF_SIZE_PRIME
  long nPrimes = nLevels;
#else
  long nPrimes = (nLevels+1)/2;
  // The first prime should be of half the size. The code below tries to find
  // a prime q0 of this size where q0-1 is divisible by 2^k * m for some k>1.
  // Then if the plaintext space is a power of two it tries to choose the
  // second prime q1 so that q0*q1 = 1 mod ptxtSpace. All the other primes are
  // chosen so that qi-1 is divisible by 2^k * m for as large k as possible.
  long twoM;
  if (ALT_CRT) 
    twoM = 2;
  else
    twoM = 2 * context.zMStar.getM();

  long bound = (1L << (context.bitsPerLevel-1));
  while (twoM < bound/(2*context.bitsPerLevel))
    twoM *= 2; // divisible by 2^k * m  for a larger k

  bound = bound - (bound % twoM) +1; // = 1 mod 2m
  long q0 = context.AddPrime(bound, twoM, false, !ALT_CRT); 
  // add next prime to chain
  
  assert(q0 != 0);
  nPrimes--;
#endif

  // Choose the next primes as large as possible
  if (nPrimes>0) AddPrimesByNumber(context, nPrimes);

  // calculate the size of the digits

  if (nDgts > nPrimes) nDgts = nPrimes; // sanity checks
  if (nDgts <= 0) nDgts = 1;
  context.digits.resize(nDgts); // allocate space

  IndexSet s1;
  double sizeSoFar = 0.0;
  double maxDigitSize = 0.0;
  if (nDgts>1) { // we break ciphetext into a few digits when key-switching
    double dsize = context.logOfProduct(context.ctxtPrimes)/nDgts; // estimate
    double target = dsize-(context.bitsPerLevel/3.0);
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
    if (!empty(s)) {
      context.digits[nDgts-1] = s;
      double thisDigitSize = context.logOfProduct(s);
      if (maxDigitSize < thisDigitSize) maxDigitSize = thisDigitSize;
    }
    else { // If last digit is empty, remove it
      nDgts--;
      context.digits.resize(nDgts);
    }
  }
  else { 
    maxDigitSize = context.logOfProduct(context.ctxtPrimes);
    context.digits[0] = context.ctxtPrimes;
  }

  // Add primes to the chain for the P factor of key-switching
  long p2r = (context.rcData.alMod)? context.rcData.alMod->getPPowR()
                                   : context.alMod.getPPowR();
  double sizeOfSpecialPrimes
    = maxDigitSize + log(nDgts/32.0)/2 + log(context.stdev *2)
      + log((double)p2r);

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

  if (rcData != other.rcData) return false;
  return true;
}

void writeContextBase(ostream& str, const FHEcontext& context)
{
  str << "[" << context.zMStar.getM()
      << " " << context.zMStar.getP()
      << " " << context.alMod.getR()
      << " [";
  for (long i=0; i<(long) context.zMStar.numOfGens(); i++) {
    str << context.zMStar.ZmStarGen(i)
	<< ((i==(long)context.zMStar.numOfGens()-1)? "]" : " ");
  }
  str << " [";
  for (long i=0; i<(long) context.zMStar.numOfGens(); i++) {
    str << context.zMStar.OrderOf(i)
	<< ((i==(long)context.zMStar.numOfGens()-1)? "]" : " ");
  }
  str << "]";
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

  str << context.rcData.mvec;
  str << " " << context.rcData.hwt;
  str << " " << context.rcData.conservative;

  str << "]\n";

  return str;
}

void readContextBase(istream& str, unsigned long& m, unsigned long& p,
		     unsigned long& r, vector<long>& gens, vector<long>& ords)
{
  // Advance str beyond first '[' 
  seekPastChar(str, '[');  // this function is defined in NumbTh.cpp

  str >> m >> p >> r;
  str >> gens;
  str >> ords;

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

    if (ALT_CRT) 
      context.moduli.push_back(Cmodulus(context.zMStar,p,1)); // a dummy object
      // FIXME: this is broken...we are not getting an FFT prime here
    else
      context.moduli.push_back(Cmodulus(context.zMStar,p,0)); // a real object

    if (s.contains(i))
      context.specialPrimes.insert(i); // special prime
    else
      context.ctxtPrimes.insert(i);    // ciphertext prime
  }

  // read in the partition to digits
  long nDigits;
  str >> nDigits;
  context.digits.resize(nDigits);
  for (long i=0; i<(long)context.digits.size(); i++)
    str >> context.digits[i];

  // Read in the partition of m into co-prime factors (if bootstrappable)
  Vec<long> mv;
  long t;
  bool consFlag;
  str >> mv;
  str >> t;
  str >> consFlag;
  if (mv.length()>0) {
    context.makeBootstrappable(mv, t, consFlag);
  }

  seekPastChar(str, ']');
  return str;
}

#include "EncryptedArray.h"
FHEcontext::~FHEcontext()
{
  delete ea;
}

// Constructors must ensure that alMod points to zMStar, and
// rcEA (if set) points to rcAlmod which points to zMStar
FHEcontext::FHEcontext(unsigned long m, unsigned long p, unsigned long r,
   const vector<long>& gens, const vector<long>& ords):
  zMStar(m, p, gens, ords), alMod(zMStar, r),
  ea(new EncryptedArray(*this, alMod))
{
  stdev=3.2;  
  bitsPerLevel = FHE_pSize;
  fftPrimeCount = 0; 

  lazy = ALT_CRT && 
    NextPowerOfTwo(zMStar.getM()) == NextPowerOfTwo(zMStar.getPhiM());
  // we only set the lazy flag if we are using ALT_CRT and if the size of
  // NTL's FFTs for m and phi(m) are the same. If NTL didn't have these
  // power-of-two jumps, we would possibly want to change this.
}
