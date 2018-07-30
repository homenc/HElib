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
 * @file primeChain.cpp
 * @brief handling the chain of moduli
 */
#include <climits>
#include <algorithm>
#include "primeChain.h"
#include "FHEContext.h"
#include "sample.h"
#include "binio.h"

NTL_CLIENT

inline bool
operator>(const ModuliSizes::Entry& a, const ModuliSizes::Entry& b)
{ return a.first>b.first; }

ostream& operator<<(ostream& s, const ModuliSizes::Entry& e)
{
  return s << '['<< e.first << ' ' << e.second << "]\n";
}
istream& operator>>(istream& s, ModuliSizes::Entry& e)
{
  seekPastChar(s,'['); // defined in NumbTh.cpp
  s >> e.first;
  s >> e.second;
  seekPastChar(s,']');
  return s;
}
void write(ostream& s, const ModuliSizes::Entry& e)
{
  write_raw_double(s, e.first);
  e.second.write(s);
}

void read(istream& s, ModuliSizes::Entry& e)
{
  e.first = read_raw_double(s);
  e.second.read(s);
}

// initialize helper table for a given chain
void ModuliSizes::init(const std::vector<Cmodulus>& chain,
                       const IndexSet& ctxtPrimes, const IndexSet& smallPrimes)
{
  long n = (1L<<smallPrimes.card()) * ctxtPrimes.card();
  sizes.reserve(n); // allocate space
  // each entry of sizes is a pair<double,IndexSet>=(size, set-of-primes)

  // Get all subsets of smallPrimes

  sizes.push_back(make_pair(0.0,IndexSet::emptySet())); // the empty set
  long idx = 1;                      // first index that's still not set

  for (long i=smallPrimes.first();   // add i to all sets upto idx-1
       i <= smallPrimes.last(); i = smallPrimes.next(i)) {
    double sizeOfQi = log(chain[i].getQ());
    for (long j=idx; j<2*idx; j++) {
      sizes.push_back(sizes[j-idx]); // make a copy
      sizes[j].first += sizeOfQi; // add sizeOfQi to size
      sizes[j].second.insert(i);  // add i to the set of primes
    }
    idx *= 2;
  }

  // For every i in ctxtPrimes, make a copy of
  // the above plus the interval [ctxt.first, i]

  IndexSet s; // empty set
  double intervalSize = 0.0;
  for (long i=ctxtPrimes.first(); // add i to all sets upto idx-1
       i <= ctxtPrimes.last(); i = ctxtPrimes.next(i)) {
    s.insert(i);                          // add prime to the interval
    intervalSize += log(chain[i].getQ()); // add its size to intervalSize
    for (long j=0; j<idx; j++) {
      sizes.push_back(sizes[j]); // make a copy
      long n = sizes.size()-1;
      sizes[n].first += intervalSize; // add size
      sizes[n].second.insert(s);      // add interval
    }
  }

  // Finally, sort the 'sizes' array
  std::sort(sizes.begin(), sizes.end());
}

// Find a suitable IndexSet of primes whose total size is in the
// target interval [low,high], trying to minimize the number of
// primes dropped from fromSet.
// If no IndexSet exsists that fits in the target interval, returns
// the IndexSet that gives the largest value smaller than low.
IndexSet ModuliSizes::getSet4Size(double low, double high,
                                  const IndexSet& fromSet)
{
  high += 0.00001; // to compensate for rounding errors
  // lower_bound returns an iterator to the first element with size>=low
  auto it = std::lower_bound(sizes.begin(), sizes.end(),
                             Entry(low, IndexSet::emptySet()));
  long idx = it - sizes.begin(); // The index of this element

  long bestOption = -1;
  long bestCost = LONG_MAX;
  for (long i=idx; i<lsize(sizes) && sizes[i].first <= high; i++) {
    long setDiffSize = empty(fromSet)?
      card(sizes[i].second) : card(fromSet / sizes[i].second);
    if (setDiffSize <= bestCost) {
      bestOption = i;
      bestCost = setDiffSize;
    }
  }

  // If nothing was found, use the closest thing below 'low'
  if (bestOption==-1 && idx>0) {
    double lowerBound = sizes[idx-1].first - 1;
    for (long i=idx-1; i>=0 && sizes[i].first >= lowerBound; --i) {
      long setDiffSize = empty(fromSet)?
        card(sizes[i].second) : card(fromSet / sizes[i].second);
      if (setDiffSize < bestCost) {
        bestOption = i;
        bestCost = setDiffSize;
      }
    }
  }
  assert(bestOption != -1); // make sure that soemthing was found

  return sizes[bestOption].second; // return the best IndexSet
}

//! Find a suitable IndexSet of primes whose total size is in the
//! target interval [low,high], trying to minimize the total number
//! of primes dropped from both from1, from2.
//! If no IndexSet exsists that fits in the target interval, returns
//! the IndexSet that gives the largest value smaller than low.
IndexSet ModuliSizes::getSet4Size(double low, double high,
                                  const IndexSet& from1, const IndexSet& from2)
{
  high += 0.00001; // to compensate for rounding errors
  // lower_bound returns an iterator to the first element with size>=low
  auto it = std::lower_bound(sizes.begin(), sizes.end(),
                             Entry(low, IndexSet::emptySet()));
  long idx = it - sizes.begin(); // The index of this element

  long bestOption = -1;
  long bestCost = LONG_MAX;
  for (long i=idx; i<lsize(sizes) && sizes[i].first <= high; i++) {
    long setDiffSize =
      card(from1 / sizes[i].second) + card(from2 / sizes[i].second);
    if (setDiffSize <= bestCost) {
      bestOption = i;
      bestCost = setDiffSize;
    }
  }

  // If nothing was found, use the closest thing below 'low'
  if (bestOption==-1 && idx>0) {
    double lowerBound = sizes[idx-1].first - 1;
    for (long i=idx-1; i>=0 && sizes[i].first >= lowerBound; --i) {
      long setDiffSize =
        card(from1 / sizes[i].second) + card(from2 / sizes[i].second);
      if (setDiffSize < bestCost) {
        bestOption = i;
        bestCost = setDiffSize;
      }
    }
  }
  assert(bestOption != -1); // make sure that soemthing was found

  return sizes[bestOption].second; // return the best IndexSet
}

ostream& operator<<(ostream& s, const ModuliSizes& szs)
{
  return s <<'['<< szs.sizes.size()<<' '<<szs.sizes<<']';
}
istream& operator>>(istream& s, ModuliSizes& szs)
{
  long n;
  seekPastChar(s,'['); // defined in NumbTh.cpp
  s >> n;
  szs.sizes.resize(n); // allocate space
  for (long i=0; i<n; i++)
    s >> szs.sizes[i];
  seekPastChar(s,']');
  return s;
}

void ModuliSizes::write(ostream& str) const
{
  write_raw_int(str, lsize(sizes));
  for (long i=0; i<lsize(sizes); i++)
    ::write(str, sizes[i]);
}

void ModuliSizes::read(istream& str)
{
  long n = read_raw_int(str);
  sizes.resize(n); // allocate space
  for (long i=0; i<n; i++)
    ::read(str, sizes[i]);
}


// Find the next prime and add it to the chain
long FHEcontext::AddPrime(long initialP, long delta, IndexSet &s)
{
  long p = initialP;
  do { p += delta; } // delta could be positive or negative
  while (p>initialP/16 && p<NTL_SP_BOUND && !(ProbPrime(p) && !inChain(p)));

  if (p<=initialP/16 || p>=NTL_SP_BOUND) return 0; // no prime found

  long i = moduli.size(); // The index of the new prime in the list
  moduli.push_back( Cmodulus(zMStar, p, 0) );
  s.insert(i);
  return p;
}

//! @brief Add small primes to get target resolution
void addSmallPrimes(FHEcontext& context, long resolution)
{
  long m = context.zMStar.getM();
  if (m<=0 || m>(1<<20))// sanity checks
    Error("AddManyPrimes: m undefined or larger than 2^20");
  // NOTE: Below we are ensured that 16m*log(m) << NTL_SP_BOUND

  if (resolution<1 || resolution>10) // set to default of 3-bit resolution
    resolution = 3;

  vector<long> sizes;
  if (NTL_SP_NBITS>=60) { // make the smallest primes 40-bit primes
    sizes.push_back(40);
    sizes.push_back(40);
  }
  else if (NTL_SP_NBITS >=50) { // make the smallest primes 35-bit primes
    sizes.push_back(35);
    sizes.push_back(35);
  }
  else { // Make the smallest ones 22-bit primes
    assert(NTL_SP_NBITS >=30);
    sizes.push_back(22);
    sizes.push_back(22);
    sizes.push_back(22);
  }

  // use sizes 60-r, 60-2r, 60-4r,... downto the sizes above
  for (long delta=resolution; NTL_SP_NBITS-delta>sizes[0]; delta*=2)
    sizes.push_back(NTL_SP_NBITS-delta);

  // Special cases: add also NTL_SP_NBITS-3*resolution,
  // and for resolution=1 also NTL_SP_NBITS-11
  if (NTL_SP_NBITS - 3*resolution > sizes[0])
    sizes.push_back(NTL_SP_NBITS- 3*resolution);
  if (resolution==1 && NTL_SP_NBITS-11 > sizes[0])
    sizes.push_back(NTL_SP_NBITS- 11);

  std::sort(sizes.begin(), sizes.end()); // order by size

  // Look for primes equal to 1 mod m*2^e for large enough e
  long e = context.zMStar.fftSizeNeeded();
  if (NTL_SP_NBITS < e + NTL::NextPowerOfTwo(m))
    e = 1;  // Sanity check: m*2^e must fit in a single-precision integer
  long m2e = m*(1L<<e);      // m times 2^e

  for (long sz : sizes) {
    long top = 1L << sz;
    long initial = top - (top % m2e) +1; // Try p= 1 mod m*2^e

    // ensure initial >= top, since addPrime will subtract m2e
    if (initial <= top) initial += m2e;
    if (context.AddPrime(initial, m2e, context.smallPrimes)==0) {
      long twoM = m*2;   // If failed, try again with e=1
      initial = top - (top % twoM) +1;
      if (initial <= top) initial += twoM;
      if (context.AddPrime(initial, twoM, context.smallPrimes)==0)
        throw(std::logic_error("addSmallPrimes: failed to find "
              +to_string(sz)+"-bit prime =1 mod "+to_string(2*m)));
    }
  }
}

// Adds several primes to the chain. If byNumber=true then totalSize specifies
// the number of primes to add. If byNumber=false then totalSize specifies the
// target natural log all the added primes.
// Returns natural log of the product of all added primes.
double AddManyPrimes(FHEcontext& context, double totalSize, 
		     bool byNumber, bool special)
{
  if (!context.zMStar.getM() || context.zMStar.getM()>(1<<20))// sanity checks
    Error("AddManyPrimes: m undefined or larger than 2^20");
  // NOTE: Below we are ensured that 16m*log(m) << NTL_SP_BOUND

  double sizeLogSoFar = 0.0; // log of added primes so far
  double addedSoFar = 0.0;   // Either size or number, depending on 'byNumber'

#ifdef NO_HALF_SIZE_PRIME
  long sizeBits = context.bitsPerLevel;
#else
  long sizeBits = 2*context.bitsPerLevel;
#endif
  if (special) { // try to use similar size for all the special primes
    double totalBits = totalSize/log(2.0);
    long numPrimes = ceil(totalBits/NTL_SP_NBITS);// how many special primes
    sizeBits = 1+ceil(totalBits/numPrimes);       // what's the size of each
    // Added one so we don't undershoot our target
  }
  if (sizeBits>NTL_SP_NBITS) sizeBits = NTL_SP_NBITS;
  long sizeBound = 1L << sizeBits;

  // Make sure that you have enough primes such that p-1 is divisible by 2m
  long twoM = 2 * context.zMStar.getM();
  if (sizeBound < twoM*log2(twoM)*8) { // bound too small to have such primes
    sizeBits = ceil(log2(twoM*log2(twoM)))+3; // increase prime size-bound
    sizeBound = 1L << sizeBits;
  }

  // make p-1 divisible by m*2^k for as large k as possible
  // (not needed when m itself a power of two)

  if (context.zMStar.getM() & 1) // m is odd, so not power of two
    while (twoM < sizeBound/(sizeBits*2)) twoM *= 2;

  long bigP = sizeBound - (sizeBound%twoM) +1; // 1 mod 2m
  while (bigP>NTL_SP_BOUND) bigP -= twoM; // sanity check

  long p = bigP+twoM; // twoM is subtracted in the AddPrime function

  // FIXME: The last prime could sometimes be slightly smaller
  while (addedSoFar < totalSize) {
    if ((p = context.AddPrime(p,-twoM,               // found a prime
                     special? context.specialPrimes: context.ctxtPrimes))) {
      sizeLogSoFar += log((double)p);
      addedSoFar = byNumber? (addedSoFar+1.0) : sizeLogSoFar;
    }
    else { // we ran out of primes, try a lower power of two
      twoM /= 2;
      assert(twoM > (long)context.zMStar.getM()); // can we go lower?
      p = bigP;
    }
  }
  return sizeLogSoFar;
}

void buildModChain(FHEcontext &context, long nLevels, long nDgts,
                   bool willBeBootstrappable)
{
  const PAlgebra& palg = context.zMStar;
  long p = palg.getP();
  long m = palg.getM();
  long p2r = context.alMod.getPPowR();

  // Ensure bitsPerLevel is large enough to surpress high-order noise terms
  { long phim = palg.getPhiM();
    double stdev = to_double(context.stdev);
    if (palg.getPow2() == 0) // not power of two
      stdev *= sqrt(m);
    long p2e = p2r;
    if (willBeBootstrappable) { // bigger p^e for bootstrapping
      double alpha; long e, ePrime;
      RecryptData::setAlphaE(alpha,e,ePrime, context);
      p2e *= NTL::power_long(p, e-ePrime);
    }
    double dBound = std::max<double>(boundFreshNoise(m, phim, stdev),
                                     boundRoundingNoise(m, phim, p2e));
    long lBound = round(log2(dBound));
    
#ifndef NO_HALF_SIZE_PRIME
    lBound = min(lBound, NTL_SP_NBITS/2);
#endif

    if (context.bitsPerLevel < lBound) {
      cerr << "buildModChain: context.bitsPerLevel upped from "
           << context.bitsPerLevel<<" to "<<lBound<< endl;
      context.bitsPerLevel = lBound;
    }
  }

#ifdef NO_HALF_SIZE_PRIME
  long nPrimes = nLevels;
#else
  long nPrimes = (nLevels+1)/2;
  // The first prime should be of half the size. The code below tries to find
  // a prime q0 of this size where q0-1 is divisible by 2^k * m for some k>1.

  long twoM = 2 * m;
  long bound = (1L << (context.bitsPerLevel-1));
  while (twoM < bound/(2*context.bitsPerLevel))
    twoM *= 2; // divisible by 2^k * m  for a larger k

  bound = bound - (bound % twoM) +1; // = 1 mod 2m
  long q0 = context.AddPrime(bound, twoM, context.ctxtPrimes); 
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

    // A hack: we break the current digit after the total size of all digits
    // so far "almost reaches" the next multiple of dsize, upto 1/3 of a level
    double target = dsize-(context.bitsPerLevel/3.0);
    long idx = context.ctxtPrimes.first();
    for (long i=0; i<nDgts-1; i++) { // set all digits but the last
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
    // The ctxt primes that are left (if any) form the last digit
    IndexSet s = context.ctxtPrimes / s1;
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
  else { // only one digit
    maxDigitSize = context.logOfProduct(context.ctxtPrimes);
    context.digits[0] = context.ctxtPrimes;
  }

  // Add special primes to the chain for the P factor of key-switching
  double sizeOfSpecialPrimes
    = maxDigitSize + log(nDgts) + log(context.stdev *2) + log((double)p2r);

  if (willBeBootstrappable)
    sizeOfSpecialPrimes += 8*log(2.0);
  // FIXME: replace 8.0 by some way of computing the real number that's needed

  AddPrimesBySize(context, sizeOfSpecialPrimes, true);
}

