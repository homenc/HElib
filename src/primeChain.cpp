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
#include <cmath>
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

  for (long i: smallPrimes) {   // add i to all sets upto idx-1
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
  for (long i: ctxtPrimes) { // add i to all sets upto idx-1
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
// If no IndexSet exists that fits in the target interval, returns
// the IndexSet that gives the largest value smaller than low, or
// else just the singleton containing the smallest prime.
IndexSet ModuliSizes::getSet4Size(double low, double high,
                                  const IndexSet& fromSet,
                                  bool reverse) const
{
  long n = sizes.size();

  // lower_bound returns an iterator to the first element with size>=low
  auto it = std::lower_bound(sizes.begin(), sizes.end(),
                             Entry(low, IndexSet::emptySet()));
  long idx = it - sizes.begin(); // The index of this element

  long bestOption = -1;
  long bestCost = LONG_MAX;
  long ii = idx;
  for (; ii < n && sizes[ii].first <= high; ii++) {
    long setDiffSize = card(fromSet / sizes[ii].second);
    if (setDiffSize <= bestCost) {
      bestOption = ii;
      bestCost = setDiffSize;
    }
  }

  // If nothing was found, use the closest set below 'low' (or
  // above 'high' if reverse).  We actually have one bit of slack, 
  // examining not the just the closest set, but those sets whose
  // size is within 1 bit of the closest.

  if (bestOption == -1) {
    if (reverse) {
      if (ii < n) {
	double upperBound = sizes[ii].first + 1.0*log(2.0);
	for (long i=ii; i < n && sizes[i].first <= upperBound; ++i) {
	  long setDiffSize = card(fromSet / sizes[i].second);
	  if (setDiffSize < bestCost) {
	    bestOption = i;
	    bestCost = setDiffSize;
	  }
	}
      }
    }
    else {
      if (idx>0) {
	double lowerBound = sizes[idx-1].first - 1.0*log(2.0);
	for (long i=idx-1; i>=0 && sizes[i].first >= lowerBound; --i) {
	  long setDiffSize = card(fromSet / sizes[i].second);
	  if (setDiffSize < bestCost) {
	    bestOption = i;
	    bestCost = setDiffSize;
	  }
	}
      }
    }
  }

  // Nothing was found. This almost surely means decryption
  // error, but we'll just display a warning and carry on.
  if (bestOption == -1) {
    Warning(__func__+std::string(": no matching IndexSet found, likely decryption error"));
    return IndexSet(0);
  }

  return sizes[bestOption].second; // return the best IndexSet
}

//! Find a suitable IndexSet of primes whose total size is in the
//! target interval [low,high], trying to minimize the total number
//! of primes dropped from both from1, from2.
//! If no IndexSet exsists that fits in the target interval, returns
//! the IndexSet that gives the largest value smaller than low.
IndexSet ModuliSizes::getSet4Size(double low, double high,
                                  const IndexSet& from1, const IndexSet& from2,
                                  bool reverse) const
{
  long n = sizes.size();

  // lower_bound returns an iterator to the first element with size>=low
  auto it = std::lower_bound(sizes.begin(), sizes.end(),
                             Entry(low, IndexSet::emptySet()));
  long idx = it - sizes.begin(); // The index of this element

  long bestOption = -1;
  long bestCost = LONG_MAX;
  long ii = idx;
  for (; ii < n && sizes[ii].first <= high; ii++) {
    long setDiffSize = card(from1 / sizes[ii].second) + 
                       card(from2 / sizes[ii].second);
    if (setDiffSize <= bestCost) {
      bestOption = ii;
      bestCost = setDiffSize;
    }
  }

  // If nothing was found, use the closest set below 'low'
  // (or above 'high' if reverse).  We actually one bit of slack,
  // examining the not just the closest set, but those sets
  // whose size is within 1 bit of the closest.
  
  if (bestOption == -1) {
    if (reverse) {
      if (ii < n) {
	double upperBound = sizes[ii].first + 1.0*log(2.0);
	for (long i=ii; i < n && sizes[i].first <= upperBound; ++i) {
	  long setDiffSize = card(from1 / sizes[i].second) + 
			     card(from2 / sizes[i].second);
	  if (setDiffSize < bestCost) {
	    bestOption = i;
	    bestCost = setDiffSize;
	  }
	}
      }
    }
    else {
      if (idx>0) {
	double lowerBound = sizes[idx-1].first - 1.0*log(2.0);
	for (long i=idx-1; i>=0 && sizes[i].first >= lowerBound; --i) {
	  long setDiffSize = card(from1 / sizes[i].second) + 
			     card(from2 / sizes[i].second);
	  if (setDiffSize < bestCost) {
	    bestOption = i;
	    bestCost = setDiffSize;
	  }
	}
      }
    }
  }

  // Nothing was found. This almost surely means decryption
  // error, but we'll just display a warning and carry on.
  if (bestOption == -1) {
    Warning(__func__+std::string(": no matching IndexSet found, likely decryption error"));
    return IndexSet(0);
  }

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


// You initialize a PrimeGenerator as follows:
//    PrimeGenerator gen(len, m);
// Each call to gen.next() generates a prime p with 
// (3/4)*2^len <= p < 2^len and p = 2^k*t*m + 1,
// where t is odd and k is as large as possible.
// If no such prime is found, then an error is raised.

struct PrimeGenerator {
  long len, m;
  long k, t;

  PrimeGenerator(long _len, long _m) : len(_len), m(_m)
  {
    //OLD: if (len > NTL_SP_NBITS || len < 2 || m >= NTL_SP_BOUND || m <= 0)
    //OLD: throw helib::InvalidArgument("PrimeGenerator: bad args");
    helib::assertInRange<helib::InvalidArgument>(len, 2l, (long)NTL_SP_NBITS, "PrimeGenerator: len is not in [2, NTL_SP_NBITS]", true);
    helib::assertInRange<helib::InvalidArgument>(m, 1l, (long)NTL_SP_BOUND, "PrimeGenerator: m is not in [1, NTL_SP_BOUND)");

    // compute k as smallest nonnegative integer such that
    // 2^{len-2} < 2^k*m
    k = 0;
    while ((m << k) <= (1L << (len-2))) k++;

    t = 8; // with above setting for k, we have 2^{len-1}/(2^k*m) < 4,
           // so setting t = 8 will trigger a new k-value with the
           // first call to next()
  }

  long next()
  {
    // we consider all odd t in the interval 
    // [ ((3/4)*2^len-1)/(2^k*m), (2^len-1)/(2^k*m) ).
    // For k satisfyng 2^{len-2} >= 2^k*m, this interval is
    // non-empty.
    // It is equivalent to consider the interval
    // of integers [tlb, tub), where tlb = ceil(((3/4)*2^len-1)/(2^k*m))
    // and tub = ceil((2^len-1)/(2^k*m)).

    long tub = divc((1L << len)-1, m << k);

    for (;;) {

      t++;

      if (t >= tub) {
	// move to smaller value of k, reset t and tub
   
	k--;

	long klb;
	if (m%2 == 0) 
	  klb = 0;
	else
	  klb = 1;

  if (k < klb) throw helib::RuntimeError("Prime generator ran out of primes");
        
	// we run k down to 0  if m is even, and down to 1
	// if m is odd.

	t = divc(3*(1L << (len-2))-1, m << k);
	tub = divc((1L << len)-1, m << k);
      }

      if (t%2 == 0) continue; // we only want to consider odd t

      long cand = ((t*m) << k) + 1; // = 2^k*t*m + 1

      // double check that cand is in the prescribed interval
      //OLD: assert(cand >= (1L << (len-2))*3 && cand < (1L << len));
      helib::assertInRange(cand, (1L << (len-2))*3, 1L << len, "Candidate cand is not in the prescribed interval");

      if (ProbPrime(cand, 60)) return cand;
      // iteration count == 60 implies 2^{-120} error probability
    }

  }

};

void FHEcontext::AddSmallPrime(long q)
{
  //OLD: assert(!inChain(q));
  helib::assertFalse(inChain(q), "Small prime q is already in the prime chain");
  long i = moduli.size(); // The index of the new prime in the list
  moduli.push_back( Cmodulus(zMStar, q, 0) );
  smallPrimes.insert(i);
}

void FHEcontext::AddCtxtPrime(long q)
{
  //OLD: assert(!inChain(q));
  helib::assertFalse(inChain(q), "Prime q is already in the prime chain");
  long i = moduli.size(); // The index of the new prime in the list
  moduli.push_back( Cmodulus(zMStar, q, 0) );
  ctxtPrimes.insert(i);
}

void FHEcontext::AddSpecialPrime(long q)
{
  //OLD: assert(!inChain(q));
  helib::assertFalse(inChain(q), "Special prime q is already in the prime chain");
  long i = moduli.size(); // The index of the new prime in the list
  moduli.push_back( Cmodulus(zMStar, q, 0) );
  specialPrimes.insert(i);
}

//! @brief Add small primes to get target resolution
void addSmallPrimes(FHEcontext& context, long resolution, long cpSize)
{
  // cpSize is the size of the ciphertext primes
  // Sanity-checks, cpSize \in [0.9*NTL_SP_NBITS, NTL_SP_NBITS]
  //OLD: assert((cpSize >= 30) && (cpSize <= NTL_SP_NBITS) && (cpSize*10 >= NTL_SP_NBITS*9));
  helib::assertTrue(cpSize >= 30, "cpSize is too small (minimum is 30)");
  helib::assertInRange(cpSize * 10, 9l * NTL_SP_NBITS, 10l * NTL_SP_NBITS, "cpSize not in [0.9*NTL_SP_NBITS, NTL_SP_NBITS]", true);
  
  long m = context.zMStar.getM();
  if (m<=0 || m>(1<<20))// sanity checks
    throw helib::RuntimeError("addSmallPrimes: m undefined or larger than 2^20");
  // NOTE: Below we are ensured that 16m*log(m) << NTL_SP_BOUND

  if (resolution<1 || resolution>10) // set to default of 3-bit resolution
    resolution = 3;

  vector<long> sizes;
  long smallest; // size of the smallest of the smallPrimes
  // We need at least two of this size, maybe three

  if (cpSize>=54)
    smallest = divc(2*cpSize,3);
  else if (cpSize >=45)
    smallest = divc(7*cpSize,10);
  else { // Make the smallest ones at least 22-bit primes
    //OLD: assert(cpSize >=30);
    // Repeted assertion
    smallest = divc(11*cpSize,15);
    sizes.push_back(smallest); // need three of them
  }
  sizes.push_back(smallest);
  sizes.push_back(smallest);

  // This ensures we can express everything to given resolution.

  // use sizes cpSize-r, cpSize-2r, cpSize-4r,... downto the sizes above
  for (long delta=resolution; cpSize-delta>smallest; delta*=2)
    sizes.push_back(cpSize-delta);

  // This helps to minimize the number of small primes needed
  // to express any particular resolution.
  // This could be removed...need to experiment.

  // Special cases: add also cpSize-3*resolution,
  // and for resolution=1 also cpSize-11
  if (cpSize - 3*resolution > smallest)
    sizes.push_back(cpSize- 3*resolution);
  if (resolution==1 && cpSize-11 > smallest)
    sizes.push_back(cpSize- 11);

  std::sort(sizes.begin(), sizes.end()); // order by size

  long last_sz = 0;
  long sz_cnt = 0;
  std::unique_ptr<PrimeGenerator> gen;
  for (long sz : sizes) {
    if (sz != last_sz) gen.reset(new PrimeGenerator(sz, m));
    long q = gen->next();
    context.AddSmallPrime(q);
    last_sz = sz;
  }
}

// Determine the target size of the ctxtPrimes. The target size is
// set at 2^n, where n is at most NTL_SP_NBITS and at least least
// ceil(0.9*NTL_SP_NBITS), so that we don't overshoot nBits by too
// much.
// The reason that we do not allow to go below 0.9*NTL_SP_NBITS is
// that we need some of the smallPrimes to be sufficiently smaller
// than the ctxtPrimes, and still we need these smallPrimes to have
// m'th roots of unity.
long ctxtPrimeSize(long nBits)
{
  // How many primes of size NTL_SP_NBITS it takes to get to nBits
  double maxPsize = NTL_SP_NBITS - 0.5;
  long nPrimes = long(ceil(nBits/maxPsize));
  long targetSize = divc(nBits, nPrimes)+1; // How large should each prime be
       // The '-0.5','+1' are because we get primes
       // that are little smaller than 2^{targetSize}

  // ensure targetSize \in [0.9*NTL_SP_NBITS, NTL_SP_NBITS]
  if (targetSize < 0.9*NTL_SP_NBITS)
    targetSize = divc(9*NTL_SP_NBITS, 10);
  if (targetSize > NTL_SP_NBITS)
    targetSize = NTL_SP_NBITS;
  if (targetSize < 30)
    targetSize = 30;

  return targetSize;
}

void addCtxtPrimes(FHEcontext& context, long nBits, long targetSize)
{
  // We add enough primes of size targetSize until their product is
  // at least 2^{nBits}

  // Sanity-checks, targetSize \in [0.9*NTL_SP_NBITS, NTL_SP_NBITS]
  //OLD: assert((targetSize >= 30) && (targetSize <= NTL_SP_NBITS) && (targetSize*10 >= NTL_SP_NBITS*9));
  helib::assertTrue(targetSize >= 30, "Target prime is too small (minimum size is 30)");
  helib::assertInRange(targetSize * 10, 9l * NTL_SP_NBITS, 10l * NTL_SP_NBITS, "targetSize not in [0.9*NTL_SP_NBITS, NTL_SP_NBITS]", true);
  const PAlgebra& palg = context.zMStar;
  long m = palg.getM();

  PrimeGenerator gen(targetSize, m);
  double bitlen = 0;     // how many bits we already have
  while (bitlen < nBits-0.5) {
    long q = gen.next();     // generate the next prime
    context.AddCtxtPrime(q); // add it to the list
    bitlen += log2(q);
  }
}


void addSpecialPrimes(FHEcontext& context, long nDgts, 
                      bool willBeBootstrappable)
{
  const PAlgebra& palg = context.zMStar;
  long p = palg.getP();
  long m = palg.getM();
  long p2r = context.alMod.getPPowR();

  long p2e = p2r;
  if (willBeBootstrappable) { // bigger p^e for bootstrapping
    long e, ePrime, a;
    RecryptData::setAE(a,e,ePrime, context);
    p2e *= NTL::power_long(p, e-ePrime);
  }

  long nCtxtPrimes = context.ctxtPrimes.card();
  if (nDgts > nCtxtPrimes) nDgts = nCtxtPrimes; // sanity checks
  if (nDgts <= 0) nDgts = 1;

  context.digits.resize(nDgts); // allocate space

  if (nDgts>1) { // we break ciphertext into a few digits when key-switching
    // NOTE: The code below assumes that all the ctxtPrimes have roughly the same size

    IndexSet remaining = context.ctxtPrimes;
    for (long dgt=0; dgt<nDgts-1; dgt++) {
      long digitCard = divc(remaining.card(), nDgts-dgt);
           // ceiling(#-of-remaining-primes, #-or-remaining-digits)

      for (long i : remaining) {
	context.digits[dgt].insert(i);
        if (context.digits[dgt].card() >= digitCard) break;
      }
      remaining.remove(context.digits[dgt]); // update the remaining set
    }
    // The last digit has everything else
    if (empty(remaining)) { // sanity chack, use one less digit
      nDgts--;
      context.digits.resize(nDgts);
    }
    else
      context.digits[nDgts-1] = remaining;
  }
  else { // only one digit
    context.digits[0] = context.ctxtPrimes;
  }

  double maxDigitLog = 0.0;
  for (auto& digit : context.digits) {
    double size = context.logOfProduct(digit);
    if (size > maxDigitLog) maxDigitLog = size;
  }

  // Add special primes to the chain for the P factor of key-switching
  double logOfSpecialPrimes
    = maxDigitLog + log(nDgts) + log(context.stdev *2) + log(p2e);

  // we now add enough special primes so that the sum of their
  // logs is at least logOfSpecial primes

  // we first calculate nbits, which is the bit length of each
  // special prime.  This is calculated so that we don't overshoot
  // logOfSpecial primes by too much because of granularity

  double totalBits = logOfSpecialPrimes/log(2.0);
  long numPrimes = ceil(totalBits/NTL_SP_NBITS);  
  // initial estimate # of special primes
  long nbits = ceil(totalBits/numPrimes);         
  // estimated size of each special prime

  nbits++;
  // add 1 so we don't undershoot 

  if (nbits > NTL_SP_NBITS) nbits = NTL_SP_NBITS;
  // make sure nbits not too large

  // now add special primes of size nbits

  PrimeGenerator gen(nbits, m);

  double logSoFar = 0.0;
  while (logSoFar < logOfSpecialPrimes) {
    long q = gen.next();

    if (context.inChain(q)) continue;
    // nbits could equal NTL_SP_BITS or the size of one 
    // of the small primes, so we have to check for duplicates here...
    // this is not the most efficient way to do this,
    // but it doesn't make sense to optimize this any further

    context.AddSpecialPrime(q);
    logSoFar += log(q);
  }

  //cerr << "****** special primes: " << logOfSpecialPrimes << " " << context.logOfProduct(context.specialPrimes) << "\n";
}

void buildModChain(FHEcontext& context, long nBits, long nDgts,
                      bool willBeBootstrappable, long resolution)
{
  long pSize = ctxtPrimeSize(nBits);
  addSmallPrimes(context, resolution, pSize);
  addCtxtPrimes(context, nBits, pSize);
  addSpecialPrimes(context, nDgts, willBeBootstrappable);
  context.setModSizeTable();
}
