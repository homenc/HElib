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
/**
 * @file primeChain.cpp
 * @brief handling the chain of moduli
 */
#include <climits>
#include <cmath>
#include <algorithm>
#include <helib/primeChain.h>
#include <helib/Context.h>
#include <helib/sample.h>
#include <helib/binio.h>
#include <helib/fhe_stats.h>
#include <helib/log.h>

namespace helib {

inline bool operator>(const ModuliSizes::Entry& a, const ModuliSizes::Entry& b)
{
  return a.first > b.first;
}

std::ostream& operator<<(std::ostream& s, const ModuliSizes::Entry& e)
{
  return s << '[' << e.first << ' ' << e.second << "]\n";
}
std::istream& operator>>(std::istream& s, ModuliSizes::Entry& e)
{
  seekPastChar(s, '['); // defined in NumbTh.cpp
  s >> e.first;
  s >> e.second;
  seekPastChar(s, ']');
  return s;
}
void write(std::ostream& s, const ModuliSizes::Entry& e)
{
  write_raw_double(s, e.first);
  e.second.write(s);
}

void read(std::istream& s, ModuliSizes::Entry& e)
{
  e.first = read_raw_double(s);
  e.second.read(s);
}

// initialize helper table for a given chain
void ModuliSizes::init(const Context& context)
{
  if (context.zMStar.getPow2())
    iFFT_cost = 0; // iFFT cost is same as FFT cost
  else
    iFFT_cost = 20; // iFFT cost is 1.20 times FFT cost.
                    // This is heuristic and derived from
                    // experimantal data.
                    // FIXME: should allow use to override
                    // this default value.

  long n = (1L << context.smallPrimes.card()) * (context.ctxtPrimes.card() + 1);
  sizes.reserve(n); // allocate space
  // each entry of sizes is a pair<double,IndexSet>=(size, set-of-primes)

  // Get all subsets of smallPrimes

  sizes.push_back(std::make_pair(0.0, IndexSet::emptySet())); // the empty set
  long idx = 1; // first index that's still not set

  for (long i : context.smallPrimes) { // add i to all sets upto idx-1
    double sizeOfQi = std::log(context.ithModulus(i).getQ());
    for (long j = idx; j < 2 * idx; j++) {
      sizes.push_back(sizes[j - idx]); // make a copy
      sizes[j].first += sizeOfQi;      // add sizeOfQi to size
      sizes[j].second.insert(i);       // add i to the set of primes
    }
    idx *= 2;
  }

  // For every i in ctxtPrimes, make a copy of
  // the above plus the interval [ctxt.first, i]

  IndexSet s; // empty set
  double intervalSize = 0.0;
  for (long i : context.ctxtPrimes) { // add i to all sets upto idx-1
    s.insert(i);                      // add prime to the interval
    intervalSize +=
        std::log(context.ithModulus(i).getQ()); // add its size to intervalSize
    for (long j = 0; j < idx; j++) {
      sizes.push_back(sizes[j]); // make a copy
      long n = sizes.size() - 1;
      sizes[n].first += intervalSize; // add size
      sizes[n].second.insert(s);      // add interval
    }
  }

  // Finally, sort the 'sizes' array
  std::sort(sizes.begin(), sizes.end());

#if 0
  // Useful for debug.
  std::cout << "\n*** sizes:" << sizes.size() << "\n";
  for (auto size : sizes) {
    std::cout << size.second << "\n";
  }
  std::cout << "\n";
#endif
}

// If the estimated cost (in terms of # of FFTs) in
// going from fromSet to toSet is card(fromSet) + C,
// this function returns the value 100*C.
//
// This uses the value iFFT_cost, which is normally set to 0
// when m is a power of two, and to 20 othewise.
// The idea is that the cost of an iFFT is estimated as
// the cost of an FFT times (1 + iFFT_cost/100).
//
// If toSet = fromSet + addSet - removeSet, then after scaling up
// to fromSet + addSet, the cost of removing removeSet via mod down is:
//    removeSet*(cost of iFFT) + toSet*(cost of FFT)
// =  (fromSet + addSet + (iFFT_cost/100)*removeSet)*(cost of FFT)

#if 1
// new method, based on above. It doesn't seem to make a big
// difference, however.
// It's not clear if this is really the best heurstic,
// but it seems to be at least as good as the old one.
// One could imagine that it may be best to measure
// just the size of addSet, as this will lead to more
// expensive operations in the future.
static inline long cost_estimate(const IndexSet& fromSet,
                                 const IndexSet& toSet,
                                 long iFFT_cost)
{
  // NOTE: addSet = toSet / fromSet, removeSet = fromSet / toSet.
  if (iFFT_cost == 0) {
    return 100 * card(toSet / fromSet);
  } else {
    return 100 * card(toSet / fromSet) + iFFT_cost * card(fromSet / toSet);
  }
}
#else
// old method, which counts just the primes in removeSet
static inline long cost_estimate(const IndexSet& fromSet,
                                 const IndexSet& toSet,
                                 long iFFT_cost)
{
  return 100 * card(fromSet / toSet);
}

#endif

// Find a suitable IndexSet of primes whose total size is in the
// target interval [low,high], trying to minimize the number of
// primes dropped from fromSet.
// If no IndexSet exists that fits in the target interval, returns
// the IndexSet that gives the largest value smaller than low, or
// else just the singleton containing the smallest prime.
IndexSet ModuliSizes::getSet4Size(double low,
                                  double high,
                                  const IndexSet& fromSet,
                                  bool reverse) const
{
  long n = sizes.size();

  // lower_bound returns an iterator to the first element with size>=low
  auto it = std::lower_bound(sizes.begin(),
                             sizes.end(),
                             Entry(low, IndexSet::emptySet()));
  long idx = it - sizes.begin(); // The index of this element

  long nchoices = 0;

  long bestOption = -1;
  long bestCost = LONG_MAX;
  long ii = idx;
  for (; ii < n && sizes[ii].first <= high; ii++) {
    nchoices++;

    long thisCost = cost_estimate(fromSet, sizes[ii].second, iFFT_cost);
    if (thisCost <= bestCost) {
      bestOption = ii;
      bestCost = thisCost;
    }
  }

  HELIB_STATS_UPDATE("window1-out", (bestOption == -1));
  HELIB_STATS_UPDATE("window1-nchoices", nchoices);

  // If nothing was found, use the closest set below 'low' (or
  // above 'high' if reverse).  We actually have one bit of slack,
  // examining not just the closest set, but those sets whose
  // size is within 1 bit of the closest.

  if (bestOption == -1) {
    if (reverse) {
      if (ii < n) {
        double upperBound = sizes[ii].first + 1.0 * std::log(2.0);
        for (long i = ii; i < n && sizes[i].first <= upperBound; ++i) {
          long thisCost = cost_estimate(fromSet, sizes[i].second, iFFT_cost);
          if (thisCost < bestCost) {
            bestOption = i;
            bestCost = thisCost;
          }
        }
      }
    } else {
      if (idx > 0) {
        double lowerBound = sizes[idx - 1].first - 1.0 * std::log(2.0);
        for (long i = idx - 1; i >= 0 && sizes[i].first >= lowerBound; --i) {
          long thisCost = cost_estimate(fromSet, sizes[i].second, iFFT_cost);
          if (thisCost < bestCost) {
            bestOption = i;
            bestCost = thisCost;
          }
        }
      }
    }
  }

  // Nothing was found. This almost surely means decryption
  // error, but we'll just display a warning and carry on.
  if (bestOption == -1) {
    Warning(__func__ + std::string(": no matching IndexSet found, "
                                   "likely decryption error"));
    return IndexSet(); // empty set
  }

  return sizes[bestOption].second; // return the best IndexSet
}

//! Find a suitable IndexSet of primes whose total size is in the
//! target interval [low,high], trying to minimize the total number
//! of primes dropped from both from1, from2.
//! If no IndexSet exists that fits in the target interval, returns
//! the IndexSet that gives the largest value smaller than low.
IndexSet ModuliSizes::getSet4Size(double low,
                                  double high,
                                  const IndexSet& from1,
                                  const IndexSet& from2,
                                  bool reverse) const
{
  long n = sizes.size();

  // lower_bound returns an iterator to the first element with size>=low
  auto it = std::lower_bound(sizes.begin(),
                             sizes.end(),
                             Entry(low, IndexSet::emptySet()));
  long idx = it - sizes.begin(); // The index of this element

  long nchoices = 0;

  long bestOption = -1;
  long bestCost = LONG_MAX;
  long ii = idx;
  for (; ii < n && sizes[ii].first <= high; ii++) {
    nchoices++;

    long thisCost = cost_estimate(from1, sizes[ii].second, iFFT_cost) +
                    cost_estimate(from2, sizes[ii].second, iFFT_cost);

    if (thisCost <= bestCost) {
      bestOption = ii;
      bestCost = thisCost;
    }
  }

  HELIB_STATS_UPDATE("window2-out", (bestOption == -1));
  HELIB_STATS_UPDATE("window2-nchoices", nchoices);

  // If nothing was found, use the closest set below 'low'
  // (or above 'high' if reverse).  We actually one bit of slack,
  // examining the not just the closest set, but those sets
  // whose size is within 1 bit of the closest.

  if (bestOption == -1) {
    if (reverse) {
      if (ii < n) {
        double upperBound = sizes[ii].first + 1.0 * std::log(2.0);
        for (long i = ii; i < n && sizes[i].first <= upperBound; ++i) {
          long thisCost = cost_estimate(from1, sizes[i].second, iFFT_cost) +
                          cost_estimate(from2, sizes[i].second, iFFT_cost);

          if (thisCost < bestCost) {
            bestOption = i;
            bestCost = thisCost;
          }
        }
      }
    } else {
      if (idx > 0) {
        double lowerBound = sizes[idx - 1].first - 1.0 * std::log(2.0);
        for (long i = idx - 1; i >= 0 && sizes[i].first >= lowerBound; --i) {
          long thisCost = cost_estimate(from1, sizes[i].second, iFFT_cost) +
                          cost_estimate(from2, sizes[i].second, iFFT_cost);

          if (thisCost < bestCost) {
            bestOption = i;
            bestCost = thisCost;
          }
        }
      }
    }
  }

  // Nothing was found. This almost surely means decryption
  // error, but we'll just display a warning and carry on.
  if (bestOption == -1) {
    Warning(__func__ + std::string(": no matching IndexSet found, "
                                   "likely decryption error"));
    return IndexSet(); // empty set
  }

  return sizes[bestOption].second; // return the best IndexSet
}

std::ostream& operator<<(std::ostream& s, const ModuliSizes& szs)
{
  return s << '[' << szs.sizes.size() << ' ' << szs.sizes << ']';
}
std::istream& operator>>(std::istream& s, ModuliSizes& szs)
{
  long n;
  seekPastChar(s, '['); // defined in NumbTh.cpp
  s >> n;
  szs.sizes.resize(n); // allocate space
  for (long i = 0; i < n; i++)
    s >> szs.sizes[i];
  seekPastChar(s, ']');
  return s;
}

void ModuliSizes::write(std::ostream& str) const
{
  write_raw_int(str, lsize(sizes));
  for (long i = 0; i < lsize(sizes); i++)
    ::helib::write(str, sizes[i]);
}

void ModuliSizes::read(std::istream& str)
{
  long n = read_raw_int(str);
  sizes.resize(n); // allocate space
  for (long i = 0; i < n; i++)
    ::helib::read(str, sizes[i]);
}

// You initialize a PrimeGenerator as follows:
//    PrimeGenerator gen(len, m);
// Each call to gen.next() generates a prime p with
// (1-1/2^B)*2^len <= p < 2^len and p = 2^k*t*m + 1,
// where t is odd and k is as large as possible
// and B is a small constant (typically, B in {2,3,4}).
// If no such prime is found, then an error is raised.

struct PrimeGenerator
{
  const static long B = 3;
  long len, m;
  long k, t;

  PrimeGenerator(long _len, long _m) : len(_len), m(_m)
  {
    assertInRange<InvalidArgument>(len,
                                   long(B),
                                   static_cast<long>(NTL_SP_NBITS),
                                   "PrimeGenerator: len is not "
                                   "in [B, NTL_SP_NBITS]",
                                   true);
    assertInRange<InvalidArgument>(m,
                                   1l,
                                   static_cast<long>(NTL_SP_BOUND),
                                   "PrimeGenerator: m is "
                                   "not in [1, NTL_SP_BOUND)");

    // compute k as smallest non-negative integer such that
    // 2^{len-B} < 2^k*m
    k = 0;
    while ((m << k) <= (1L << (len - B)))
      k++;

    t = divc((1L << len) - 1, m << k);
    // this ensures the fist call to next will trigger a new k-value
  }

  long next()
  {
    // we consider all odd t in the interval
    // [ (1-1/2^B)*2^len-1)/(2^k*m), (2^len-1)/(2^k*m) ).
    // For k satisfying 2^{len-B} >= 2^k*m, this interval is
    // contains at least one integer.
    // It is equivalent to consider the interval
    // of integers [tlb, tub), where tlb = ceil(((1-1/2^B)*2^len-1)/(2^k*m))
    // and tub = ceil((2^len-1)/(2^k*m)).

    long tub = divc((1L << len) - 1, m << k);

    for (;;) {

      t++;

      if (t >= tub) {
        // move to smaller value of k, reset t and tub

        k--;

        long klb;
        if (m % 2 == 0)
          klb = 0;
        else
          klb = 1;

        // we run k down to 0  if m is even, and down to 1
        // if m is odd.

        if (k < klb)
          throw RuntimeError("Prime generator ran out of primes");

        t = divc((1L << len) - (1L << (len - B)) - 1, m << k);
        tub = divc((1L << len) - 1, m << k);
      }

      if (t % 2 == 0)
        continue; // we only want to consider odd t

      long cand = ((t * m) << k) + 1; // = 2^k*t*m + 1

      // double check that cand is in the prescribed interval
      assertInRange(cand,
                    (1L << len) - (1L << (len - B)),
                    1L << len,
                    "Candidate cand is not in the prescribed interval");

      if (NTL::ProbPrime(cand, 60))
        return cand;
      // iteration count == 60 implies 2^{-120} error probability
    }
  }
};

void Context::AddSmallPrime(long q)
{
  assertFalse(inChain(q), "Small prime q is already in the prime chain");
  long i = moduli.size(); // The index of the new prime in the list
  moduli.push_back(Cmodulus(zMStar, q, 0));
  smallPrimes.insert(i);
}

void Context::AddCtxtPrime(long q)
{
  assertFalse(inChain(q), "Prime q is already in the prime chain");
  long i = moduli.size(); // The index of the new prime in the list
  moduli.push_back(Cmodulus(zMStar, q, 0));
  ctxtPrimes.insert(i);
}

void Context::AddSpecialPrime(long q)
{
  assertFalse(inChain(q), "Special prime q is already in the prime chain");
  long i = moduli.size(); // The index of the new prime in the list
  moduli.push_back(Cmodulus(zMStar, q, 0));
  specialPrimes.insert(i);
}

//! @brief Add small primes to get target resolution
// FIXME: there is some black magic here.
// we need to better document the strategy.
static void addSmallPrimes(Context& context, long resolution, long cpSize)
{
  // cpSize is the size of the ciphertext primes
  // Sanity-checks, cpSize \in [0.9*NTL_SP_NBITS, NTL_SP_NBITS]
  assertTrue(cpSize >= 30, "cpSize is too small (minimum is 30)");
  assertInRange(cpSize * 10,
                9l * NTL_SP_NBITS,
                10l * NTL_SP_NBITS,
                "cpSize not in [0.9*NTL_SP_NBITS, NTL_SP_NBITS]",
                true);

  long m = context.zMStar.getM();
  if (m <= 0 || m > (1 << 20)) // sanity checks
    throw RuntimeError("addSmallPrimes: m undefined or larger than 2^20");
  // NOTE: Below we are ensured that 16m*log(m) << NTL_SP_BOUND

  if (resolution < 1 || resolution > 10) // set to default of 3-bit resolution
    resolution = 3;

  std::vector<long> sizes;
  long smallest; // size of the smallest of the smallPrimes
  // We need at least two of this size, maybe three

  if (cpSize >= 54)
    smallest = divc(2 * cpSize, 3);
  else if (cpSize >= 45)
    smallest = divc(7 * cpSize, 10);
  else { // Make the smallest ones at least 22-bit primes
    smallest = divc(11 * cpSize, 15);
    sizes.push_back(smallest); // need three of them
  }
  sizes.push_back(smallest);
  sizes.push_back(smallest);

  // This ensures we can express everything to given resolution.

  // use sizes cpSize-r, cpSize-2r, cpSize-4r,... down to the sizes above
  for (long delta = resolution; cpSize - delta > smallest; delta *= 2)
    sizes.push_back(cpSize - delta);

  // This helps to minimize the number of small primes needed
  // to express any particular resolution.
  // This could be removed...need to experiment.

  // Special cases: add also cpSize-3*resolution,
  // and for resolution=1 also cpSize-11
  if (cpSize - 3 * resolution > smallest)
    sizes.push_back(cpSize - 3 * resolution);
  if (resolution == 1 && cpSize - 11 > smallest)
    sizes.push_back(cpSize - 11);

  std::sort(sizes.begin(), sizes.end()); // order by size

  long last_sz = 0;
  std::unique_ptr<PrimeGenerator> gen;
  for (long sz : sizes) {
    if (sz != last_sz)
      gen.reset(new PrimeGenerator(sz, m));
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
static long ctxtPrimeSize(long nBits)
{
  double bit_loss =
      -std::log1p(-1.0 / double(1L << PrimeGenerator::B)) / std::log(2.0);
  // std::cerr << "*** bit_loss=" << bit_loss;

  // How many primes of size NTL_SP_NBITS it takes to get to nBits
  double maxPsize = NTL_SP_NBITS - bit_loss;
  // primes of length len are guaranteed to be at least (1-1/2^B)*2^len,

  long nPrimes = long(ceil(nBits / maxPsize));
  // this is sufficiently many primes

  // now we want to trim the size to avoid unnecssary overshooting
  // so we decrease targetSize, while guaranteeing that
  // nPrimes primes of length targetSize multiply out to
  // at least nBits bits.

  long targetSize = NTL_SP_NBITS;
  while (10 * (targetSize - 1) >= 9 * NTL_SP_NBITS && (targetSize - 1) >= 30 &&
         ((targetSize - 1) - bit_loss) * nPrimes >= nBits)
    targetSize--;

  if (((targetSize - 1) - bit_loss) * nPrimes >= nBits)
    Warning(__func__ + std::string(": non-optimal targetSize"));

  return targetSize;
}

static void addCtxtPrimes(Context& context, long nBits, long targetSize)
{
  // We add enough primes of size targetSize until their product is
  // at least 2^{nBits}

  // Sanity-checks, targetSize \in [0.9*NTL_SP_NBITS, NTL_SP_NBITS]
  assertTrue(targetSize >= 30,
             "Target prime is too small (minimum size is 30)");
  assertInRange(targetSize * 10,
                9l * NTL_SP_NBITS,
                10l * NTL_SP_NBITS,
                "targetSize not in [0.9*NTL_SP_NBITS, NTL_SP_NBITS]",
                true);
  const PAlgebra& palg = context.zMStar;
  long m = palg.getM();

  PrimeGenerator gen(targetSize, m);
  double bitlen = 0; // how many bits we already have
  while (bitlen < nBits - 0.5) {
    long q = gen.next();     // generate the next prime
    context.AddCtxtPrime(q); // add it to the list
    bitlen += std::log2(q);
  }

  // std::cerr << "*** ctxtPrimes excess: " << (bitlen - nBits) << "\n";
  HELIB_STATS_UPDATE("excess-ctxtPrimes", bitlen - nBits);
}

static void addSpecialPrimes(Context& context,
                             long nDgts,
                             bool willBeBootstrappable,
                             UNUSED long skHwt,
                             long bitsInSpecialPrimes)
{
  const PAlgebra& palg = context.zMStar;
  long p = std::abs(palg.getP()); // for CKKS, palg.getP() == -1
  long m = palg.getM();
  long phim = palg.getPhiM();
  long p2r = context.isCKKS() ? 1 : context.alMod.getPPowR();

  long p2e = p2r;
  if (willBeBootstrappable && !context.isCKKS()) {
    // bigger p^e for bootstrapping
    long e, ePrime;
    RecryptData::setAE(e, ePrime, context);
    p2e *= NTL::power_long(p, e - ePrime);

    // initialize e and ePrime parameters in the context
    context.e_param = e;
    context.ePrime_param = ePrime;
  }

  long nCtxtPrimes = context.ctxtPrimes.card();
  if (nDgts > nCtxtPrimes)
    nDgts = nCtxtPrimes; // sanity checks
  if (nDgts <= 0)
    nDgts = 1;

  context.digits.resize(nDgts); // allocate space

  if (nDgts > 1) { // we break ciphertext into a few digits when key-switching
    // NOTE: The code below assumes that all the ctxtPrimes have roughly the
    // same size

    IndexSet remaining = context.ctxtPrimes;
    for (long dgt = 0; dgt < nDgts - 1; dgt++) {
      long digitCard = divc(remaining.card(), nDgts - dgt);
      // ceiling(#-of-remaining-primes, #-or-remaining-digits)

      for (long i : remaining) {
        context.digits[dgt].insert(i);
        if (context.digits[dgt].card() >= digitCard)
          break;
      }
      remaining.remove(context.digits[dgt]); // update the remaining set
    }
    // The last digit has everything else
    if (empty(remaining)) { // sanity check, use one less digit
      nDgts--;
      context.digits.resize(nDgts);
    } else
      context.digits[nDgts - 1] = remaining;
  } else { // only one digit
    context.digits[0] = context.ctxtPrimes;
  }

  double maxDigitLog = 0.0;
  for (auto& digit : context.digits) {
    double size = context.logOfProduct(digit);
    if (size > maxDigitLog)
      maxDigitLog = size;
  }

  // Add special primes to the chain for the P factor of key-switching
  double nBits;

  if (bitsInSpecialPrimes)
    nBits = bitsInSpecialPrimes;
  else {
#if 0
    nBits = (maxDigitLog + std::log(nDgts) + NTL::log(context.stdev * 2) +
             std::log(p2e)) /
            std::log(2.0);
    // FIXME: Victor says: the above calculation does not make much sense to me
#else
    double h;
    if (context.hwt_param == 0)
      h = phim / 2.0;
    else
      h = context.hwt_param;

    double log_phim = std::log(phim);
    if (log_phim < 1)
      log_phim = 1;

    if (context.isCKKS()) {
      // This is based on a smaller noise estimate so as
      // to better protext precision...this is based on
      // a noise level equal to the mod switch added noise.
      // Note that the relin_CKKS_adjust function in Ctxt.cpp
      // depends on this estimate.
      nBits = (maxDigitLog + NTL::log(context.stdev) + std::log(nDgts) -
               0.5 * std::log(h)) /
              std::log(2.0);
    } else if (palg.getPow2()) {
      nBits = (maxDigitLog + std::log(p2e) + NTL::log(context.stdev) +
               0.5 * std::log(12.0) + std::log(nDgts) -
               0.5 * std::log(log_phim) - 2 * std::log(p) - std::log(h)) /
              std::log(2.0);
    } else {
      nBits =
          (maxDigitLog + std::log(m) + std::log(p2e) + NTL::log(context.stdev) +
           0.5 * std::log(12.0) + std::log(nDgts) - 0.5 * log_phim -
           0.5 * std::log(log_phim) - 2 * std::log(p) - std::log(h)) /
          std::log(2.0);
    }

    // Both of the above over-estimate nBits by a factor of log2(context.scale).
    // That should provide a sufficient safety margin.
    // See design document

#endif
  }

  if (nBits < 1)
    nBits = 1;

  double bit_loss =
      -std::log1p(-1.0 / double(1L << PrimeGenerator::B)) / std::log(2.0);

  // How many primes of size NTL_SP_NBITS it takes to get to nBits
  double maxPsize = NTL_SP_NBITS - bit_loss;
  // primes of length len are guaranteed to be at least (1-1/2^B)*2^len,

  long nPrimes = long(ceil(nBits / maxPsize));
  // this is sufficiently many prime

  // now we want to trim the size to avoid unnecssary overshooting
  // so we decrease targetSize, while guaranteeing that
  // nPrimes primes of length targetSize multiply out to
  // at least nBits bits.

  long targetSize = NTL_SP_NBITS;
  while ((targetSize - 1) >= 0.55 * NTL_SP_NBITS && (targetSize - 1) >= 30 &&
         ((targetSize - 1) - bit_loss) * nPrimes >= nBits)
    targetSize--;

  if (((targetSize - 1) - bit_loss) * nPrimes >= nBits)
    Warning(__func__ + std::string(": non-optimal targetSize"));

  PrimeGenerator gen(targetSize, m);

  while (nPrimes > 0) {
    long q = gen.next();

    if (context.inChain(q))
      continue;
    // nbits could equal NTL_SP_BITS or the size of one
    // of the small primes, so we have to check for duplicates here...
    // this is not the most efficient way to do this,
    // but it doesn't make sense to optimize this any further

    context.AddSpecialPrime(q);
    nPrimes--;
  }

  // std::cerr << "*** specialPrimes excess: " <<
  // (context.logOfProduct(context.specialPrimes)/std::log(2.0) - nBits) <<
  // "\n";
  HELIB_STATS_UPDATE(
      "excess-specialPrimes",
      context.logOfProduct(context.specialPrimes) / std::log(2.0) - nBits);
}

void endBuildModChain(Context& context)
{
  context.setModSizeTable();

  long m = context.zMStar.getM();
  std::vector<long> mvec;
  pp_factorize(mvec, m);
  NTL::Vec<long> mmvec;
  convert(mmvec, mvec);
  context.pwfl_converter = std::make_shared<PowerfulDCRT>(context, mmvec);
}

static void CheckPrimes(const Context& context,
                        const IndexSet& s,
                        const char* name)
{
  for (long i : s) {
    NTL::zz_pPush push;
    context.ithModulus(i).restoreModulus();
    if (!NTL::zz_p::IsFFTPrime()) {
      Warning(__func__ + std::string(": non-FFT prime in ") + name);
    }
  }
}

void buildModChain(Context& context,
                   long nBits,
                   long nDgts,
                   bool willBeBootstrappable,
                   long skHwt,
                   long resolution,
                   long bitsInSpecialPrimes)
{
  // Cannot build modulus chain with nBits < 0
  assertTrue<InvalidArgument>(nBits > 0,
                              "Cannot initialise modulus chain with nBits < 1");

  assertTrue(skHwt >= 0, "invalid skHwt parameter");

  if (context.isCKKS())
    willBeBootstrappable = false;
  // ignore for CKKS

  if (skHwt == 0) {
    // default skHwt: if bootstrapping, set to BOOT_DFLT_SK_HWT
    if (willBeBootstrappable)
      skHwt = BOOT_DFLT_SK_HWT;
  }

  // initialize hwt param in context
  context.hwt_param = skHwt;

  long pSize = ctxtPrimeSize(nBits);
  addSmallPrimes(context, resolution, pSize);
  addCtxtPrimes(context, nBits, pSize);
  addSpecialPrimes(context,
                   nDgts,
                   willBeBootstrappable,
                   skHwt,
                   bitsInSpecialPrimes);

  CheckPrimes(context, context.smallPrimes, "smallPrimes");
  CheckPrimes(context, context.ctxtPrimes, "ctxtPrimes");
  CheckPrimes(context, context.specialPrimes, "specialPrimes");

  endBuildModChain(context);
}

} // namespace helib
