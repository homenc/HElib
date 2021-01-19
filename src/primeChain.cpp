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
#include "binio.h"
#include <helib/fhe_stats.h>
#include <helib/log.h>

#include "io.h"

namespace helib {

inline bool operator>(const ModuliSizes::Entry& a, const ModuliSizes::Entry& b)
{
  return a.first > b.first;
}

std::ostream& operator<<(std::ostream& s, const ModuliSizes::Entry& e)
{
  executeRedirectJsonError<void>([&]() {
    json j{{"first", e.first}, {"second", unwrap(e.second.writeToJSON())}};
    s << j.dump() << std::endl;
  });
  return s;
}

std::istream& operator>>(std::istream& s, ModuliSizes::Entry& e)
{
  executeRedirectJsonError<void>([&]() {
    json j;
    s >> j;
    e.first = j.at("first").get<double>();
    e.second = IndexSet::readFromJSON(wrap(j.at("second")));
  });
  return s;
}

void write(std::ostream& s, const ModuliSizes::Entry& e)
{
  write_raw_double(s, e.first);
  e.second.writeTo(s);
}

void read(std::istream& s, ModuliSizes::Entry& e)
{
  e.first = read_raw_double(s);
  e.second = IndexSet::readFrom(s);
}

// initialize helper table for a given chain
void ModuliSizes::init(const Context& context)
{
  if (context.getZMStar().getPow2())
    iFFT_cost = 0; // iFFT cost is same as FFT cost
  else
    iFFT_cost = 20; // iFFT cost is 1.20 times FFT cost.
                    // This is heuristic and derived from
                    // experimantal data.
                    // FIXME: should allow use to override
                    // this default value.

  long n = (1L << context.getSmallPrimes().card()) *
           (context.getCtxtPrimes().card() + 1);
  sizes.reserve(n); // allocate space
  // each entry of sizes is a pair<double,IndexSet>=(size, set-of-primes)

  // Get all subsets of smallPrimes

  sizes.push_back(std::make_pair(0.0, IndexSet::emptySet())); // the empty set
  long idx = 1; // first index that's still not set

  for (long i : context.getSmallPrimes()) { // add i to all sets upto idx-1
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
  for (long i : context.getCtxtPrimes()) { // add i to all sets upto idx-1
    s.insert(i);                           // add prime to the interval
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

} // namespace helib
