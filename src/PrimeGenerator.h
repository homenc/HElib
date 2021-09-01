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

/* Intel HEXL integration and code reorg.
 * Copyright (C) 2021 Intel Corporation
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

// You initialize a PrimeGenerator as follows:
//    PrimeGenerator gen(len, m);
// Each call to gen.next() generates a prime p with
// (1-1/2^B)*2^len <= p < 2^len and p = 2^k*t*m + 1,
// where t is odd and k is as large as possible
// and B is a small constant (typically, B in {2,3,4}).
// If no such prime is found, then an error is raised.

#ifndef HELIB_PRIME_GENERATOR_H
#define HELIB_PRIME_GENERATOR_H

#include "macro.h" // Private Header

namespace helib {

class PrimeGenerator
{
private:
  long len;
  long m;
  long k;
  long t;

public:
  const static long B = 3;

  PrimeGenerator(long _len, long _m) : len(_len), m(_m)
  {
    assertInRange<InvalidArgument>(len,
                                   long(B),
                                   static_cast<long>(HELIB_SP_NBITS),
                                   "PrimeGenerator: len is not "
                                   "in [B, HELIB_SP_NBITS]",
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
    // It is equivalent to consider the interval of integers
    // [t_lower_bound, t_upper_bound),
    // where t_lower_bound = ceil(((1-1/2^B)*2^len-1)/(2^k*m))
    // and t_upper_bound = ceil((2^len-1)/(2^k*m)).

    long t_upper_bound = divc((1L << len) - 1, m << k);

    for (;;) {

      t++;

      if (t >= t_upper_bound) {
        // move to smaller value of k, reset t and t_upper_bound

        k--;

        long k_lower_bound = (m % 2 == 0) ? 0 : 1;

        // we run k down to 0  if m is even, and down to 1
        // if m is odd.

        if (k < k_lower_bound)
          throw RuntimeError("Prime generator ran out of primes");

        t = divc((1L << len) - (1L << (len - B)) - 1, m << k);
        t_upper_bound = divc((1L << len) - 1, m << k);
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

} // namespace helib
#endif // HELIB_PRIME_GENERATOR_H
