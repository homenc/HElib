// TODO add copyright 

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
    // It is equivalent to consider the interval
    // of integers [tlb, tub), where tlb = ceil(((1-1/2^B)*2^len-1)/(2^k*m))
    // and tub = ceil((2^len-1)/(2^k*m)).

    long tub = divc((1L << len) - 1, m << k);

    for (;;) {

      t++;

      if (t >= tub) {
        // move to smaller value of k, reset t and tub

        k--;

        long klb = (m % 2 == 0) ? 0 : 1;

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

} // namespace helib
#endif // HELIB_PRIME_GENERATOR_H
