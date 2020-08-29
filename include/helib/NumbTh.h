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
#ifndef HELIB_NUMBTH_H
#define HELIB_NUMBTH_H
/**
 * @file NumbTh.h
 * @brief Miscellaneous utility functions.
 **/
#include <vector>
#include <set>
#include <cmath>
#include <complex>
#include <string>
#include <climits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ctime>
#include <memory>

#include <NTL/version.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/xdouble.h>

#include <NTL/mat_GF2.h>
#include <NTL/mat_GF2E.h>
#include <NTL/GF2XFactoring.h>

#include <NTL/mat_lzz_p.h>
#include <NTL/mat_lzz_pE.h>
#include <NTL/lzz_pXFactoring.h>

#include <NTL/GF2EX.h>
#include <NTL/lzz_pEX.h>

#include <NTL/FFT.h>

// Test for the "right version" of NTL (currently 11.0.0)
#if (NTL_MAJOR_VERSION < 11)
#error "This version of HElib requires NTL version 11.0.0 or above"
#endif

#include <helib/assertions.h>
#include <helib/apiAttributes.h>

namespace helib {

extern const long double PI;

namespace FHEglobals {
//! @brief A dry-run flag
//! The dry-run option disables most operations, to save time. This lets
//! us quickly go over the evaluation of a circuit and estimate the
//! resulting noise magnitude, without having to actually compute anything.
extern bool dryRun;

//! @brief A list of required automorphisms
//! When non-nullptr, causes Ctxt::smartAutomorphism to just record the
//! requested automorphism rather than actually performing it. This can
//! be used to get a list of needed automorphisms for certain operations
//! and then generate all these key-switching matrices. Should only be
//! used in conjunction with dryRun=true
extern std::set<long>* automorphVals;
extern std::set<long>* automorphVals2;

} // namespace FHEglobals

inline bool setDryRun(bool toWhat = true)
{
  return (FHEglobals::dryRun = toWhat);
}

inline bool isDryRun() { return FHEglobals::dryRun; }

inline void setAutomorphVals(std::set<long>* aVals)
{
  FHEglobals::automorphVals = aVals;
}

inline bool isSetAutomorphVals()
{
  return FHEglobals::automorphVals != nullptr;
}

inline void recordAutomorphVal(long k) { FHEglobals::automorphVals->insert(k); }

inline void setAutomorphVals2(std::set<long>* aVals)
{
  FHEglobals::automorphVals2 = aVals;
}

inline bool isSetAutomorphVals2()
{
  return FHEglobals::automorphVals2 != nullptr;
}

inline void recordAutomorphVal2(long k)
{
  FHEglobals::automorphVals2->insert(k);
}

typedef long LONG; // using this to identify casts that we should
                   // really get rid of at some point in the future

/**
 * @brief Considers `bits` as a vector of bits and returns the value it
 * represents when interpreted as a n-bit 2's complement number, where n is
 * given by `bitSize`.
 * @param bits The value containing the bits to be reinterpreted.
 * @param bitSize The number of bits to use, taken from the least significant
 * end of `bits`.
 * @return The value of the reinterpreted number as a long.
 **/
long bitSetToLong(long bits, long bitSize);

//! @brief Routines for computing mathematically correct mod and div.
//!
//! mcDiv(a, b) = floor(a / b), mcMod(a, b) = a - b*mcDiv(a, b);
//! in particular, mcMod(a, b) is 0 or has the same sign as b

long mcMod(long a, long b);
long mcDiv(long a, long b);

//! Return balanced remainder. Assumes a in [0, q) and returns
//! balanced remainder in (-q/2, q/2]
inline long balRem(long a, long q)
{
  if (a > q / 2)
    return a - q;
  else
    return a;
}

//! Return the square of a number as a double
inline double fsquare(double x) { return x * x; }

//! Return multiplicative order of p modulo m, or 0 if GCD(p, m) != 1
long multOrd(long p, long m);

//! @brief Prime power solver.
//!
//! A is an n x n matrix, b is a length n (row) vector, this function finds a
//! solution for the matrix-vector equation x A = b. An error is raised if A
//! is not invertible mod p.
//!
//! NTL's current smallint modulus, zz_p::modulus(), is assumed to be p^r,
//! for p prime, r >= 1 integer.
void ppsolve(NTL::vec_zz_pE& x,
             const NTL::mat_zz_pE& A,
             const NTL::vec_zz_pE& b,
             long p,
             long r);

//! @brief A version for GF2: must have p == 2 and r == 1
void ppsolve(NTL::vec_GF2E& x,
             const NTL::mat_GF2E& A,
             const NTL::vec_GF2E& b,
             long p,
             long r);

//! @brief Compute the inverse mod p^r of an n x n matrix.
//!
//! NTL's current smallint modulus zz_p::modulus() is assumed to be p^r for
//! p prime, r >= 1 integer. For the zz_pE variant also zz_pE::modulus() must
//! be initialized. An error is raised if A is not invertible mod p.
void ppInvert(NTL::mat_zz_p& X, const NTL::mat_zz_p& A, long p, long r);
void ppInvert(NTL::mat_zz_pE& X, const NTL::mat_zz_pE& A, long p, long r);

// variants for GF2/GF2E to help with template code
inline void ppInvert(NTL::mat_GF2& X,
                     const NTL::mat_GF2& A,
                     UNUSED long p,
                     UNUSED long r)
{
  NTL::inv(X, A);
}

inline void ppInvert(NTL::mat_GF2E& X,
                     const NTL::mat_GF2E& A,
                     UNUSED long p,
                     UNUSED long r)
{
  NTL::inv(X, A);
}

void buildLinPolyMatrix(NTL::mat_zz_pE& M, long p);
void buildLinPolyMatrix(NTL::mat_GF2E& M, long p);

//! @brief Combination of buildLinPolyMatrix and ppsolve.
//!
//! Obtain the linearized polynomial coefficients from a vector L representing
//! the action of a linear map on the standard basis for zz_pE over zz_p.
//!
//! NTL's current smallint modulus, zz_p::modulus(), is assumed to be p^r,
//! for p prime, r >= 1 integer.
void buildLinPolyCoeffs(NTL::vec_zz_pE& C,
                        const NTL::vec_zz_pE& L,
                        long p,
                        long r);

//! @brief A version for GF2: must be called with p == 2 and r == 1
void buildLinPolyCoeffs(NTL::vec_GF2E& C,
                        const NTL::vec_GF2E& L,
                        long p,
                        long r);

//! @brief Apply a linearized polynomial with coefficient vector C.
//!
//! NTL's current smallint modulus, zz_p::modulus(), is assumed to be p^r,
//! for p prime, r >= 1 integer.
void applyLinPoly(NTL::zz_pE& beta,
                  const NTL::vec_zz_pE& C,
                  const NTL::zz_pE& alpha,
                  long p);

//! @brief A version for GF2: must be called with p == 2 and r == 1
void applyLinPoly(NTL::GF2E& beta,
                  const NTL::vec_GF2E& C,
                  const NTL::GF2E& alpha,
                  long p);

//! Base-2 logarithm
inline double log2(const NTL::xdouble& x) { return log(x) * 1.442695040889; }

//! @brief Factoring by trial division, only works for N<2^{60}, only the
//! primes are recorded, not their multiplicity.
void factorize(std::vector<long>& factors, long N);
void factorize(std::vector<NTL::ZZ>& factors, const NTL::ZZ& N);

//! @brief Factoring by trial division, only works for N<2^{60}
//! primes and multiplicities are recorded
void factorize(NTL::Vec<NTL::Pair<long, long>>& factors, long N);

//! @brief Prime-power factorization
void pp_factorize(std::vector<long>& factors, long N);

//! Compute Phi(N) and also factorize N.
void phiN(long& phiN, std::vector<long>& facts, long N);
void phiN(NTL::ZZ& phiN, std::vector<NTL::ZZ>& facts, const NTL::ZZ& N);

//! Compute Phi(N).
long phi_N(long N);

//! Returns in gens a generating set for Zm* /<p>, and in ords the
//! order of these generators. Return value is the order of p in Zm*.
long findGenerators(std::vector<long>& gens,
                    std::vector<long>& ords,
                    long m,
                    long p,
                    const std::vector<long>& candidates = std::vector<long>());

//! Find e-th root of unity modulo the current modulus.
void FindPrimitiveRoot(NTL::zz_p& r, unsigned long e);
void FindPrimitiveRoot(NTL::ZZ_p& r, unsigned long e);

//! Compute mobius function (naive method as n is small).
long mobius(long n);

//! Compute cyclotomic polynomial.
NTL::ZZX Cyclotomic(long N);

//! Return a degree-d irreducible polynomial mod p
NTL::ZZX makeIrredPoly(long p, long d);

//! Find a primitive root modulo N.
long primroot(long N, long phiN);

//! Compute the highest power of p that divides N.
long ord(long N, long p);

inline bool is2power(long m)
{
  long k = NTL::NextPowerOfTwo(m);
  return (((unsigned long)m) == (1UL << k));
}

// Returns a random mod p polynomial of degree < n
NTL::ZZX RandPoly(long n, const NTL::ZZ& p);

///@{
/**
 * @brief Reduce all the coefficients of a polynomial modulo q.
 *
 * When abs=false reduce to interval (-q/2,...,q/2), when abs=true reduce
 * to [0,q). When abs=false and q=2, maintains the same sign as the input.
 */
void PolyRed(NTL::ZZX& out, const NTL::ZZX& in, long q, bool abs = false);
void PolyRed(NTL::ZZX& out,
             const NTL::ZZX& in,
             const NTL::ZZ& q,
             bool abs = false);

inline void PolyRed(NTL::ZZX& F, long q, bool abs = false)
{
  PolyRed(F, F, q, abs);
}

inline void PolyRed(NTL::ZZX& F, const NTL::ZZ& q, bool abs = false)
{
  PolyRed(F, F, q, abs);
}

void vecRed(NTL::Vec<NTL::ZZ>& out,
            const NTL::Vec<NTL::ZZ>& in,
            long q,
            bool abs);

void vecRed(NTL::Vec<NTL::ZZ>& out,
            const NTL::Vec<NTL::ZZ>& in,
            const NTL::ZZ& q,
            bool abs);
///@}

// The interface has changed so that abs defaults to false,
// which is more consistent with the other interfaces.
// Calls without any explicit value for abs should generate a
// "deprecated" warning.

void MulMod(NTL::ZZX& out, const NTL::ZZX& f, long a, long q, bool abs);

[[deprecated("Please use MulMod with explicit abs argument.")]] inline void
MulMod(NTL::ZZX& out, const NTL::ZZX& f, long a, long q)
{
  MulMod(out, f, a, q, false);
}

inline NTL::ZZX MulMod(const NTL::ZZX& f, long a, long q, bool abs)
{
  NTL::ZZX res;
  MulMod(res, f, a, q, abs);
  return res;
}

[[deprecated("Please use MulMod with explicit abs argument.")]] inline NTL::ZZX
MulMod(const NTL::ZZX& f, long a, long q)
{
  NTL::ZZX res;
  MulMod(res, f, a, q, false);
  return res;
}

//! Multiply the polynomial f by the integer a modulo q
//! output coefficients are balanced (appropriately randomized for even q)
void balanced_MulMod(NTL::ZZX& out, const NTL::ZZX& f, long a, long q);

///@{
//! @name Some enhanced conversion routines
inline void convert(long& x1, const NTL::GF2X& x2) { x1 = rep(ConstTerm(x2)); }
inline void convert(long& x1, const NTL::zz_pX& x2) { x1 = rep(ConstTerm(x2)); }
void convert(NTL::vec_zz_pE& X, const std::vector<NTL::ZZX>& A);
void convert(NTL::mat_zz_pE& X, const std::vector<std::vector<NTL::ZZX>>& A);
void convert(std::vector<NTL::ZZX>& X, const NTL::vec_zz_pE& A);
void convert(std::vector<std::vector<NTL::ZZX>>& X, const NTL::mat_zz_pE& A);
void convert(NTL::Vec<long>& out, const NTL::ZZX& in);
void convert(NTL::Vec<long>& out, const NTL::zz_pX& in, bool symmetric = true);
void convert(NTL::Vec<long>& out, const NTL::GF2X& in);
void convert(NTL::ZZX& out, const NTL::Vec<long>& in);
void convert(NTL::GF2X& out, const NTL::Vec<long>& in);
// right now, this is just a place-holder...it may or may not
// eventually be further fleshed out

///@}

//! A generic template that resolves to NTL's conv routine
template <typename T1, typename T2>
void convert(T1& x1, const T2& x2)
{
  NTL::conv(x1, x2);
}

//! generic vector conversion routines
template <typename T1, typename T2>
void convert(std::vector<T1>& v1, const std::vector<T2>& v2)
{
  long n = v2.size();
  v1.resize(n);
  for (long i = 0; i < n; i++)
    convert(v1[i], v2[i]);
}

template <typename T1, typename T2>
void convert(std::vector<T1>& v1, const NTL::Vec<T2>& v2)
{
  long n = v2.length();
  v1.resize(n);
  for (long i = 0; i < n; i++)
    convert(v1[i], v2[i]);
}

template <typename T1, typename T2>
void convert(NTL::Vec<T1>& v1, const std::vector<T2>& v2)
{
  long n = v2.size();
  v1.SetLength(n);
  for (long i = 0; i < n; i++)
    convert(v1[i], v2[i]);
}

//! Trivial type conversion, useful for generic code
template <typename T>
void convert(std::vector<T>& v1, const std::vector<T>& v2)
{
  v1 = v2;
}

template <typename T1, typename T2>
T1 convert(const T2& v2)
{
  T1 v1;
  convert(v1, v2);
  return v1;
}

template <typename T>
std::vector<T> vector_replicate(const T& a, long n)
{
  std::vector<T> res;
  res.resize(n);
  for (long i = 0; i < n; i++)
    res[i] = a;
  return res;
}

template <typename T>
std::vector<T> Vec_replicate(const T& a, long n)
{
  NTL::Vec<T> res;
  res.SetLength(n);
  for (long i = 0; i < n; i++)
    res[i] = a;
  return res;
}

//! returns \prod_d vec[d]
long computeProd(const NTL::Vec<long>& vec);
long computeProd(const std::vector<long>& vec);

// some useful operations
void mul(std::vector<NTL::ZZX>& x, const std::vector<NTL::ZZX>& a, long b);
void div(std::vector<NTL::ZZX>& x, const std::vector<NTL::ZZX>& a, long b);
void add(std::vector<NTL::ZZX>& x,
         const std::vector<NTL::ZZX>& a,
         const std::vector<NTL::ZZX>& b);

//! @brief Finds whether x is an element of the set X of size sz,
//! Returns -1 it not and the location if true
long is_in(long x, int* X, long sz);

//! @brief Returns a CRT coefficient: x = (0 mod p, 1 mod q).
//! If symmetric is set then x \in [-pq/2, pq/2), else x \in [0,pq)
inline long CRTcoeff(long p, long q, bool symmetric = false)
{
  long pInv = NTL::InvMod(p, q); // p^-1 mod q \in [0,q)
  if (symmetric && 2 * pInv >= q)
    return p * (pInv - q);
  else
    return p * pInv;
}

/**
 * @brief Incremental integer CRT for vectors.
 *
 * Expects co-primes p,q with q odd, and such that all the entries in v1 are
 * in [-p/2,p/2). Returns in v1 the CRT of vp mod p and vq mod q, as integers
 * in [-pq/2, pq/2). Uses the formula:
 * \f[               CRT(vp,p,vq,q) = vp + [(vq-vp) * p^{-1}]_q * p, \f]
 * where [...]_q means reduction to the interval [-q/2,q/2). Notice that if
 * q is odd then this is the same as reducing to [-(q-1)/2,(q-1)/2], which
 * means that [...]_q * p is in [-p(q-1)/2, p(q-1)/2], and since vp is in
 * [-p/2,p/2) then the sum is indeed in [-pq/2,pq/2).
 *
 * Return true is both vectors are of the same length, false otherwise
 */
template <class zzvec> // zzvec can be vec_NTL::ZZ, vec_long, or Vec<zz_p>
bool intVecCRT(NTL::vec_ZZ& vp, const NTL::ZZ& p, const zzvec& vq, long q);

/**
 * @brief Find the index of the (first) largest/smallest element.
 *
 * These procedures are roughly just simpler variants of std::max_element and
 * std::min_element. argmin/argmax are implemented as a template, so the code
 * must be placed in the header file for the compiler to find it. The class T
 * must have an implementation of operator> and operator< for this template to
 * work.
 * @tparam maxFlag A boolean value: true - argmax, false - argmin
 **/
template <typename T, bool maxFlag>
long argminmax(std::vector<T>& v)
{
  if (v.size() < 1)
    return -1; // error: this is an empty array
  unsigned long idx = 0;
  T target = v[0];
  for (unsigned long i = 1; i < v.size(); i++)
    if (maxFlag) {
      if (v[i] > target) {
        target = v[i];
        idx = i;
      }
    } else {
      if (v[i] < target) {
        target = v[i];
        idx = i;
      }
    }
  return (long)idx;
}

template <typename T>
long argmax(std::vector<T>& v)
{
  return argminmax<T, true>(v);
}

template <typename T>
long argmin(std::vector<T>& v)
{
  return argminmax<T, false>(v);
}

//! @brief A variant with a specialized comparison function
//! (*moreThan)(a,b) returns the comparison a>b
inline long argmax(std::vector<long>& v, bool (*moreThan)(long, long))
{
  if (v.size() < 1)
    return -INT_MAX; // error: this is an empty array
  unsigned long idx = 0;
  long target = v[0];
  for (unsigned long i = 1; i < v.size(); i++)
    if ((*moreThan)(v[i], target)) {
      target = v[i];
      idx = i;
    }
  return (long)idx;
}

// Check that x is in 1 += epsilon
inline bool closeToOne(const NTL::xdouble& x, long p)
{
  double pinv = 1.0 / p;
  return (x < (1.0 + pinv) && x > (1 - pinv));
}

// Use continued fractions to approximate a float x as x ~ a/b
std::pair<long, long> rationalApprox(double x, long denomBound = 0);
std::pair<NTL::ZZ, NTL::ZZ> rationalApprox(
    NTL::xdouble x,
    NTL::xdouble denomBound = NTL::xdouble(0.0));
/**
 * @brief Facility for "restoring" the NTL PRG state.
 *
 * NTL's random number generation facility is pretty limited, and does not
 * provide a way to save/restore the state of a pseudo-random stream. This
 * class gives us that ability: Constructing a RandomState object uses the PRG
 * to generate 512 bits and stores them. Upon destruction (or an explicit call
 * to restore()), these bits are used to re-set the seed of the PRG. A typical
 * usage of the class is as follows:
 * \code
 *   {
 *     RandomState r;      // save the random state
 *
 *     SetSeed(something); // set the PRG seed to something
 *     ...                 // more code that uses the new PRG seed
 *
 *   } // The destructor is called implicitly, PRG state is restored
 * \endcode
 **/
class RandomState
{
private:
  NTL::ZZ state;
  bool restored;

public:
  RandomState()
  {
    RandomBits(state, 512);
    restored = false;
  }

  //! Restore the PRG state of NTL
  void restore()
  {
    if (!restored) {
      SetSeed(state);
      restored = true;
    }
  }

  ~RandomState() { restore(); }

private:
  RandomState(const RandomState&);            // disable copy constructor
  RandomState& operator=(const RandomState&); // disable assignment
};

//! @brief Advance the input stream beyond white spaces and a single instance of
//! the char cc
void seekPastChar(std::istream& str, int cc);

/**
 * @brief Advance the input stream `str` beyond white spaces and a single
 * `separator` in the region-of-interest delimited by `begin_char` and
 * `end_char`.
 * @param str The stream to be advanced.
 * @param begin_char The character determining the beginning of the
 * region-of-interest (to advance beyond of).
 * @param separator The separator character to advance beyond of.
 * @param end_char The character determining the end of the region-of-interest
 * (to advance beyond of).
 * @return `true` if the region-of-interest is not completed (i.e.: `end_char`
 * is not reached). `false` otherwise.
 * @note Throws `helib::RuntimeError` if after spaces there is a character
 * different from `begin_char`, `beyond`, or `end_char`.
 */
bool iterateInterestRegion(std::istream& str,
                           int begin_char,
                           int separator,
                           int end_char);

/**
 * @brief Advance the input stream `istr` beyond white spaces.  Then split the
 * region delimited by `begin_char` and `end_char` at each occurrence of
 * `separator` that is not contained in an inner `begin_char` - `end_char`
 * section.  The function returns a `std::vector<std::stringstream>` with the
 * stream of every section of the input region.
 * @param istr The stream to be advanced.
 * @param begin_char The character determining the beginning of the
 * region-of-interest.
 * @param end_char The character determining the end of the
 * region-of-interest
 * @param separator The separator character to split at.
 * @param skip_space Boolean value determining whether to skip spaces when
 * extracting the sub-streams (default = `true`).
 * @return A `std::vector<std::stringstream>` with the stream of every section
 * of the input region.
 * @throws IOError If the stream is badly formatted (i.e. it does not start with
 * `begin_char` or it does not end with `end_char`).
 * @note Requires `begin_char`, `end_char` and `separator` to be distinct and
 * different from ` ` (space).
 */
std::vector<std::stringstream> extractTokenizeRegion(std::istream& istr,
                                                     char begin_char,
                                                     char end_char,
                                                     char separator,
                                                     bool skip_space = true);

//! @brief Reverse a vector in place
template <typename T>
void reverse(NTL::Vec<T>& v, long lo, long hi)
{
  long n = v.length();
  assertInRange(lo, 0l, hi, "Invalid argument: Bad interval", true);
  assertTrue(hi < n, "Invalid argument: Interval exceeds vector size");

  if (lo >= hi)
    return;

  for (long i = lo, j = hi; i < j; i++, j--)
    swap(v[i], v[j]);
}

//! @brief Rotate a vector in place using swaps
// Example: rotate by 1 means [0 1 2 3] -> [3 0 1 2]
//          rotate by -1 means [0 1 2 3] -> [1 2 3 0]
template <typename T>
void rotate(NTL::Vec<T>& v, long k)
{
  long n = v.length();
  if (n <= 1)
    return;

  k %= n;
  if (k < 0)
    k += n;

  if (k == 0)
    return;

  reverse(v, 0, n - 1);
  reverse(v, 0, k - 1);
  reverse(v, k, n - 1);
}

// An experimental facility as it is annoying that vector::size() is an
// unsigned quantity. This leads to all kinds of annoying warning messages.
//! @brief Size of STL vector as a long (rather than unsigned long)
template <typename T>
inline long lsize(const std::vector<T>& v)
{
  return (long)v.size();
}

//! NTL/std compatibility

// Utility functions, release memory of std::vector and NTL::Vec
template <typename T>
void killVec(std::vector<T>& vec)
{
  std::vector<T>().swap(vec);
}

template <typename T>
void killVec(NTL::Vec<T>& vec)
{
  vec.kill();
}

// Set length to zero, but don't necessarily release memory
template <typename T>
void setLengthZero(std::vector<T>& vec)
{
  if (vec.size() > 0)
    vec.resize(0, vec[0]);
}

template <typename T>
void setLengthZero(NTL::Vec<T>& vec)
{
  if (vec.length() > 0)
    vec.SetLength(0, vec[0]);
}

template <typename T>
inline long lsize(const NTL::Vec<T>& v)
{
  return v.length();
}

template <typename T>
void resize(NTL::Vec<T>& v, long sz, const T& val)
{
  return v.SetLength(sz, val);
}

template <typename T>
void resize(std::vector<T>& v, long sz, const T& val)
{
  return v.resize(sz, val);
}

template <typename T>
void resize(NTL::Vec<T>& v, long sz)
{
  return v.SetLength(sz);
}

template <typename T>
void resize(std::vector<T>& v, long sz)
{
  return v.resize(sz);
}

//! @brief Testing if two vectors point to the same object
// Believe it or not, this is really the way to do it...
template <typename T1, typename T2>
bool sameObject(const T1* p1, const T2* p2)
{
  return dynamic_cast<const void*>(p1) == dynamic_cast<const void*>(p2);
}

//! @brief Modular composition of polynomials: res = g(h) mod f
void ModComp(NTL::ZZX& res,
             const NTL::ZZX& g,
             const NTL::ZZX& h,
             const NTL::ZZX& f);

//! @brief Evaluates a modular integer polynomial, returns poly(x) mod p
long polyEvalMod(const NTL::ZZX& poly, long x, long p);

//! @brief Interpolate polynomial such that poly(x[i] mod p)=y[i] (mod p^e)
//! It is assumed that the points x[i] are all distinct modulo p
void interpolateMod(NTL::ZZX& poly,
                    const NTL::vec_long& x,
                    const NTL::vec_long& y,
                    long p,
                    long e = 1);

//! @brief returns ceiling(a/b); assumes a >=0, b>0, a+b <= MAX_LONG
inline long divc(long a, long b) { return (a + b - 1) / b; }

//! @class zz_pXModulus1
//! @brief Auxiliary classes to facilitate faster reduction mod Phi_m(X)
//!        when the input has degree less than m
class zz_pXModulus1
{
public:
  long m;
  NTL::zz_pX f;
  long n;

  bool specialLogic;

  long k, k1;
  NTL::fftRep R0, R1;

  NTL::zz_pXModulus fm; // just in case...

  zz_pXModulus1(long _m, const NTL::zz_pX& _f);

  const NTL::zz_pXModulus& upcast() const { return fm; }
};

void rem(NTL::zz_pX& r, const NTL::zz_pX& a, const zz_pXModulus1& ff);

//! placeholder for pXModulus ...no optimizations
class ZZ_pXModulus1 : public NTL::ZZ_pXModulus
{
public:
  ZZ_pXModulus1(UNUSED long _m, const NTL::ZZ_pX& _f) : NTL::ZZ_pXModulus(_f) {}
  const NTL::ZZ_pXModulus& upcast() const { return *this; }
};

template <typename T>
std::ostream& operator<<(std::ostream& s, std::vector<T> v)
{
  if (v.size() == 0)
    return (s << "[]");

  s << '[';
  for (long i = 0; i < (long)v.size() - 1; i++)
    s << v[i] << ' ';
  return (s << v[v.size() - 1] << ']');
}

template <typename T>
std::istream& operator>>(std::istream& s, std::vector<T>& v)
{
  NTL::Vec<T> vv; // read into an NTL vector, then convert
  s >> vv;
  convert(v, vv);
  return s;
}

template <typename T>
std::string vecToStr(const std::vector<T>& v)
{
  std::stringstream ss;
  ss << v;
  return ss.str();
}

template <typename T>
NTL::Vec<T> atoVec(const char* a)
{
  NTL::Vec<T> v;
  std::string s(a);
  std::stringstream ss(s);
  ss >> v;
  return v;
}

template <typename T>
std::vector<T> atovector(const char* a)
{
  NTL::Vec<T> v1 = atoVec<T>(a);
  std::vector<T> v2;
  convert(v2, v1);
  return v2;
}

#ifndef NTL_PROVIDES_TRUNC_FFT
// Define truncated FFT routines if not provided by NTL
inline void TofftRep_trunc(NTL::fftRep& y,
                           const NTL::zz_pX& x,
                           long k,
                           UNUSED long len,
                           long lo,
                           long hi)
{
  TofftRep(y, x, k, lo, hi);
}

inline void TofftRep_trunc(NTL::fftRep& y,
                           const NTL::zz_pX& x,
                           long k,
                           long len)
{
  TofftRep_trunc(y, x, k, len, 0, deg(x));
}
#endif

#if 0
//! @brief stand-in for make_unique, which is C++14, not C++11
template<typename T, typename... Args>
std::unique_ptr<T> build_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
#endif

//! Simple routine for computing the max-abs of a vector
//! of complex numbers and real numbers

inline double max_abs(const std::vector<std::complex<double>>& vec)
{
  double res = 0;
  for (auto x : vec) {
    double t = std::abs(x);
    if (res < t)
      res = t;
  }
  return res;
}

inline double max_abs(const std::vector<double>& vec)
{
  double res = 0;
  for (auto x : vec) {
    double t = std::abs(x);
    if (res < t)
      res = t;
  }
  return res;
}

//! This should go in NTL some day...
//! Just call as make_lazy(obj, ...) to initialize a lazy object
//! via a call to a constructor T(...)
template <typename T, typename P, typename... Args>
void make_lazy(const NTL::Lazy<T, P>& obj, Args&&... args)
{
  typename NTL::Lazy<T, P>::Builder builder(obj);
  if (!builder())
    return;
  NTL::UniquePtr<T, P> ptr;
  ptr.make(std::forward<Args>(args)...);
  builder.move(ptr);
}

//! This should go in NTL some day...
//! Just call as make_lazy(obj, f, ....) to initialize a lazy object
//! via a call to f(*obj, ...)
template <typename T, typename P, typename F, typename... Args>
void make_lazy_with_fun(const NTL::Lazy<T, P>& obj, F f, Args&&... args)
{
  typename NTL::Lazy<T, P>::Builder builder(obj);
  if (!builder())
    return;
  NTL::UniquePtr<T, P> ptr;
  ptr.make();
  f(*ptr, std::forward<Args>(args)...);
  builder.move(ptr);
}

// An array of inverse erfc values.
// erfc_inverse[i] = x means 2^{-i} = erfc(x/sqrt(2))

const double erfc_inverse[] = {0,
                               0.6744897501960817432,
                               1.1503493803760081782,
                               1.5341205443525463117,
                               1.8627318674216514554,
                               2.1538746940614562129,
                               2.4175590162365050618,
                               2.6600674686174596585,
                               2.8856349124267571473,
                               3.0972690781987844623,
                               3.2971933456919633418,
                               3.4871041041144311068,
                               3.6683292851213230192,
                               3.8419306855019108708,
                               4.0087725941685849622,
                               4.1695693233491057549,
                               4.3249190408260462571,
                               4.4753284246542033544,
                               4.6212310014992471565,
                               4.7630010342678139569,
                               4.9009642079631930118};

#define ERFC_INVERSE_SIZE (long(sizeof(erfc_inverse) / sizeof(erfc_inverse[0])))

} // namespace helib

#endif // ifndef HELIB_NUMBTH_H
