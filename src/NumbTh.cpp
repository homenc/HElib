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
#include <helib/NumbTh.h>
#include <helib/timing.h>
#include <helib/range.h>
#include <helib/log.h>

#include <fstream>
#include <cctype>
#include <algorithm> // defines count(...), min(...)

namespace helib {

// First 70 odd decimal places of Pi == 4.0 * atan(1).
const long double PI =
    3.1415926535897932384626433832795028841971693993751058209749445923078164L;

bool FHEglobals::dryRun = false;

// Considering bits as a vector of bits, return the value it represents when
// interpreted as a bitSize-bit 2's complement number.
// For example, bitSetToLong(0b10111, 5) = -9.
long bitSetToLong(long bits, long bitSize)
{
  assertTrue<InvalidArgument>(bitSize >= 0, "bitSize must be non-negative.");
  long result = 0;
  for (long multiplier = 1; bitSize > 0; bits >>= 1, multiplier <<= 1)
    if (--bitSize != 0) // NB: The decrement for the loop is done here.
      result += (bits & 1) * multiplier;
    else
      result -= (bits & 1) * multiplier;
  return result;
}

// Mathematically correct mod and div, avoids overflow
long mcMod(long a, long b)
{
  long r = a % b;

  if (r != 0 && (b < 0) != (r < 0))
    return r + b;
  else
    return r;
}

long mcDiv(long a, long b)
{

  long r = a % b;
  long q = a / b;

  if (r != 0 && (b < 0) != (r < 0))
    return q + 1;
  else
    return q;
}

// return multiplicative order of p modulo m, or 0 if GCD(p, m) != 1
long multOrd(long p, long m)
{
  if (NTL::GCD(p, m) != 1)
    return 0;

  p = p % m;
  long ord = 1;
  long val = p;
  while (val != 1) {
    ord++;
    val = NTL::MulMod(val, p, m);
  }
  return ord;
}

// returns \prod_d vec[d]
template <typename T>
static inline long computeProd(const T& vec, long k)
{
  long prod = 1;
  for (long d = 0; d < k; d++)
    prod = prod * vec[d];
  return prod;
}
long computeProd(const NTL::Vec<long>& vec)
{
  return computeProd(vec, vec.length());
}
long computeProd(const std::vector<long>& vec)
{
  return computeProd(vec, vec.size());
}

// return a degree-d irreducible polynomial mod p
NTL::ZZX makeIrredPoly(long p, long d)
{
  assertTrue<InvalidArgument>(d >= 1l, "polynomial degree is less than 1");
  assertTrue<InvalidArgument>((bool)NTL::ProbPrime(p),
                              "modulus p is not prime");

  if (d == 1)
    return NTL::ZZX(1, 1); // the monomial X

  NTL::zz_pBak bak;
  bak.save();
  NTL::zz_p::init(p);
  return NTL::to_ZZX(NTL::BuildIrred_zz_pX(d));
}

// Factoring by trial division, only works for N<2^{60}.
// Only the primes are recorded, not their multiplicity
template <typename zz>
static void factorT(std::vector<zz>& factors, const zz& N)
{
  HELIB_TIMER_START;

  factors.resize(0); // reset the factors

  if (N < 2)
    return; // sanity check

  NTL::PrimeSeq s;
  zz n = N;
  while (n > 1) {

    if (NTL::ProbPrime(n)) { // we are left with just a single prime
      factors.push_back(n);
      return;
    }

    // n is a composite, so find the next prime that divides it
    long p = s.next();
    while (p && n % p != 0)
      p = s.next();

    if (!p)
      throw RuntimeError("ran out out small primes");

    zz pp;
    NTL::conv(pp, p);
    factors.push_back(pp);
    do {
      n /= p;
    } while (n % p == 0);
  }
}

void factorize(std::vector<long>& factors, long N)
{
  factorT<long>(factors, N);
}

void factorize(std::vector<NTL::ZZ>& factors, const NTL::ZZ& N)
{
  factorT<NTL::ZZ>(factors, N);
}

// Returns a list of prime factors and their multiplicity,
// N = \prod_i factors[i].first^{factors[i].second}
void factorize(NTL::Vec<NTL::Pair<long, long>>& factors, long N)
{
  factors.SetLength(0);

  if (N < 2)
    return;

  NTL::PrimeSeq s;
  long n = N;
  while (n > 1) {

    if (NTL::ProbPrime(n)) { // n itself is a prime, add (n,1) to the list
      append(factors, NTL::cons(n, 1L));
      return;
    }

    // n is composite, so find next prime that divides it

    long p = s.next();
    while (p && n % p != 0)
      p = s.next();

    if (!p)
      throw RuntimeError("ran out out small primes");

    long e = 1;
    n = n / p;
    while ((n % p) == 0) {
      n = n / p;
      e++;
    }
    append(factors, NTL::cons(p, e)); // add (p,e) to the list
  }
}

// Prime-power factorization
void pp_factorize(std::vector<long>& factors, long N)
{
  NTL::Vec<NTL::Pair<long, long>> pf;
  factorize(pf, N); // prime factors, N = \prod_i pf[i].first^{pf[i].second}
  factors.resize(pf.length());
  for (long i = 0; i < pf.length(); i++)
    factors[i] = NTL::power_long(pf[i].a, pf[i].b); // p_i^e_i
}

template <typename zz>
static void phiNT(zz& phin, std::vector<zz>& facts, const zz& N)
{
  if (facts.size() == 0)
    factorize(facts, N);

  zz n = N;
  NTL::conv(phin, 1); // initialize phiN=1
  for (unsigned long i = 0; i < facts.size(); i++) {
    zz p = facts[i];
    phin *= (p - 1); // first factor of p
    for (n /= p; (n % p) == 0; n /= p)
      phin *= p; // multiple factors of p
  }
}
// Specific template instantiations for long and ZZ
void phiN(long& pN, std::vector<long>& fs, long N) { phiNT<long>(pN, fs, N); }
void phiN(NTL::ZZ& pN, std::vector<NTL::ZZ>& fs, const NTL::ZZ& N)
{
  phiNT<NTL::ZZ>(pN, fs, N);
}

/* Compute Phi(N) */
long phi_N(long N)
{
  long phiN = 1, p, e;
  NTL::PrimeSeq s;
  while (N != 1) {
    p = s.next();
    e = 0;
    while ((N % p) == 0) {
      N = N / p;
      e++;
    }
    if (e != 0) {
      phiN = phiN * (p - 1) * NTL::power_long(p, e - 1);
    }
  }
  return phiN;
}

/* While generating the representation of (Z/mZ)^*, we keep the elements in
 * equivalence classes, and each class has a representative element (called
 * a pivot), which is the smallest element in the class. Initially each element
 * is in its own class. When we add a new generator g we unify classes if
 * their members are a factor of g from each other, repeating this process
 * until no further unification is possible.
 *
 * We begin by adding p as a generator, thus computing the equivalence
 * classes of (Z/mZ)^* /<p>. Then we repeatedly compute the orders of
 * all elements in the current quotient group, choose the highest-order
 * element and add it as a generator, then recompute the new quotient
 * group and so on, until the remaining quotient group is the trivial
 * one, containing just a a single element. A twist is that when choosing
 * the highest-order generator, we try to find one whose order in the
 * current quotient group is the same as in the original group (Z/mZ)^*,
 * and failing that, we find one whose order in the current quotient group
 * is the same as in (Z/mZ)^* /<p>.  This last type of generator is
 * always guaranteed.
 **/

// The function conjClasses(classes,g,m) unifies equivalence classes that
// have elements which are a factor of g apart, the pivot of the unified
// class is the smallest element in that class.
static void conjClasses(std::vector<long>& classes, long g, long m)
{
  NTL::mulmod_t minv = NTL::PrepMulMod(m);

  for (long i = 0; i < m; i++) {
    if (classes[i] == 0)
      continue; // i \notin (Z/mZ)^*

    if (classes[i] < i) { // i is not a pivot, updated its pivot
      classes[i] = classes[classes[i]];
      continue;
    }

    // If i is a pivot, update other pivots to point to it
    NTL::mulmod_precon_t gminv = NTL::PrepMulModPrecon(g, m, minv);
    long j = NTL::MulModPrecon(i, g, m, gminv);
    while (classes[j] != i) {
      classes[classes[j]] = i; // Merge the equivalence classes of j and i

      // Note: if classes[j]!=j then classes[j] will be updated later,
      //       when we get to i=j and use the code for "i not pivot".

      j = NTL::MulModPrecon(j, g, m, gminv);
    }
  }
}

// The function compOrder(orders, classes,flag,m) computes the order of elements
// of the quotient group, relative to current equivalent classes. If flag==1
// then also check if the order is the same as in (Z/mZ)^* and store the order
// with negative sign if not.
static void compOrder(std::vector<long>& orders,
                      std::vector<long>& classes,
                      long m)
{
  NTL::mulmod_t minv = NTL::PrepMulMod(m);

  orders[0] = 0;
  orders[1] = 1;
  for (long i = 2; i < m; i++) {
    if (classes[i] <= 1) { // ignore i not in Z_m^* and order-0 elements
      orders[i] = (classes[i] == 1) ? 1 : 0;
      continue;
    }

    if (classes[i] < i) { // not a pivot
      orders[i] = orders[classes[i]];
      continue;
    }

    NTL::mulmod_precon_t iminv = NTL::PrepMulModPrecon(i, m, minv);

    // For an element i>1, the order is at least 2
    long j = NTL::MulModPrecon(i, i, m, iminv);
    long ord = 2;
    while (classes[j] != 1) {
      j = NTL::MulModPrecon(j, i, m, iminv); // next element in <i>
      ord++; // count how many steps until we reach 1
    }

    orders[i] = ord;
  }
}

// Returns in gens a generating set for Zm* /<p>, and in ords the
// order of these generators. Return value is the order of p in Zm*.
long findGenerators(std::vector<long>& gens,
                    std::vector<long>& ords,
                    long m,
                    long p,
                    const std::vector<long>& candidates)
{
  gens.clear();
  ords.clear();
  // Compute the generators for (Z/mZ)^*
  std::vector<long> classes(m);
  std::vector<long> orders(m);

  for (long i = 0; i < m; i++) { // initially each element in its own class
    if (NTL::GCD(i, m) != 1)
      classes[i] = 0; // i is not in (Z/mZ)^*
    else
      classes[i] = i;
  }

  // Start building a representation of (Z/mZ)^*, first use the generator p
  conjClasses(classes, p % m, m); // merge classes that have a factor of p

  std::vector<int> p_subgp(m);
  for (long i : range(m))
    p_subgp[i] = 0;

  // The order of p is the size of the equivalence class of 1
  long ordP = 0;
  for (long i : range(m)) {
    if (classes[i] == 1) {
      ordP++;
      p_subgp[i] = 1;
    }
  }

  long candIdx = 0;
  while (true) {
    compOrder(orders, classes, m);

    long idx = 0;
    if (candIdx < lsize(candidates)) { // try next candidate
      idx = candidates[candIdx++];
      if (orders[idx] <= 1)
        idx = 0;
    }
    if (idx == 0) { // no viable candidates supplied externally
      long largest_ord = 1;

      for (long i : range(m)) {
        if (largest_ord < orders[i])
          largest_ord = orders[i];
      }

      if (largest_ord > 1) {
        long best_quality = 0;
        long best_idx = -1;

        for (long i = 0; i < m && best_quality < 2; i++) {
          if (orders[i] == largest_ord) {
            long j = NTL::PowerMod(i, largest_ord, m);
            if (j == 1) {
              best_idx = i;
              best_quality = 2;
            } else if (best_quality < 1 && p_subgp[j]) {
              best_idx = i;
              best_quality = 1;
            }
          }
        }

        if (!best_quality)
          Warning("low quality generator");
        idx = best_idx;
      }
    }

    if (!idx)
      break; // we are done

    // store generator with same order as in (Z/mZ)^*
    gens.push_back(idx);
    ords.push_back(orders[idx]);
    conjClasses(classes, idx, m); // merge classes that have a factor of idx
  }
  return ordP;
}

// finding e-th root of unity modulo the current modulus
// VJS: rewritten to be both faster and deterministic,
//  and assumes that current modulus is prime

template <typename zp, typename zz>
void FindPrimRootT(zp& root, unsigned long e)
{
  zz qm1 = zp::modulus() - 1;

  assertEq(static_cast<long>(qm1 % e), 0l, "e does not divide zp::modulus()-1");

  std::vector<long> facts;
  factorize(facts, e); // factorization of e

  root = 1;

  for (unsigned long i = 0; i < facts.size(); i++) {
    long p = facts[i];
    long pp = p;
    long ee = e / p;
    while (ee % p == 0) {
      ee = ee / p;
      pp = pp * p;
    }
    // so now we have e = pp * ee, where pp is
    // the power of p that divides e.
    // Our goal is to find an element of order pp

    NTL::PrimeSeq s;
    long q;
    zp qq, qq1;
    long iter = 0;
    do {
      iter++;
      if (iter > 1000000)
        throw RuntimeError("FindPrimitiveRoot: possible infinite loop?");
      q = s.next();
      NTL::conv(qq, q);
      power(qq1, qq, qm1 / p);
    } while (qq1 == 1);
    power(qq1, qq, qm1 / pp); // qq1 has order pp

    mul(root, root, qq1);
  }

  // independent check that we have an e-th root of unity
  {
    zp s;

    power(s, root, e);
    if (s != 1)
      throw RuntimeError("FindPrimitiveRoot: internal error (1)");

    // check that s^{e/p} != 1 for any prime divisor p of e
    for (unsigned long i = 0; i < facts.size(); i++) {
      long e2 = e / facts[i];
      power(s, root, e2); // s = root^{e/p}
      if (s == 1)
        throw RuntimeError("FindPrimitiveRoot: internal error (2)");
    }
  }
}
// instantiations of the template
void FindPrimitiveRoot(NTL::zz_p& r, unsigned long e)
{
  FindPrimRootT<NTL::zz_p, long>(r, e);
}
void FindPrimitiveRoot(NTL::ZZ_p& r, unsigned long e)
{
  FindPrimRootT<NTL::ZZ_p, NTL::ZZ>(r, e);
}

/* Compute mobius function (naive method as n is small) */
long mobius(long n)
{
  long p, e, arity = 0;
  NTL::PrimeSeq s;
  while (n != 1) {
    p = s.next();
    e = 0;
    while ((n % p == 0)) {
      n = n / p;
      e++;
    }
    if (e > 1) {
      return 0;
    }
    if (e != 0) {
      arity ^= 1;
    }
  }
  if (arity == 0) {
    return 1;
  }
  return -1;
}

// Based on Algorithm 4 (the Sparse Power Series Algorithm) from
// ANDREW ARNOLD AND MICHAEL MONAGAN, CALCULATING CYCLOTOMIC POLYNOMIALS,
// MATHEMATICS OF COMPUTATION, Volume 80, Number 276, October 2011,
// Pages 2359-2379

NTL::ZZX Cyclotomic(long n)
{
  assertEq(n >= 1, true, "n >= 1");

  // remove 2's

  long num_twos = 0;
  while (n % 2 == 0) {
    n /= 2;
    num_twos++;
  }

  if (n == 1) {
    if (num_twos == 0) {
      NTL::ZZX res; // X-1
      SetCoeff(res, 1);
      SetCoeff(res, 0, -1);
      return res;
    } else {
      NTL::ZZX res; // X^{2^{num_twos-1}}+1
      SetCoeff(res, 1L << (num_twos - 1));
      SetCoeff(res, 0);
      return res;
    }
  }

  std::vector<long> facs;
  factorize(facs, n);

  long k = facs.size();

  long radn = 1;
  long phi_radn = 1;
  for (long i : range(k)) {
    radn *= facs[i];
    phi_radn *= (facs[i] - 1);
  }

  long D = phi_radn / 2;

  if (radn <= 10000000L) {
    // for n <= 10^6, results in Arnold and Monogan
    // imply that all coefficients of Phi_n(X) are
    // less than 2^27 in absolute value.
    // So we compute the coefficients using 32-bit arithmetic.

    // NOTE: _ntl_uint32 is either unsigned int or unsigned long.

    NTL::Vec<_ntl_uint32> A;
    A.SetLength(D + 1);
    A[0] = 1;
    for (long i = 1; i <= D; i++)
      A[i] = 0;

    for (long bits : range(1L << k)) {
      long d = 1;
      long parity = k & 1;
      for (long pos : range(k)) {
        if ((1L << pos) & bits) {
          d *= facs[pos];
          parity = 1 - parity;
        }
      }

      if (parity == 0) {
        for (long i = D; i >= d; i--)
          A[i] -= A[i - d];
      } else {
        for (long i = d; i <= D; i++)
          A[i] += A[i - d];
      }
    }

    if (num_twos > 0) {
      for (long i = 1; i <= D; i += 2)
        A[i] = -A[i];
    }

    long q = n / radn;
    if (num_twos > 0)
      q = q << (num_twos - 1);
    long phi_n = phi_radn * q;

    NTL::ZZX res;
    res.rep.SetLength(phi_n + 1);
    for (long i = 0; i <= D; i++)
      conv(res.rep[i * q], NTL::cast_signed(A[i]));
    for (long i = D + 1; i <= phi_radn; i++)
      conv(res.rep[i * q], NTL::cast_signed(A[phi_radn - i]));

    return res;
  } else {
    // Exactly the same logic, but with bigint arithmetic.
    // This is pretty academic...

    NTL::Vec<NTL::ZZ> A;
    A.SetLength(D + 1);
    A[0] = 1;
    for (long i = 1; i <= D; i++)
      A[i] = 0;

    for (long bits : range(1L << k)) {
      long d = 1;
      long parity = k & 1;
      for (long pos : range(k)) {
        if ((1L << pos) & bits) {
          d *= facs[pos];
          parity = 1 - parity;
        }
      }

      if (parity == 0) {
        for (long i = D; i >= d; i--)
          A[i] -= A[i - d];
      } else {
        for (long i = d; i <= D; i++)
          A[i] += A[i - d];
      }
    }

    if (num_twos > 0) {
      for (long i = 1; i <= D; i += 2)
        A[i] = -A[i];
    }

    long q = n / radn;
    if (num_twos > 0)
      q = q << (num_twos - 1);
    long phi_n = phi_radn * q;

    NTL::ZZX res;
    res.rep.SetLength(phi_n + 1);
    for (long i = 0; i <= D; i++)
      res.rep[i * q] = A[i];
    for (long i = D + 1; i <= phi_radn; i++)
      res.rep[i * q] = A[phi_radn - i];

    return res;
  }
}

/* Find a primitive root modulo N */
long primroot(long N, long phiN)
{
  long g = 2, p;
  NTL::PrimeSeq s;
  bool flag = false;

  while (flag == false) {
    flag = true;
    s.reset(1);
    do {
      p = s.next();
      if ((phiN % p) == 0) {
        if (NTL::PowerMod(g, phiN / p, N) == 1) {
          flag = false;
        }
      }
    } while (p < phiN && flag);
    if (flag == false) {
      g++;
    }
  }
  return g;
}

long ord(long N, long p)
{
  long o = 0;
  while ((N % p) == 0) {
    o++;
    N /= p;
  }
  return o;
}

NTL::ZZX RandPoly(long n, const NTL::ZZ& p)
{
  NTL::ZZX F;
  F.SetMaxLength(n);
  NTL::ZZ p2;
  p2 = p >> 1;
  for (long i = 0; i < n; i++) {
    SetCoeff(F, i, RandomBnd(p) - p2);
  }
  return F;
}

/* When q=2 maintains the same sign as the input */
void PolyRed(NTL::ZZX& out, const NTL::ZZX& in, const NTL::ZZ& q, bool abs)
{
  // ensure that out has the same degree as in
  out.SetMaxLength(deg(in) + 1); // allocate space if needed
  if (deg(out) > deg(in))
    trunc(out, out, deg(in) + 1); // remove high degrees

  NTL::ZZ q2;
  q2 = q >> 1;
  for (long i = 0; i <= deg(in); i++) {
    NTL::ZZ c = coeff(in, i);
    c %= q;
    if (abs) {
      if (c < 0)
        c += q;
    } else if (q != 2) {
      if (c > q2) {
        c = c - q;
      } else if (c < -q2) {
        c = c + q;
      }
    } else // q=2
    {
      if (sign(coeff(in, i)) != sign(c)) {
        c = -c;
      }
    }
    SetCoeff(out, i, c);
  }
}

void PolyRed(NTL::ZZX& out, const NTL::ZZX& in, long q, bool abs)
{
  // ensure that out has the same degree as in
  out.SetMaxLength(deg(in) + 1); // allocate space if needed
  if (deg(out) > deg(in))
    trunc(out, out, deg(in) + 1); // remove high degrees

  long q2;
  q2 = q >> 1;
  for (long i = 0; i <= deg(in); i++) {
    long c = coeff(in, i) % q;
    if (abs) {
      if (c < 0) {
        c = c + q;
      }
    } else if (q == 2) {
      if (coeff(in, i) < 0) {
        c = -c;
      }
    } else {
      if (c >= q2) {
        c = c - q;
      } else if (c < -q2) {
        c = c + q;
      }
    }
    SetCoeff(out, i, c);
  }
}

void vecRed(NTL::Vec<NTL::ZZ>& out,
            const NTL::Vec<NTL::ZZ>& in,
            long q,
            bool abs)
{
  out.SetLength(in.length()); // allocate space if needed

  for (long i = 0; i < in.length(); i++) {
    long c = in[i] % q;
    if (abs) {
      if (c < 0)
        c += q;
    } else if (q == 2) {
      if (in[i] < 0)
        c = -c;
    } else {
      if (c >= q / 2)
        c -= q;
      else if (c < -(q / 2))
        c += q;
    }
    out[i] = c;
  }
}

void vecRed(NTL::Vec<NTL::ZZ>& out,
            const NTL::Vec<NTL::ZZ>& in,
            const NTL::ZZ& q,
            bool abs)
{
  out.SetLength(in.length()); // allocate space if needed

  for (long i = 0; i < in.length(); i++) {
    NTL::ZZ c = in[i] % q;
    if (abs) {
      if (c < 0)
        c += q;
    } else if (q == 2) {
      if (in[i] < 0)
        c = -c;
    } else {
      if (c >= q / 2)
        c -= q;
      else if (c < -(q / 2))
        c += q;
    }
    out[i] = c;
  }
}

void MulMod(NTL::ZZX& out,
            const NTL::ZZX& f,
            long a,
            long q,
            bool abs /*default=false*/)
{
  long n = f.rep.length();
  out.rep.SetLength(n);

  NTL::mulmod_precon_t aqinv = NTL::PrepMulModPrecon(a, q);
  for (long i : range(n)) {
    long c = rem(f.rep[i], q);
    c = NTL::MulModPrecon(c, a, q, aqinv); // returns c \in [0,q-1]
    if (!abs && c >= q / 2)
      c -= q;
    out.rep[i] = c;
  }

  out.normalize();
}

void balanced_MulMod(NTL::ZZX& out, const NTL::ZZX& f, long a, long q)
{
  long n = f.rep.length();
  out.rep.SetLength(n);

  NTL::mulmod_precon_t aqinv = NTL::PrepMulModPrecon(a, q);
  for (long i : range(n)) {
    long c = rem(f.rep[i], q);
    c = NTL::MulModPrecon(c, a, q, aqinv); // returns c \in [0,q-1]
    if (c > q / 2 || (q % 2 == 0 && c == q / 2 && NTL::RandomBnd(2)))
      c -= q;
    out.rep[i] = c;
  }

  out.normalize();
}

long is_in(long x, int* X, long sz)
{
  for (long i = 0; i < sz; i++) {
    if (x == X[i]) {
      return i;
    }
  }
  return -1;
}

/* Incremental integer CRT for vectors. Expects co-primes p>0,q>0 with q odd,
 * and such that all the entries in vp are in [-p/2,p/2) and all entries in
 * vq are in [0,q-1). Returns in vp the CRT of vp mod p and vq mod q, as
 * integers in [-pq/2, pq/2). Uses the formula:
 *
 *   CRT(vp,p,vq,q) = vp + p*[ (vq-vp)*p^{-1} ]_q
 *
 * where [...]_q means reduction to the interval [-q/2,q/2). As q is odd then
 * this is the same as reducing to [-(q-1)/2,(q-1)/2], hence [...]_q * p is
 * in [-p(q-1)/2, p(q-1)/2], and since vp is in [-p/2,p/2) then the sum is
 * indeed in [-pq/2,pq/2).
 *
 * Returns true if both vectors are of the same length, false otherwise
 */
template <typename zzvec>
bool intVecCRT(NTL::vec_ZZ& vp, const NTL::ZZ& p, const zzvec& vq, long q)
{
  long pInv = NTL::InvMod(rem(p, q), q); // p^{-1} mod q
  long n = std::min(vp.length(), vq.length());
  long q_over_2 = q / 2;
  NTL::ZZ tmp;
  long vqi;
  NTL::mulmod_precon_t pqInv = NTL::PrepMulModPrecon(pInv, q);
  for (long i = 0; i < n; i++) {
    NTL::conv(vqi, vq[i]); // convert to single precision
    long vq_minus_vp_mod_q = NTL::SubMod(vqi, rem(vp[i], q), q);

    long delta_times_pInv =
        NTL::MulModPrecon(vq_minus_vp_mod_q, pInv, q, pqInv);
    if (delta_times_pInv > q_over_2)
      delta_times_pInv -= q;

    mul(tmp, delta_times_pInv, p); // tmp = [(vq_i-vp_i)*p^{-1}]_q * p
    vp[i] += tmp;
  }
  // other entries (if any) are 0 mod q
  for (long i = vq.length(); i < vp.length(); i++) {
    long minus_vp_mod_q = NTL::NegateMod(rem(vp[i], q), q);

    long delta_times_pInv = NTL::MulModPrecon(minus_vp_mod_q, pInv, q, pqInv);
    if (delta_times_pInv > q_over_2)
      delta_times_pInv -= q;

    mul(tmp, delta_times_pInv, p); // tmp = [(vq_i-vp_i)*p^{-1}]_q * p
    vp[i] += tmp;
  }
  return (vp.length() == vq.length());
}
// specific instantiations: vq can be vec_long, vec_ZZ, or Vec<zz_p>
template bool intVecCRT(NTL::vec_ZZ&, const NTL::ZZ&, const NTL::vec_ZZ&, long);
template bool intVecCRT(NTL::vec_ZZ&,
                        const NTL::ZZ&,
                        const NTL::vec_long&,
                        long);
template bool intVecCRT(NTL::vec_ZZ&,
                        const NTL::ZZ&,
                        const NTL::Vec<NTL::zz_p>&,
                        long);

// ModComp: a pretty lame implementation

void ModComp(NTL::ZZX& res,
             const NTL::ZZX& g,
             const NTL::ZZX& h,
             const NTL::ZZX& f)
{
  assertEq<InvalidArgument>(LeadCoeff(f),
                            NTL::ZZ(1l),
                            "polynomial is not monic");

  NTL::ZZX hh = h % f;
  NTL::ZZX r = NTL::to_ZZX(0);

  for (long i = deg(g); i >= 0; i--)
    r = (r * hh + coeff(g, i)) % f;

  res = r;
}

long polyEvalMod(const NTL::ZZX& poly, long x, long p)
{
  long ret = 0;
  x %= p;
  if (x < 0)
    x += p;
  NTL::mulmod_precon_t xpinv = NTL::PrepMulModPrecon(x, p);
  for (long i = deg(poly); i >= 0; i--) {
    long coeff = rem(poly[i], p);
    ret = NTL::AddMod(ret, coeff, p); // Add the coefficient of x^i
    if (i > 0)
      ret = NTL::MulModPrecon(ret, x, p, xpinv); // then mult by x
  }
  return ret;
}

static void recursiveInterpolateMod(NTL::ZZX& poly,
                                    const NTL::vec_long& x,
                                    NTL::vec_long& y,
                                    const NTL::vec_zz_p& xmod,
                                    NTL::vec_zz_p& ymod,
                                    long p,
                                    long p2e)
{
  if (p2e <= 1) { // recursion edge condition, mod-1 poly = 0
    clear(poly);
    return;
  }

  // convert y input to zz_p
  for (long j = 0; j < y.length(); j++)
    ymod[j] = NTL::to_zz_p(y[j] % p);

  // a polynomial p_i s.t. p_i(x[j]) = i'th p-base digit of poly(x[j])
  NTL::zz_pX polyMod;
  interpolate(polyMod, xmod, ymod); // interpolation modulo p
  NTL::ZZX polyTmp;
  NTL::conv(polyTmp, polyMod); // convert to ZZX

  // update ytmp by subtracting the new digit, then dividing by p
  for (long j = 0; j < y.length(); j++) {
    y[j] -= polyEvalMod(polyTmp, x[j], p2e); // mod p^e
    if (y[j] < 0)
      y[j] += p2e;
    // if (y[j] % p != 0) {
    //   cerr << "@@error (p2^e="<<p2e<<"): y["<<j<<"] not divisible by "<<p<<
    //   endl; exit(0);
    // }
    y[j] /= p;
  } // maybe it's worth optimizing above by using multi-point evaluation

  // recursive call to get the solution of poly'(x)=y mod p^{e-1}
  recursiveInterpolateMod(poly, x, y, xmod, ymod, p, p2e / p);

  // return poly = p*poly' + polyTmp
  poly *= p;
  poly += polyTmp;
}

// Interpolate the integer polynomial such that poly(x[i] mod p)=y[i] (mod p^e)
// It is assumed that the points x[i] are all distinct modulo p
void interpolateMod(NTL::ZZX& poly,
                    const NTL::vec_long& x,
                    const NTL::vec_long& y,
                    long p,
                    long e)
{
  poly = NTL::ZZX::zero();          // initialize to zero
  long p2e = NTL::power_long(p, e); // p^e

  NTL::vec_long ytmp(NTL::INIT_SIZE, y.length()); // A temporary writable copy
  for (long j = 0; j < y.length(); j++) {
    ytmp[j] = y[j] % p2e;
    if (ytmp[j] < 0)
      ytmp[j] += p2e;
  }

  NTL::zz_pBak bak;
  bak.save(); // Set the current modulus to p
  NTL::zz_p::init(p);

  NTL::vec_zz_p xmod(NTL::INIT_SIZE, x.length()); // convert to zz_p
  for (long j = 0; j < x.length(); j++)
    xmod[j] = NTL::to_zz_p(x[j] % p);

  NTL::vec_zz_p ymod(NTL::INIT_SIZE, y.length()); // scratch space
  recursiveInterpolateMod(poly, x, ytmp, xmod, ymod, p, p2e);
}

// advance the input stream beyond white spaces and a single instance of cc
void seekPastChar(std::istream& str, int cc)
{
  int c = str.get();
  while (isspace(c))
    c = str.get();
  if (c != cc) {
    std::stringstream ss;
    ss << "Seeking past character='" << static_cast<char>(cc) << "' (ascii "
       << cc << ")"
       << ", found an unknown character='" << static_cast<char>(c)
       << "' (ascii " << c << ")";
    throw IOError(ss.str());
  }
}

// Advance the input stream `str` beyond white spaces and a single
// `separator` in the region-of-interest delimited by `begin_char` and
// `end_char`.
bool iterateInterestRegion(std::istream& str,
                           int begin_char,
                           int separator,
                           int end_char)
{
  int c = str.get();
  while (isspace(c)) {
    c = str.get();
  }
  if (c == begin_char || c == separator) {
    // Reached beginning of region or reached a separator. Return true
    return true;
  } else if (c == end_char) {
    // Reached end_char. Return false
    return false;
  } else {
    // Reached something different. Throw
    std::stringstream ss;
    ss << "Iterating on region found a non-delimiting "
       << "character='" << static_cast<char>(c) << "' (ascii " << c << "). "
       << "Delimiters: "
       << "begin_char='" << static_cast<char>(begin_char) << "' (ascii "
       << begin_char << "), "
       << "separator='" << static_cast<char>(separator) << "' (ascii "
       << separator << "), or "
       << "end_char='" << static_cast<char>(end_char) << "' (ascii " << end_char
       << ")";
    throw IOError(ss.str());
  }
}

// Advance the input stream `str` beyond white spaces and then split the
// region-of-interest delimited by `begin_char` and `end_char` at every
// top-level occurrence of `separator`.
std::vector<std::stringstream> extractTokenizeRegion(std::istream& istr,
                                                     char begin_char,
                                                     char end_char,
                                                     char separator,
                                                     bool skip_space)
{
  // Check if the arguments are valid.
  assertNeq<InvalidArgument>(begin_char,
                             ' ',
                             "Invalid begin_char. "
                             "Should be different from ' ' (space).");
  assertNeq<InvalidArgument>(begin_char,
                             end_char,
                             "Invalid begin_char. "
                             "Should be different from end_char.");
  assertNeq<InvalidArgument>(begin_char,
                             separator,
                             "Invalid begin_char. "
                             "Should be different from separator.");
  assertNeq<InvalidArgument>(end_char,
                             ' ',
                             "Invalid end_char. "
                             "Should be different from ' ' (space).");
  assertNeq<InvalidArgument>(end_char,
                             separator,
                             "Invalid end_char. "
                             "Should be different from separator.");
  assertNeq<InvalidArgument>(separator,
                             ' ',
                             "Invalid separator. "
                             "Should be different from ' ' (space).");
  int ch = istr.get();
  // Skip leading whitespaces
  while (isspace(ch)) {
    ch = istr.get();
  }
  // Fail if the input is not starting with a '['
  if (ch != begin_char) {
    std::stringstream ss;
    ss << "Extract and tokenization of stream failed with: Region beginning "
          "with character='"
       << static_cast<char>(ch) << "' (ascii " << ch << "). "
       << "Delimiters: "
       << "begin_char='" << static_cast<char>(begin_char) << "' (ascii "
       << begin_char << "), "
       << "separator='" << static_cast<char>(separator) << "' (ascii "
       << separator << "), or "
       << "end_char='" << static_cast<char>(end_char) << "' (ascii " << end_char
       << ")";
    throw IOError(ss.str());
  }
  std::vector<std::stringstream> res;
  std::stringstream current_stream;
  int depth = 0;

  while (istr.peek() != EOF) {
    ch = istr.get();
    if (ch == ' ' && skip_space) {
      // Skip whitespaces
      continue;
    } else if (ch == end_char && depth == 0) {
      // We found the section closure
      if (!current_stream.str().empty()) {
        // If something is in the stream add it to the returns
        res.emplace_back(std::move(current_stream));
      }
      return res;
    } else if (ch == end_char && depth > 0) {
      // We found the end of an inner section
      depth--;
      current_stream << static_cast<char>(ch);
    } else if (ch == begin_char) {
      // We found the begin of an inner section
      depth++;
      current_stream << static_cast<char>(ch);
    } else if (ch == separator && depth == 0) {
      res.emplace_back(std::move(current_stream));
      current_stream = std::stringstream();
    } else {
      current_stream << static_cast<char>(ch);
    }
  }
  // Throw as the section was not closed
  std::stringstream ss;
  ss << "Extract and tokenization of stream failed with: Region not closed. "
     << "Delimiters: "
     << "begin_char='" << static_cast<char>(begin_char) << "' (ascii "
     << begin_char << "), "
     << "separator='" << static_cast<char>(separator) << "' (ascii "
     << separator << "), or "
     << "end_char='" << static_cast<char>(end_char) << "' (ascii " << end_char
     << ")";
  throw IOError(ss.str());
}

// stuff added relating to linearized polynomials and support routines

// Builds the matrix defining the linearized polynomial transformation.
//
// NTL's current smallint modulus, zz_p::modulus(), is assumed to be p^r,
// for p prime, r >= 1 integer.
//
// After calling this function, one can call ppsolve(C, L, M, p, r) to get
// the coefficients C for the linearized polynomial represented the linear
// map defined by its action on the standard basis for zz_pE over zz_p:
// for i = 0..zz_pE::degree()-1: x^i -> L[i], where x = (X mod zz_pE::modulus())

void buildLinPolyMatrix(NTL::mat_zz_pE& M, long p)
{
  long d = NTL::zz_pE::degree();

  M.SetDims(d, d);

  for (long j = 0; j < d; j++)
    NTL::conv(M[0][j], NTL::zz_pX(j, 1));

  for (long i = 1; i < d; i++)
    for (long j = 0; j < d; j++)
      M[i][j] = power(M[i - 1][j], p);
}

void buildLinPolyMatrix(NTL::mat_GF2E& M, long p)
{
  assertEq<InvalidArgument>(p,
                            2l,
                            "p is not 2 when building "
                            "a mat_GF2E (Galois field 2)");

  long d = NTL::GF2E::degree();

  M.SetDims(d, d);

  for (long j = 0; j < d; j++)
    NTL::conv(M[0][j], NTL::GF2X(j, 1));

  for (long i = 1; i < d; i++)
    for (long j = 0; j < d; j++)
      M[i][j] = power(M[i - 1][j], p);
}

// some auxiliary conversion routines

void convert(NTL::vec_zz_pE& X, const std::vector<NTL::ZZX>& A)
{
  long n = A.size();
  NTL::zz_pX tmp;
  X.SetLength(n);
  for (long i = 0; i < n; i++) {
    NTL::conv(tmp, A[i]);
    NTL::conv(X[i], tmp);
  }
}

void convert(NTL::mat_zz_pE& X, const std::vector<std::vector<NTL::ZZX>>& A)
{
  long n = A.size();

  if (n == 0) {
    long m = X.NumCols();
    X.SetDims(0, m);
    return;
  }

  long m = A[0].size();
  X.SetDims(n, m);

  for (long i = 0; i < n; i++)
    convert(X[i], A[i]);
}

void convert(std::vector<NTL::ZZX>& X, const NTL::vec_zz_pE& A)
{
  long n = A.length();
  X.resize(n);
  for (long i = 0; i < n; i++)
    NTL::conv(X[i], rep(A[i]));
}

void convert(std::vector<std::vector<NTL::ZZX>>& X, const NTL::mat_zz_pE& A)
{
  long n = A.NumRows();
  X.resize(n);
  for (long i = 0; i < n; i++)
    convert(X[i], A[i]);
}

void convert(NTL::Vec<long>& out, const NTL::ZZX& in)
{
  out.SetLength(in.rep.length());
  for (long i = 0; i < out.length(); i++)
    out[i] = NTL::conv<long>(in[i]);
}

void convert(NTL::Vec<long>& out, const NTL::zz_pX& in, bool symmetric)
{
  out.SetLength(in.rep.length());
  for (long i = 0; i < out.length(); i++)
    out[i] = NTL::conv<long>(in[i]);

  if (symmetric) { // convert to representation symmetric around 0
    long p = NTL::zz_p::modulus();
    for (long i = 0; i < out.length(); i++)
      if (out[i] > p / 2)
        out[i] -= p;
  }
}

void convert(NTL::Vec<long>& out, const NTL::GF2X& in)
{
  out.SetLength(1 + deg(in));
  for (long i = 0; i < out.length(); i++)
    out[i] = NTL::conv<long>(in[i]);
}

void convert(NTL::ZZX& out, const NTL::Vec<long>& in)
{
  out.SetLength(in.length());
  for (long i = 0; i < in.length(); i++)
    out[i] = NTL::conv<NTL::ZZ>(in[i]);
  out.normalize();
}

void convert(NTL::GF2X& out, const NTL::Vec<long>& in)
{
  out.SetLength(in.length());
  for (long i = 0; i < in.length(); i++)
    out[i] = NTL::conv<NTL::GF2>(in[i]);
  out.normalize();
}

void mul(std::vector<NTL::ZZX>& x, const std::vector<NTL::ZZX>& a, long b)
{
  long n = a.size();
  x.resize(n);
  for (long i = 0; i < n; i++)
    mul(x[i], a[i], b);
}

void div(std::vector<NTL::ZZX>& x, const std::vector<NTL::ZZX>& a, long b)
{
  long n = a.size();
  x.resize(n);
  for (long i = 0; i < n; i++)
    div(x[i], a[i], b);
}

void add(std::vector<NTL::ZZX>& x,
         const std::vector<NTL::ZZX>& a,
         const std::vector<NTL::ZZX>& b)
{
  long n = a.size();
  if (n != (long)b.size())
    throw InvalidArgument("add: a and b dimension differ");
  for (long i = 0; i < n; i++)
    add(x[i], a[i], b[i]);
}

// prime power solver
// zz_p::modulus() is assumed to be p^r, for p prime, r >= 1
// A is an n x n matrix, b is a length n (row) vector,
// and a solution for the matrix-vector equation x A = b is found.
// If A is not invertible mod p, then error is raised.
void ppsolve(NTL::vec_zz_pE& x,
             const NTL::mat_zz_pE& A,
             const NTL::vec_zz_pE& b,
             long p,
             long r)
{

  if (r == 1) {
    NTL::zz_pE det;
    solve(det, x, A, b);
    if (det == 0)
      throw InvalidArgument("ppsolve: matrix not invertible");
    return;
  }

  long n = A.NumRows();
  if (n != A.NumCols())
    throw InvalidArgument("ppsolve: matrix not square");
  if (n == 0)
    throw InvalidArgument("ppsolve: matrix of dimension 0");

  NTL::zz_pContext pr_context;
  pr_context.save();

  NTL::zz_pEContext prE_context;
  prE_context.save();

  NTL::zz_pX G = NTL::zz_pE::modulus();

  NTL::ZZX GG = to_ZZX(G);

  std::vector<std::vector<NTL::ZZX>> AA;
  convert(AA, A);

  std::vector<NTL::ZZX> bb;
  convert(bb, b);

  NTL::zz_pContext p_context(p);
  p_context.restore();

  NTL::zz_pX G1 = NTL::to_zz_pX(GG);
  NTL::zz_pEContext pE_context(G1);
  pE_context.restore();

  // we are now working mod p...

  // invert A mod p

  NTL::mat_zz_pE A1;
  convert(A1, AA);

  NTL::mat_zz_pE I1;
  NTL::zz_pE det;

  inv(det, I1, A1);
  if (det == 0) {
    throw LogicError("ppsolve: matrix not invertible");
  }

  NTL::vec_zz_pE b1;
  convert(b1, bb);

  NTL::vec_zz_pE y1;
  y1 = b1 * I1;

  std::vector<NTL::ZZX> yy;
  convert(yy, y1);

  // yy is a solution mod p

  for (long k = 1; k < r; k++) {
    // lift solution yy mod p^k to a solution mod p^{k+1}

    pr_context.restore();
    prE_context.restore();
    // we are now working mod p^r

    NTL::vec_zz_pE d, y;
    convert(y, yy);

    d = b - y * A;

    std::vector<NTL::ZZX> dd;
    convert(dd, d);

    long pk = NTL::power_long(p, k);
    std::vector<NTL::ZZX> ee;
    div(ee, dd, pk);

    p_context.restore();
    pE_context.restore();

    // we are now working mod p

    NTL::vec_zz_pE e1;
    convert(e1, ee);
    NTL::vec_zz_pE z1;
    z1 = e1 * I1;

    std::vector<NTL::ZZX> zz, ww;
    convert(zz, z1);

    mul(ww, zz, pk);
    add(yy, yy, ww);
  }

  pr_context.restore();
  prE_context.restore();

  convert(x, yy);

  assertEq(x * A, b, "Failed to found solution x to matrix equation x*A == b");
}

void ppsolve(NTL::vec_GF2E& x,
             const NTL::mat_GF2E& A,
             const NTL::vec_GF2E& b,
             long p,
             long r)
{
  assertEq<InvalidArgument>(p,
                            2l,
                            "modulus p is not 2 with GF2E (Galois field 2)");
  assertEq<InvalidArgument>(r,
                            1l,
                            "Hensel lifting r is not 2 with"
                            " GF2E (Galois field 2)");

  NTL::GF2E det;
  solve(det, x, A, b);
  if (det == 0)
    throw InvalidArgument("ppsolve: matrix not invertible");
}

// prime power solver
// A is an n x n matrix, we compute its inverse mod p^r. An error is raised
// if A is not invertible mod p. zz_p::modulus() is assumed to be p^r, for
// p prime, r >= 1. Also zz_pE::modulus() is assumed to be initialized.
void ppInvert(NTL::mat_zz_pE& X, const NTL::mat_zz_pE& A, long p, long r)
{
  if (r == 1) { // use native inversion from NTL
    inv(X, A);  // X = A^{-1}
    return;
  }

  // begin by inverting A modulo p

  // convert to ZZX for a safe translation to mod-p objects
  std::vector<std::vector<NTL::ZZX>> tmp;
  convert(tmp, A);
  { // open a new block for mod-p computation
    NTL::ZZX G;
    convert(G, NTL::zz_pE::modulus());
    NTL::zz_pBak bak_pr;
    bak_pr.save(); // backup the mod-p^r moduli
    NTL::zz_pEBak bak_prE;
    bak_prE.save();
    NTL::zz_p::init(p); // Set the mod-p moduli
    NTL::zz_pE::init(NTL::conv<NTL::zz_pX>(G));

    NTL::mat_zz_pE A1, Inv1;
    convert(A1, tmp);   // Recover A as a mat_zz_pE object modulo p
    inv(Inv1, A1);      // Inv1 = A^{-1} (mod p)
    convert(tmp, Inv1); // convert to ZZX for translation to a mod-p^r object
  } // mod-p^r moduli restored on destruction of bak_pr and bak_prE
  NTL::mat_zz_pE XX;
  convert(XX, tmp); // XX = A^{-1} (mod p)

  // Now lift the solution modulo p^r

  // Compute the "correction factor" Z, s.t. XX*A = I - p*Z (mod p^r)
  long n = A.NumRows();
  const NTL::mat_zz_pE I = NTL::ident_mat_zz_pE(n); // identity matrix
  NTL::mat_zz_pE Z = I - XX * A;

  convert(tmp, Z); // Convert to ZZX to divide by p
  for (long i = 0; i < n; i++)
    for (long j = 0; j < n; j++)
      tmp[i][j] /= p;
  convert(Z, tmp); // convert back to a mod-p^r object

  // The inverse of A is ( I+(pZ)+(pZ)^2+...+(pZ)^{r-1} )*XX (mod p^r). We use
  // O(log r) products to compute it as (I+pZ)* (I+(pZ)^2)* (I+(pZ)^4)*...* XX

  long e = NTL::NextPowerOfTwo(r); // 2^e is smallest power of two >= r

  Z *= p;                      // = pZ
  NTL::mat_zz_pE prod = I + Z; // = I + pZ
  for (long i = 1; i < e; i++) {
    sqr(Z, Z);       // = (pZ)^{2^i}
    prod *= (I + Z); // = sum_{j=0}^{2^{i+1}-1} (pZ)^j
  }
  mul(X, prod, XX); // X = A^{-1} mod p^r
  assertEq(X * A,
           I,
           "Failed to found solution X to matrix equation X*A == I "
           "where I is the identity matrix");
}

// FIXME: at some point need to make a template for these two functions
// prime power solver
// A is an n x n matrix, we compute its inverse mod p^r. An error is raised
// if A is not invertible mod p. zz_p::modulus() is assumed to be p^r, for
// p prime, r >= 1.
void ppInvert(NTL::mat_zz_p& X, const NTL::mat_zz_p& A, long p, long r)
{
  if (r == 1) { // use native inversion from NTL
    inv(X, A);  // X = A^{-1}
    return;
  }

  // begin by inverting A modulo p

  // convert to long for a safe translation to mod-p objects
  NTL::Mat<long> tmp;
  conv(tmp, A);
  { // open a new block for mod-p computation
    NTL::zz_pBak bak_pr;
    bak_pr.save();      // backup the mod-p^r moduli
    NTL::zz_p::init(p); // Set the mod-p moduli

    NTL::mat_zz_p A1, Inv1;
    conv(A1, tmp);   // Recover A as a mat_zz_pE object modulo p
    inv(Inv1, A1);   // Inv1 = A^{-1} (mod p)
    conv(tmp, Inv1); // convert to long for translation to a mod-p^r object
  } // mod-p^r moduli restored on destruction of bak_pr and bak_prE
  NTL::mat_zz_p XX;
  conv(XX, tmp); // XX = A^{-1} (mod p)

  // Now lift the solution modulo p^r

  // Compute the "correction factor" Z, s.t. XX*A = I - p*Z (mod p^r)
  long n = A.NumRows();
  const NTL::mat_zz_p I = NTL::ident_mat_zz_p(n); // identity matrix
  NTL::mat_zz_p Z = I - XX * A;

  conv(tmp, Z); // Convert to long to divide by p
  for (long i = 0; i < n; i++)
    for (long j = 0; j < n; j++)
      tmp[i][j] /= p;
  conv(Z, tmp); // convert back to a mod-p^r object

  // The inverse of A is ( I+(pZ)+(pZ)^2+...+(pZ)^{r-1} )*XX (mod p^r). We use
  // O(log r) products to compute it as (I+pZ)* (I+(pZ)^2)* (I+(pZ)^4)*...* XX

  long e = NTL::NextPowerOfTwo(r); // 2^e is smallest power of two >= r

  Z *= p;                     // = pZ
  NTL::mat_zz_p prod = I + Z; // = I + pZ
  for (long i = 1; i < e; i++) {
    sqr(Z, Z);       // = (pZ)^{2^i}
    prod *= (I + Z); // = sum_{j=0}^{2^{i+1}-1} (pZ)^j
  }
  mul(X, prod, XX); // X = A^{-1} mod p^r
  assertEq(X * A,
           I,
           "Failed to found solution X to matrix equation X*A == I "
           "where I is the identity matrix");
}

void buildLinPolyCoeffs(NTL::vec_zz_pE& C_out,
                        const NTL::vec_zz_pE& L,
                        long p,
                        long r)
{
  HELIB_TIMER_START;
  NTL::mat_zz_pE M;
  buildLinPolyMatrix(M, p);

  NTL::vec_zz_pE C;
  ppsolve(C, M, L, p, r);

  C_out = C;
  HELIB_TIMER_STOP;
}

void buildLinPolyCoeffs(NTL::vec_GF2E& C_out,
                        const NTL::vec_GF2E& L,
                        long p,
                        long r)
{
  HELIB_TIMER_START;
  assertEq<InvalidArgument>(p,
                            2l,
                            "modulus p is not 2 with GF2E (Galois field 2)");
  assertEq<InvalidArgument>(r,
                            1l,
                            "Hensel lifting r is not 2 "
                            "with GF2E (Galois field 2)");

  NTL::mat_GF2E M;
  buildLinPolyMatrix(M, p);

  NTL::vec_GF2E C;
  ppsolve(C, M, L, p, r);

  C_out = C;
  HELIB_TIMER_STOP;
}

void applyLinPoly(NTL::zz_pE& beta,
                  const NTL::vec_zz_pE& C,
                  const NTL::zz_pE& alpha,
                  long p)
{
  long d = NTL::zz_pE::degree();
  assertEq<InvalidArgument>(d,
                            C.length(),
                            "C length is not equal to NTL::zz_pE::degree()");

  NTL::zz_pE gamma, res;

  gamma = NTL::to_zz_pE(NTL::zz_pX(1, 1));
  res = C[0] * alpha;
  for (long i = 1; i < d; i++) {
    gamma = power(gamma, p);
    res += C[i] * NTL::to_zz_pE(
                      CompMod(rep(alpha), rep(gamma), NTL::zz_pE::modulus()));
  }

  beta = res;
}

void applyLinPoly(NTL::GF2E& beta,
                  const NTL::vec_GF2E& C,
                  const NTL::GF2E& alpha,
                  long p)
{
  long d = NTL::GF2E::degree();
  assertEq<InvalidArgument>(d,
                            C.length(),
                            "C length is not equal to GF2E::degree()");

  NTL::GF2E gamma, res;

  gamma = NTL::to_GF2E(NTL::GF2X(1, 1));
  res = C[0] * alpha;
  for (long i = 1; i < d; i++) {
    gamma = power(gamma, p);
    res +=
        C[i] * to_GF2E(CompMod(rep(alpha), rep(gamma), NTL::GF2E::modulus()));
  }

  beta = res;
}

// use continued fractions to get "best" rational approximation
std::pair<long, long> rationalApprox(double x, long denomBound)
{
  int sign = 1;
  if (x < 0) {
    sign = -1;
    x = -x;
  }
  if (denomBound <= 0)
    denomBound = 1L << (NTL_SP_NBITS / 2);
  double epsilon = 1.0 / (denomBound * 8.0); // "smudge factor"
  double a = floor(x + epsilon);
  double xi = x - a;
  long prevDenom = 0;
  long denom = 1;

  // Continued fractions: a_{i+1}=floor(1/xi), x_{i+1} = 1/xi - a_{i+1}
  while (xi > 0) {
    xi = 1 / xi;
    double ai = floor(
        xi + epsilon); // NOTE: epsilon is meant to counter rounding errors
    xi = xi - ai;

    double tmpDenom = denom * ai + prevDenom;
    if (tmpDenom > denomBound) // bound exceeded: return previous denominator
      break;
    // update denominator
    prevDenom = denom;
    denom = tmpDenom;
    //    cout << "  ai="<<ai<<", xi="<<xi<<", denominator="<<denom<<endl;
  }
  assertTrue<RuntimeError>(denom * x < NTL_SP_BOUND,
                           "Single-precision bound exceeded");
  long numer = long(round(denom * x)) * sign;

  return std::make_pair(numer, denom);
}

// use continued fractions to get "best" rational approximation
std::pair<NTL::ZZ, NTL::ZZ> rationalApprox(NTL::xdouble x,
                                           NTL::xdouble denomBound)
{
  int sign = 1;
  if (x < 0) {
    sign = -1;
    x = -x;
  }
  if (denomBound <= 0)
    denomBound = NTL::conv<NTL::xdouble>(1L << (NTL_SP_NBITS / 2));

  NTL::xdouble epsilon = 0.125 / denomBound; // "smudge factor"
  NTL::xdouble a = floor(x + epsilon);

  NTL::xdouble xi = x - a;
  NTL::xdouble prevDenom(0.0);
  NTL::xdouble xdenom(1.0);

  // Continued fractions: a_{i+1}=floor(1/xi), x_{i+1} = 1/xi - a_{i+1}
  while (xi > 0) {
    xi = 1 / xi;
    NTL::xdouble ai = floor(
        xi + epsilon); // NOTE: epsilon is meant to counter rounding errors
    xi = xi - ai;

    NTL::xdouble tmpDenom = xdenom * ai + prevDenom;
    if (tmpDenom > denomBound) // bound exceeded: return previous denominator
      break;
    // update denominator
    prevDenom = xdenom;
    xdenom = tmpDenom;
  }

  NTL::ZZ numer = NTL::conv<NTL::ZZ>(xdenom * x + 0.5) * sign;
  NTL::ZZ denom = NTL::conv<NTL::ZZ>(xdenom);

  return std::make_pair(numer, denom);
}

// Auxilliary classes to facilitate faster reduction mod Phi_m(X)
// when the input has degree less than m

static void LocalCopyReverse(NTL::zz_pX& x,
                             const NTL::zz_pX& a,
                             long lo,
                             long hi)

// x[0..hi-lo] = reverse(a[lo..hi]), with zero fill
// input may not alias output

{
  long i, j, n, m;

  n = hi - lo + 1;
  m = a.rep.length();

  x.rep.SetLength(n);

  const NTL::zz_p* ap = a.rep.elts();
  NTL::zz_p* xp = x.rep.elts();

  for (i = 0; i < n; i++) {
    j = hi - i;
    if (j < 0 || j >= m)
      clear(xp[i]);
    else
      xp[i] = ap[j];
  }

  x.normalize();
}

static void LocalCyclicReduce(NTL::zz_pX& x, const NTL::zz_pX& a, long m)

// computes x = a mod X^m-1

{
  long n = deg(a);
  long i, j;
  NTL::zz_p accum;

  if (n < m) {
    x = a;
    return;
  }

  if (&x != &a)
    x.rep.SetLength(m);

  for (i = 0; i < m; i++) {
    accum = a.rep[i];
    for (j = i + m; j <= n; j += m)
      add(accum, accum, a.rep[j]);
    x.rep[i] = accum;
  }

  if (&x == &a)
    x.rep.SetLength(m);

  x.normalize();
}

zz_pXModulus1::zz_pXModulus1(long _m, const NTL::zz_pX& _f) :
    m(_m), f(_f), n(deg(f))
{
  assertTrue<InvalidArgument>(m > n, "_m is less or equal than _f's degree");

  specialLogic = (m - n > 10 && m < 2 * n);

  build(fm, f);

  if (specialLogic) {
    // std::cout << "*** special\n";
    NTL::zz_pX P1, P2, P3;

    LocalCopyReverse(P3, f, 0, n);
    InvTrunc(P2, P3, m - n);
    LocalCopyReverse(P1, P2, 0, m - n - 1);

    k = NTL::NextPowerOfTwo(2 * (m - 1 - n) + 1);
    k1 = NTL::NextPowerOfTwo(n);

    TofftRep(R0, P1, k);
    TofftRep(R1, f, k1);
  }
  // else { std::cout << "=== non-special\n"; }
}

void rem(NTL::zz_pX& r, const NTL::zz_pX& a, const zz_pXModulus1& ff)
{
  if (!ff.specialLogic) {
    rem(r, a, ff.fm);
    return;
  }

  long m = ff.m;
  long n = ff.n;
  long k = ff.k;
  long k1 = ff.k1;
  const NTL::fftRep& R0 = ff.R0;
  const NTL::fftRep& R1 = ff.R1;

  if (deg(a) < n) {
    r = a;
    return;
  }

  NTL::zz_pX P2, P3;

  NTL::fftRep R2, R3;

  // TofftRep(R2, a, k, n, m-1);
  TofftRep_trunc(R2, a, k, 2 * (m - 1 - n) + 1, n, m - 1);
  mul(R2, R2, R0);
  FromfftRep(P3, R2, m - 1 - n, 2 * (m - 1 - n));

  long l = 1L << k1;

  TofftRep(R3, P3, k1);
  mul(R3, R3, R1);
  FromfftRep(P3, R3, 0, n - 1);
  LocalCyclicReduce(P2, a, l);
  trunc(P2, P2, n);
  sub(P2, P2, P3);
  r = P2;
}

} // namespace helib
