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
/* EncryptedArray.cpp - Data-movement operations on arrays of slots
 */
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <helib/EncryptedArray.h>
#include <helib/polyEval.h>
#include <helib/debugging.h>

namespace helib {

// Compute a degree-p polynomial poly(x) s.t. for any t<e and integer z of the
// form z = z0 + p^t*z1 (with 0<=z0<p), we have poly(z) = z0 (mod p^{t+1}).
//
// We get poly(x) by interpolating a degree-(p-1) polynomial poly'(x)
// s.t. poly'(z0)=z0 - z0^p (mod p^e) for all 0<=z0<p, and then setting
// poly(x) = x^p + poly'(x).
static void buildDigitPolynomial(NTL::ZZX& result, long p, long e)
{
  if (p < 2 || e <= 1)
    return; // nothing to do
  HELIB_TIMER_START;
  long p2e = NTL::power_long(p, e); // the integer p^e

  // Compute x - x^p (mod p^e), for x=0,1,...,p-1
  NTL::vec_long x(NTL::INIT_SIZE, p);
  NTL::vec_long y(NTL::INIT_SIZE, p);
  long bottom = -(p / 2);
  for (long j = 0; j < p; j++) {
    long z = bottom + j;
    x[j] = z;
    y[j] =
        z - NTL::PowerMod((z < 0 ? z + p2e : z), p, p2e); // x - x^p (mod p^e)

    while (y[j] > p2e / 2)
      y[j] -= p2e;
    while (y[j] < -(p2e / 2))
      y[j] += p2e;
  }
  interpolateMod(result, x, y, p, e);
  // interpolating p points, should get deg<=p-1
  assertTrue(deg(result) < p, "Interpolation error.  Degree too high.");
  SetCoeff(result, p); // return result = x^p + poly'(x)
  //  cerr << "# digitExt mod "<<p<<"^"<<e<<"="<<result<<endl;
  HELIB_TIMER_STOP;
}

// extractDigits assumes that the slots of *this contains integers mod p^r
// i.e., that only the free terms are nonzero. (If that assumptions does
// not hold then the result will not be a valid ciphertext anymore.)
//
// It returns in the slots of digits[j] the j'th-lowest digits from the
// integers in the slots of the input. Namely, the i'th slot of digits[j]
// contains the j'th digit in the p-base expansion of the integer in the
// i'th slot of the *this. The plaintext space of digits[j] is mod p^{r-j},
// and all the digits are at the same level.

int fhe_watcher = 0;

void extractDigits(std::vector<Ctxt>& digits, const Ctxt& c, long r)
{
  const Context& context = c.getContext();
  long rr = c.effectiveR();
  if (r <= 0 || r > rr)
    r = rr; // how many digits to extract

  long p = context.zMStar.getP();

  NTL::ZZX x2p;
  if (p > 3) {
    buildDigitPolynomial(x2p, p, r);
  }

  Ctxt tmp(c.getPubKey(), c.getPtxtSpace());
  digits.resize(r, tmp); // allocate space

#ifdef HELIB_DEBUG
  fprintf(stderr, "***\n");
#endif
  for (long i = 0; i < r; i++) {
    tmp = c;
    for (long j = 0; j < i; j++) {

      if (p == 2)
        digits[j].square();
      else if (p == 3)
        digits[j].cube();
      else
        polyEval(digits[j], x2p, digits[j]);
        // "in spirit" digits[j] = digits[j]^p

#ifdef HELIB_DEBUG
      fprintf(stderr, "%5ld", digits[j].bitCapacity());
#endif

      tmp -= digits[j];
      tmp.divideByP();
    }
    digits[i] = tmp; // needed in the next round

#ifdef HELIB_DEBUG
    if (dbgKey) {
      double ratio = log(embeddingLargestCoeff(digits[i], *dbgKey) /
                         digits[i].getNoiseBound()) /
                     log(2.0);
      fprintf(stderr, "%5ld [%f]", digits[i].bitCapacity(), ratio);
      if (ratio > 0)
        fprintf(stderr, " BAD-BOUND");
      fprintf(stderr, "\n");
    } else {
      fprintf(stderr, "%5ld\n", digits[i].bitCapacity());
    }
#endif
  }

#ifdef HELIB_DEBUG
  fprintf(stderr, "***\n");
#endif
}

static void compute_a_vals(NTL::Vec<NTL::ZZ>& a, long p, long e)
// computes a[m] = a(m)/m! for m = p..(e-1)(p-1)+1,
// as defined by Chen and Han.
// a.length() is set to (e-1)(p-1)+2

{
  NTL::ZZ p_to_e = NTL::power_ZZ(p, e);
  NTL::ZZ p_to_2e = NTL::power_ZZ(p, 2 * e);

  long len = (e - 1) * (p - 1) + 2;

  NTL::ZZ_pPush push(p_to_2e);

  NTL::ZZ_pX x_plus_1_to_p = power(NTL::ZZ_pX(NTL::INIT_MONO, 1) + 1, p);
  NTL::ZZ_pX denom =
      InvTrunc(x_plus_1_to_p - NTL::ZZ_pX(NTL::INIT_MONO, p), len);
  NTL::ZZ_pX poly = MulTrunc(x_plus_1_to_p, denom, len);
  poly *= p;

  a.SetLength(len);

  NTL::ZZ m_fac(1);
  for (long m = 2; m < p; m++) {
    m_fac = MulMod(m_fac, m, p_to_2e);
  }

  for (long m = p; m < len; m++) {
    m_fac = MulMod(m_fac, m, p_to_2e);
    NTL::ZZ c = rep(coeff(poly, m));
    NTL::ZZ d = GCD(m_fac, p_to_2e);
    if (d == 0 || d > p_to_e || c % d != 0)
      throw RuntimeError("cannot divide");
    NTL::ZZ m_fac_deflated = (m_fac / d) % p_to_e;
    NTL::ZZ c_deflated = (c / d) % p_to_e;
    a[m] = MulMod(c_deflated, InvMod(m_fac_deflated, p_to_e), p_to_e);
  }
}

// This computes Chen and Han's magic polynomial G, which
// has the property that G(x) = (x mod p) (mod p^e).
// Here, (x mod p) is in the interval [0,1] if p == 2,
// and otherwise, is in the interval (-p/2, p/2).
static void compute_magic_poly(NTL::ZZX& poly1, long p, long e)
{
  HELIB_TIMER_START;

  NTL::Vec<NTL::ZZ> a;

  compute_a_vals(a, p, e);

  NTL::ZZ p_to_e = NTL::power_ZZ(p, e);
  long len = (e - 1) * (p - 1) + 2;

  NTL::ZZ_pPush push(p_to_e);

  NTL::ZZ_pX poly(0);
  NTL::ZZ_pX term(1);
  NTL::ZZ_pX X(NTL::INIT_MONO, 1);

  poly = 0;
  term = 1;

  for (long m = 0; m < p; m++) {
    term *= (X - m);
  }

  for (long m = p; m < len; m++) {
    poly += term * NTL::conv<NTL::ZZ_p>(a[m]);
    term *= (X - m);
  }

  // replace poly by poly(X+(p-1)/2) for odd p
  if (p % 2 == 1) {
    NTL::ZZ_pX poly2(0);

    for (long i = deg(poly); i >= 0; i--)
      poly2 = poly2 * (X + (p - 1) / 2) + poly[i];

    poly = poly2;
  }

  poly = X - poly;
  poly1 = NTL::conv<NTL::ZZX>(poly);
}

// extendExtractDigits assumes that the slots of *this contains integers mod
// p^{r+e} i.e., that only the free terms are nonzero. (If that assumptions
// does not hold then the result will not be a valid ciphertext anymore.)
//
// It returns in the slots of digits[j] the j'th-lowest digits from the
// integers in the slots of the input. Namely, the i'th slot of digits[j]
// contains the j'th digit in the p-base expansion of the integer in the i'th
// slot of the *this.  The plaintext space of digits[j] is mod p^{e+r-j}.

void extendExtractDigits(std::vector<Ctxt>& digits,
                         const Ctxt& c,
                         long r,
                         long e)
{
  const Context& context = c.getContext();

  long p = context.zMStar.getP();
  NTL::ZZX x2p;
  if (p > 3) {
    buildDigitPolynomial(x2p, p, r);
  }

  // we should pre-compute this table
  // for i = 0..r-1, entry i is G_{e+r-i} in Chen and Han
  NTL::Vec<NTL::ZZX> G;
  G.SetLength(r);
  for (long i : range(r)) {
    compute_magic_poly(G[i], p, e + r - i);
  }

  std::vector<Ctxt> digits0;

  Ctxt tmp(c.getPubKey(), c.getPtxtSpace());

  digits.resize(r, tmp); // allocate space
  digits0.resize(r, tmp);

#ifdef HELIB_DEBUG
  fprintf(stderr, "***\n");
#endif
  for (long i : range(r)) {
    tmp = c;
    for (long j : range(i)) {
      if (digits[j].capacity() >= digits0[j].capacity()) {
        // optimization: digits[j] is better than digits0[j],
        // so just use it

        tmp -= digits[j];
#ifdef HELIB_DEBUG
        fprintf(stderr, "%5ld*", digits[j].bitCapacity());
#endif
      } else {
        if (p == 2)
          digits0[j].square();
        else if (p == 3)
          digits0[j].cube();
        else
          polyEval(digits0[j],
                   x2p,
                   digits0[j]); // "in spirit" digits0[j] = digits0[j]^p

        tmp -= digits0[j];
#ifdef HELIB_DEBUG
        fprintf(stderr, "%5ld ", digits0[j].bitCapacity());
#endif
      }
      tmp.divideByP();
    }
    digits0[i] = tmp; // needed in the next round
    polyEval(digits[i], G[i], tmp);

#ifdef HELIB_DEBUG
    if (dbgKey) {
      double ratio = log(embeddingLargestCoeff(digits[i], *dbgKey) /
                         digits[i].getNoiseBound()) /
                     log(2.0);
      fprintf(stderr,
              "%5ld  --- %5ld",
              digits0[i].bitCapacity(),
              digits[i].bitCapacity());
      fprintf(stderr, " [%f]", ratio);
      if (ratio > 0)
        fprintf(stderr, " BAD-BOUND");
      fprintf(stderr, "\n");
    } else {
      fprintf(stderr,
              "%5ld  --- %5ld\n",
              digits0[i].bitCapacity(),
              digits[i].bitCapacity());
    }
#endif
  }
}

} // namespace helib
