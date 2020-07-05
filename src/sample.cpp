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
/* sample.cpp - implementing various sampling routines */
#include <vector>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/BasicThreadPool.h>
#include <helib/NumbTh.h>
#include <helib/Context.h>
#include <helib/sample.h>
#include <helib/norms.h>
#include <helib/apiAttributes.h>

#include <helib/powerful.h>
// only used in experimental Hwt sampler

namespace helib {

// Sample a degree-(n-1) poly, with only Hwt nonzero coefficients
void sampleHWt(zzX& poly, long n, long Hwt)
{
  if (n <= 0)
    n = lsize(poly);
  if (n <= 0)
    return;
  if (Hwt >= n) {
#ifdef HELIB_DEBUG
    std::cerr << "Hwt=" << Hwt << ">=n=" << n << ", is this ok?\n";
#endif
    Hwt = n - 1;
  }
  poly.SetLength(n); // allocate space
  for (long i = 0; i < n; i++)
    poly[i] = 0;

  long i = 0;
  while (i < Hwt) { // continue until exactly Hwt nonzero coefficients
    long u = NTL::RandomBnd(n);             // The next coefficient to choose
    if (poly[u] == 0) {                     // if we didn't choose it already
      long b = NTL::RandomBits_long(2) & 2; // b random in {0,2}
      poly[u] = b - 1;                      //   random in {-1,1}

      i++; // count another nonzero coefficient
    }
  }
}
// Sample a degree-(n-1) NTL::ZZX, with only Hwt nonzero coefficients
void sampleHWt(NTL::ZZX& poly, long n, long Hwt)
{
  zzX pp;
  sampleHWt(pp, n, Hwt);
  convert(poly, pp);
}

// Sample a degree-(n-1) poly, with -1/0/+1 coefficients.
// Each coefficients is +-1 with probability prob/2 each,
// and 0 with probability 1-prob. By default, pr[nonzero]=1/2.
void sampleSmall(zzX& poly, long n, double prob)
{
  if (n <= 0)
    n = lsize(poly);
  if (n <= 0)
    return;
  // prob must be in [2^{-15}, 1]
  assertInRange<InvalidArgument>(prob,
                                 3.05e-5,
                                 1.0,
                                 "prob must be between 2^{-15}"
                                 " and 1 inclusive",
                                 true);
  poly.SetLength(n);

  constexpr long bitSize = 16;
  constexpr long hiMask = (1 << (bitSize - 1)); // top bit = 2^15
  constexpr long loMask = hiMask - 1;           // bottom 15 bits

  long threshold = round(hiMask * prob); // threshold/2^15 = Pr[nonzero]

  NTL_EXEC_RANGE(n, first, last)
  for (long i = first; i < last; i++) {
    long u = NTL::RandomBits_long(bitSize); // a random 16-bit number
    long uLo = u & loMask;                  // bottom 15 bits
    long uHi = u & hiMask;                  // top bit

    // with probability threshold/2^15, choose between +-1
    if (uLo < threshold) {                  // compare low 15 bits to threshold
      poly[i] = (uHi >> (bitSize - 2)) - 1; // topBit*2 - 1 \in {+-1}
    }

    // with probability 1-prob, set to zero
    else
      poly[i] = 0;
  }
  NTL_EXEC_RANGE_END
}
void sampleSmall(NTL::ZZX& poly, long n, double prob)
{
  zzX pp;
  sampleSmall(pp, n, prob);
  convert(poly.rep, pp);
  poly.normalize();
}

// Choose a vector of continuous Gaussians
void sampleGaussian(std::vector<double>& dvec, long n, double stdev)
{
  static double const bignum = static_cast<double>(LONG_MAX);
  // THREADS: C++11 guarantees these are initialized only once

  if (n <= 0)
    n = lsize(dvec);
  if (n <= 0)
    return;
  dvec.resize(n, 0.0); // allocate space for n variables

  // Uses the Box-Muller method to get two Normal(0,stdev^2) variables
  for (long i = 0; i < n; i += 2) {
    // r1, r2 are "uniform in (0,1)"
    double r1 = (1 + NTL::RandomBnd(LONG_MAX)) / (bignum + 1);
    double r2 = (1 + NTL::RandomBnd(LONG_MAX)) / (bignum + 1);
    double theta = 2.0L * PI * r1;
    double rr = sqrt(-2.0 * log(r2)) * stdev;
    if (rr > 8 * stdev) // sanity-check, truncate at 8 standard deviations
      rr = 8 * stdev;

    // Generate two Gaussians RV's
    dvec[i] = rr * cos(theta);
    if (i + 1 < n)
      dvec[i + 1] = rr * sin(theta);
  }
}

// Sample a degree-(n-1) NTL::ZZX, with rounded Gaussian coefficients
void sampleGaussian(zzX& poly, long n, double stdev)
{
  if (n <= 0)
    n = lsize(poly);
  if (n <= 0)
    return;
  std::vector<double> dvec;
  sampleGaussian(dvec, n, stdev); // sample continuous Gaussians

  // round and copy to coefficients of poly
  clear(poly);
  poly.SetLength(n); // allocate space for degree-(n-1) polynomial
  for (long i = 0; i < n; i++)
    poly[i] = long(round(dvec[i])); // round to nearest integer
}
// Sample a degree-(n-1) NTL::ZZX, with rounded Gaussian coefficients
void sampleGaussian(NTL::ZZX& poly, long n, double stdev)
{
  zzX pp;
  sampleGaussian(pp, n, stdev);
  convert(poly.rep, pp);
  poly.normalize();
}
#if 0
void sampleGaussian(NTL::ZZX &poly, long n, double stdev)
{
  static long const bignum = 0xfffffff;
  // THREADS: C++11 guarantees these are initialized only once

  if (n<=0) n=deg(poly)+1; if (n<=0) return;
  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial
  for (long i=0; i<n; i++) SetCoeff(poly, i, NTL::ZZ::zero());

  // Uses the Box-Muller method to get two Normal(0,stdev^2) variables
  for (long i=0; i<n; i+=2) {
    double r1 = (1+RandomBnd(bignum))/((double)bignum+1);
    double r2 = (1+RandomBnd(bignum))/((double)bignum+1);
    double theta = 2.0L * PI * r1;
    double rr= sqrt(-2.0*log(r2))*stdev;

    // sanity-check, no more than 8 standard deviations
    assertTrue(rr < 8*stdev, "no more than 8 standard deviations");

    // Generate two Gaussians RV's, rounded to integers
    long x = (long) floor(rr*cos(theta) +0.5);
    SetCoeff(poly, i, x);
    if (i+1 < n) {
      x = (long) floor(rr*sin(theta) +0.5);
      SetCoeff(poly, i+1, x);
    }
  }
  poly.normalize(); // need to call this after we work on the coeffs
}
#endif

// Sample a degree-(n-1) zzX, with coefficients uniform in [-B,B]
void sampleUniform(zzX& poly, long n, long B)
{
  assertTrue<InvalidArgument>(B > 0l, "Invalid coefficient interval");
  if (n <= 0)
    n = lsize(poly);
  if (n <= 0)
    return;
  poly.SetLength(n); // allocate space for degree-(n-1) polynomial

  for (long i = 0; i < n; i++)
    poly[i] = NTL::RandomBnd(2 * B + 1) - B;
}

// Sample a degree-(n-1) NTL::ZZX, with coefficients uniform in [-B,B]
void sampleUniform(NTL::ZZX& poly, long n, const NTL::ZZ& B)
{
  assertTrue<InvalidArgument>(static_cast<bool>(B > 0l),
                              "Invalid coefficient interval");
  if (n <= 0)
    n = deg(poly) + 1;
  if (n <= 0)
    return;
  clear(poly);
  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial

  NTL::ZZ UB = 2 * B + 1;
  for (long i = n - 1; i >= 0; i--) {
    NTL::ZZ tmp = RandomBnd(UB) - B;
    SetCoeff(poly, i, tmp);
  }
}

/********************************************************************
 * Below are versions of the sampling routines that sample modulo
 * X^m-1 and then reduce mod Phi_m(X). The exception is when m is
 * a power of two, where we still sample directly mod Phi_m(X).
 ********************************************************************/
double sampleHWt(zzX& poly, const Context& context, long Hwt)
{
  const PAlgebra& palg = context.zMStar;
  double retval;

  if (palg.getPow2() == 0) { // not power of two
    long m = palg.getM();
    sampleHWt(poly, m, Hwt);
    reduceModPhimX(poly, palg);
    retval = context.noiseBoundForHWt(Hwt, m);
  } else { // power of two
    long phim = palg.getPhiM();
    sampleHWt(poly, phim, Hwt);
    retval = context.noiseBoundForHWt(Hwt, phim);
  }

  return retval;
}

double sampleHWtBoundedEffectiveBound(const Context& context, long Hwt)
{
  const PAlgebra& palg = context.zMStar;
  long phim = palg.getPhiM();

  double bound = sqrt(Hwt * log(phim));
  // should be good with probability at least 1/2
  // NOTE: the general formula is sigma*sqrt(log(phim)),
  // assuming we are sampling from a zero mean complex Gaussian
  // with std deviation sigma

  return bound;
}

double sampleHWtBounded(zzX& poly, const Context& context, long Hwt)
{
  double bound = sampleHWtBoundedEffectiveBound(context, Hwt);
  const PAlgebra& palg = context.zMStar;

#if 1
  double val;
  long count = 0;
  do {
    sampleHWt(poly, context, Hwt);
    val = embeddingLargestCoeff(poly, palg);
  } while (++count < 1000 && val > bound); // repeat until <= bound
#else
  double val, val1;
  zzX poly1;
  long succ = 0;

  long count = 100;

  val = 1e9;

  cout << "sampleHWtBounded:\n";

  for (long i = 0; i < count; i++) {
    sampleHWt(poly1, context, Hwt);
    val1 = embeddingLargestCoeff(poly1, palg);
    if (val1 <= bound) {
      succ++;
      val = val1;
      poly = poly1;
      cout << "*";
    } else {
      cout << ".";
    }
  }

  cout << "\nsucc%=" << ((double(succ) / count) * 100) << "\n";
  cout << "bound=" << bound << "\n";
  cout << "Hwt=" << Hwt << "\n";

#endif

  if (val > bound) {
    std::stringstream ss;
    ss << "Error: sampleHWtBounded, after " << count
       << " trials, still val=" << val << '>' << "bound=" << bound;
    throw RuntimeError(ss.str());
  }
  return bound;
}

double sampleSmall(zzX& poly, const Context& context)
{
  const PAlgebra& palg = context.zMStar;
  double retval;

  if (palg.getPow2() == 0) { // not power of two
    long m = palg.getM();
    long phim = palg.getPhiM();
    sampleSmall(poly, m, phim / (2.0 * m)); // nonzero with prob phi(m)/2m
    // FIXME: does this probability make sense?  What is the goal there??
    reduceModPhimX(poly, palg);
    retval = context.noiseBoundForSmall(phim / (2.0 * m), m);
  } else { // power of two
    long phim = palg.getPhiM();
    sampleSmall(poly, phim);
    retval = context.noiseBoundForSmall(0.5, phim);
  }

  return retval;
}

// Same as above, but ensure the result is not too much larger than typical
double sampleSmallBounded(zzX& poly, const Context& context)
{
  const PAlgebra& palg = context.zMStar;
  long phim = palg.getPhiM();

  double bound = sqrt(phim * log(phim) / 2.0);
  // should be good with probability at least 1/2
  // NOTE: the general formula is sigma*sqrt(log(phim)),
  // assuming we are sampling from a zero mean complex Gaussian
  // with std deviation sigma

#if 1
  double val;
  long count = 0;
  do {
    sampleSmall(poly, context);
    val = embeddingLargestCoeff(poly, palg);
  } while (++count < 1000 && val > bound); // repeat until <= bound
#else
  double val, val1;
  zzX poly1;
  long succ = 0;

  long count = 100;

  val = 1e9;

  cout << "sampleSmallBounded:\n";

  for (long i = 0; i < count; i++) {
    sampleSmall(poly1, context);
    val1 = embeddingLargestCoeff(poly1, palg);
    if (val1 <= bound) {
      succ++;
      val = val1;
      poly = poly1;
      cout << "*";
    } else {
      cout << ".";
    }
  }

  cout << "\nsucc%=" << ((double(succ) / count) * 100) << "\n";
  cout << "bound=" << bound << "\n";

#endif

  if (val > bound) {
    std::stringstream ss;
    ss << "Error: sampleSmallBounded, after " << count
       << " trials, still val=" << val << '>' << "bound=" << bound;
    throw RuntimeError(ss.str());
  }
  return bound;
}

double sampleGaussian(zzX& poly, const Context& context, double stdev)
{
  const PAlgebra& palg = context.zMStar;
  double retval;

  if (palg.getPow2() == 0) { // not power of two
    long m = palg.getM();
    sampleGaussian(poly, m, stdev);
    reduceModPhimX(poly, palg);
    retval = context.noiseBoundForGaussian(stdev, m);
  } else { // power of two
    long phim = palg.getPhiM();
    sampleGaussian(poly, phim, stdev);
    retval = context.noiseBoundForGaussian(stdev, phim);
  }

  return retval;
}
// Same as above, but ensure the result is not too much larger than typical
double sampleGaussianBounded(zzX& poly, const Context& context, double stdev)
{
  const PAlgebra& palg = context.zMStar;
  long m = palg.getM();
  long phim = palg.getPhiM();

  double bound = stdev * (((palg.getPow2() == 0) ? sqrt(m * log(phim))
                                                 : sqrt(phim * log(phim))));
  // should be good with probability at least 1/2
  // NOTE: the general formula is sigma*sqrt(log(phim)),
  // assuming we are sampling from a zero mean complex Gaussian
  // with std deviation sigma

#if 1
  double val;
  long count = 0;
  do {
    sampleGaussian(poly, context, stdev);
    val = embeddingLargestCoeff(poly, palg);
  } while (++count < 1000 && val > bound); // repeat until <=bound
#else
  double val, val1;
  zzX poly1;
  long succ = 0;

  long count = 100;

  val = 1e9;

  std::cout << "sampleGaussianBounded:\n";

  for (long i = 0; i < count; i++) {
    sampleGaussian(poly1, context, stdev);
    val1 = embeddingLargestCoeff(poly1, palg);
    if (val1 <= bound) {
      succ++;
      val = val1;
      poly = poly1;
      std::cout << "*";
    } else {
      std::cout << ".";
    }
  }

  std::cout << "\nsucc%=" << ((double(succ) / count) * 100) << "\n";
  std::cout << "bound=" << bound << "\n";
  std::cout << "val=" << val << "\n";

#endif

  if (val > bound) {
    std::stringstream ss;
    ss << "Error: sampleGaussianBounded, after " << count
       << " trials, still val=" << val << '>' << "bound=" << bound;
    throw RuntimeError(ss.str());
  }
  return bound;
}

double sampleUniform(zzX& poly, const Context& context, long B)
{
  const PAlgebra& palg = context.zMStar;
  double retval;

  if (palg.getPow2() == 0) { // not power of two
    long m = palg.getM();
    sampleUniform(poly, m, B);
    reduceModPhimX(poly, palg);
    retval = context.noiseBoundForUniform(B, m);
  } else { // power of two
    long phim = palg.getPhiM();
    sampleUniform(poly, phim, B);
    retval = context.noiseBoundForUniform(B, phim);
  }

  return retval;
}

NTL::xdouble sampleUniform(NTL::ZZX& poly,
                           const Context& context,
                           const NTL::ZZ& B)
{
  const PAlgebra& palg = context.zMStar;
  NTL::xdouble retval;

  if (palg.getPow2() == 0) { // not power of two
    long m = palg.getM();
    sampleUniform(poly, m, B);
    NTL::rem(poly, poly, palg.getPhimX());
    retval = context.noiseBoundForUniform(NTL::conv<NTL::xdouble>(B), m);
  } else { // power of two
    long phim = palg.getPhiM();
    sampleUniform(poly, phim, B);
    retval = context.noiseBoundForUniform(NTL::conv<NTL::xdouble>(B), phim);
  }

  return retval;
}

// Helper functions, return a bound B such that for random noise
// terms we have Pr[|canonicalEmbed(noise)|_{\infty} > B] < epsilon.
// (The default is epsilon = 2^{-40}.)
double boundFreshNoise(long m, long phim, double sigma, double epsilon)
{
  // The various constants in this function were determined experimentally.

  /* We begin by computing the standard deviation of the magnitude of
   * f(zeta), where f = sampleSmallBounded * sampleGaussian(sigma)
   *                 +  sampleSmall * sampleGaussianBounded(sigma)
   *                 + sampleGaussian(sigma)
   * and zeta is an m'th root-of-unity, which we approximate as:
   */
  double stdev =
      (sigma + 0.1) * 0.54 * (1 + (is2power(m) ? sqrt(phim * m) : phim));

  /* Then we use the following rules:
   *      Pr[|f(zeta)| > stdev] = 0.644
   *      Pr[|f(zeta)| > 2*stdev] = 0.266
   *      Pr[|f(zeta)| > 3*stdev] = 9.04e-2
   *      Pr[|f(zeta)| > 4*stdev] = 2.76e-2
   *      Pr[|f(zeta)| > 5*stdev] = 7.89e-3
   *      Pr[|f(zeta)| > 6*stdev] = 2.16e-3
   *      Pr[|f(zeta)| > n*stdev] = 2.16e-3 * 4^{6-n} for n>6
   *
   * We return the smallest number of standard deviations n satisfying
   *      Pr[|f(zeta)|>(n stdev)] = epsilon / phi(m)
   */
  epsilon /= phim;

  if (epsilon >= 1.87e-3) { // use the values from above
    if (epsilon >= 0.64) {
      return stdev;
    } else if (epsilon >= 0.26) {
      return 2 * stdev;
    } else if (epsilon >= 8.52e-2) {
      return 3 * stdev;
    } else if (epsilon >= 2.54e-2) {
      return 4 * stdev;
    } else if (epsilon >= 7.06e-3) {
      return 5 * stdev;
    } else {
      return 6 * stdev;
    }
  }
  long num = 7;
  for (double prob = (1.87e-3) / 4; prob > epsilon; prob /= 4)
    num++;

  return stdev * num;
}
double boundRoundingNoise(UNUSED long m, long phim, long p2r, double epsilon)
{
  // The various constants in this function were determined experimentally.

  /* We begin by computing the standard deviation of the magnitude of
   * f(zeta), where
   *          f= sampleSmallBounded*sampleUniform(p) +sampleUniform(p),
   * and zeta is an m'th root-of-unity, which we approximate as:
   */
  double stdev = (2 * p2r + 1) * (phim - 2) / 8.0;

  /* Then we use the following rules:
   *      Pr[|f(zeta)| > stdev] = 0.514
   *      Pr[|f(zeta)| > 2*stdev] = 0.194
   *      Pr[|f(zeta)| > 3*stdev] = 6.7e-2
   *      Pr[|f(zeta)| > 4*stdev] = 2.23e-2
   *      Pr[|f(zeta)| > 5*stdev] = 7.21e-3
   *      Pr[|f(zeta)| > 6*stdev] = 2.31e-3
   *      Pr[|f(zeta)| > 7*stdev] = 7.25e-4
   *      Pr[|f(zeta)| > n*stdev] = 7.25e-4 * 3.3^{7-n} for n>5
   *
   * We return the smallest number of standard deviations n satisfying
   *      Pr[|f(zeta)|>(n stdev)] = epsilon / phi(m)
   */
  epsilon /= phim;

  if (epsilon >= 7.25e-4) { // use the values from above
    if (epsilon >= 0.514) {
      return stdev;
    } else if (epsilon >= 0.194) {
      return 2 * stdev;
    } else if (epsilon >= 6.7e-2) {
      return 3 * stdev;
    } else if (epsilon >= 2.23e-2) {
      return 4 * stdev;
    } else if (epsilon >= 7.21e-3) {
      return 5 * stdev;
    } else if (epsilon >= 2.31e-3) {
      return 6 * stdev;
    } else {
      return 7 * stdev;
    }
  }
  long num = 8;
  for (double prob = (7.25e-4) / 3.2; prob > epsilon; prob /= 3.2)
    num++;

  return stdev * num;
}

} // namespace helib
