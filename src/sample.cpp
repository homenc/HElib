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
/* sample.cpp - implementing various sampling routines */
#include <vector>
#include <cassert>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/BasicThreadPool.h>
#include "NumbTh.h"
#include "PAlgebra.h"
#include "sample.h"
NTL_CLIENT

// Sample a degree-(n-1) poly, with only Hwt nonzero coefficients
void sampleHWt(zzX &poly, long n, long Hwt)
{
  if (n<=0) n=lsize(poly); if (n<=0) return;
  if (Hwt>=n) {
#ifdef DEBUG_PRINTOUT
    std::cerr << "Hwt="<<Hwt<<">=n="<<n<<", is this ok?\n";
#endif
    Hwt = n-1;
  }
  poly.SetLength(n); // allocate space
  for (long i=0; i<n; i++) poly[i] = 0;

  long i=0;
  while (i<Hwt) {  // continue until exactly Hwt nonzero coefficients
    long u = NTL::RandomBnd(n);  // The next coefficient to choose
    if (poly[u]==0) { // if we didn't choose it already
      long b = NTL::RandomBits_long(2)&2; // b random in {0,2}
      poly[u] = b-1;                      //   random in {-1,1}

      i++; // count another nonzero coefficient
    }
  }
}
// Sample a degree-(n-1) ZZX, with only Hwt nonzero coefficients
void sampleHWt(ZZX &poly, long n, long Hwt)
{
  zzX pp;
  sampleHWt(pp, n, Hwt);
  convert(poly, pp);
}

// Sample a degree-(n-1) poly, with -1/0/+1 coefficients.
// Each coefficients is +-1 with probability prob/2 each,
// and 0 with probability 1-prob. By default, pr[nonzero]=1/2.
void sampleSmall(zzX &poly, long n, double prob)
{
  if (n<=0) n=lsize(poly); if (n<=0) return;
  assert(prob>3.05e-5 && prob<=1); // prob must be in [2^{-15},1/2]
  poly.SetLength(n);

  constexpr long bitSize=16;
  constexpr long hiMask = (1<<(bitSize-1)); // top bit = 2^15
  constexpr long loMask = hiMask-1;         // bottom 15 bits

  long threshold = round(hiMask*prob); // threshold/2^15 = Pr[nonzero]

  NTL_EXEC_RANGE(n, first, last)
  for (long i=first; i<last; i++) {
    long u = NTL::RandomBits_long(bitSize); // a random 16-bit number
    long uLo = u & loMask; // bottom 15 bits
    long uHi = u & hiMask; // top bit

    // with probability threshold/2^15, choose between +-1
    if (uLo<threshold) { // compare low 15 bits to threshold
      poly[i] = (uHi>>(bitSize-2))-1; // topBit*2 - 1 \in {+-1}
    }

    // with probability 1-prob, set to zero
    else poly[i] = 0;
  }
  NTL_EXEC_RANGE_END
}
void sampleSmall(ZZX &poly, long n, double prob)
{  
  zzX pp;
  sampleSmall(pp, n, prob);
  convert(poly.rep, pp);
  poly.normalize();
}

// Choose a vector of continuous Gaussians
void sampleGaussian(std::vector<double> &dvec, long n, double stdev)
{
  static double const Pi=4.0*atan(1.0);  // Pi=3.1415..
  static double const bignum = LONG_MAX; // convert to double
  // THREADS: C++11 guarantees these are initialized only once

  if (n<=0) n=lsize(dvec); if (n<=0) return;
  dvec.resize(n, 0.0);        // allocate space for n variables

  // Uses the Box-Muller method to get two Normal(0,stdev^2) variables
  for (long i=0; i<n; i+=2) {
    // r1, r2 are "uniform in (0,1)"
    double r1 = (1+NTL::RandomBnd(LONG_MAX))/(bignum+1);
    double r2 = (1+NTL::RandomBnd(LONG_MAX))/(bignum+1);
    double theta=2*Pi*r1;
    double rr= sqrt(-2.0*log(r2))*stdev;
    if (rr > 8*stdev) // sanity-check, trancate at 8 standard deviations
      rr = 8*stdev;

    // Generate two Gaussians RV's
    dvec[i] = rr*cos(theta);
    if (i+1 < n)
      dvec[i+1] = rr*sin(theta);
  }
}

// Sample a degree-(n-1) ZZX, with rounded Gaussian coefficients
void sampleGaussian(zzX &poly, long n, double stdev)
{
  if (n<=0) n=lsize(poly); if (n<=0) return;
  std::vector<double> dvec;
  sampleGaussian(dvec, n, stdev); // sample continuous Gaussians

  // round and copy to coefficients of poly
  clear(poly);
  poly.SetLength(n); // allocate space for degree-(n-1) polynomial
  for (long i=0; i<n; i++)
    poly[i] = long(round(dvec[i])); // round to nearest integer
}
// Sample a degree-(n-1) ZZX, with rounded Gaussian coefficients
void sampleGaussian(ZZX &poly, long n, double stdev)
{
  zzX pp;
  sampleGaussian(pp, n, stdev);
  convert(poly.rep, pp);
  poly.normalize();
}
#if 0
void sampleGaussian(ZZX &poly, long n, double stdev)
{
  static double const Pi=4.0*atan(1.0); // Pi=3.1415..
  static long const bignum = 0xfffffff;
  // THREADS: C++11 guarantees these are initialized only once

  if (n<=0) n=deg(poly)+1; if (n<=0) return;
  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial
  for (long i=0; i<n; i++) SetCoeff(poly, i, ZZ::zero());

  // Uses the Box-Muller method to get two Normal(0,stdev^2) variables
  for (long i=0; i<n; i+=2) {
    double r1 = (1+RandomBnd(bignum))/((double)bignum+1);
    double r2 = (1+RandomBnd(bignum))/((double)bignum+1);
    double theta=2*Pi*r1;
    double rr= sqrt(-2.0*log(r2))*stdev;

    assert(rr < 8*stdev); // sanity-check, no more than 8 standard deviations

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
  assert (B>0);
  if (n<=0) n=lsize(poly); if (n<=0) return;
  poly.SetLength(n); // allocate space for degree-(n-1) polynomial

  for (long i = 0; i < n; i++)
    poly[i] = NTL::RandomBnd(2*B +1) - B;
}

// Sample a degree-(n-1) ZZX, with coefficients uniform in [-B,B]
void sampleUniform(ZZX& poly, long n, const ZZ& B)
{
  assert (B>0);
  if (n<=0) n=deg(poly)+1; if (n<=0) return;
  clear(poly);
  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial

  ZZ UB = 2*B +1;
  for (long i = n-1; i >= 0; i--) {
    ZZ tmp = RandomBnd(UB) - B;
    SetCoeff(poly, i, tmp);
  }
}


/********************************************************************
 * Below are versions of the sampling routines that sample modulo
 * X^m-1 and then reduce mod Phi_m(X). The exception is when m is
 * a power of two, where we still sample directly mod Phi_m(X).
 ********************************************************************/
void sampleHWt(zzX &poly, const PAlgebra& palg, long Hwt)
{
  if (palg.getPow2() == 0) { // not power of two
    sampleHWt(poly, palg.getM(), Hwt);
    reduceModPhimX(poly, palg);
  }
  else // power of two
    sampleHWt(poly, palg.getPhiM(), Hwt);
}
void sampleSmall(zzX &poly, const PAlgebra& palg)
{
  if (palg.getPow2() == 0) { // not power of two
    long m = palg.getM();
    long phim = palg.getPhiM();
    sampleSmall(poly, m, phim/(2.0*m)); // nonzero with prob phi(m)/2m
    reduceModPhimX(poly, palg);
  }
  else // power of two
    sampleSmall(poly, palg.getPhiM());
}
void sampleGaussian(zzX &poly, const PAlgebra& palg, double stdev)
{
  if (palg.getPow2() == 0) { // not power of two
    sampleGaussian(poly, palg.getM(), stdev);
    reduceModPhimX(poly, palg);
  }
  else // power of two
    sampleGaussian(poly, palg.getPhiM(), stdev);
}
void sampleUniform(zzX &poly, const PAlgebra& palg, long B)
{
  if (palg.getPow2() == 0) { // not power of two
    sampleUniform(poly, palg.getM(), B);
    reduceModPhimX(poly, palg);
  }
  else // power of two
    sampleUniform(poly, palg.getPhiM(), B);
}
void sampleUniform(ZZX &poly, const PAlgebra& palg, const ZZ& B)
{
  if (palg.getPow2() == 0) { // not power of two
    sampleUniform(poly, palg.getM(), B);
    NTL::rem(poly, poly, palg.getPhimX());
  }
  else // power of two
    sampleUniform(poly, palg.getPhiM(), B);
}


// Helper functions, return a bound B such that for random noise
// terms we have Pr[|canonicalEmbed(noise)|_{\infty} > B] < epsilon.
// (The default is epsilon = 2^{-40}.)
double boundFreshNoise(long m, long phim, double sigma, double epsilon)
{
  // The various constants in this function were determined experimentally.

  /* We begin by computing the standard deviation of the magnitude of
   * f(zeta), where f = sampleSmall*sampleGaussian(sigma)
   *                 +  sampleSmall*sampleGaussian(sigma) + sampleGaussian(sigma)
   * and zeta is an m'th root-of-unity, which we approximate as:
   */
  double stdev = 0.56*sigma*((m&1)? sqrt(m*phim) : phim);

  /* Then we use the following rules:
   *      Pr[|f(zeta)| > stdev] = 0.63
   *      Pr[|f(zeta)| > 2*stdev] = 0.25
   *      Pr[|f(zeta)| > 3*stdev] = 8.14e-2
   *      Pr[|f(zeta)| > 4*stdev] = 2.38e-2
   *      Pr[|f(zeta)| > 5*stdev] = 6.5e-3
   *      Pr[|f(zeta)| > n*stdev] = 6.5e-3 * 4^{5-n} for n>5
   *
   * We return the smallest number of standard deviations n satifying
   *      Pr[|f(zeta)|>(n stdev)] = epsilon / phi(m)
   */
  epsilon /= phim;
}
double boundRoundingNoise(long m, long phim, long p2r, double epsilon)
{
  // The various constants in this function were determined experimentally.

  /* We begin by computing the standard deviation of the magnitude of
   * f(zeta), where f= sampleSmall*sampleUniform(p) +sampleUniform(p),
   * and zeta is an m'th root-of-unity, which we approximate as:
   */
  double stdev = (2*p2r+1)*phim/8.0;

  /* Then we use the following rules:
   *      Pr[|f(zeta)| > stdev] = 0.508
   *      Pr[|f(zeta)| > 2*stdev] = 0.189
   *      Pr[|f(zeta)| > 3*stdev] = 6.46e-2
   *      Pr[|f(zeta)| > 4*stdev] = 2.11e-2
   *      Pr[|f(zeta)| > 5*stdev] = 6.75e-3
   *      Pr[|f(zeta)| > n*stdev] = 6.75e-3 * 3.2^{5-n} for n>5
   *
   * We return the smallest number of standard deviations n satifying
   *      Pr[|f(zeta)|>(n stdev)] = epsilon / phi(m)
   */
  epsilon /= phim;

  if (epsilon >= 6.75e-3) { // use the values from above
    if (epsilon >= 0.508)        { return stdev; }
    else if (epsilon >= 0.189)   { return 2*stdev; }
    else if (epsilon >= 6.46e-2) { return 3*stdev; }
    else if (epsilon >= 2.11e-2) { return 4*stdev; }
    else                         { return 5*stdev; }
  }
  long num = 6;
  for (double prob=(6.75e-3)/(3.2); prob>epsilon; prob /= 3.2)
    num++;

  return stdev * num;
}
