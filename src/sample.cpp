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
#include <cassert>
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include "sample.h"
NTL_CLIENT

// MinGW hack
#ifndef lrand48
#if defined(__MINGW32__) || defined(WIN32)
#define drand48() (((double)rand()) / RAND_MAX)
#define lrand48() rand()
#endif
#endif

void sampleHWt(ZZX &poly, long Hwt, long n)
{
  if (n<=0) n=deg(poly)+1; if (n<=0) return;
  clear(poly);          // initialize to zero
  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial

  long b,u,i=0;
  if (Hwt>n) Hwt=n;
  while (i<Hwt) {  // continue until exactly Hwt nonzero coefficients
    u=lrand48()%n; // The next coefficient to choose
    if (IsZero(coeff(poly,u))) { // if we didn't choose it already
      b = lrand48()&2; // b random in {0,2}
      b--;             //   random in {-1,1}
      SetCoeff(poly,u,b);

      i++; // count another nonzero coefficient
    }
  }
  poly.normalize(); // need to call this after we work on the coeffs
}

void sampleSmall(ZZX &poly, long n)
{
  if (n<=0) n=deg(poly)+1; if (n<=0) return;
  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial

  for (long i=0; i<n; i++) {    // Chosse coefficients, one by one
    long u = lrand48();
    if (u&1) {                 // with prob. 1/2 choose between -1 and +1
      u = (u & 2) -1;
      SetCoeff(poly, i, u);
    }
    else SetCoeff(poly, i, 0); // with ptob. 1/2 set to 0
  }
  poly.normalize(); // need to call this after we work on the coeffs
}

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

void sampleUniform(ZZX& poly, const ZZ& B, long n)
{
  if (n<=0) n=deg(poly)+1; if (n<=0) return;
  if (B <= 0) {
    clear(poly);
    return;
  }

  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial

  ZZ UB, tmp;

  UB =  2*B + 1;
  for (long i = 0; i < n; i++) {
    RandomBnd(tmp, UB);
    tmp -= B; 
    poly.rep[i] = tmp;
  }

  poly.normalize();
}
