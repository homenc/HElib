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
#ifndef HELIB_NORMS_H
#define HELIB_NORMS_H
/**
 * @file norms.h - computing various norms of ring elements
 **/
#include <vector>
#include <complex>
#include <NTL/ZZX.h>
#include <NTL/xdouble.h>
//#include <NTL/RR.h>
#include "zzX.h"

class DoubleCRT;

long sumOfCoeffs(const zzX& f);         // = f(1)
NTL::ZZ sumOfCoeffs(const NTL::ZZX& f); // = f(1)
NTL::ZZ sumOfCoeffs(const DoubleCRT& f); // somewhat lame implementation

//! The L-infinity norm of an element (in coefficient representation)
template<class T>
double largestCoeff(const NTL::Vec<T>& f)
{
  double mx = 0;
  for (auto& x : f) {
    auto sz = abs(x);
    if (mx < sz)
      mx = NTL::conv<double>(sz);
  }
  return mx;
}
template<class T>
double largestCoeff(const std::vector<T>& f)
{
  double mx = 0;
  for (auto& x : f) {
    auto sz = abs(x);
    if (mx < sz)
      mx = NTL::conv<double>(sz);
  }
  return mx;
}
NTL::ZZ largestCoeff(const NTL::ZZX& f);
NTL::ZZ largestCoeff(const DoubleCRT& f);

//! The L2-norm of an element (in coefficient representation)
double coeffsL2NormSquared(const zzX& f); // l2 norm^2
NTL::xdouble coeffsL2NormSquared(const NTL::ZZX& f); // l2 norm^2
NTL::xdouble coeffsL2NormSquared(const DoubleCRT& f); // l2 norm^2

inline double coeffsL2Norm(const zzX& f) // l2 norm
{ return sqrt(coeffsL2NormSquared(f)); }
inline NTL::xdouble coeffsL2Norm(const NTL::ZZX& f) // l2 norm
{ return sqrt(coeffsL2NormSquared(f)); }
inline NTL::xdouble coeffsL2Norm(const DoubleCRT& f) // l2 norm
{ return sqrt(coeffsL2NormSquared(f)); }

// Choosing between implementations 
#if defined(FFT_NATIVE) || defined(FFT_ARMA)
#define FFT_IMPL 1
#else
#define FFT_IMPL 0 // no implementation of FFT
#endif

typedef std::complex<double> cx_double;
typedef std::complex<NTL::xdouble> cx_xdouble;
typedef std::complex<long double> cx_ldouble;

#if FFT_IMPL
//! Computing the L2 norm of the canonical embedding
double embeddingL2NormSquared(const zzX& f, const PAlgebra& palg);
inline double embeddingL2Norm(const zzX& f, const PAlgebra& palg)
{ return sqrt(embeddingL2NormSquared(f,palg)); }

//! Computing the L-infinity norm of the canonical embedding
double embeddingLargestCoeff(const zzX& f, const PAlgebra& palg);

// Same as above, for ZZX
NTL::xdouble embeddingL2NormSquared(const NTL::ZZX& f, const PAlgebra& palg);
inline NTL::xdouble embeddingL2Norm(const NTL::ZZX& f, const PAlgebra& palg)
{ return sqrt(embeddingL2NormSquared(f,palg)); }
NTL::xdouble embeddingLargestCoeff(const NTL::ZZX& f, const PAlgebra& palg);

//! Computing the canonical embedding (in fft.cpp). This function
//! returns in v only the first half of the entries, the others are
//! v[phi(m)-i] = conj(v[i])
void canonicalEmbedding(std::vector<cx_double>& v,
                        const zzX& f, const PAlgebra& palg);

void canonicalEmbedding(std::vector<cx_double>& v,
                        const NTL::ZZX& f, const PAlgebra& palg);

//! Roughly the inverse of canonicalEmbedding, except for scaling and rounding issues
void embedInSlots(zzX& f, const std::vector<cx_double>& v,
                  const PAlgebra& palg,
                  double scaling=1.0, bool strictInverse=false);
// When m is a power of two, the strictInverse flag has no effect.
// Otherwise, calling embedInSlots(f,v,palg,1.0,strictInverse=true) after
// setting canonicalEmbedding(v, f, palg), is sure to recover the same f,
// but embedInSlots(f,v,palg,1.0,strictInverse=false) may fail to recover
// the same f.

#endif // FFT_IMPL

#endif // ifndef HELIB_NORMS_H
