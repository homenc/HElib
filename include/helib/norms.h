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
#include <helib/zzX.h>

namespace helib {

class DoubleCRT;

long sumOfCoeffs(const zzX& f);          // = f(1)
NTL::ZZ sumOfCoeffs(const NTL::ZZX& f);  // = f(1)
NTL::ZZ sumOfCoeffs(const DoubleCRT& f); // somewhat lame implementation

//! The L-infinity norm of an element (in coefficient representation)
template <typename T>
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
template <typename T>
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

NTL::ZZ largestCoeff(const NTL::Vec<NTL::ZZ>& f);
// somebody eliminated this...please leave it here!

NTL::ZZ largestCoeff(const DoubleCRT& f);

//! The L2-norm of an element (in coefficient representation)
double coeffsL2NormSquared(const zzX& f);             // l2 norm^2
NTL::xdouble coeffsL2NormSquared(const NTL::ZZX& f);  // l2 norm^2
NTL::xdouble coeffsL2NormSquared(const DoubleCRT& f); // l2 norm^2

inline double coeffsL2Norm(const zzX& f) // l2 norm
{
  return sqrt(coeffsL2NormSquared(f));
}
inline NTL::xdouble coeffsL2Norm(const NTL::ZZX& f) // l2 norm
{
  return sqrt(coeffsL2NormSquared(f));
}
inline NTL::xdouble coeffsL2Norm(const DoubleCRT& f) // l2 norm
{
  return sqrt(coeffsL2NormSquared(f));
}

typedef std::complex<double> cx_double;

//! Computing the L-infinity norm of the canonical embedding
//! Assumed: deg(f) < phi(m).
double embeddingLargestCoeff(const zzX& f, const PAlgebra& palg);

double embeddingLargestCoeff(const std::vector<double>& f,
                             const PAlgebra& palg);

// computes two for the price of one
void embeddingLargestCoeff_x2(double& norm1,
                              double& norm2,
                              const std::vector<double>& f1,
                              const std::vector<double>& f2,
                              const PAlgebra& palg);

NTL::xdouble embeddingLargestCoeff(const NTL::ZZX& f, const PAlgebra& palg);

//! Computes canonical embedding.
//! Requires p==-1 and m==2^k where k >=2 and f.length() < m/2.
//! Sets v[m/4-1-i] = DFT[palg.ith_rep(i)] for i in range(m/4),
//! where DFT[j] = f(W^j) for j in range(m), and W = exp(-2*pi*I/m).
// FIXME: need to to understand and document why te v array
// gets initialized in the order that it does...what else
// in the library depends on this particular order.
void CKKS_canonicalEmbedding(std::vector<cx_double>& v,
                             const zzX& f,
                             const PAlgebra& palg);

void CKKS_canonicalEmbedding(std::vector<cx_double>& v,
                             const NTL::ZZX& f,
                             const PAlgebra& palg);

void CKKS_canonicalEmbedding(std::vector<cx_double>& v,
                             const std::vector<double>& f,
                             const PAlgebra& palg);

//! Requires p==-1 and m==2^k where k >=2.
//! Computes the inverse of canonical embedding, scaled by scaling
//! and then rounded to nearest integer.
void CKKS_embedInSlots(zzX& f,
                       const std::vector<cx_double>& v,
                       const PAlgebra& palg,
                       double scaling);

} // namespace helib

#endif // ifndef HELIB_NORMS_H
