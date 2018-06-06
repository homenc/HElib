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
/**
 * @file norms.cpp - computing various norms of ring elements
 **/
#include <numeric>
#include <NTL/BasicThreadPool.h>
#include "NumbTh.h"
#include "DoubleCRT.h"
#include "norms.h"
NTL_CLIENT

long sumOfCoeffs(const zzX& f) // = f(1)
{
  long sum = 0;
  for (long i=0; i<lsize(f); i++) sum += f[i];
  return sum;
}
ZZ sumOfCoeffs(const ZZX& f) // = f(1)
{
  ZZ sum = ZZ::zero();
  for (long i=0; i<=deg(f); i++) sum += coeff(f,i);
  return sum;
}
NTL::ZZ sumOfCoeffs(const DoubleCRT& f)
{
  ZZX poly;
  f.toPoly(poly);
  return sumOfCoeffs(poly);
}


long largestCoeff(const zzX& f) // l_infty norm
{
  long mx = 0;
  for (long i=0; i<lsize(f); i++) {
    if (mx < abs(f[i]))
      mx = abs(f[i]);
  }
  return mx;
}
ZZ largestCoeff(const ZZX& f)
{
  ZZ mx = ZZ::zero();
  for (long i=0; i<=deg(f); i++) {
    if (mx < abs(coeff(f,i)))
      mx = abs(coeff(f,i));
  }
  return mx;
}
NTL::ZZ largestCoeff(const DoubleCRT& f)
{
  ZZX poly;
  f.toPoly(poly);
  return largestCoeff(poly);
}


double coeffsL2NormSquared(const zzX& f) // l_2 norm square
{
  double s = 0.0;
  for (long i=0; i<lsize(f); i++) {
    double coef = f[i];
    s += coef * coef;
  }
  return s;
}
xdouble coeffsL2NormSquared(const ZZX& f) // l_2 norm square
{
  xdouble s(0.0);
  for (long i=0; i<=deg(f); i++) {
    xdouble coef(conv<xdouble>(coeff(f,i)));
    s += coef * coef;
  }
  return s;
}
xdouble coeffsL2NormSquared(const DoubleCRT& f) // l2 norm^2
{
  ZZX poly;
  f.toPoly(poly);
  return coeffsL2NormSquared(poly);
}

// An extremely lame implementation of approximating the canonical
// embedding norm of a polynomial
#include <complex>
#include <cmath>

// evaluate poly(x) using Horner's rule
complex<double> complexEvalPoly(const zzX& poly, const complex<double>& x)
{
  if (lsize(poly)<=0) return complex<double>(0.0,0.0);
  complex<double> res(double(poly[0]), 0.0);
  for (long i=1; i<lsize(poly); i++) {
    res *= x;
    res += complex<double>(double(poly[i]));
  }
  return res;
}
complex<xdouble> complexEvalPoly(const ZZX& poly, const complex<double>& x)
{
  if (deg(poly)<0) return complex<xdouble>(xdouble(0), xdouble(0));
  complex<xdouble> res(conv<xdouble>(poly[0]), xdouble(0));
  for (long i=1; i<=deg(poly); i++) {
    res *= x;
    res += complex<xdouble>(conv<xdouble>(poly[i]), xdouble(0));
  }
  return res;
}
complex<xdouble> complexEvalPoly(const DoubleCRT& f,
                                 const complex<double>& x)
{
  ZZX tmp;
  f.toPoly(tmp);
  return complexEvalPoly(tmp, x);
}

// l_2 norm square of canonical embedding
constexpr double pi = 4 * std::atan(1);
double embeddingL2NormSquared(const zzX& f, long m)
{
  vector<double> acc(m, 0.0);
  NTL_EXEC_RANGE(m/2, first, last)// no need to compute the second half
  for (long i=1+first; i<=last; ++i)
    if (NTL::GCD(i,m)==1) {
      auto rou = std::polar<double>(1.0, (2*pi*i)/m); // root of unity
      acc[i] += std::norm(complexEvalPoly(f,rou));
    }
  NTL_EXEC_RANGE_END
  return 2*std::accumulate(acc.begin(), acc.end(), 0.0);
}
xdouble embeddingL2NormSquared(const ZZX& f, long m)
{
  vector<xdouble> acc(m, xdouble(0.0));
  NTL_EXEC_RANGE(m/2, first, last)// no need to compute the second half
  for (long i=1+first; i <=last; ++i) // no need to compute the second half
    if (NTL::GCD(i,m)==1) {
      auto rou = std::polar<double>(1.0, (2*pi*i)/m); // root of unity
      acc[i] += std::norm(complexEvalPoly(f,rou));
    }
  NTL_EXEC_RANGE_END
  return 2*std::accumulate(acc.begin(), acc.end(), xdouble(0.0));
}
xdouble embeddingL2NormSquared(const DoubleCRT& f)
{
  ZZX tmp;
  f.toPoly(tmp);
  return embeddingL2NormSquared(tmp, f.getContext().zMStar.getM());
}


void canonicalEmbedding(std::vector<dcomplex>& v, const zzX& f, long m)
{
  cout << "canonicalEmbedding: m="<<m<<", phi(m)="<<phi_N(m)
       << ", m'th ROU="<<std::polar<double>(1.0, (2*pi)/m)<<endl;
  v.resize(phi_N(m)/2);
  long idx = 0;
  for (long i=1; i <= m/2; ++i) // no need to compute the second half
    if (NTL::GCD(i,m)==1) {
      auto rou = std::polar<double>(1.0, (2*pi*i)/m); // root of unity
      v[idx++] = complexEvalPoly(f,rou);
    }
}
void canonicalEmbedding(std::vector<xcomplex>& v, const NTL::ZZX& f, long m)
{
  v.resize(phi_N(m)/2);
  long idx = 0;
  for (long i=1; i <= m/2; ++i) // no need to compute the second half
    if (NTL::GCD(i,m)==1) {
      auto rou = std::polar<double>(1.0, (2*pi*i)/m); // root of unity
      v[idx++] = complexEvalPoly(f,rou);
    }
}
void canonicalEmbedding(std::vector<xcomplex>& v, const DoubleCRT& f)
{
  ZZX tmp;
  f.toPoly(tmp);
  canonicalEmbedding(v, tmp, f.getContext().zMStar.getM());
}
