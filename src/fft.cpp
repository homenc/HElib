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
 * @file fft.cpp - computing the canonical embedding and related norms
 **/
#include <complex>
#include <cmath>
#include <numeric> // std::accumulate
#include <NTL/BasicThreadPool.h>
#include "NumbTh.h"
#include "timing.h"
#include "norms.h"
NTL_CLIENT

constexpr double pi = 4 * std::atan(1);

#define FFT_ARMA

#ifdef FFT_NATIVE
// An extremely lame implementation of the canonical embedding

// evaluate poly(x) using Horner's rule
cx_double complexEvalPoly(const zzX& poly, const cx_double& x)
{
  if (lsize(poly)<=0) return cx_double(0.0,0.0);
  cx_double res(double(poly[0]), 0.0);
  for (long i=1; i<lsize(poly); i++) {
    res *= x;
    res += cx_double(double(poly[i]));
  }
  return res;
}

void canonicalEmbedding(std::vector<cx_double>& v, const zzX& f, long m)
{
  FHE_TIMER_START;
  vector<long> zmstar(phi_N(m)/2); // the first half of Zm*
  for (long i=1, idx=0; i<=m/2; i++)
    if (NTL::GCD(i,m)==1) zmstar[idx++] = i;

  v.resize(zmstar.size());
  NTL_EXEC_RANGE(zmstar.size(), first, last)
  for (long i=first; i < last; ++i) {
    auto rou = std::polar<double>(1.0, -(2*pi*zmstar[i])/m); // root of unity
    v[i] = complexEvalPoly(f,rou);
  }
  NTL_EXEC_RANGE_END
}
void canonicalUnembedding(zzX& f, const std::vector<cx_double>& v, long m)
{
  NTL::Error("canonicalUnembedding not implemented\n");
}
#else // ifdef FFT_NATIVE
#ifdef FFT_ARMA
#warning "canonicalEmbedding implemented via Armadillo"

#include <armadillo>
void convert(zzX& to, const arma::vec& from)
{
  to.SetLength(from.size());
  NTL_EXEC_RANGE(to.length(), first, last)
  for (long i=first; i<last; i++)
    to[i] = std::round(from[i]);
  NTL_EXEC_RANGE_END
}
void convert(arma::vec& to, const zzX& from)
{
  to.resize(from.length());
  NTL_EXEC_RANGE(from.length(), first, last)
  for (long i=first; i<last; i++)
    to[i] = from[from.length()-i-1];
  NTL_EXEC_RANGE_END
}

void canonicalEmbedding(std::vector<cx_double>& v, const zzX& f, long m)
{
  FHE_TIMER_START;
  arma::vec av; // convert to vector of doubles
  convert(av, f);
  arma::cx_vec avv = arma::fft(av,m); // compute the full FFT

  v.resize(phi_N(m)/2); // the first half of Zm*
  for (long i=1, idx=0; i<=m/2; i++)
    if (NTL::GCD(i,m)==1) v[idx++] = avv[i];
}
#else // ifdef FFT_ARMA
#error "No implementation found for canonicalEmbedding"
#endif // ifdef FFT_ARMA
#endif // ifdef FFT_NATIVE

// l_2 norm square of canonical embedding
double embeddingL2NormSquared(const zzX& f, long m)
{
  std::vector<cx_double> emb;
  canonicalEmbedding(emb, f, m);
  double acc = 0.0;
  for (auto& x : emb)
    acc += std::norm(x);
  return 2*acc;
}
