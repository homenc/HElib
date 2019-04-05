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
#include "PAlgebra.h"
NTL_CLIENT

const double pi = 4 * std::atan(1);

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
    to[i] = from[i];
  NTL_EXEC_RANGE_END
}

void convert(arma::vec& to, const ZZX& from)
{
  to.resize(from.rep.length());
  NTL_EXEC_RANGE(from.rep.length(), first, last)
  for (long i=first; i<last; i++) {
    double x = conv<double>(from[i]);
    to[i] = x;
  }
  NTL_EXEC_RANGE_END
}

#if 0

//======================

namespace arma {

  template<> 
  struct is_supported_elem_type<RR> {
    enum {value = 1};
  }; 

  template<> 
  struct is_supported_elem_type<cx_RR> {
    enum {value = 1};
  }; 

  template<> 
  struct is_real<RR> {
    enum {value = 1};
  }; 

}


void convert(zzX& to, const arma::Col<RR>& from)
{
  to.SetLength(from.size());
  for (long i: range(to.length()))
    to[i] = conv<long>(RoundToZZ(from[i]));
}

void convert(arma::Col<RR>& to, const zzX& from)
{
  to.resize(from.length());
  for (long i: range(from.length()))
    to[i] = from[i];
}

void convert(arma::Col<RR>& to, const ZZX& from)
{
  to.resize(from.rep.length());
  for (long i: range(from.rep.length())) {
    RR x = conv<RR>(from[i]);
    to[i] = x;
  }
}

#endif

// Computing the canonical embedding. This function returns in v only
// the first half of the entries, the others are v[phi(m)-i]=conj(v[i])
void canonicalEmbedding(std::vector<cx_double>& v,
                        const zzX& f, const PAlgebra& palg)
{
  FHE_TIMER_START;
  long m = palg.getM();
  long phimBy2 = divc(palg.getPhiM(),2);
  arma::vec av; // convert to vector of doubles
  convert(av, f);
  arma::cx_vec avv = arma::fft(av,m); // compute the full FFT

  v.resize(phimBy2); // the first half of Zm*

  if (palg.getNSlots()==phimBy2) // order roots by the palg order
    for (long i=0; i<phimBy2; i++)
      v[phimBy2-i-1] = avv[palg.ith_rep(i)];
  else                           // order roots sequentially
    for (long i=1, idx=0; i<=m/2; i++)
      if (palg.inZmStar(i)) v[idx++] = avv[i];
}

void canonicalEmbedding(std::vector<cx_double>& v,
                        const ZZX& f, const PAlgebra& palg)
{
  FHE_TIMER_START;
  long m = palg.getM();
  long phimBy2 = divc(palg.getPhiM(),2);
  arma::vec av; // convert to vector of doubles
  convert(av, f);
  arma::cx_vec avv = arma::fft(av,m); // compute the full FFT

  v.resize(phimBy2); // the first half of Zm*

  if (palg.getNSlots()==phimBy2) // order roots by the palg order
    for (long i=0; i<phimBy2; i++)
      v[phimBy2-i-1] = avv[palg.ith_rep(i)];
  else                           // order roots sequentially
    for (long i=1, idx=0; i<=m/2; i++)
      if (palg.inZmStar(i)) v[idx++] = avv[i];
}


#if 0
void canonicalEmbedding(std::vector<cx_RR>& v,
                        const ZZX& f, const PAlgebra& palg)
{
  FHE_TIMER_START;
  long m = palg.getM();
  long phimBy2 = divc(palg.getPhiM(),2);
  arma::Col<RR> av; // convert to vector of doubles
  convert(av, f);
  arma::Col<cx_RR> avv = arma::fft(av,m); // compute the full FFT

  v.resize(phimBy2); // the first half of Zm*

  if (palg.getNSlots()==phimBy2) // order roots by the palg order
    for (long i=0; i<phimBy2; i++)
      v[phimBy2-i-1] = avv[palg.ith_rep(i)];
  else                           // order roots sequentially
    for (long i=1, idx=0; i<=m/2; i++)
      if (palg.inZmStar(i)) v[idx++] = avv[i];
}
#endif

// Roughly the inverse of canonicalEmbedding, except for scaling and
// rounding issues. Calling embedInSlots(f,v,palg,1.0,strictInverse=true)
// after setting canonicalEmbedding(v, f, palg), is sure to recover the
// same f, but embedInSlots(f,v,palg,1.0,strictInverse=false) may return
// a different "nearby" f.
void embedInSlots(zzX& f, const std::vector<cx_double>& v,
                  const PAlgebra& palg, double scaling, bool strictInverse)
{
  FHE_TIMER_START;
  long m = palg.getM();
  long phimBy2 = divc(palg.getPhiM(),2);
  arma::cx_vec avv(m);
  for (auto& x: avv) x = 0.0;

  if (palg.getNSlots()==phimBy2) // roots ordered by the palg order
    for (long i=0; i<palg.getNSlots(); i++) {
      long j = palg.ith_rep(i);
      long ii = palg.getNSlots()-i-1;
      if (ii < lsize(v)) {
        avv[j] = scaling*v[ii];
        avv[m-j] = std::conj(avv[j]);
      }
    }
  else                           // roots ordered sequentially
    for (long i=1, idx=0; i<=m/2 && idx<lsize(v); i++) {
      if (palg.inZmStar(i)) {
        avv[i] = scaling*v[idx++];
        avv[m-i] = std::conj(avv[i]);
      }
    }
  arma::vec av = arma::real(arma::ifft(avv,m)); // compute the inverse FFT

  // If v was obtained by canonicalEmbedding(v,f,palg,1.0) then we have
  // the guarantee that m*av is an integral polynomial, and moreover
  // m*av mod Phi_m(x) is in m*Z[X].
  if (strictInverse) av *= m; // scale up by m
  convert(f, av);    // round to an integer polynomial
  reduceModPhimX(f, palg);
  if (strictInverse) f /= m;  // scale down by m
  normalize(f);
}
#else
#ifdef FFT_NATIVE
#warning "canonicalEmbedding implemented via slow DFT, expect very slow key-generation"
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

void canonicalEmbedding(std::vector<cx_double>& v, const zzX& f, const PAlgebra& palg)
{
  FHE_TIMER_START;
  long m = palg.getM();
  long phimBy2 = divc(palg.getPhiM(),2);
  vector<long> zmstar(phimBy2); // the first half of Zm*

  if (palg.getNSlots()==phimBy2) // order roots by the palg order
    for (long i=0; i<phimBy2; i++)
      zmstar[phimBy2-i-1] = palg.ith_rep(i);
  else                           // order roots sequentially
    for (long i=1, idx=0; i<=m/2; i++)
      if (palg.inZmStar(i)) zmstar[idx++] = i;

  v.resize(phimBy2);
  NTL_EXEC_RANGE(phimBy2, first, last)
  for (long i=first; i < last; ++i) {
    auto rou = std::polar<double>(1.0, -(2*pi*zmstar[i])/m); // root of unity
    v[i] = complexEvalPoly(f,rou);
  }
  NTL_EXEC_RANGE_END
  FHE_TIMER_STOP;
}

// evaluate poly(x) using Horner's rule
// FIXME: this is actually evaluating the reverse polynomial
cx_double complexEvalPoly(const Vec<double>& poly, const cx_double& x)
{
  if (poly.length()<=0) return cx_double(0.0,0.0);
  cx_double res(poly[0], 0.0);
  for (long i: range(1, poly.length())) {
    res *= x;
    res += cx_double(poly[i]);
  }
  return res;
}

void canonicalEmbedding(std::vector<cx_double>& v, const ZZX& f, const PAlgebra& palg)
{
  FHE_TIMER_START;
  long m = palg.getM();
  long phimBy2 = divc(palg.getPhiM(),2);
  vector<long> zmstar(phimBy2); // the first half of Zm*

  Vec<double> ff;
  conv(ff, f.rep);

  if (palg.getNSlots()==phimBy2) // order roots by the palg order
    for (long i=0; i<phimBy2; i++)
      zmstar[phimBy2-i-1] = palg.ith_rep(i);
  else                           // order roots sequentially
    for (long i=1, idx=0; i<=m/2; i++)
      if (palg.inZmStar(i)) zmstar[idx++] = i;

  v.resize(phimBy2);
  NTL_EXEC_RANGE(phimBy2, first, last)
  for (long i=first; i < last; ++i) {
    auto rou = std::polar<double>(1.0, -(2*pi*zmstar[i])/m); // root of unity
    v[i] = complexEvalPoly(ff,rou);
  }
  NTL_EXEC_RANGE_END
  FHE_TIMER_STOP;
}

void embedInSlots(zzX& f, const std::vector<cx_double>& v,
                  const PAlgebra& palg, double scaling, bool strictInverse)
{
  throw helib::LogicError("embedInSlots not implemented with FFT_NATIVE");
}
#endif // ifdef FFT_NATIVE
#endif // ifdef FFT_ARMA
