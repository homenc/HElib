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
#include <algorithm>
#include <NTL/BasicThreadPool.h>
#include "NumbTh.h"
#include "timing.h"
#include "norms.h"
#include "PAlgebra.h"
NTL_CLIENT


static void
basicCanonicalEmbedding(std::vector<cx_double>& v, 
                        const std::vector<double>& in, 
                        const PAlgebra& palg)
{
  long m = palg.getM();
  long phimBy2 = divc(palg.getPhiM(),2);

  vector<cx_double> buf(m);
  for (long i: range(in.size())) buf[i] = in[i];
  for (long i: range(in.size(), m)) buf[i] = 0;
  palg.getFFTInfo().apply(&buf[0]);

  v.resize(phimBy2); // the first half of Zm*

  // FIXME: need to document these two different strategies
  if (palg.getNSlots()==phimBy2) // order roots by the palg order
    for (long i=0; i<phimBy2; i++)
      v[phimBy2-i-1] = buf[palg.ith_rep(i)];
  else                           // order roots sequentially
    for (long i=1, idx=0; i<=m/2; i++)
      if (palg.inZmStar(i)) v[idx++] = buf[i];
}



// Computing the canonical embedding. This function returns in v only
// the first half of the entries, the others are v[phi(m)-i]=conj(v[i])
void canonicalEmbedding(std::vector<cx_double>& v,
                        const zzX& f, const PAlgebra& palg)
{
  FHE_TIMER_START;

  vector<double> x;
  convert(x, f);

  basicCanonicalEmbedding(v, x, palg);
}
   


void canonicalEmbedding(std::vector<cx_double>& v,
                        const ZZX& f, const PAlgebra& palg)
{
  FHE_TIMER_START;

  vector<double> x;
  convert(x, f.rep);

  basicCanonicalEmbedding(v, x, palg);
}

void canonicalEmbedding(std::vector<cx_double>& v,
                        const std::vector<double>& f, const PAlgebra& palg)
{
  FHE_TIMER_START;

  basicCanonicalEmbedding(v, f, palg);
}

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
  vector<cx_double> avv(m);
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


  // Compute the inverse FFT and extract the real part.

  // NOTES:
  // For a polynomial f with complex coeffs, and w a root of unity,
  // we have f(conj(w)) = conj(conj(f)(w)).  So we can compute
  // an inverse FFT of f as conj(FFT(conj(f)))/m.  Since we only
  // extract the real part, we can skip the outer conj.


  for (long i: range(m)) avv[i] = conj(avv[i]);
  palg.getFFTInfo().apply(&avv[0]);
  vector<double> av(m);

  // if strictInverse we need to scale up by m, so we just skip
  // the division by m step required by the inverse fft
  if (!strictInverse) {
    double m_inv = 1/double(m);
    for (long i: range(m)) av[i] = avv[i].real() * m_inv;
  }

  // If v was obtained by canonicalEmbedding(v,f,palg,1.0) then we have
  // the guarantee that m*av is an integral polynomial, and moreover
  // m*av mod Phi_m(x) is in m*Z[X].

  // round to an integer polynomial
  f.SetLength(m);
  for (long i: range(m)) f[i] = std::round(av[i]);

  reduceModPhimX(f, palg);

  if (strictInverse) f /= m;  // scale down by m
  normalize(f);
}


