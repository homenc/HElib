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
/**
 * @file norms.cpp - computing various norms of ring elements
 **/
#include <complex>
#include <cmath>
#include <algorithm>

#include <helib/NumbTh.h>
#include <helib/DoubleCRT.h>
#include <helib/norms.h>
#include <helib/PAlgebra.h>
#include <helib/fhe_stats.h>
#include <helib/range.h>

namespace helib {

#define USE_HALF_FFT (1)
// setting this to one makes use if a faster "half FFT"
// strategy when m is even.  Experimentally, setting to 1 seems best.

#define USE_QUARTER_FFT (1)
// setting this to one makes use if a faster "quarter FFT"
// strategy when m  = 0 (mod 4).
// Experimentally, setting to 1 seems best.

#define USE_TWO_QUARTERS (0)
// if m = 0 (4), then one can do two quarter FFT's of one half FFT.
// This is only used in embeddingLargestCoeff_x2.
// Experimentally, one half FFT seems a bit better, so setting to 0 is best.

long sumOfCoeffs(const zzX& f) // = f(1)
{
  long sum = 0;
  for (long i = 0; i < lsize(f); i++)
    sum += f[i];
  return sum;
}

NTL::ZZ sumOfCoeffs(const NTL::ZZX& f) // = f(1)
{
  NTL::ZZ sum = NTL::ZZ::zero();
  for (long i = 0; i <= deg(f); i++)
    sum += coeff(f, i);
  return sum;
}

NTL::ZZ sumOfCoeffs(const DoubleCRT& f)
{
  NTL::ZZX poly;
  f.toPoly(poly);
  return sumOfCoeffs(poly);
}

NTL::ZZ largestCoeff(const NTL::ZZX& f)
{
  NTL::ZZ mx = NTL::ZZ::zero();
  for (long i = 0; i <= deg(f); i++) {
    if (mx < abs(coeff(f, i)))
      mx = abs(coeff(f, i));
  }
  return mx;
}

NTL::ZZ largestCoeff(const NTL::Vec<NTL::ZZ>& f)
{
  NTL::ZZ mx = NTL::ZZ::zero();
  for (auto& x : f) {
    if (mx < abs(x))
      mx = abs(x);
  }
  return mx;
}

NTL::ZZ largestCoeff(const DoubleCRT& f)
{
  NTL::ZZX poly;
  f.toPoly(poly);
  return largestCoeff(poly);
}

double coeffsL2NormSquared(const zzX& f) // l_2 norm square
{
  double s = 0.0;
  for (long i = 0; i < lsize(f); i++) {
    double coef = f[i];
    s += coef * coef;
  }
  return s;
}

NTL::xdouble coeffsL2NormSquared(const NTL::ZZX& f) // l_2 norm square
{
  NTL::xdouble s(0.0);
  for (long i = 0; i <= deg(f); i++) {
    NTL::xdouble coef(NTL::conv<NTL::xdouble>(coeff(f, i)));
    s += coef * coef;
  }
  return s;
}

NTL::xdouble coeffsL2NormSquared(const DoubleCRT& f) // l2 norm^2
{
  NTL::ZZX poly;
  f.toPoly(poly);
  return coeffsL2NormSquared(poly);
}

// ===  Computing the L-infinity norm of the canonical embedding ===

// The built-in complex mul has a lot of extra tests related to
// infinites and NaNs that slows it down considerably
static inline cx_double MUL(cx_double a, cx_double b)
{
  double x = a.real(), y = a.imag(), u = b.real(), v = b.imag();
  return cx_double(x * u - y * v, x * v + y * u);
}

static double basic_embeddingLargestCoeff(const std::vector<double>& f,
                                          const PAlgebra& palg)
{

  long m = palg.getM();
  long sz = f.size();

  if (sz > m)
    throw LogicError("vector too big f canonicalEmbedding");

  std::vector<cx_double> buf(m);
  for (long i : range(0, sz))
    buf[i] = f[i];
  for (long i : range(sz, m))
    buf[i] = 0;
  palg.getFFTInfo().apply(&buf[0]);

  double mx = 0;

  for (long i = 1; i <= m / 2; i++) {
    if (palg.inZmStar(i)) {
      double n = std::norm(buf[i]);
      if (mx < n)
        mx = n;
    }
  }

  return sqrt(mx);
}

// Odd-Power Trick.  Suppose m is even and deg(f) < m/2,
// and we want to compute f(W^(2j+1)) =  \sum_{i < m/2} f_i W^{i(2j+1)})
// for j in range(m/2), where W is a primitive m-th root of unity.
// Observe that W^{i(2j+1)} = W^i V^{ij}, where V = W^2 is a primitive
// (m/2)-th root of unity.  So if we define g_i = f_i W^i,
// then f(W^(2j+1)) = \sum_{i < m/2} g_i V^{ij}.  So this
// reduces the original problem to a single (m/2)-point DFT.

static double half_embeddingLargestCoeff(const std::vector<double>& f,
                                         const PAlgebra& palg)
{

  long m = palg.getM();
  long sz = f.size();

  if (sz > m / 2)
    throw LogicError("vector too big f canonicalEmbedding");

  const half_FFT& hfft = palg.getHalfFFTInfo();
  const cx_double* pow = &hfft.pow[0];

  std::vector<cx_double> buf(m / 2);
  for (long i : range(0, sz))
    buf[i] = f[i] * pow[i];
  for (long i : range(sz, m / 2))
    buf[i] = 0;
  hfft.fft.apply(&buf[0]);

  double mx = 0;

  for (long i = 1; i <= m / 2; i += 2) {
    if (palg.inZmStar(i)) {
      double n = std::norm(buf[i >> 1]);
      if (mx < n)
        mx = n;
    }
  }

  return sqrt(mx);
}

// This combines the Odd-Power Trick (see above) with the standard technique
// for cutting the cost of a real DFT in half.  Provided m is a multiple of 4,
// we can reduce to a single complex (m/4)-point DFT.

static double quarter_embeddingLargestCoeff(const std::vector<double>& f,
                                            const PAlgebra& palg)
{

  long m = palg.getM();
  long sz = f.size();

  if (sz > m / 2)
    throw LogicError("vector too big f canonicalEmbedding");

  const quarter_FFT& qfft = palg.getQuarterFFTInfo();
  const cx_double* pow1 = &qfft.pow1[0];
  const cx_double* pow2 = &qfft.pow2[0];

  std::vector<cx_double> buf(m / 4);
  for (long i : range(0, sz / 2))
    buf[i] = MUL(cx_double(f[2 * i], f[2 * i + 1]), pow2[i]);
  for (long i : range(sz / 2, m / 4))
    buf[i] = 0;
  if (sz % 2)
    buf[sz / 2] = MUL(cx_double(f[sz - 1], 0), pow2[sz / 2]);
  qfft.fft.apply(&buf[0]);

  double mx = 0;

  for (long i = 1; i <= m / 2; i += 2) {
    if (palg.inZmStar(i)) {
      cx_double z1 = buf[i >> 1];
      cx_double z2 = std::conj(buf[(m / 2 - i) >> 1]);
      cx_double Xe = 0.5 * (z1 + z2);
      cx_double tmp = 0.5 * (z1 - z2);
      cx_double Xo = cx_double(tmp.imag(), -tmp.real());
      // Xo = -I*tmp

      cx_double Y = Xe + MUL(Xo, pow1[i >> 1]);
      double n = std::norm(Y);
      if (mx < n)
        mx = n;
    }
  }

  return sqrt(mx);
}

double embeddingLargestCoeff(const std::vector<double>& f, const PAlgebra& palg)
{
  HELIB_NTIMER_START(AAA_embeddingLargest);

  long m = palg.getM();
  if (USE_HALF_FFT && m % 2 == 0) {
    if (USE_QUARTER_FFT && m % 4 == 0)
      return quarter_embeddingLargestCoeff(f, palg);
    else
      return half_embeddingLargestCoeff(f, palg);
  } else
    return basic_embeddingLargestCoeff(f, palg);
}

// Computes two for the price of one!
// This is a standard technique for computing the two real DFT's
// using one complex DFT.

static void basic_embeddingLargestCoeff_x2(double& norm1,
                                           double& norm2,
                                           const std::vector<double>& f1,
                                           const std::vector<double>& f2,
                                           const PAlgebra& palg)
{

  long m = palg.getM();
  long sz1 = f1.size();
  long sz2 = f2.size();

  if (sz1 > m || sz2 > m)
    throw LogicError("vector too big in canonicalEmbedding");

  long sz_max = std::max(sz1, sz2);
  long sz_min = std::min(sz1, sz2);

  std::vector<cx_double> buf(m);
  for (long i : range(0, sz_min))
    buf[i] = cx_double(f1[i], f2[i]);
  for (long i : range(sz_min, sz1))
    buf[i] = cx_double(f1[i], 0);
  for (long i : range(sz_min, sz2))
    buf[i] = cx_double(0, f2[i]);
  for (long i : range(sz_max, m))
    buf[i] = 0;

  palg.getFFTInfo().apply(&buf[0]);

  double mx1 = 0, mx2 = 0;

  for (long i = 1; i <= m / 2; i++) {
    if (palg.inZmStar(i)) {
      cx_double x1 = 0.5 * (buf[i] + std::conj(buf[m - i]));
      cx_double x2 = 0.5 * (buf[i] - std::conj(buf[m - i]));
      // should multiply by -sqrt(-1) to get the right value
      // of x2, but we only need its norm, so we don't bother

      double n1 = std::norm(x1);
      double n2 = std::norm(x2);

      if (mx1 < n1)
        mx1 = n1;
      if (mx2 < n2)
        mx2 = n2;
    }
  }

  norm1 = sqrt(mx1);
  norm2 = sqrt(mx2);

#if 0
// debugging code
  double xx[2], yy[2];
  xx[0] = basic_embeddingLargestCoeff(f1, palg);
  xx[1] = basic_embeddingLargestCoeff(f2, palg);
  yy[0] = norm1;
  yy[1] = norm2;

  for (long i: range(2)) {
    double relerr = 0;
    if (xx[i] == 0) {
      if (yy[i] > 0) relerr = 1;
    }
    else {
      relerr = abs(xx[i]-yy[i])/xx[i];
    }

    HELIB_STATS_UPDATE("embeddingLargestCoeff_x2", relerr);
  }
#endif
}

static void half_embeddingLargestCoeff_x2(double& norm1,
                                          double& norm2,
                                          const std::vector<double>& f1,
                                          const std::vector<double>& f2,
                                          const PAlgebra& palg)
{

  long m = palg.getM();
  long sz1 = f1.size();
  long sz2 = f2.size();

  if (sz1 > m / 2 || sz2 > m / 2)
    throw LogicError("vector too big in canonicalEmbedding");

  long sz_max = std::max(sz1, sz2);
  long sz_min = std::min(sz1, sz2);

  const half_FFT& hfft = palg.getHalfFFTInfo();
  const cx_double* pow = &hfft.pow[0];

  // Odd-Power Trick.  See above.

  std::vector<cx_double> buf(m / 2);
  for (long i : range(0, sz_min))
    buf[i] = MUL(cx_double(f1[i], f2[i]), pow[i]);
  for (long i : range(sz_min, sz1))
    buf[i] = MUL(cx_double(f1[i], 0), pow[i]);
  for (long i : range(sz_min, sz2))
    buf[i] = MUL(cx_double(0, f2[i]), pow[i]);
  for (long i : range(sz_max, m / 2))
    buf[i] = 0;

  hfft.fft.apply(&buf[0]);

  double mx1 = 0, mx2 = 0;

  for (long i = 1; i <= m / 2; i += 2) {
    if (palg.inZmStar(i)) {
      cx_double x1 = 0.5 * (buf[i >> 1] + std::conj(buf[(m - i) >> 1]));
      cx_double x2 = 0.5 * (buf[i >> 1] - std::conj(buf[(m - i) >> 1]));
      // should multiply by -I to get the right value
      // of x2, but we only need its norm, so we don't bother

      double n1 = std::norm(x1);
      double n2 = std::norm(x2);

      if (mx1 < n1)
        mx1 = n1;
      if (mx2 < n2)
        mx2 = n2;
    }
  }

  norm1 = sqrt(mx1);
  norm2 = sqrt(mx2);

#if 0
// debugging code
  double xx[2], yy[2];
  xx[0] = basic_embeddingLargestCoeff(f1, palg);
  xx[1] = basic_embeddingLargestCoeff(f2, palg);
  yy[0] = norm1;
  yy[1] = norm2;

  for (long i: range(2)) {
    double relerr = 0;
    if (xx[i] == 0) {
      if (yy[i] > 0) relerr = 1;
    }
    else {
      relerr = abs(xx[i]-yy[i])/xx[i];
    }

    HELIB_STATS_UPDATE("embeddingLargestCoeff_x2", relerr);
  }
#endif
}

void embeddingLargestCoeff_x2(double& norm1,
                              double& norm2,
                              const std::vector<double>& f1,
                              const std::vector<double>& f2,
                              const PAlgebra& palg)
{
  HELIB_NTIMER_START(AAA_embeddingLargest_x2);

  long m = palg.getM();
  if (USE_HALF_FFT && m % 2 == 0) {
    if (USE_QUARTER_FFT && USE_TWO_QUARTERS && m % 4 == 0) {
      norm1 = quarter_embeddingLargestCoeff(f1, palg);
      norm2 = quarter_embeddingLargestCoeff(f2, palg);
    } else
      half_embeddingLargestCoeff_x2(norm1, norm2, f1, f2, palg);
  } else
    basic_embeddingLargestCoeff_x2(norm1, norm2, f1, f2, palg);
}

double embeddingLargestCoeff(const zzX& f, const PAlgebra& palg)
{
  std::vector<double> ff;
  convert(ff, f);
  return embeddingLargestCoeff(ff, palg);
}

static NTL::xdouble scale(std::vector<double>& ff, const NTL::ZZX& f)
{
  // max allowed bits to avoid double overflow in computations.
  // Note that the implementation of the L-infty norms naively
  // computes the complex norm and then takes a sqrt of that.
  // Given that an IEEE double can handle numbers only less than 1024
  // bits in length, this limits the bit lengths that are input
  // to the FFT to 512 - log_2(m), where m = size of the transform
  // (in the worst case, the FFT will add log_2(m) bits).
  // Using the extremely pessimistic bound of log_2(m) <= 64,
  // we get MAX_BITS = 448.  We set it to 400 just to be even
  // more on the safe side.

  const long MAX_BITS = 400;
  NTL::xdouble retval;
  long size = NTL::MaxBits(f);

  if (size > MAX_BITS) {
    long m = f.rep.length();
    ff.resize(m);

    NTL::ZZ scaled;
    for (long i : range(m)) {
      // NOTE: too bad NTL does not provide a scale and convert routine
      RightShift(scaled, f.rep[i], size - MAX_BITS);
      ff[i] = NTL::conv<double>(scaled);
    }

    retval = NTL::power2_xdouble(size - MAX_BITS);
  } else {
    convert(ff, f.rep);
    retval = 1.0;
  }

  return retval;
}

NTL::xdouble embeddingLargestCoeff(const NTL::ZZX& f, const PAlgebra& palg)
{
  std::vector<double> ff; // to hold a scaled-down version of ff;
  NTL::xdouble factor = scale(ff, f);
  return embeddingLargestCoeff(ff, palg) * factor;
}

// === Computing the canonical embedding and inverse ===

// Odd-Power Trick. See above.

// NOTE: since m is a multiple of 4, we could use the "quarter FFT"
// logic.  We can revisit this if this code ever becomes a bottleneck,
// but this does not seem to be a significant issue at the moment.

void CKKS_canonicalEmbedding(std::vector<cx_double>& v,
                             const std::vector<double>& in,
                             const PAlgebra& palg)
{
  HELIB_TIMER_START;

  long sz = in.size();
  long m = palg.getM();

  if (!(palg.getP() == -1 && palg.getPow2() >= 2 && sz <= m / 2))
    throw LogicError("bad args to CKKS_canonicalEmbedding");

  const half_FFT& hfft = palg.getHalfFFTInfo();
  const cx_double* pow = &hfft.pow[0];

  std::vector<cx_double> buf(m / 2);
  for (long i : range(0, sz))
    buf[i] = in[i] * pow[i];
  for (long i : range(sz, m / 2))
    buf[i] = 0;
  hfft.fft.apply(&buf[0]);

  v.resize(m / 4);
  for (long i : range(m / 4))
    v[m / 4 - i - 1] = buf[palg.ith_rep(i) >> 1];
}

void CKKS_canonicalEmbedding(std::vector<cx_double>& v,
                             const zzX& f,
                             const PAlgebra& palg)
{
  HELIB_TIMER_START;

  std::vector<double> x;
  convert(x, f);

  CKKS_canonicalEmbedding(v, x, palg);
}

void CKKS_canonicalEmbedding(std::vector<cx_double>& v,
                             const NTL::ZZX& f,
                             const PAlgebra& palg)
{
  HELIB_TIMER_START;

  std::vector<double> x;
  convert(x, f.rep);

  CKKS_canonicalEmbedding(v, x, palg);
}

// The forward transform (as computed by CKKS_canonicalEmbedding)
// is computing the linear transformation:
//    DFT * D
// where D is a diagonal matrix with W^i on the diagonal in the i-th row,
// and DFT is an (m/2)-point DFT matrix for an evaluation at powers of V=W^2.
// D is applied first, and DFT is applied second.
//
// As described above in the odd-power trick, the j-th entry is
// f(W^{2j+1}).
//
// Here, we want to compute the inverse transformation:
//   D^{-1} * DFT^{-1}
// Now, DFT^{-1} = 1/(m/2) times a DFT matrix for conj(V).
// So this computation effectively maps a polynomial g of degree < m/2
// to a vector whose j-th component is
//   g(conj(V)^j) * conj(V)^j = conj(conj(g)(V^j)) * conj(V^j)
//                            = conj(conj(g)(V^j) * V^j)
// Moreover, this value is a real number, so we can drop the
// outer conj, and just compute conj(g)(V^j) * V^j.
//
// Note also that the forward transform takes a real (m/2)-vector
// and expands it to a complex (m/2)-vector, then applies DFT*D,
// and then contracts this to a complex (m/4)-vector by keeping
// only one out of each conjugate pair.  The inverse transform
// reverses the contracting step by reinserting the missing conjugate,
// then applies D^{-1} * DFT^{-1}, and then reverses the expanding
// step by dropping the complex part.

void CKKS_embedInSlots(zzX& f,
                       const std::vector<cx_double>& v,
                       const PAlgebra& palg,
                       double scaling)

{
  HELIB_TIMER_START;

  long v_sz = v.size();
  long m = palg.getM();

  if (!(palg.getP() == -1 && palg.getPow2() >= 2))
    throw LogicError("bad args to CKKS_canonicalEmbedding");

  std::vector<cx_double> buf(m / 2, cx_double(0));
  for (long i : range(m / 4)) {
    long j = palg.ith_rep(i);
    long ii = m / 4 - i - 1;
    if (ii < v_sz) {
      buf[j >> 1] = std::conj(v[ii]);
      buf[(m - j) >> 1] = v[ii];
    }
  }

  const half_FFT& hfft = palg.getHalfFFTInfo();
  const cx_double* pow = &hfft.pow[0];

  scaling /= (m / 2);
  // This is becuase DFT^{-1} = 1/(m/2) times a DFT matrix for conj(V)

  hfft.fft.apply(&buf[0]);
  f.SetLength(m / 2);
  for (long i : range(m / 2))
    f[i] = std::round(MUL(buf[i], pow[i]).real() * scaling);

  normalize(f);
}

// === obsolete versions of canonical embedding and inverse ===

// These are less efficient, and seem to have some logic errors.
// At the very least, it is not clear what, exactly, these functions
// are meant to compute.

#if 0

static void
dft_to_canonical(std::vector<cx_double>& v,
                 const std::vector<cx_double>& buf,
                 const PAlgebra& palg)
{
  long m = palg.getM();
  long phimBy2 = divc(palg.getPhiM(),2);
  v.resize(phimBy2); // the first half of Zm*

  // FIXME: need to document these two different strategies
  if (palg.getNSlots()==phimBy2) // order roots by the palg order
    for (long i=0; i<phimBy2; i++)
      v[phimBy2-i-1] = buf[palg.ith_rep(i)];
  else                           // order roots sequentially
    for (long i=1, idx=0; i<=m/2; i++)
      if (palg.inZmStar(i)) v[idx++] = buf[i];
}


void
canonicalEmbedding(std::vector<cx_double>& v,
		   const std::vector<double>& in,
		   const PAlgebra& palg)
{
  HELIB_TIMER_START;
  long m = palg.getM();

  if (long(in.size()) > m)
    throw LogicError("std::vector too big in canonicalEmbedding");

  vector<cx_double> buf(m);
  for (long i: range(in.size())) buf[i] = in[i];
  for (long i: range(in.size(), m)) buf[i] = 0;
  palg.getFFTInfo().apply(&buf[0]);

  dft_to_canonical(v, buf, palg);
}



// Computing the canonical embedding. This function returns in v only
// the first half of the entries, the others are v[phi(m)-i]=conj(v[i])
void canonicalEmbedding(std::vector<cx_double>& v,
                        const zzX& f, const PAlgebra& palg)
{
  HELIB_TIMER_START;

  vector<double> x;
  convert(x, f);

  canonicalEmbedding(v, x, palg);
}



void canonicalEmbedding(std::vector<cx_double>& v,
                        const ZZX& f, const PAlgebra& palg)
{
  HELIB_TIMER_START;

  vector<double> x;
  convert(x, f.rep);

  canonicalEmbedding(v, x, palg);
}




// Roughly the inverse of canonicalEmbedding, except for scaling and
// rounding issues. Calling embedInSlots(f,v,palg,1.0,strictInverse=true)
// after setting canonicalEmbedding(v, f, palg), is sure to recover the
// same f, but embedInSlots(f,v,palg,1.0,strictInverse=false) may return
// a different "nearby" f.
void embedInSlots(zzX& f, const std::vector<cx_double>& v,
                  const PAlgebra& palg, double scaling, bool strictInverse)
{
  HELIB_TIMER_START;
  long m = palg.getM();
  long phimBy2 = divc(palg.getPhiM(),2);
  vector<cx_double> avv(m);
  for (auto& x: avv) x = 0.0;

  if (palg.getNSlots()==phimBy2) { // roots ordered by the palg order
    for (long i=0; i<palg.getNSlots(); i++) {
      long j = palg.ith_rep(i);
      long ii = palg.getNSlots()-i-1;
      if (ii < lsize(v)) {
        avv[j] = scaling*v[ii];
        avv[m-j] = std::conj(avv[j]);
      }
    }
  }
  else {                           // roots ordered sequentially
    for (long i=1, idx=0; i<=m/2 && idx<lsize(v); i++) {
      if (palg.inZmStar(i)) {
        avv[i] = scaling*v[idx++];
        avv[m-i] = std::conj(avv[i]);
      }
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

#if 0
  for (long i: range(m)) {
    if (abs(f[i]-av[i]) > 0.25) {
       cerr << "*** too much rounding: " << abs(f[i]-av[i]) << "\n";
    }
  }
#endif


  reduceModPhimX(f, palg);

  if (strictInverse) f /= m;  // scale down by m
  normalize(f);
}

#endif

} // namespace helib
