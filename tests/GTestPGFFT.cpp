#include <helib/PGFFT.h>

#include <iostream>
#include <cstdlib>
#include <ctime>

// these are just for the Fft stuff
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <limits>
#include <gtest/gtest.h>

namespace {

// RandomBnd(n) returns a random number in [0..n).
// Assumes n > 0.
// FIXME: uses brain-dead rand() function
static long RandomBnd(long n)
{
  const int BPL = std::numeric_limits<unsigned long>::digits;
  const int ROTAMT = 7;
  unsigned long x = 0;
  for (long i = 0; i < 12; i++) {
    unsigned long x1 = std::rand();
    // rotate x ROTAMT bits
    x = (x << ROTAMT) | (x >> (BPL - ROTAMT));
    x = x ^ x1;
  }

  return long(x % ((unsigned long)n));
}

static void SetSeed() { srand(time(0)); }

//================== Fft ====================

// I've modified this a bit from code I got here:
// https://www.nayuki.io/page/free-small-fft-in-multiple-languages

// Specifically, I modified it to use long doubles instead of doubles

// Here is the original copyright notice:

/*
 * Free FFT and convolution (C++)
 *
 * Copyright (c) 2017 Project Nayuki. (MIT License)
 * https://www.nayuki.io/page/free-small-fft-in-multiple-languages
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * - The above copyright notice and this permission notice shall be included in
 *   all copies or substantial portions of the Software.
 * - The Software is provided "as is", without warranty of any kind, express or
 *   implied, including but not limited to the warranties of merchantability,
 *   fitness for a particular purpose and noninfringement. In no event shall the
 *   authors or copyright holders be liable for any claim, damages or other
 *   liability, whether in an action of contract, tort or otherwise, arising
 * from, out of or in connection with the Software or the use or other dealings
 * in the Software.
 */

typedef long double ldbl;
typedef std::complex<ldbl> lcx;

namespace Fft {

/*
 * Computes the discrete Fourier transform (DFT) of the given complex vector,
 * storing the result back into the vector. The vector can have any length. This
 * is a wrapper function.
 */
void transform(std::vector<lcx>& vec);

/*
 * Computes the inverse discrete Fourier transform (IDFT) of the given complex
 * vector, storing the result back into the vector. The vector can have any
 * length. This is a wrapper function. This transform does not perform scaling,
 * so the inverse is not a true inverse.
 */
void inverseTransform(std::vector<lcx>& vec);

/*
 * Computes the discrete Fourier transform (DFT) of the given complex vector,
 * storing the result back into the vector. The vector's length must be a power
 * of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
 */
void transformRadix2(std::vector<lcx>& vec);

/*
 * Computes the discrete Fourier transform (DFT) of the given complex vector,
 * storing the result back into the vector. The vector can have any length. This
 * requires the convolution function, which in turn requires the radix-2 FFT
 * function. Uses Bluestein's chirp z-transform algorithm.
 */
void transformBluestein(std::vector<lcx>& vec);

/*
 * Computes the circular convolution of the given complex vectors. Each vector's
 * length must be the same.
 */
void convolve(const std::vector<lcx>& vecx,
              const std::vector<lcx>& vecy,
              std::vector<lcx>& vecout);

} // namespace Fft

using std::complex;
using std::size_t;
using std::vector;

// Private function prototypes
static size_t reverseBits(size_t x, int n);

void Fft::transform(vector<lcx>& vec)
{
  size_t n = vec.size();
  if (n == 0)
    return;
  else if ((n & (n - 1)) == 0) // Is power of 2
    transformRadix2(vec);
  else // More complicated algorithm for arbitrary sizes
    transformBluestein(vec);
}

void Fft::inverseTransform(vector<lcx>& vec)
{
  std::transform(vec.cbegin(),
                 vec.cend(),
                 vec.begin(),
                 static_cast<lcx (*)(const lcx&)>(std::conj));
  transform(vec);
  std::transform(vec.cbegin(),
                 vec.cend(),
                 vec.begin(),
                 static_cast<lcx (*)(const lcx&)>(std::conj));
}

void Fft::transformRadix2(vector<lcx>& vec)
{
  // Length variables
  size_t n = vec.size();
  int levels = 0; // Compute levels = floor(log2(n))
  for (size_t temp = n; temp > 1U; temp >>= 1)
    levels++;
  if (static_cast<size_t>(1U) << levels != n)
    throw std::domain_error("Length is not a power of 2");

  const ldbl pi = atan(ldbl(1)) * 4.0;

  // Trignometric table
  vector<lcx> expTable(n / 2);
  for (size_t i = 0; i < n / 2; i++) {
    // expTable[i] = std::exp(lcx(0, -2 * M_PI * i / n));
    ldbl angle = -2 * pi * i / n;
    expTable[i] = lcx(std::cos(angle), std::sin(angle));
  }

  // Bit-reversed addressing permutation
  for (size_t i = 0; i < n; i++) {
    size_t j = reverseBits(i, levels);
    if (j > i)
      std::swap(vec[i], vec[j]);
  }

  // Cooley-Tukey decimation-in-time radix-2 FFT
  for (size_t size = 2; size <= n; size *= 2) {
    size_t halfsize = size / 2;
    size_t tablestep = n / size;
    for (size_t i = 0; i < n; i += size) {
      for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
        lcx temp = vec[j + halfsize] * expTable[k];
        vec[j + halfsize] = vec[j] - temp;
        vec[j] += temp;
      }
    }
    if (size == n) // Prevent overflow in 'size *= 2'
      break;
  }
}

void Fft::transformBluestein(vector<lcx>& vec)
{
  // Find a power-of-2 convolution length m such that m >= n * 2 + 1
  size_t n = vec.size();
  size_t m = 1;
  while (m / 2 <= n) {
    if (m > SIZE_MAX / 2)
      throw std::length_error("Vector too large");
    m *= 2;
  }

  const ldbl pi = atan(ldbl(1)) * 4.0;

  // Trignometric table
  vector<lcx> expTable(n);
  for (size_t i = 0; i < n; i++) {
    unsigned long long temp = static_cast<unsigned long long>(i) * i;
    temp %= static_cast<unsigned long long>(n) * 2;
    ldbl angle = pi * temp / n;
    // Less accurate alternative if long long is unavailable: double angle =
    // M_PI * i * i / n;
    expTable[i] = lcx(std::cos(-angle), std::sin(-angle));
  }

  // Temporary vectors and preprocessing
  vector<lcx> av(m);
  for (size_t i = 0; i < n; i++)
    av[i] = vec[i] * expTable[i];
  vector<lcx> bv(m);
  bv[0] = expTable[0];
  for (size_t i = 1; i < n; i++)
    bv[i] = bv[m - i] = std::conj(expTable[i]);

  // Convolution
  vector<lcx> cv(m);
  convolve(av, bv, cv);

  // Postprocessing
  for (size_t i = 0; i < n; i++)
    vec[i] = cv[i] * expTable[i];
}

void Fft::convolve(const vector<lcx>& xvec,
                   const vector<lcx>& yvec,
                   vector<lcx>& outvec)
{

  size_t n = xvec.size();
  if (n != yvec.size() || n != outvec.size())
    throw std::domain_error("Mismatched lengths");
  vector<lcx> xv = xvec;
  vector<lcx> yv = yvec;
  transform(xv);
  transform(yv);
  for (size_t i = 0; i < n; i++)
    xv[i] *= yv[i];
  inverseTransform(xv);
  for (size_t i = 0; i < n;
       i++) // Scaling (because this FFT implementation omits it)
    outvec[i] = xv[i] / static_cast<ldbl>(n);
}

static size_t reverseBits(size_t x, int n)
{
  size_t result = 0;
  for (int i = 0; i < n; i++, x >>= 1)
    result = (result << 1) | (x & 1U);
  return result;
}

//===========================================

typedef complex<double> cmplx_t;

static void TestIt(long n)
{
  helib::PGFFT pgfft(n);

  for (long j = 0; j < 10; j++) {

    vector<cmplx_t> v(n);
    for (int i = 0; i < n; i++) {
      v[i] = RandomBnd(20) - 10;
    }

    vector<cmplx_t> v0(v);

    pgfft.apply(v.data());

    vector<lcx> vv(n);
    for (int i = 0; i < n; i++)
      vv[i] = v0[i];

    Fft::transform(vv);

    ldbl vv_norm = 0;
    for (int i = 0; i < n; i++) {
      vv_norm += std::norm(vv[i]);
    }

    vv_norm = sqrt(vv_norm);

    ldbl diff_norm = 0;
    for (int i = 0; i < n; i++) {
      lcx val = v[i];
      lcx diff = val - vv[i];
      diff_norm += std::norm(diff);
    }
    diff_norm = sqrt(diff_norm);

    if (vv_norm == 0) {
      // vv has norm = 0. Cheching if the fft is correct looking only the
      // enumerator.
      EXPECT_EQ(diff_norm, 0);
    } else {
      // Check if the fft relative error is smaller than the treshold.
      ldbl rel_err = diff_norm / vv_norm;
      EXPECT_LE(rel_err, 1e-9);
    }
  }
}

TEST(GTestPGFFT, PGFFTWorksInRange1to100Points)
{
  SetSeed();

  for (long n = 1; n <= 100; n++)
    TestIt(n);
}

TEST(GTestPGFFT, PGFFTWorksInRange256to32768PowerOfTwoPoints)
{
  SetSeed();

  for (long n = 256; n <= 32 * 1024; n *= 2)
    TestIt(n);
}

TEST(GTestPGFFT, PGFFTWorksInRange10000to20000Points)
{
  SetSeed();

  for (long i = 0; i < 100; i++) {
    long n = RandomBnd(10000) + 10000;
    TestIt(n);
  }
}
} // namespace
