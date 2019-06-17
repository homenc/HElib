#include <NTL/ZZ.h>

#include <iostream>


// these are just for the Fft stuff
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>






#include <vector>
#include <complex>
#include <cassert>




class PGFFT {
public:
   explicit
   PGFFT(long n_);

   void apply(std::vector<std::complex<double>>& v) const;


private:
   long n;
   long k;

   long strategy;

   // holds all of the twiddle factors
   std::vector<std::vector<std::complex<double>>> tab;

   // additional data structures needed for Bluestein
   std::vector<std::complex<double>> powers;
   std::vector<std::complex<double>> Rb;

   // additonal data structures needed for 2^k-point FFT
   std::vector<long> rev, rev1;


};



// set to 0 to disable the truncated Bluestein
#define PGFFT_USE_TRUNCATED_BLUE (1)

using std::vector;
using std::complex;

typedef complex<double> cmplx_t;

typedef long double ldbl;


#if (defined(__GNUC__) && (__GNUC__ >= 4))

// on relative modern versions of gcc, we can
// decalare "restricted" pointers in C++

#define PGFFT_RESTRICT __restrict

#else

#define PGFFT_RESTRICT

#endif



/***************************************************************

The code here is translated from NTL's small FFT prime, which in turn
is translated from some older software by David Harvey.

For completeness, here is David Harvey's original copyright notice:

==============================================================================

fft62: a library for number-theoretic transforms

Copyright (C) 2013, David Harvey

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


==============================================================================


TRUNCATED FFT

This code is derived from code originally developed
by David Harvey.  I include his original documentation,
annotated appropriately to highlight differences in
the implemebtation (see NOTEs).

The DFT is defined as follows.

Let the input sequence be a_0, ..., a_{N-1}.

Let w = standard primitive N-th root of 1, i.e. w = g^(2^FFT62_MAX_LGN / N),
where g = some fixed element of Z/pZ of order 2^FFT62_MAX_LGN.

Let Z = an element of (Z/pZ)^* (twisting parameter).

Then the output sequence is
  b_j = \sum_{0 <= i < N} Z^i a_i w^(ij'), for 0 <= j < N,
where j' is the length-lgN bit-reversal of j.

Some of the FFT routines can operate on truncated sequences of certain
"admissible" sizes. A size parameter n is admissible if 1 <= n <= N, and n is
divisible by a certain power of 2. The precise power depends on the recursive
array decomposition of the FFT. The smallest admissible n' >= n can be
obtained via fft62_next_size().

NOTE: the twising parameter is not implemented.
NOTE: the next admissible size function is called FFTRoundUp,


Truncated FFT interface is as follows:

xn and yn must be admissible sizes for N.

Input in xp[] is a_0, a_1, ..., a_{xn-1}. Assumes a_i = 0 for xn <= i < N.

Output in yp[] is b_0, ..., b_{yn-1}, i.e. only first yn outputs are computed.

Twisting parameter Z is described by z and lgH. If z == 0, then Z = basic
2^lgH-th root of 1, and must have lgH >= lgN + 1. If z != 0, then Z = z
(and lgH is ignored).

The buffers {xp,xn} and {yp,yn} may overlap, but only if xp == yp.

Inputs are in [0, 2p), outputs are in [0, 2p).

threads = number of OpenMP threads to use.



Inverse truncated FFT interface is as follows.

xn and yn must be admissible sizes for N, with yn <= xn.

Input in xp[] is b_0, b_1, ..., b_{yn-1}, N*a_{yn}, ..., N*a_{xn-1}.

Assumes a_i = 0 for xn <= i < N.

Output in yp[] is N*a_0, ..., N*a_{yn-1}.

Twisting parameter Z is described by z and lgH. If z == 0, then Z = basic
2^lgH-th root of 1, and must have lgH >= lgN + 1. If z != 0, then Z = z^(-1)
(and lgH is ignored).

The buffers {xp,xn} and {yp,yn} may overlap, but only if xp == yp.

Inputs are in [0, 4p), outputs are in [0, 4p).

threads = number of OpenMP threads to use.

(note: no function actually implements this interface in full generality!
This is because it is tricky (and not that useful) to implement the twisting
parameter when xn != yn.)

NOTE: threads and twisting parameter are not used here.
NOTE: the code has been re-written and simplified so that
  everything is done in place, so xp == yp.


***************************************************************/



#define PGFFT_FFT_RDUP (4)
// Currently, this should be at least 2 to support
// loop unrolling in the FFT implementation


static inline long 
FFTRoundUp(long xn, long k)
{
   long n = 1L << k;
   if (xn <= 0) return n;
   // default truncation value of 0 gets converted to n

   xn = ((xn+((1L << PGFFT_FFT_RDUP)-1)) >> PGFFT_FFT_RDUP) << PGFFT_FFT_RDUP;

   if (k >= 10) {
      if (xn > n - (n >> 4)) xn = n;
   }
   else {
      if (xn > n - (n >> 3)) xn = n;
   }
   // truncation just a bit below n does not really help
   // at all, and can sometimes slow things down slightly, so round up
   // to n.  This also takes care of cases where xn > n.
   // Actually, for smallish n, we should round up sooner,
   // at n-n/8, and for larger n, we should round up later,
   // at n-m/16.  At least, experimentally, this is what I see.

   return xn;
}




#define fwd_butterfly(xx0, xx1, w)  \
do \
{ \
   cmplx_t x0_ = xx0; \
   cmplx_t x1_ = xx1; \
   cmplx_t t_  = x0_ -  x1_; \
   xx0 = x0_ + x1_; \
   xx1 = t_ * w; \
}  \
while (0)



#define fwd_butterfly0(xx0, xx1) \
do   \
{  \
   cmplx_t x0_ = xx0;  \
   cmplx_t x1_ = xx1;  \
   xx0 = x0_ + x1_; \
   xx1 = x0_ - x1_; \
}  \
while (0)




// requires size divisible by 8
static void
new_fft_layer(cmplx_t* xp, long blocks, long size,
              const cmplx_t* PGFFT_RESTRICT wtab)
{
  size /= 2;

  do
    {
      cmplx_t* PGFFT_RESTRICT xp0 = xp;
      cmplx_t* PGFFT_RESTRICT xp1 = xp + size;

      // first 4 butterflies
      fwd_butterfly0(xp0[0+0], xp1[0+0]);
      fwd_butterfly(xp0[0+1], xp1[0+1], wtab[0+1]);
      fwd_butterfly(xp0[0+2], xp1[0+2], wtab[0+2]);
      fwd_butterfly(xp0[0+3], xp1[0+3], wtab[0+3]);

      // 4-way unroll
      for (long j = 4; j < size; j += 4) {
        fwd_butterfly(xp0[j+0], xp1[j+0], wtab[j+0]);
        fwd_butterfly(xp0[j+1], xp1[j+1], wtab[j+1]);
        fwd_butterfly(xp0[j+2], xp1[j+2], wtab[j+2]);
        fwd_butterfly(xp0[j+3], xp1[j+3], wtab[j+3]);
      }

      xp += 2 * size;
    }
  while (--blocks != 0);
}



static void
new_fft_last_two_layers(cmplx_t* xp, long blocks, const cmplx_t* wtab)
{
  // 4th root of unity
  cmplx_t w = wtab[1];

  do
    {
      cmplx_t u0 = xp[0];
      cmplx_t u1 = xp[1];
      cmplx_t u2 = xp[2];
      cmplx_t u3 = xp[3];

      cmplx_t v0 = u0 + u2;
      cmplx_t v2 = u0 - u2;
      cmplx_t v1 = u1 + u3; 
      cmplx_t t  = u1 - u3; 
      cmplx_t v3 = t * w;

      xp[0] = v0 + v1;
      xp[1] = v0 - v1;
      xp[2] = v2 + v3;
      xp[3] = v2 - v3; 

      xp += 4;
    }
  while (--blocks != 0);
}


static void 
new_fft_base(cmplx_t* xp, long lgN, const vector<vector<cmplx_t>>& tab)
{
  if (lgN == 0) return;

  if (lgN == 1)
    {
      cmplx_t x0 = xp[0];
      cmplx_t x1 = xp[1];
      xp[0] = x0 + x1;
      xp[1] = x0 - x1;
      return;
    }


  long N = 1L << lgN;

  for (long j = lgN, size = N, blocks = 1;
       j > 2; j--, blocks <<= 1, size >>= 1)
    new_fft_layer(xp, blocks, size, &tab[j][0]);

  new_fft_last_two_layers(xp, N/4, &tab[2][0]);
}


// Implements the truncated FFT interface, described above.
// All computations done in place, and xp should point to 
// an array of size N, all of which may be overwitten
// during the computation.

#define PGFFT_NEW_FFT_THRESH (10)

static
void new_fft_short(cmplx_t* xp, long yn, long xn, long lgN, 
                   const vector<vector<cmplx_t>>& tab)
{
  long N = 1L << lgN;

  if (yn == N)
    {
      if (xn == N && lgN <= PGFFT_NEW_FFT_THRESH)
	{
	  // no truncation
	  new_fft_base(xp, lgN, tab);
	  return;
	}
    }

  // divide-and-conquer algorithm

  long half = N >> 1;

  if (yn <= half)
    {
      if (xn <= half)
	{
	  new_fft_short(xp, yn, xn, lgN - 1, tab);
	}
      else
	{
	  xn -= half;

	  // (X, Y) -> X + Y
	  for (long j = 0; j < xn; j++)
	    xp[j] = xp[j] + xp[j + half];

	  new_fft_short(xp, yn, half, lgN - 1, tab);
	}
    }
  else
    {
      yn -= half;
      
      cmplx_t* PGFFT_RESTRICT xp0 = xp;
      cmplx_t* PGFFT_RESTRICT xp1 = xp + half;
      const cmplx_t* PGFFT_RESTRICT wtab = &tab[lgN][0];

      if (xn <= half)
	{
	  // X -> (X, w*X)
	  for (long j = 0; j < xn; j++)
	    xp1[j] = xp0[j] * wtab[j];

	  new_fft_short(xp0, half, xn, lgN - 1, tab);
	  new_fft_short(xp1, yn, xn, lgN - 1, tab);
	}
      else
	{
	  xn -= half;

	  // (X, Y) -> (X + Y, w*(X - Y))
          // DIRT: assumes xn is a multiple of 4
          fwd_butterfly0(xp0[0], xp1[0]);
          fwd_butterfly(xp0[1], xp1[1], wtab[1]);
          fwd_butterfly(xp0[2], xp1[2], wtab[2]);
          fwd_butterfly(xp0[3], xp1[3], wtab[3]);
	  for (long j = 4; j < xn; j+=4) {
            fwd_butterfly(xp0[j+0], xp1[j+0], wtab[j+0]);
            fwd_butterfly(xp0[j+1], xp1[j+1], wtab[j+1]);
            fwd_butterfly(xp0[j+2], xp1[j+2], wtab[j+2]);
            fwd_butterfly(xp0[j+3], xp1[j+3], wtab[j+3]);
          }

	  // X -> (X, w*X)
	  for (long j = xn; j < half; j++)
	    xp1[j] = xp0[j] * wtab[j];

	  new_fft_short(xp0, half, half, lgN - 1, tab);
	  new_fft_short(xp1, yn, half, lgN - 1, tab);
	}
    }
}

static void new_fft(cmplx_t* xp, long lgN, const vector<vector<cmplx_t>>& tab)
{
   long N = 1L << lgN;
   new_fft_short(xp, N, N, lgN, tab);
}





#define inv_butterfly0(xx0, xx1)  \
do   \
{  \
   cmplx_t x0_ = xx0;  \
   cmplx_t x1_ = xx1;  \
   xx0 = x0_ + x1_;  \
   xx1 = x0_ - x1_;  \
} while (0)


#define inv_butterfly_neg(xx0, xx1, w)  \
do  \
{  \
   cmplx_t x0_ = xx0;  \
   cmplx_t x1_ = xx1;  \
   cmplx_t t_ = x1_ * w;  \
   xx0 = x0_ - t_; /* NEG */   \
   xx1 = x0_ + t_; /* NEG */   \
} while (0)


// requires size divisible by 8
static void
new_ifft_layer(cmplx_t* xp, long blocks, long size, const cmplx_t* wtab)
{

  size /= 2;
  const cmplx_t* PGFFT_RESTRICT wtab1 = wtab + size;

  do
    {

      cmplx_t* PGFFT_RESTRICT xp0 = xp;
      cmplx_t* PGFFT_RESTRICT xp1 = xp + size;


      // first 4 butterflies
      inv_butterfly0(xp0[0], xp1[0]);
      inv_butterfly_neg(xp0[1], xp1[1], wtab1[-1]);
      inv_butterfly_neg(xp0[2], xp1[2], wtab1[-2]);
      inv_butterfly_neg(xp0[3], xp1[3], wtab1[-3]);

      // 4-way unroll
      for (long j = 4; j < size; j+= 4) {
         inv_butterfly_neg(xp0[j+0], xp1[j+0], wtab1[-(j+0)]);
         inv_butterfly_neg(xp0[j+1], xp1[j+1], wtab1[-(j+1)]);
         inv_butterfly_neg(xp0[j+2], xp1[j+2], wtab1[-(j+2)]);
         inv_butterfly_neg(xp0[j+3], xp1[j+3], wtab1[-(j+3)]);
      }

      xp += 2 * size;
    }
  while (--blocks != 0);
}

static void
new_ifft_first_two_layers(cmplx_t* xp, long blocks, const cmplx_t* wtab)
{
  // 4th root of unity
  cmplx_t w = wtab[1];

  do
    {
      cmplx_t u0 = xp[0];
      cmplx_t u1 = xp[1];
      cmplx_t u2 = xp[2];
      cmplx_t u3 = xp[3];

      cmplx_t v0 = u0 + u1;
      cmplx_t v1 = u0 - u1;
      cmplx_t v2 = u2 + u3;
      cmplx_t t  = u2 - u3;
      cmplx_t v3 = t * w;

      xp[0] = v0 + v2;
      xp[2] = v0 - v2;
      xp[1] = v1 - v3;  // NEG
      xp[3] = v1 + v3;  // NEG

      xp += 4;
    }
  while (--blocks != 0);
}


static void
new_ifft_base(cmplx_t* xp, long lgN, const vector<vector<cmplx_t>>& tab)
{
  if (lgN == 0) return;


  if (lgN == 1)
    {
      cmplx_t x0 = xp[0];
      cmplx_t x1 = xp[1];
      xp[0] = x0 + x1;
      xp[1] = x0 - x1;
      return;
    }


  long blocks = 1L << (lgN - 2);
  new_ifft_first_two_layers(xp, blocks, &tab[2][0]);
  blocks >>= 1;

  long size = 8;
  for (long j = 3; j <= lgN; j++, blocks >>= 1, size <<= 1)
    new_ifft_layer(xp, blocks, size, &tab[j][0]);
}

static
void new_ifft_short2(cmplx_t* yp, long yn, long lgN, const vector<vector<cmplx_t>>& tab);



static
void new_ifft_short1(cmplx_t* xp, long yn, long lgN, const vector<vector<cmplx_t>>& tab)

// Implements truncated inverse FFT interface, but with xn==yn.
// All computations are done in place.

{
  long N = 1L << lgN;

  if (yn == N && lgN <= PGFFT_NEW_FFT_THRESH)
    {
      // no truncation
      new_ifft_base(xp, lgN, tab);
      return;
    }

  // divide-and-conquer algorithm

  long half = N >> 1;

  if (yn <= half)
    {
      // X -> 2X
      for (long j = 0; j < yn; j++)
      	xp[j] = 2.0 * xp[j];

      new_ifft_short1(xp, yn, lgN - 1, tab);
    }
  else
    {
      cmplx_t* PGFFT_RESTRICT xp0 = xp;
      cmplx_t* PGFFT_RESTRICT xp1 = xp + half;
      const cmplx_t* PGFFT_RESTRICT wtab = &tab[lgN][0];

      new_ifft_short1(xp0, half, lgN - 1, tab);

      yn -= half;

      // X -> (2X, w*X)
      for (long j = yn; j < half; j++)
	{
	  cmplx_t x0 = xp0[j];
	  xp0[j] = 2.0 * x0;
	  xp1[j] = x0 * wtab[j];
	}

      new_ifft_short2(xp1, yn, lgN - 1, tab);

      // (X, Y) -> (X + Y/w, X - Y/w)
      {
	const cmplx_t* PGFFT_RESTRICT wtab1 = wtab + half;

	// DIRT: assumes yn is a multiple of 4
	inv_butterfly0(xp0[0], xp1[0]);
	inv_butterfly_neg(xp0[1], xp1[1], wtab1[-1]);
	inv_butterfly_neg(xp0[2], xp1[2], wtab1[-2]);
	inv_butterfly_neg(xp0[3], xp1[3], wtab1[-3]);
	for (long j = 4; j < yn; j+=4) {
	  inv_butterfly_neg(xp0[j+0], xp1[j+0], wtab1[-(j+0)]);
	  inv_butterfly_neg(xp0[j+1], xp1[j+1], wtab1[-(j+1)]);
	  inv_butterfly_neg(xp0[j+2], xp1[j+2], wtab1[-(j+2)]);
	  inv_butterfly_neg(xp0[j+3], xp1[j+3], wtab1[-(j+3)]);
	}
      }
    }
}



static
void new_ifft_short2(cmplx_t* xp, long yn, long lgN, const vector<vector<cmplx_t>>& tab)

// Implements truncated inverse FFT interface, but with xn==N.
// All computations are done in place.

{
  long N = 1L << lgN;

  if (yn == N && lgN <= PGFFT_NEW_FFT_THRESH)
    {
      // no truncation
      new_ifft_base(xp, lgN, tab);
      return;
    }

  // divide-and-conquer algorithm

  long half = N >> 1;

  if (yn <= half)
    {
      // X -> 2X
      for (long j = 0; j < yn; j++)
     	xp[j] = 2.0 * xp[j];
      // (X, Y) -> X + Y
      for (long j = yn; j < half; j++)
	xp[j] = xp[j] + xp[j + half];

      new_ifft_short2(xp, yn, lgN - 1, tab);

      // (X, Y) -> X - Y
      for (long j = 0; j < yn; j++)
	xp[j] = xp[j] - xp[j + half];
    }
  else
    {
      cmplx_t* PGFFT_RESTRICT xp0 = xp;
      cmplx_t* PGFFT_RESTRICT xp1 = xp + half;
      const cmplx_t* PGFFT_RESTRICT wtab = &tab[lgN][0];

      new_ifft_short1(xp0, half, lgN - 1, tab);

      yn -= half;


      // (X, Y) -> (2X - Y, w*(X - Y))
      for (long j = yn; j < half; j++)
	{
	  cmplx_t x0 = xp0[j];
	  cmplx_t x1 = xp1[j];
	  cmplx_t u = x0 - x1;
	  xp0[j] = x0 + u;
	  xp1[j] = u * wtab[j];
	}

      new_ifft_short2(xp1, yn, lgN - 1, tab);

      // (X, Y) -> (X + Y/w, X - Y/w)
      {
	const cmplx_t* PGFFT_RESTRICT wtab1 = wtab + half;

	// DIRT: assumes yn is a multiple of 4
	inv_butterfly0(xp0[0], xp1[0]);
	inv_butterfly_neg(xp0[1], xp1[1], wtab1[-1]);
	inv_butterfly_neg(xp0[2], xp1[2], wtab1[-2]);
	inv_butterfly_neg(xp0[3], xp1[3], wtab1[-3]);
	for (long j = 4; j < yn; j+=4) {
	  inv_butterfly_neg(xp0[j+0], xp1[j+0], wtab1[-(j+0)]);
	  inv_butterfly_neg(xp0[j+1], xp1[j+1], wtab1[-(j+1)]);
	  inv_butterfly_neg(xp0[j+2], xp1[j+2], wtab1[-(j+2)]);
	  inv_butterfly_neg(xp0[j+3], xp1[j+3], wtab1[-(j+3)]);
	}
      }
    }
}



static void 
new_ifft(cmplx_t* xp, long lgN, const vector<vector<cmplx_t>>& tab)
{
   long N = 1L << lgN;
   new_ifft_short1(xp, N, lgN, tab);
}


static void
compute_table(vector<vector<cmplx_t>>& tab, long k)
{
  if (k < 2) return;

  const ldbl pi = atan(ldbl(1)) * 4.0;

  tab.resize(k+1);
  for (long s = 2; s <= k; s++) {
    long m = 1L << s;
    tab[s].resize(m/2);
    for (long j = 0; j < m/2; j++) {
      ldbl angle = -((2 * pi) * (ldbl(j)/ldbl(m)));
      tab[s][j] = cmplx_t(cos(angle), sin(angle));
    }
  }
}

static long 
RevInc(long a, long k)
{
   long j, m;

   j = k;
   m = 1L << (k-1);

   while (j && (m & a)) {
      a ^= m;
      m >>= 1;
      j--;
   }
   if (j) a ^= m;
   return a;
}

static void
BRC_init(long k, vector<long>& rev)
{
   long n = (1L << k);
   rev.resize(n);
   long i, j;
   for (i = 0, j = 0; i < n; i++, j = RevInc(j, k))
      rev[i] = j;
}


#define PGFFT_BRC_THRESH (11)
#define PGFFT_BRC_Q (5)

// Must have PGFFT_BRC_THRESH >= 2*PGFFT_BRC_Q
// Should also have (1L << (2*PGFFT_BRC_Q)) small enough
// so that we can fit that many long's into the cache

static
void BasicBitReverseCopy(cmplx_t *B, 
                         const cmplx_t *A, long k, const vector<long>& rev)
{
   long n = 1L << k;
   long i, j;

   for (i = 0; i < n; i++)
      B[rev[i]] = A[i];
}

static void 
COBRA(cmplx_t * PGFFT_RESTRICT B, const cmplx_t * PGFFT_RESTRICT A, long k,
      const vector<long>& rev, const vector<long> rev1)
{
   long q = PGFFT_BRC_Q;
   long k1 = k - 2*q;

   vector<cmplx_t> BRC_temp(1L << (2*q));

   cmplx_t * PGFFT_RESTRICT T = &BRC_temp[0];
   const long * PGFFT_RESTRICT rev_k1 = &rev[0];
   const long * PGFFT_RESTRICT rev_q = &rev1[0];
   

   for (long b = 0; b < (1L << k1); b++) {
      long b1 = rev_k1[b]; 
      for (long a = 0; a < (1L << q); a++) {
         long a1 = rev_q[a]; 
         for (long c = 0; c < (1L << q); c++) 
            T[(a1 << q) + c] = A[(a << (k1+q)) + (b << q) + c]; 
      }

      for (long c = 0; c < (1L << q); c++) {
         long c1 = rev_q[c];
         for (long a1 = 0; a1 < (1L << q); a1++) 
            B[(c1 << (k1+q)) + (b1 << q) + a1] = T[(a1 << q) + c];
      }
   }
}


static long 
pow2_precomp(long n, vector<long>& rev, vector<long>& rev1, vector<vector<cmplx_t>>& tab)
{
   // k = least k such that 2^k >= n
   long k = 0;
   while ((1L << k) < n) k++;

   compute_table(tab, k);

   if (k <= PGFFT_BRC_THRESH) {
      BRC_init(k, rev);
   }
   else {
      long q = PGFFT_BRC_Q;
      long k1 = k - 2*q;
      BRC_init(k1, rev);
      BRC_init(q, rev1);
   }
   

   return k;
}

static void
pow2_comp(vector<cmplx_t>& a, 
                  long n, long k, const vector<long>& rev, const vector<long>& rev1,
                  const vector<vector<cmplx_t>>& tab)
{
   vector<cmplx_t> x(a);

   new_fft(&x[0], k, tab);
   if (k <= PGFFT_BRC_THRESH)
      BasicBitReverseCopy(&a[0], &x[0], k, rev);
   else
      COBRA(&a[0], &x[0], k, rev, rev1);
}

static long
bluestein_precomp(long n, vector<cmplx_t>& powers, vector<cmplx_t>& Rb, 
                  vector<vector<cmplx_t>>& tab)
{
   // k = least k such that 2^k >= 2*n-1
   long k = 0;
   while ((1L << k) < 2*n-1) k++;

   compute_table(tab, k);

   const ldbl pi = atan(ldbl(1)) * 4.0;

   powers.resize(n);
   powers[0] = 1;
   long i_sqr = 0;
   for (long i = 1; i < n; i++) {
      // i^2 = (i-1)^2 + 2*i-1
      i_sqr = (i_sqr + 2*i - 1) % (2*n);
      ldbl angle = -((2 * pi) * (ldbl(i_sqr)/ldbl(2*n)));
      powers[i] = cmplx_t(cos(angle), sin(angle));
   }

   long N = 1L << k;
   Rb.resize(N);
   for (long i = 0; i < N; i++) Rb[i] = 0;

   Rb[n-1] = 1;
   i_sqr = 0;
   for (long i = 1; i < n; i++) {
      // i^2 = (i-1)^2 + 2*i-1
      i_sqr = (i_sqr + 2*i - 1) % (2*n);
      ldbl angle = (2 * pi) * (ldbl(i_sqr)/ldbl(2*n));
      Rb[n-1+i] = Rb[n-1-i] = cmplx_t(cos(angle), sin(angle));
   }
  
   new_fft(&Rb[0], k, tab);

   return k;

}


static void
bluestein_comp(vector<cmplx_t>& a, 
                  long n, long k, const vector<cmplx_t>& powers, const vector<cmplx_t>& Rb, 
                  const vector<vector<cmplx_t>>& tab)
{
   long N = 1L << k;

   vector<cmplx_t> x(N);

   for (long i = 0; i < n; i++)
      x[i] = a[i] * powers[i];

   for (long i = n; i < N; i++)
      x[i] = 0;

   new_fft(&x[0], k, tab);

   for (long i = 0; i < N; i++)
      x[i] *= Rb[i];

   new_ifft(&x[0], k, tab);

   double Ninv = 1/double(N);
   
   for (long i = 0; i < n; i++) 
      a[i] = x[n-1+i] * powers[i] * Ninv; 

}




static long
bluestein_precomp1(long n, vector<cmplx_t>& powers, vector<cmplx_t>& Rb, 
                  vector<vector<cmplx_t>>& tab)
{
   // k = least k such that 2^k >= 2*n-1
   long k = 0;
   while ((1L << k) < 2*n-1) k++;

   compute_table(tab, k);

   const ldbl pi = atan(ldbl(1)) * 4.0;

   powers.resize(n);
   powers[0] = 1;
   long i_sqr = 0;

   if (n % 2 == 0) {
      for (long i = 1; i < n; i++) {
	 // i^2 = (i-1)^2 + 2*i-1
	 i_sqr = (i_sqr + 2*i - 1) % (2*n);
	 ldbl angle = -((2 * pi) * (ldbl(i_sqr)/ldbl(2*n)));
	 powers[i] = cmplx_t(cos(angle), sin(angle));
      }
   }
   else {
      for (long i = 1; i < n; i++) {
	 // i^2*((n+1)/2) = (i-1)^2*((n+1)/2) + i + ((n-1)/2) (mod n)
	 i_sqr = (i_sqr + i + (n-1)/2) % n;
	 ldbl angle = -((2 * pi) * (ldbl(i_sqr)/ldbl(n)));
	 powers[i] = cmplx_t(cos(angle), sin(angle));
      }
   }

   long N = 1L << k;
   Rb.resize(N);
   for (long i = 0; i < N; i++) Rb[i] = 0;

   Rb[0] = 1;
   i_sqr = 0;

   if (n % 2 == 0) {
      for (long i = 1; i < n; i++) {
	 // i^2 = (i-1)^2 + 2*i-1
	 i_sqr = (i_sqr + 2*i - 1) % (2*n);
	 ldbl angle = (2 * pi) * (ldbl(i_sqr)/ldbl(2*n));
	 Rb[i] = cmplx_t(cos(angle), sin(angle));
      }
   }
   else {
      for (long i = 1; i < n; i++) {
	 // i^2*((n+1)/2) = (i-1)^2*((n+1)/2) + i + ((n-1)/2) (mod n)
	 i_sqr = (i_sqr + i + (n-1)/2) % n;
	 ldbl angle = (2 * pi) * (ldbl(i_sqr)/ldbl(n));
	 Rb[i] = cmplx_t(cos(angle), sin(angle));
      }
   }
  
   new_fft(&Rb[0], k, tab);

   return k;

}


static void
bluestein_comp1(vector<cmplx_t>& a, 
                  long n, long k, const vector<cmplx_t>& powers, const vector<cmplx_t>& Rb, 
                  const vector<vector<cmplx_t>>& tab)
{
   long N = 1L << k;

   vector<cmplx_t> x(N);

   for (long i = 0; i < n; i++)
      x[i] = a[i] * powers[i];

   long len = FFTRoundUp(2*n-1, k);
   long ilen = FFTRoundUp(n, k);

   for (long i = n; i < ilen; i++)
      x[i] = 0;

   new_fft_short(&x[0], len, ilen, k, tab);

   for (long i = 0; i < len; i++)
      x[i] *= Rb[i];

   new_ifft_short1(&x[0], len, k, tab);

   double Ninv = 1/double(N);
   
   for (long i = 0; i < n-1; i++) 
      a[i] = (x[i] + x[n+i]) * powers[i] * Ninv; 

   a[n-1] = x[n-1] * powers[n-1] * Ninv;

}

#define PGFFT_STRATEGY_NULL  (0)
#define PGFFT_STRATEGY_POW2  (1)
#define PGFFT_STRATEGY_BLUE  (2)
#define PGFFT_STRATEGY_TBLUE (3)

static long choose_strategy(long n)
{
   if (n == 1) return PGFFT_STRATEGY_NULL; 

   if ((n & (n - 1)) == 0) return PGFFT_STRATEGY_POW2; 

   if (!PGFFT_USE_TRUNCATED_BLUE) return PGFFT_STRATEGY_BLUE;

   // choose between Bluestein and truncated Bluestein

   // k = least k such that 2^k >= 2*n-1
   long k = 0;
   while ((1L << k) < 2*n-1) k++;

   long rdup = FFTRoundUp(2*n-1, k);
   if (rdup == (1L << k)) return PGFFT_STRATEGY_BLUE;

   return PGFFT_STRATEGY_TBLUE;
}

PGFFT::PGFFT(long n_)
{
   assert(n_ > 0);
   n = n_;

   strategy = choose_strategy(n);

   switch (strategy) {

   case PGFFT_STRATEGY_NULL: 
      break;

   case PGFFT_STRATEGY_POW2:
      k = pow2_precomp(n, rev, rev1, tab);
      break;

   case PGFFT_STRATEGY_BLUE:
      k = bluestein_precomp(n, powers, Rb, tab);
      break;

   case PGFFT_STRATEGY_TBLUE:
      k = bluestein_precomp1(n, powers, Rb, tab);
      break;

   default: ;

   }
}

void PGFFT::apply(vector<cmplx_t>& v) const
{
   assert(long(v.size()) == n);

   switch (strategy) {

   case PGFFT_STRATEGY_NULL: 
      break;

   case PGFFT_STRATEGY_POW2:
      pow2_comp(v, n, k, rev, rev1, tab);
      break;

   case PGFFT_STRATEGY_BLUE:
      bluestein_comp(v, n, k, powers, Rb, tab);
      break;

   case PGFFT_STRATEGY_TBLUE:
      bluestein_comp1(v, n, k, powers, Rb, tab);
      break;

   default: ;

   }
}

//================== Fft ====================

typedef complex<ldbl> lcx;

namespace Fft {
   
   /* 
    * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
    * The vector can have any length. This is a wrapper function.
    */
   void transform(std::vector<lcx> &vec);
   
   
   /* 
    * Computes the inverse discrete Fourier transform (IDFT) of the given complex vector, storing the result back into the vector.
    * The vector can have any length. This is a wrapper function. This transform does not perform scaling, so the inverse is not a true inverse.
    */
   void inverseTransform(std::vector<lcx> &vec);
   
   
   /* 
    * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
    * The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
    */
   void transformRadix2(std::vector<lcx> &vec);
   
   
   /* 
    * Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
    * The vector can have any length. This requires the convolution function, which in turn requires the radix-2 FFT function.
    * Uses Bluestein's chirp z-transform algorithm.
    */
   void transformBluestein(std::vector<lcx> &vec);
   
   
   /* 
    * Computes the circular convolution of the given complex vectors. Each vector's length must be the same.
    */
   void convolve(
      const std::vector<lcx> &vecx,
      const std::vector<lcx> &vecy,
      std::vector<lcx> &vecout);
   
}



using std::complex;
using std::size_t;
using std::vector;


// Private function prototypes
static size_t reverseBits(size_t x, int n);


void Fft::transform(vector<lcx> &vec) {
   size_t n = vec.size();
   if (n == 0)
      return;
   else if ((n & (n - 1)) == 0)  // Is power of 2
      transformRadix2(vec);
   else  // More complicated algorithm for arbitrary sizes
      transformBluestein(vec);
}


void Fft::inverseTransform(vector<lcx> &vec) {
   std::transform(vec.cbegin(), vec.cend(), vec.begin(),
      static_cast<lcx (*)(const lcx &)>(std::conj));
   transform(vec);
   std::transform(vec.cbegin(), vec.cend(), vec.begin(),
      static_cast<lcx (*)(const lcx &)>(std::conj));
}


void Fft::transformRadix2(vector<lcx> &vec) {
   // Length variables
   size_t n = vec.size();
   int levels = 0;  // Compute levels = floor(log2(n))
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
      expTable[i] = lcx(cos(angle), sin(angle));
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
      if (size == n)  // Prevent overflow in 'size *= 2'
         break;
   }
}


void Fft::transformBluestein(vector<lcx> &vec) {
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
      // Less accurate alternative if long long is unavailable: double angle = M_PI * i * i / n;
      expTable[i] = lcx(cos(-angle), sin(-angle));
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


void Fft::convolve(
      const vector<lcx> &xvec,
      const vector<lcx> &yvec,
      vector<lcx> &outvec) {
   
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
   for (size_t i = 0; i < n; i++)  // Scaling (because this FFT implementation omits it)
      outvec[i] = xv[i] / static_cast<ldbl>(n);
}


static size_t reverseBits(size_t x, int n) {
   size_t result = 0;
   for (int i = 0; i < n; i++, x >>= 1)
      result = (result << 1) | (x & 1U);
   return result;
}



//===========================================


#define TIME_IT(t, action) \
do { \
   double _t0, _t1; \
   long _iter = 1; \
   long _cnt = 0; \
   do { \
      _t0 = NTL::GetTime(); \
      for (long _i = 0; _i < _iter; _i++) { action; _cnt++; } \
      _t1 = NTL::GetTime(); \
   } while ( _t1 - _t0 < 3 && (_iter *= 2)); \
   t = (_t1 - _t0)/_iter; \
} while(0)


static cmplx_t check_sum = 0;


static void
TestIt(long n)
{
   PGFFT pgfft(n);

   for (long i = 0; i < 10; i++) {
      vector<cmplx_t> v(n);
      for (int i = 0; i < n; i++)
	 v[i] = NTL::RandomBnd(20)-10;

      vector<cmplx_t> v0(v);

      pgfft.apply(v);

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

      ldbl rel_err = diff_norm/vv_norm;

      if (rel_err > 1e-9) {
         std::cout << "BAD: " << n << " " << rel_err << " " << choose_strategy(n) <<  "\n";
      }
   }
}

int main()
{
#if 0
   long k = 6;
   long n = 64;

   vector<vector<cmplx_t>> tab;
   compute_table(tab, k);

   vector<cmplx_t> v(n);
   for (int i = 0; i < n; i++)
      v[i] = ((i+1)%6)*((i+1)%6); 


   //new_fft_base(&v[0], k, tab);
   new_fft(&v[0], k, tab);

   //new_ifft_base(&v[0], k, tab);
   new_ifft(&v[0], k, tab);

   for (long i = 0; i < n; i++) v[i] /= double(n);

#if 0
   vector<long> rev;
   BRC_init(k, rev);

   vector<cmplx_t> w(n);
   BasicBitReverseCopy(&w[0], &v[0], k, rev);
#endif

   for (int i = 0; i < n; i++)
      std::cout << v[i] << "\n";
#elif 0
   for (long n = 1; n <= 100; n++)
      TestIt(n);

   for (long n = 256; n <= 32*1024; n *= 2)
      TestIt(n);

   for (long i = 0; i < 100; i++) {
      long n = NTL::RandomBnd(10000) + 10000;
      TestIt(n);
   }
#elif 1
   //long n = 28679;
   long n = 10000;
   // long n = 45551;
   // long n = 17;

   vector<cmplx_t> v(n);
   for (int i = 0; i < n; i++)
      v[i] = NTL::RandomBnd(20)-10;

   vector<cmplx_t> v0(v);



   PGFFT pgfft(n);
   pgfft.apply(v);


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

   ldbl rel_err = diff_norm/vv_norm;

   std::cout << log(rel_err)/log(2.0) << "\n";
   
   
#else
   //long n = 28679;
   long n = 1000;
   //long n = 45551;
   //long n = 90001;

   vector<cmplx_t> v(n);
   for (long i = 0; i < n; i++)
      v[i] = NTL::RandomBnd(10)-5;

   vector<cmplx_t> w = v;

   PGFFT pgfft(n);

   double t;
   TIME_IT(t, (v=w, pgfft.apply(v), check_sum += v[0])); 

   std::cout << t << "\n";
   std::cout << check_sum << "\n";


#endif
}
