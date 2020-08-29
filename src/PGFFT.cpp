/****************************************************************************

PGFFT: Pretty Good FFT (v1.8)

Copyright (C) 2019, victor Shoup

See below for more details.

****************************************************************************/



#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#ifdef __GNUC__
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif
#endif

#define PGFFT_USE_TRUNCATED_BLUE (1)
// set to 0 to disable the truncated Bluestein

#define PGFFT_USE_EXPLICIT_MUL (1)
// Set to 0 to disable explict complex multiplication.
// The built-in complex multiplication routines are
// incredibly slow, because the standard requires special handling
// of non-finite complex values.
// To fix the problem, by default, PGFFT will override these routines
// with explicitly defined multiplication functions.
// Another way to solve this problem is to compile with the -ffast-math
// option (at least, on gcc, that's the right flag).


//============================================

#ifndef PGFFT_DISABLE_SIMD

#ifdef __AVX__
#define HAVE_AVX
//#warning "HAVE_AVX"
#endif

#ifdef __AVX2__
#define HAVE_AVX2
//#warning "HAVE_AVX2"
#endif

#if defined(HAVE_AVX) || defined(HAVE_AVX2)
#define USE_PD4
#endif

#endif




#include <helib/PGFFT.h>
#include <cassert>
#include <cstdlib>
#include <limits>

#ifdef USE_PD4
#include <immintrin.h>
#endif

namespace helib {

using std::vector;
using std::complex;

template<class T>
using aligned_vector = PGFFT::aligned_vector<T>;

typedef complex<double> cmplx_t;
typedef long double ldbl;
//typedef double ldbl;

#ifdef USE_PD4
bool PGFFT::simd_enabled() { return true; }
#else
bool PGFFT::simd_enabled() { return false; }
#endif


#if (PGFFT_USE_EXPLICIT_MUL)

static inline cmplx_t
MUL(cmplx_t a, cmplx_t b)
{
   double x = a.real(), y = a.imag(), u = b.real(), v = b.imag();
   return cmplx_t(x*u-y*v, x*v+y*u);
}

static inline cmplx_t
CMUL(cmplx_t a, cmplx_t b)
{
   double x = a.real(), y = a.imag(), u = b.real(), v = b.imag();
   return cmplx_t(x*u+y*v, y*u-x*v);
}

#else

#define MUL(a, b) ((a) * (b))
#define CMUL(a, b) ((a) * std::conj(b))

#endif


#if (defined(__GNUC__) && (__GNUC__ >= 4))

// on relative modern versions of gcc, we can
// decalare "restricted" pointers in C++

#define RESTRICT __restrict

#else

#define RESTRICT

#endif

/**************************************************************

   Aligned allocation

**************************************************************/

#ifdef USE_PD4

#define PGFFT_ALIGN (64)

void *
PGFFT::aligned_allocate(std::size_t n, std::size_t nelts)
{
   if (n > std::numeric_limits<std::size_t>::max() / nelts) return 0;
   std::size_t sz = n * nelts;
   std::size_t alignment = PGFFT_ALIGN;
   if (sz > std::numeric_limits<std::size_t>::max() - alignment) return 0;

   sz += alignment;
   char* buf = (char*) std::malloc(sz);

   if (!buf) return 0;

   int remainder = ((unsigned long long)buf) % alignment;
   int offset = alignment - remainder;
   char* ret = buf + offset;

   ret[-1] = offset;

   return ret;
}


void
PGFFT::aligned_deallocate(void *p)
{
   if (!p) return;
   char *cp = (char *) p;
   int offset = cp[-1];
   std::free(cp - offset);
}

#else

void *
PGFFT::aligned_allocate(std::size_t n, std::size_t sz)
{
   if (n > std::numeric_limits<std::size_t>::max() / sz) return 0;
   std::size_t size = n * sz;
   return std::malloc(size);
}

void
PGFFT::aligned_deallocate(void *p)
{
   if (!p) return;
   std::free(p);
}

#endif

/**************************************************************

   Packed Double abstraction layer

**************************************************************/



namespace {



//=================== PD4 implementation ===============

#if defined(USE_PD4)

struct PD4 {
   __m256d data;


   PD4() = default;
   PD4(double x) : data(_mm256_set1_pd(x)) { }
   PD4(__m256d _data) : data(_data) { }
   PD4(double d0, double d1, double d2, double d3)
      : data(_mm256_set_pd(d3, d2, d1, d0)) { }

   static PD4 load(const double *p) { return _mm256_load_pd(p); }

   // load from unaligned address
   static PD4 loadu(const double *p) { return _mm256_loadu_pd(p); }
};

inline void
load(PD4& x, const double *p)
{ x = PD4::load(p); }

// load from unaligned address
inline void
loadu(PD4& x, const double *p)
{ x = PD4::loadu(p); }

inline void
store(double *p, PD4 a)
{ _mm256_store_pd(p, a.data); }

// store to unaligned address
inline void
storeu(double *p, PD4 a)
{ _mm256_storeu_pd(p, a.data); }


// swap even/odd slots
// e.g., 0123 -> 1032
inline PD4
swap2(PD4 a)
{ return _mm256_permute_pd(a.data, 0x5); }

// 0123 -> 0022
inline PD4
dup2even(PD4 a)
{ return _mm256_permute_pd(a.data, 0);   }

// 0123 -> 1133
inline PD4
dup2odd(PD4 a)
{ return _mm256_permute_pd(a.data, 0xf);   }

// blend even/odd slots
// 0123, 4567 -> 0527
inline PD4
blend2(PD4 a, PD4 b)
{ return _mm256_blend_pd(a.data, b.data, 0xa); }

// 0123, 4567 -> 0426
inline PD4
blend_even(PD4 a, PD4 b)
{ return _mm256_unpacklo_pd(a.data, b.data); }


// 0123, 4567 -> 1537
inline PD4
blend_odd(PD4 a, PD4 b)
{ return _mm256_unpackhi_pd(a.data, b.data); }


inline void
clear(PD4& x)
{ x.data = _mm256_setzero_pd(); }

inline PD4
operator+(PD4 a, PD4 b)
{ return _mm256_add_pd(a.data, b.data); }

inline PD4
operator-(PD4 a, PD4 b)
{ return _mm256_sub_pd(a.data, b.data); }

inline PD4
operator*(PD4 a, PD4 b)
{ return _mm256_mul_pd(a.data, b.data); }

inline PD4
operator/(PD4 a, PD4 b)
{ return _mm256_div_pd(a.data, b.data); }

inline PD4&
operator+=(PD4& a, PD4 b)
{ a = a + b; return a; }

inline PD4&
operator-=(PD4& a, PD4 b)
{ a = a - b; return a; }

inline PD4&
operator*=(PD4& a, PD4 b)
{ a = a * b; return a; }

inline PD4&
operator/=(PD4& a, PD4 b)
{ a = a / b; return a; }

#ifdef HAVE_AVX2

// a*b+c (fused)
inline PD4
fused_muladd(PD4 a, PD4 b, PD4 c)
{ return _mm256_fmadd_pd(a.data, b.data, c.data); }
// NEEDS: FMA

// a*b-c (fused)
inline PD4
fused_mulsub(PD4 a, PD4 b, PD4 c)
{ return _mm256_fmsub_pd(a.data, b.data, c.data); }
// NEEDS: FMA

// -a*b+c (fused)
inline PD4
fused_negmuladd(PD4 a, PD4 b, PD4 c)
{ return _mm256_fnmadd_pd(a.data, b.data, c.data); }
// NEEDS: FMA

// (a0,a1,a2,a3), (b0,b1,b2,b3), (c0,c1,c2,c3) ->
// (a0*b0-c0, a1*b1+c1, a2*b2-c2, a3*b3+c3)
inline PD4
fmaddsub(PD4 a, PD4 b, PD4 c)
{ return _mm256_fmaddsub_pd(a.data, b.data, c.data); }
// NEEDS: FMA
// (plain addsub only needs AVX)

// (a0,a1,a2,a3), (b0,b1,b2,b3), (c0,c1,c2,c3) ->
// (a0*b0+c0, a1*b1-c1, a2*b2+c2, a3*b3-c3)
inline PD4
fmsubadd(PD4 a, PD4 b, PD4 c)
{ return _mm256_fmsubadd_pd(a.data, b.data, c.data); }
// NEEDS: FMA
// (there is no plain subadd)
#endif


#endif

}


/***************************************************************


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
   xx1 = MUL(t_, w); \
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


#define inv_butterfly0(xx0, xx1)  \
do   \
{  \
   cmplx_t x0_ = xx0;  \
   cmplx_t x1_ = xx1;  \
   xx0 = x0_ + x1_;  \
   xx1 = x0_ - x1_;  \
} while (0)


#define inv_butterfly(xx0, xx1, w)  \
do  \
{  \
   cmplx_t x0_ = xx0;  \
   cmplx_t x1_ = xx1;  \
   cmplx_t t_ = CMUL(x1_, w);  \
   xx0 = x0_ + t_;  \
   xx1 = x0_ - t_;  \
} while (0)



#ifdef USE_PD4

#ifdef HAVE_AVX2

static inline PD4
complex_mul(PD4 ab, PD4 cd)
{
   PD4 cc = dup2even(cd);
   PD4 dd = dup2odd(cd);
   PD4 ba = swap2(ab);
   return fmaddsub(ab, cc, ba*dd);
}

static inline PD4
complex_conj_mul(PD4 ab, PD4 cd)
// (ac+bd,bc-ad)
{
   PD4 cc = dup2even(cd);
   PD4 dd = dup2odd(cd);
   PD4 ba = swap2(ab);
   return fmsubadd(ab, cc, ba*dd);
}


#define MUL2(x_0, x_1, a_0, a_1, b_0, b_1) \
do { \
   x_0 = complex_mul(a_0, b_0); \
   x_1 = complex_mul(a_1, b_1); \
} while (0)

#define CMUL2(x_0, x_1, a_0, a_1, b_0, b_1) \
do { \
   x_0 = complex_conj_mul(a_0, b_0); \
   x_1 = complex_conj_mul(a_1, b_1); \
} while (0)

#else
// This code sequence works without FMA
#define MUL2(x_0, x_1, a_0, a_1, b_0, b_1) \
do { \
    PD4 a_re_ = blend_even(a_0, a_1); \
    PD4 a_im_ = blend_odd(a_0, a_1); \
 \
    PD4 b_re_ = blend_even(b_0, b_1); \
    PD4 b_im_ = blend_odd(b_0, b_1); \
 \
    PD4 x_re_ = a_re_*b_re_ - a_im_*b_im_; \
    PD4 x_im_ = a_re_*b_im_ + a_im_*b_re_; \
 \
    x_0 = blend_even(x_re_, x_im_); \
    x_1 = blend_odd(x_re_, x_im_); \
} while (0)

#define CMUL2(x_0, x_1, a_0, a_1, b_0, b_1) \
do { \
    PD4 a_re_ = blend_even(a_0, a_1); \
    PD4 a_im_ = blend_odd(a_0, a_1); \
 \
    PD4 b_re_ = blend_even(b_0, b_1); \
    PD4 b_im_ = blend_odd(b_0, b_1); \
 \
    PD4 x_re_ = a_re_*b_re_ + a_im_*b_im_; \
    PD4 x_im_ = a_im_*b_re_ - a_re_*b_im_; \
 \
    x_0 = blend_even(x_re_, x_im_); \
    x_1 = blend_odd(x_re_, x_im_); \
} while (0)

#endif



static inline void
fwd_butterfly_loop_simd(
   long size,
   double * RESTRICT xp0,
   double * RESTRICT xp1,
   const double * RESTRICT wtab)
{
  for (long j = 0; j < size; j += 4) {
    PD4 x0_0 = PD4::load(xp0+2*(j+0));
    PD4 x0_1 = PD4::load(xp0+2*(j+2));
    PD4 x1_0 = PD4::load(xp1+2*(j+0));
    PD4 x1_1 = PD4::load(xp1+2*(j+2));
    PD4 w_0  = PD4::load(wtab+2*(j+0));
    PD4 w_1  = PD4::load(wtab+2*(j+2));

    PD4 xx0_0 = x0_0 + x1_0;
    PD4 xx0_1 = x0_1 + x1_1;

    PD4 diff_0 = x0_0 - x1_0;
    PD4 diff_1 = x0_1 - x1_1;

    PD4 xx1_0, xx1_1;
    MUL2(xx1_0, xx1_1, diff_0, diff_1, w_0, w_1);

    store(xp0+2*(j+0), xx0_0);
    store(xp0+2*(j+2), xx0_1);
    store(xp1+2*(j+0), xx1_0);
    store(xp1+2*(j+2), xx1_1);
  }
}

static inline void
fwd_butterfly_loop(
   long size,
   cmplx_t * RESTRICT xp0,
   cmplx_t * RESTRICT xp1,
   const cmplx_t * RESTRICT wtab)
{
   // NOTE: C++11 guarantees that these reinterpret_cast's work as expected
   fwd_butterfly_loop_simd(
      size,
      reinterpret_cast<double*>(xp0),
      reinterpret_cast<double*>(xp1),
      reinterpret_cast<const double*>(wtab));
}

static inline void
inv_butterfly_loop_simd(
   long size,
   double * RESTRICT xp0,
   double * RESTRICT xp1,
   const double * RESTRICT wtab)
{
  for (long j = 0; j < size; j += 4) {
    PD4 x0_0 = PD4::load(xp0+2*(j+0));
    PD4 x0_1 = PD4::load(xp0+2*(j+2));
    PD4 x1_0 = PD4::load(xp1+2*(j+0));
    PD4 x1_1 = PD4::load(xp1+2*(j+2));
    PD4 w_0  = PD4::load(wtab+2*(j+0));
    PD4 w_1  = PD4::load(wtab+2*(j+2));

    PD4 t_0, t_1;
    CMUL2(t_0, t_1, x1_0, x1_1, w_0, w_1);

    PD4 xx0_0 = x0_0 + t_0;
    PD4 xx0_1 = x0_1 + t_1;

    PD4 xx1_0 = x0_0 - t_0;
    PD4 xx1_1 = x0_1 - t_1;

    store(xp0+2*(j+0), xx0_0);
    store(xp0+2*(j+2), xx0_1);
    store(xp1+2*(j+0), xx1_0);
    store(xp1+2*(j+2), xx1_1);
  }
}

static inline void
inv_butterfly_loop(
   long size,
   cmplx_t * RESTRICT xp0,
   cmplx_t * RESTRICT xp1,
   const cmplx_t * RESTRICT wtab)
{
   // NOTE: C++11 guarantees that these reinterpret_cast's work as expected
   inv_butterfly_loop_simd(
      size,
      reinterpret_cast<double*>(xp0),
      reinterpret_cast<double*>(xp1),
      reinterpret_cast<const double*>(wtab));
}

#else

static inline void
fwd_butterfly_loop(
   long size,
   cmplx_t * RESTRICT xp0,
   cmplx_t * RESTRICT xp1,
   const cmplx_t * RESTRICT wtab)
{
   fwd_butterfly0(xp0[0+0], xp1[0+0]);
   fwd_butterfly(xp0[0+1], xp1[0+1], wtab[0+1]);
   fwd_butterfly(xp0[0+2], xp1[0+2], wtab[0+2]);
   fwd_butterfly(xp0[0+3], xp1[0+3], wtab[0+3]);
   for (long j = 4; j < size; j += 4) {
     fwd_butterfly(xp0[j+0], xp1[j+0], wtab[j+0]);
     fwd_butterfly(xp0[j+1], xp1[j+1], wtab[j+1]);
     fwd_butterfly(xp0[j+2], xp1[j+2], wtab[j+2]);
     fwd_butterfly(xp0[j+3], xp1[j+3], wtab[j+3]);
   }
}

static inline void
inv_butterfly_loop(
   long size,
   cmplx_t * RESTRICT xp0,
   cmplx_t * RESTRICT xp1,
   const cmplx_t * RESTRICT wtab)
{
   inv_butterfly0(xp0[0+0], xp1[0+0]);
   inv_butterfly(xp0[0+1], xp1[0+1], wtab[0+1]);
   inv_butterfly(xp0[0+2], xp1[0+2], wtab[0+2]);
   inv_butterfly(xp0[0+3], xp1[0+3], wtab[0+3]);
   for (long j = 4; j < size; j += 4) {
     inv_butterfly(xp0[j+0], xp1[j+0], wtab[j+0]);
     inv_butterfly(xp0[j+1], xp1[j+1], wtab[j+1]);
     inv_butterfly(xp0[j+2], xp1[j+2], wtab[j+2]);
     inv_butterfly(xp0[j+3], xp1[j+3], wtab[j+3]);
   }
}

#endif


#if (defined(USE_PD4))


static inline void
mul_loop_simd(
   long size,
   double * RESTRICT xp,
   const double * yp)
{
  long j;
  for (j = 0; j < size; j += 4) {
    PD4 x_0 = PD4::load(xp+2*(j+0));
    PD4 x_1 = PD4::load(xp+2*(j+2));
    PD4 y_0 = PD4::load(yp+2*(j+0));
    PD4 y_1 = PD4::load(yp+2*(j+2));

    PD4 z_0, z_1;
    MUL2(z_0, z_1, x_0, x_1, y_0, y_1);

    store(xp+2*(j+0), z_0);
    store(xp+2*(j+2), z_1);
  }
}


static inline void
mul_loop(
   long size,
   cmplx_t * xp,
   const cmplx_t * yp)
{
   // NOTE: C++11 guarantees that these reinterpret_cast's work as expected
   mul_loop_simd(
      size,
      reinterpret_cast<double*>(xp),
      reinterpret_cast<const double*>(yp));
}


#else


static inline void
mul_loop(
   long size,
   cmplx_t * xp,
   const cmplx_t * yp)
{
  for (long j = 0; j < size; j++)
    xp[j] = MUL(xp[j], yp[j]);
}

#endif


// requires size divisible by 8
static void
new_fft_layer(cmplx_t* xp, long blocks, long size,
              const cmplx_t* RESTRICT wtab)
{
  size /= 2;

  do
    {
      cmplx_t* RESTRICT xp0 = xp;
      cmplx_t* RESTRICT xp1 = xp + size;

      fwd_butterfly_loop(size, xp0, xp1, wtab);

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

      //cmplx_t v3 = MUL(t, w);
      // DIRT: relies on w == (0,-1)
      cmplx_t v3(t.imag(), -t.real());


      xp[0] = v0 + v1;
      xp[1] = v0 - v1;
      xp[2] = v2 + v3;
      xp[3] = v2 - v3;

      xp += 4;
    }
  while (--blocks != 0);
}


static void
new_fft_base(cmplx_t* xp, long lgN, const vector<aligned_vector<cmplx_t>>& tab)
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
                   const vector<aligned_vector<cmplx_t>>& tab)
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

      cmplx_t* RESTRICT xp0 = xp;
      cmplx_t* RESTRICT xp1 = xp + half;
      const cmplx_t* RESTRICT wtab = &tab[lgN][0];

      if (xn <= half)
	{
	  // X -> (X, w*X)
	  for (long j = 0; j < xn; j++)
	    xp1[j] = MUL(xp0[j], wtab[j]);

	  new_fft_short(xp0, half, xn, lgN - 1, tab);
	  new_fft_short(xp1, yn, xn, lgN - 1, tab);
	}
      else
	{
	  xn -= half;

	  // (X, Y) -> (X + Y, w*(X - Y))
          // DIRT: assumes xn is a multiple of 4
          fwd_butterfly_loop(xn, xp0, xp1, wtab);

	  // X -> (X, w*X)
	  for (long j = xn; j < half; j++)
	    xp1[j] = MUL(xp0[j], wtab[j]);

	  new_fft_short(xp0, half, half, lgN - 1, tab);
	  new_fft_short(xp1, yn, half, lgN - 1, tab);
	}
    }
}

static void new_fft(cmplx_t* xp, long lgN, const vector<aligned_vector<cmplx_t>>& tab)
{
   long N = 1L << lgN;
   new_fft_short(xp, N, N, lgN, tab);
}





// requires size divisible by 8
static void
new_ifft_layer(cmplx_t* xp, long blocks, long size,
               const cmplx_t* RESTRICT wtab)
{

  size /= 2;

  do
    {

      cmplx_t* RESTRICT xp0 = xp;
      cmplx_t* RESTRICT xp1 = xp + size;


      inv_butterfly_loop(size, xp0, xp1, wtab);

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

      //cmplx_t v3 = CMUL(t, w);
      // DIRT: relies on w == (0,1)
      cmplx_t v3(-t.imag(), t.real());

      xp[0] = v0 + v2;
      xp[2] = v0 - v2;
      xp[1] = v1 + v3;
      xp[3] = v1 - v3;

      xp += 4;
    }
  while (--blocks != 0);
}


static void
new_ifft_base(cmplx_t* xp, long lgN, const vector<aligned_vector<cmplx_t>>& tab)
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
void new_ifft_short2(cmplx_t* yp, long yn, long lgN, const vector<aligned_vector<cmplx_t>>& tab);



static
void new_ifft_short1(cmplx_t* xp, long yn, long lgN, const vector<aligned_vector<cmplx_t>>& tab)

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
      cmplx_t* RESTRICT xp0 = xp;
      cmplx_t* RESTRICT xp1 = xp + half;
      const cmplx_t* RESTRICT wtab = &tab[lgN][0];

      new_ifft_short1(xp0, half, lgN - 1, tab);

      yn -= half;

      // X -> (2X, w*X)
      for (long j = yn; j < half; j++)
	{
	  cmplx_t x0 = xp0[j];
	  xp0[j] = 2.0 * x0;
	  xp1[j] = MUL(x0, wtab[j]);
	}

      new_ifft_short2(xp1, yn, lgN - 1, tab);

      // (X, Y) -> (X + Y/w, X - Y/w)
      {
        inv_butterfly_loop(yn, xp0, xp1, wtab);
      }
    }
}



static
void new_ifft_short2(cmplx_t* xp, long yn, long lgN, const vector<aligned_vector<cmplx_t>>& tab)

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
      cmplx_t* RESTRICT xp0 = xp;
      cmplx_t* RESTRICT xp1 = xp + half;
      const cmplx_t* RESTRICT wtab = &tab[lgN][0];

      new_ifft_short1(xp0, half, lgN - 1, tab);

      yn -= half;


      // (X, Y) -> (2X - Y, w*(X - Y))
      for (long j = yn; j < half; j++)
	{
	  cmplx_t x0 = xp0[j];
	  cmplx_t x1 = xp1[j];
	  cmplx_t u = x0 - x1;
	  xp0[j] = x0 + u;
	  xp1[j] = MUL(u, wtab[j]);
	}

      new_ifft_short2(xp1, yn, lgN - 1, tab);

      // (X, Y) -> (X + Y/w, X - Y/w)
      {
        inv_butterfly_loop(yn, xp0, xp1, wtab);
      }
    }
}



static void
new_ifft(cmplx_t* xp, long lgN, const vector<aligned_vector<cmplx_t>>& tab)
{
   long N = 1L << lgN;
   new_ifft_short1(xp, N, lgN, tab);
}


static void
compute_table(vector<aligned_vector<cmplx_t>>& tab, long k)
{
  if (k < 2) return;

  const ldbl pi = std::atan(ldbl(1)) * 4.0;

  tab.resize(k+1);
  for (long s = 2; s <= k; s++) {
    long m = 1L << s;
    tab[s].resize(m/2);
    for (long j = 0; j < m/2; j++) {
      ldbl angle = -((2 * pi) * (ldbl(j)/ldbl(m)));
      tab[s][j] = cmplx_t(std::cos(angle), std::sin(angle));
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
// so that we can fit that many cmplx_t's into the cache


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
COBRA(cmplx_t * RESTRICT B, const cmplx_t * RESTRICT A, long k,
      const vector<long>& rev, const vector<long> rev1)
{
   constexpr long q = PGFFT_BRC_Q;
   long k1 = k - 2*q;

   aligned_vector<cmplx_t> BRC_temp(1L << (2*q));

   cmplx_t * RESTRICT T = &BRC_temp[0];
   const long * RESTRICT rev_k1 = &rev[0];
   const long * RESTRICT rev_q = &rev1[0];


   for (long b = 0; b < (1L << k1); b++) {
      long b1 = rev_k1[b];
      for (long a = 0; a < (1L << q); a++) {
         long a1 = rev_q[a];
         cmplx_t *T_p = &T[a1 << q];
         const cmplx_t *A_p = &A[(a << (k1+q)) + (b << q)];
#ifdef USE_PD4
         for (long c = 0; c < (1 << q); c += 4) {
            PD4 x0 = PD4::load(reinterpret_cast<const double*>(&A_p[c+0]));
            PD4 x1 = PD4::load(reinterpret_cast<const double*>(&A_p[c+2]));
            store(reinterpret_cast<double*>(&T_p[c+0]), x0);
            store(reinterpret_cast<double*>(&T_p[c+2]), x1);
         }
#else
         for (long c = 0; c < (1L << q); c++) T_p[c] = A_p[c];
#endif
      }

      for (long c = 0; c < (1L << q); c++) {
         long c1 = rev_q[c];
         cmplx_t *B_p = &B[(c1 << (k1+q)) + (b1 << q)];
         cmplx_t *T_p = &T[c];
         for (long a1 = 0; a1 < (1l << q); a1++)
            B_p[a1] = T_p[a1 << q];
      }
   }
}


static long
pow2_precomp(long n, vector<long>& rev, vector<long>& rev1, vector<aligned_vector<cmplx_t>>& tab)
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
pow2_comp(const cmplx_t* src, cmplx_t* dst,
                  long n, long k, const vector<long>& rev, const vector<long>& rev1,
                  const vector<aligned_vector<cmplx_t>>& tab)
{
   aligned_vector<cmplx_t> x;
   x.assign(src, src+n);

   new_fft(&x[0], k, tab);
#if 0
   for (long i = 0; i < n; i++) dst[i] = x[i];
#else
   if (k <= PGFFT_BRC_THRESH)
      BasicBitReverseCopy(&dst[0], &x[0], k, rev);
   else
      COBRA(&dst[0], &x[0], k, rev, rev1);
#endif
}

static long
bluestein_precomp(long n, aligned_vector<cmplx_t>& powers,
                  aligned_vector<cmplx_t>& Rb,
                  vector<aligned_vector<cmplx_t>>& tab)
{
   // k = least k such that 2^k >= 2*n-1
   long k = 0;
   while ((1L << k) < 2*n-1) k++;

   compute_table(tab, k);

   const ldbl pi = std::atan(ldbl(1)) * 4.0;

   powers.resize(n);
   powers[0] = 1;
   long i_sqr = 0;
   for (long i = 1; i < n; i++) {
      // i^2 = (i-1)^2 + 2*i-1
      i_sqr = (i_sqr + 2*i - 1) % (2*n);
      ldbl angle = -((2 * pi) * (ldbl(i_sqr)/ldbl(2*n)));
      powers[i] = cmplx_t(std::cos(angle), std::sin(angle));
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
      Rb[n-1+i] = Rb[n-1-i] = cmplx_t(std::cos(angle), std::sin(angle));
   }

   new_fft(&Rb[0], k, tab);

   double Ninv = 1/double(N);
   for (long i = 0; i < N; i++)
      Rb[i] *= Ninv;

   return k;

}


static void
bluestein_comp(const cmplx_t* src, cmplx_t* dst,
                  long n, long k, const aligned_vector<cmplx_t>& powers,
                  const aligned_vector<cmplx_t>& Rb,
                  const vector<aligned_vector<cmplx_t>>& tab)
{
   long N = 1L << k;

   aligned_vector<cmplx_t> x(N);

   for (long i = 0; i < n; i++)
      x[i] = MUL(src[i], powers[i]);

   for (long i = n; i < N; i++)
      x[i] = 0;

   new_fft(&x[0], k, tab);

   // for (long i = 0; i < N; i++) x[i] = MUL(x[i], Rb[i]);
   mul_loop(N, &x[0], &Rb[0]);

   new_ifft(&x[0], k, tab);

   double Ninv = 1/double(N);

   for (long i = 0; i < n; i++)
      dst[i] = MUL(x[n-1+i], powers[i]);

}




static long
bluestein_precomp1(long n, aligned_vector<cmplx_t>& powers,
                  aligned_vector<cmplx_t>& Rb,
                  vector<aligned_vector<cmplx_t>>& tab)
{
   // k = least k such that 2^k >= 2*n-1
   long k = 0;
   while ((1L << k) < 2*n-1) k++;

   compute_table(tab, k);

   const ldbl pi = std::atan(ldbl(1)) * 4.0;

   powers.resize(n);
   powers[0] = 1;
   long i_sqr = 0;

   if (n % 2 == 0) {
      for (long i = 1; i < n; i++) {
	 // i^2 = (i-1)^2 + 2*i-1
	 i_sqr = (i_sqr + 2*i - 1) % (2*n);
	 ldbl angle = -((2 * pi) * (ldbl(i_sqr)/ldbl(2*n)));
	 powers[i] = cmplx_t(std::cos(angle), std::sin(angle));
      }
   }
   else {
      for (long i = 1; i < n; i++) {
	 // i^2*((n+1)/2) = (i-1)^2*((n+1)/2) + i + ((n-1)/2) (mod n)
	 i_sqr = (i_sqr + i + (n-1)/2) % n;
	 ldbl angle = -((2 * pi) * (ldbl(i_sqr)/ldbl(n)));
	 powers[i] = cmplx_t(std::cos(angle), std::sin(angle));
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
	 Rb[i] = cmplx_t(std::cos(angle), std::sin(angle));
      }
   }
   else {
      for (long i = 1; i < n; i++) {
	 // i^2*((n+1)/2) = (i-1)^2*((n+1)/2) + i + ((n-1)/2) (mod n)
	 i_sqr = (i_sqr + i + (n-1)/2) % n;
	 ldbl angle = (2 * pi) * (ldbl(i_sqr)/ldbl(n));
	 Rb[i] = cmplx_t(std::cos(angle), std::sin(angle));
      }
   }

   new_fft(&Rb[0], k, tab);

   double Ninv = 1/double(N);
   for (long i = 0; i < N; i++)
      Rb[i] *= Ninv;

   return k;

}


static void
bluestein_comp1(const cmplx_t* src, cmplx_t* dst,
                  long n, long k, const aligned_vector<cmplx_t>& powers,
                  const aligned_vector<cmplx_t>& Rb,
                  const vector<aligned_vector<cmplx_t>>& tab)
{
   long N = 1L << k;

   aligned_vector<cmplx_t> x(N);

   for (long i = 0; i < n; i++)
      x[i] = MUL(src[i], powers[i]);

   long len = FFTRoundUp(2*n-1, k);
   long ilen = FFTRoundUp(n, k);

   for (long i = n; i < ilen; i++)
      x[i] = 0;

   new_fft_short(&x[0], len, ilen, k, tab);

   // for (long i = 0; i < len; i++) x[i] = MUL(x[i], Rb[i]);
   mul_loop(len, &x[0], &Rb[0]);

   new_ifft_short1(&x[0], len, k, tab);

   double Ninv = 1/double(N);

   for (long i = 0; i < n-1; i++)
      dst[i] = MUL(x[i] + x[n+i], powers[i]);

   dst[n-1] = MUL(x[n-1], powers[n-1]);

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

   //std::cout << strategy << "\n";

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

void PGFFT::apply(const cmplx_t* src, cmplx_t* dst) const
{
   switch (strategy) {

   case PGFFT_STRATEGY_NULL:
      break;

   case PGFFT_STRATEGY_POW2:
      pow2_comp(src, dst, n, k, rev, rev1, tab);
      break;

   case PGFFT_STRATEGY_BLUE:
      bluestein_comp(src, dst, n, k, powers, Rb, tab);
      break;

   case PGFFT_STRATEGY_TBLUE:
      bluestein_comp1(src, dst, n, k, powers, Rb, tab);
      break;

   default: ;

   }
}

}


/****************************************************************************

PGFFT: Pretty Good FFT (v1.8)

Copyright (C) 2019, victor Shoup

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

****************************************************************************

The logic of this code is derived from code originally developed by David Harvey,
even though the code itself has been essentially rewritten from scratch.
Here is David Harvey's original copyright notice.

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

****************************************************************************/

#pragma GCC diagnostic pop
