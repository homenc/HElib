/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
/* bluestein.cpp -
 * An implementation of non-power-of-two FFT using Bluestein's trick
 *
 * Entry points: 
 *
 *   void DFT(ZZ_pX& x, const ZZ_pX& a, long n, const ZZ_p& root);
 *   void DFT(zz_pX& x, const zz_pX& a, long n, const zz_p& root);
 *   void BluesteinFFT(ZZ_pX& x, const ZZ_pX& a, long n, const ZZ_p& root,
 *                     ZZ_pX& powers, FFTRep& Rb);
 *   void BluesteinFFT(zz_pX& x, const zz_pX& a, long n, const zz_p& root,
 *                     zz_pX& powers, fftRep& Rb);
 *
 */
#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ_pX.h>

#include "bluestein.h"

#include "timing.h"

NTL_CLIENT

/* This module builds on Shoup's NTL, and contains both a bigint version
 * with types ZZ_p and ZZ_pX and a smallint version with types zz_p and zz_pX.
 *
 * These two versions are otherwise identical, hence a function template
 * is used to describe both. (NTL itself usually uses macros rather than
 * templates for this purpose, but in 2011 compiler support for templates
 * is already stable enough to use.)  We use the name tXYZ(...) for the
 * template and XYZ(...) for the instantiated function.
 */


// NTL uses different function names for conversion between standard and
// FFT representations for ZZ_pX vs zz_pX. Hence we need a wrapper so that
// we can use the same name in the template below.
inline void ToFFTRep(fftRep& y, const zz_pX& x, long k)
   { TofftRep(y, x, k, 0, deg(x));  }
inline void FromFFTRep(zz_pX& x, fftRep& y, long lo, long hi)
   { FromfftRep(x, y, lo, hi); }


/* BluesteinFFT(x, a, n, root, powers, Rb):
 * ----------------------------------------
 * Compute length-n FFT of the coefficient-vector of a and put the result in
 * x. If the degree of a is less than n then it treats the top coefficients
 * as 0, if the degree of a is more than n then the extra coefficients are
 * ignored. Similarly, if the top entries in x are zeros then x will have
 * degree smaller than n. The argument root is a 2n-th root of unity, namely
 * BluesteinFFT(...,root,...)=DFT(...,root^2,...).
 *
 * The inverse-FFT is obtained just by calling BluesteinFFT(... root^{-1}),
 * but this procedure is *NOT SCALED*, so BluesteinFFT(x,a,n,root,...) and
 * then BluesteinFFT(a2,x,n,root^{-1},...) will result in a2=a*n.
 *
 * In addition to the size-n FFT of a which is returned in x, this procedure
 * also returns the powers of root in the powers argument:
 *    powers = [1, root, root^4, root^9, ..., root^{(n-1)^2}]
 * and in Rb it returns the size-N FFT representation of the negative
 * powers (with N>=2n-1, N a power of two):
 *    b = [0,...,0, root^{-(n-1)^2},...,root^{-4}, root^{-1}, 1, 
 *                  root^{-1},root^{-4},...,root^{-(n-1)^2}, 0...,0]
 * On subsequent calls with these 'powers' and 'Rb', these arrays are
 * not computed again but taken from these pre-comuted variables.
 *
 * If the powers and Rb arguments are initialized, then it is assumed that
 * they were computed correctly from root. The bahavior is undefined when
 * calling with initialized powers and Rb but a different root. In particular,
 * to compute the inverse-FFT (using root^{-1}), one must provide different
 * powers and Rb arguments than those that were given when computing in the
 * forward direction using root. To reset these arguments between calls
 * with different root values, use clear(powers); Rb.SetSize(0);
 *
 */

void tBluesteinFFT(zz_pX& x, long n, const zz_p& root,
		  zz_pX& powers, Vec<mulmod_precon_t>& powers_aux, 
                  fftRep& Rb, fftrep_aux& Rb_aux, fftRep& Ra)
{
  if (IsZero(x)) return;
  if (n<=0) {
    clear(x);
    return;
  }

  long p = zz_p::modulus();

  zz_p one; one=1;
  if (powers.rep.length()<n) powers.SetMaxLength(n);

  if (deg(powers)!=n-1) { // If the degree does not match, re-compute powers
    SetCoeff(powers,0,one);
    for (long i=1; i<n; i++) {
      long iSqr = MulMod(i, i, 2*n); // i^2 mod 2n
      SetCoeff(powers,i, power(root,iSqr)); // powers[i] = root^{i^2}
    }

    // powers_aux tracks powers
    powers_aux.SetLength(n);
    for (long i = 0; i < n; i++)
      powers_aux[i] = PrepMulModPrecon(rep(powers[i]), p, 1/((double) p));

  } // if deg(powers)==n, assume that it already includes powers of root

  long dx = deg(x);
  for (long i=0; i<=dx; i++) {
    x[i].LoopHole() = MulModPrecon(rep(x[i]), rep(powers[i]), p, powers_aux[i]);
  }
  x.normalize();

  long k = NextPowerOfTwo(2*n-1);
  long k2 = 1L << k; // k2 = 2^k
  Ra.SetSize(k);
  ToFFTRep(Ra, x, k);

  if (Rb.k!=k) { // If the transform size doesn't match, re-do the transform
    Rb.SetSize(k);
    zz_pX b(INIT_SIZE, k2);

    zz_p rInv = inv(root);
    SetCoeff(b,n-1,one); // b[n-1] = 1
    for (long i=1; i<n; i++) {
      long iSqr = MulMod(i, i, 2*n); // i^2 mod 2n
      zz_p bi = power(rInv,iSqr);
      SetCoeff(b,n-1+i, bi); // b[n-1+i] = b[n-1-i] = root^{-i^2}
      SetCoeff(b,n-1-i,bi);              
    }
    //    cout << "  b="<<b<<"\n";
    //    cout << "  (a.powers)*b="<<b*x<<"\n";
    ToFFTRep(Rb, b, k);

    // Rb_aux tracks Rb...
    // We are experimentally using some undocumented 
    // features of NTL.

    Rb_aux.SetLength(Rb.NumPrimes);
    for (long i = 0; i < Rb.NumPrimes; i++) {
      long q = FFTPrime[i];
      double qinv = 1/((double) q);
      long len = 1L << k;
      Rb_aux[i].SetLength(len);
      for (long j = 0; j < len; j++)
        Rb_aux[i][j] = PrepMulModPrecon(Rb.tbl[i][j], q, qinv);
    }
    

  } // if Rb.k==k, assume that Rb already contains a transform of b

  // mul(Ra,Ra,Rb);           // multiply in FFT representation

#if 1
  for (long i = 0; i < Ra.NumPrimes; i++) {
    long len = 1L << Ra.k;
    long q = FFTPrime[i];
    long *ap = &Ra.tbl[i][0];
    long *bp = &Rb.tbl[i][0];
    mulmod_precon_t *auxp = &Rb_aux[i][0]; 

    for (long j = 0; j < len; j++)
      ap[j] = MulModPrecon(ap[j], bp[j], q, auxp[j]);
  }
#endif


  FromFFTRep(x, Ra, n-1, 2*(n-1)); // then convert back
  dx = deg(x); 
  for (long i=0; i<=dx; i++) {
    x[i].LoopHole() = MulModPrecon(rep(x[i]), rep(powers[i]), p, powers_aux[i]);
  }
  x.normalize();
}


// NOTE: in the above routine, I tried a couple of optimizations.
// Namely, for the powers vector and the fftRep Rb, I tried
// storing all the MulModPrecon multipliers, to see if that would
// speed things up...unfortunately, it did not... :-(
// It had very little impact (it seems to help...but just slightly).
// Restructuring the powers multiplication seemed to actually help
// slightly more than the Rb preconditioning.
//
// I should probably clean out the code a bit, and get rid of this stuff...
// but for now, I'm leaving it in...
//    -- VJS



void tBluesteinFFT(ZZ_pX& x, long n, const ZZ_p& root,
		  ZZ_pX& powers, FFTRep& Rb, FFTRep& Ra)
{
  if (IsZero(x)) return;
  if (n<=0) {
    clear(x);
    return;
  }

  ZZ_p one; one=1;
  if (powers.rep.length()<n) powers.SetMaxLength(n);

  if (deg(powers)!=n-1) { // If the degree does not match, re-compute powers
    SetCoeff(powers,0,one);
    for (long i=1; i<n; i++) {
      long iSqr = MulMod(i, i, 2*n); // i^2 mod 2n
      SetCoeff(powers,i, power(root,iSqr)); // powers[i] = root^{i^2}
    }
  } // if deg(powers)==n, assume that it already includes powers of root

  long dx = deg(x);
  for (long i=0; i<=dx; i++) {
    x.rep[i] *= coeff(powers, i);
  }
  x.normalize();

  long k = NextPowerOfTwo(2*n-1);
  long k2 = 1L << k; // k2 = 2^k
  Ra.SetSize(k);
  ToFFTRep(Ra, x, k);

  if (Rb.k!=k) { // If the transform size doesn't match, re-do the transform
    Rb.SetSize(k);
    ZZ_pX b(INIT_SIZE, k2);

    ZZ_p rInv = inv(root);
    SetCoeff(b,n-1,one); // b[n-1] = 1
    for (long i=1; i<n; i++) {
      long iSqr = MulMod(i, i, 2*n); // i^2 mod 2n
      ZZ_p bi = power(rInv,iSqr);
      SetCoeff(b,n-1+i, bi); // b[n-1+i] = b[n-1-i] = root^{-i^2}
      SetCoeff(b,n-1-i,bi);              
    }
    //    cout << "  b="<<b<<"\n";
    //    cout << "  (a.powers)*b="<<b*x<<"\n";
    ToFFTRep(Rb, b, k);
  } // if Rb.k==k, assume that Rb already contains a transform of b

  mul(Ra,Ra,Rb);           // multiply in FFT representation
  FromFFTRep(x, Ra, n-1, 2*(n-1)); // then convert back
  dx = deg(x); 
  for (long i=0; i<=dx; i++) {
    x.rep[i] *= coeff(powers, i);
  }
  x.normalize();
}


// For degugging purposes: a slow DFT procedure
// x[k] = \sum_{i=0}^{n-1} a[i] root^{ki}
template <class zp, class zpX>
void tDFT(zpX& x, const zpX& a, long n, const zp& root)
{
  clear(x);
  if (IsZero(a) || n<=0) return;

  x.SetMaxLength(n);

  zp sum = coeff(a,0);
  for (long i=1; i<n; i++) sum += coeff(a,i);
  SetCoeff(x,0,sum);

  zp term, base = root;  // base = root^k
  for (long k=1; k<n; k++) {
    sum = coeff(a,0);
    term = base;    // term = root^{ki}
    for (long i=1; i<n; i++) {
      sum += coeff(a,i) * term; // add a[i] * root^{ki}
      term *= base; // term = root^{k(i+1)}
    }
    SetCoeff(x, k, sum);
    base *= root;          // base = root^{k+1}
  }
  x.normalize();
}

// Instantiations of the templates above for ZZ_p/ZZ_pX/FFTRep
// and for zz_p/zz_pX/fftrep

void BluesteinFFT(ZZ_pX& x, long n, const ZZ_p& root,
		  ZZ_pX& powers, Vec<mulmod_precon_t>& powers_aux, 
                  FFTRep& Rb, fftrep_aux& Rb_aux, FFTRep& Ra)
{  tBluesteinFFT(x,n,root,powers,Rb,Ra); }

void BluesteinFFT(zz_pX& x, long n, const zz_p& root,
		  zz_pX& powers, Vec<mulmod_precon_t>& powers_aux, 
                  fftRep& Rb, fftrep_aux& Rb_aux, fftRep& Ra)
{  
FHE_TIMER_START
   tBluesteinFFT(x,n,root,powers,powers_aux,Rb,Rb_aux,Ra); 
FHE_TIMER_STOP
}

void DFT(ZZ_pX& x, const ZZ_pX& a, long n, const ZZ_p& root)
{ tDFT<ZZ_p,ZZ_pX>(x,a,n,root); }

void DFT(zz_pX& x, const zz_pX& a, long n, const zz_p& root)
{ tDFT<zz_p,zz_pX>(x,a,n,root); }
