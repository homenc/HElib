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
 */

#include "bluestein.h"
#include "timing.h"
#include "CModulus.h"



void BluesteinInit(long n, const zz_p& root, zz_pX& powers, 
                   Vec<mulmod_precon_t>& powers_aux, fftRep& Rb)
{
  long p = zz_p::modulus();

  zz_p one; one=1;
  powers.SetMaxLength(n);

  SetCoeff(powers,0,one);
  for (long i=1; i<n; i++) {
    long iSqr = MulMod(i, i, 2*n); // i^2 mod 2n
    SetCoeff(powers,i, power(root,iSqr)); // powers[i] = root^{i^2}
  }

  // powers_aux tracks powers
  powers_aux.SetLength(n);
  for (long i = 0; i < n; i++)
    powers_aux[i] = PrepMulModPrecon(rep(powers[i]), p);


  long k = NextPowerOfTwo(2*n-1);
  long k2 = 1L << k; // k2 = 2^k

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

  TofftRep(Rb, b, k);
}

void BluesteinFFT(zz_pX& x, long n, const zz_p& root,
		  const zz_pX& powers, const Vec<mulmod_precon_t>& powers_aux, 
                  const fftRep& Rb)
{
  // FHE_TIMER_START;

  if (IsZero(x)) return;
  if (n<=0) {
    clear(x);
    return;
  }

  long p = zz_p::modulus();

  long dx = deg(x);
  for (long i=0; i<=dx; i++) {
    x[i].LoopHole() = MulModPrecon(rep(x[i]), rep(powers[i]), p, powers_aux[i]);
  }
  x.normalize();

  long k = NextPowerOfTwo(2*n-1);
  fftRep& Ra = Cmodulus::getScratch_fftRep(k);
  TofftRep(Ra, x, k);

  mul(Ra,Ra,Rb);           // multiply in FFT representation

  FromfftRep(x, Ra, n-1, 2*(n-1)); // then convert back
  dx = deg(x); 
  for (long i=0; i<=dx; i++) {
    x[i].LoopHole() = MulModPrecon(rep(x[i]), rep(powers[i]), p, powers_aux[i]);
  }
  x.normalize();
}




