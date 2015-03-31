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
 *
 *
 */
#ifndef _Bluestein
#define _Bluestein

/**
* @file bluestein.h
* @brief declaration of BluesteinFFT(x, n, root, powers, powers_aux, Rb):
*
* Compute length-n FFT of the coefficient-vector of x (in place) 
* If the degree of x is less than n then it treats the top coefficients
* as 0, if the degree of x is more than n then the extra coefficients are
* ignored. Similarly, if the top entries in x are zeros then x will have
* degree smaller than n. The argument root is a 2n-th root of unity, namely
* BluesteinFFT(...,root,...)=DFT(...,root^2,...).
*
* The inverse-FFT is obtained just by calling BluesteinFFT(... root^{-1}),
* but this procedure is *NOT SCALED*, so BluesteinFFT(x,n,root,...) and
* then BluesteinFFT(x,n,root^{-1},...) will result in x = n * x_original
*
* The values powers, powers_aux, and Rb must be precomputed by first
* calling BluesteinInit(n, root, powers, powers_aux, Rb).
*
**/


#include "NumbTh.h"



//! @brief initialize bluestein
void BluesteinInit(long n, const zz_p& root, zz_pX& powers, 
                   Vec<mulmod_precon_t>& powers_aux, fftRep& Rb);


//! @brief apply bluestein
void BluesteinFFT(zz_pX& x, long n, const zz_p& root, 
                  const zz_pX& powers, const Vec<mulmod_precon_t>& powers_aux, 
                  const fftRep& Rb);

#endif
