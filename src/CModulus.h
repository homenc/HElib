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
#ifndef _CModulus_H_
#define _CModulus_H_
/**
 * @file CModulus.h
 * @brief Supports forward and backward length-m FFT transformations
 *
 * This is a wrapper around the bluesteinFFT routines, for one modulus q.
 **/
#include "NumbTh.h"
#include "PAlgebra.h"
#include "bluestein.h"
#include "cloned_ptr.h"



/**
* @class Cmodulus
*
* On initialization, it initizlies NTL's zz_pContext for this q
* and computes a 2m-th root of unity r mod q and also r^{-1} mod q.
* Thereafter this class provides FFT and iFFT routines that converts between
* time & frequency domains. Some tables are computed the first time that
* each directions is called, which are then used in subsequent computations.
* 
* The "time domain" polynomials are represented as ZZX, which are reduced
* modulo Phi_m(X). The "frequency domain" are just vectors of integers
* (vec_long), that store only the evaluation in primitive m-th
* roots of unity.
**/

class Cmodulus {
  long          q;       // the modulus
  zz_pContext   context; // NTL's tables for this modulus

  const PAlgebra* zMStar;  // points to the Zm* structure, m is FFT size

  long        m_inv;   // m^{-1} mod q

  long        root;    // 2m-th root of unity modulo q
  long        rInv;    // root^{-1} mod q

  copied_ptr<zz_pX>    powers;  // tables for forward FFT
  Vec<mulmod_precon_t> powers_aux;
  copied_ptr<fftRep>   Rb;

  copied_ptr<zz_pX>    ipowers; // tables for backward FFT
  Vec<mulmod_precon_t> ipowers_aux;
  copied_ptr<fftRep>   iRb;

  copied_ptr<zz_pXModulus1> phimx; // PhimX modulo q, for faster division w/ remainder


  // Allocate memory and compute roots
  void privateInit(const PAlgebra&, long rt);

 public:

  // Destructor and constructors

  // Default constructor
  Cmodulus() {}

  Cmodulus(const Cmodulus &other) { *this = other; }

  // Specify m and q, and optionally also the root
  // if q == 0, then the current context is used
  Cmodulus(const PAlgebra &zms, long qq, long rt);

  // Copy operator
  Cmodulus& operator=(const Cmodulus &other);

  // utility methods

  const PAlgebra &getZMStar() const { return *zMStar; }
  unsigned long getM() const    { return zMStar->getM(); }
  unsigned long getPhiM() const { return zMStar->getPhiM(); }
  long getQ() const          { return q; }
  long getRoot() const       { return root; }
  const zz_pXModulus1& getPhimX() const  { return *phimx; }

  //! @brief Restore NTL's current modulus
  void restoreModulus() const {context.restore();}

  // FFT routines

  // sets zp context internally
  void FFT(vec_long &y, const ZZX& x) const;  // y = FFT(x)

  // expects zp context to be set externally
  void iFFT(zz_pX &x, const vec_long& y) const; // x = FFT^{-1}(y)

  // returns thread-local scratch space
  // DIRT: this zz_pX is used for several zz_p moduli,
  // which is not officially sanctioned by NTL, but should be OK.
  static zz_pX& getScratch_zz_pX();

  static Vec<long>& getScratch_vec_long();

  // returns thread-local scratch space
  // DIRT: this use a couple of internal, undocumented
  // NTL interfaces
  static fftRep& getScratch_fftRep(long k);
};


#endif // ifdef _CModulus_H_
