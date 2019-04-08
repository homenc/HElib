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
#ifndef HELIB_CMODULUS_H
#define HELIB_CMODULUS_H
/**
 * @file CModulus.h
 * @brief Supports forward and backward length-m FFT transformations
 *
 * This is a wrapper around the bluesteinFFT routines, for one modulus q.
 **/
#include "NumbTh.h"
#include "PAlgebra.h"
#include "bluestein.h"
#include "clonedPtr.h"

/**
* @class Cmodulus
* @brief Provides FFT and iFFT routines modulo a single-precision prime
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
  NTL::mulmod_t      qinv;    // PrepMulMod(q);

  NTL::zz_pContext   context; // NTL's tables for this modulus

  const PAlgebra* zMStar;  // points to the Zm* structure, m is FFT size

  long        m_inv;   // m^{-1} mod q

  long        root;    // 2m-th root of unity modulo q
  long        rInv;    // root^{-1} mod q

  copied_ptr<NTL::zz_pX>    powers;  // tables for forward FFT
  NTL::Vec<NTL::mulmod_precon_t> powers_aux;
  copied_ptr<NTL::fftRep>   Rb;

  copied_ptr<NTL::zz_pX>    ipowers; // tables for backward FFT
  NTL::Vec<NTL::mulmod_precon_t> ipowers_aux;
  copied_ptr<NTL::fftRep>   iRb;

  copied_ptr<zz_pXModulus1> phimx; // PhimX modulo q, for faster division w/ remainder


  // Allocate memory and compute roots
  void privateInit(const PAlgebra&, long rt);

 public:

#ifdef FHE_OPENCL
  SmartPtr<AltFFTPrimeInfo> altFFTInfo;
  // We need to allow copying...the underlying object
  // is immutable
#endif

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
  NTL::mulmod_t getQInv() const          { return qinv; }
  long getRoot() const       { return root; }
  const zz_pXModulus1& getPhimX() const  { return *phimx; }

  //! @brief Restore NTL's current modulus
  void restoreModulus() const {context.restore();}

  // FFT routines

  // sets zp context internally
  void FFT(NTL::vec_long &y, const NTL::ZZX& x) const;  // y = FFT(x)
  void FFT(NTL::vec_long &y, const zzX& x) const;  // y = FFT(x)

  // auxilliary routine used by above two routines
  void FFT_aux(NTL::vec_long &y, NTL::zz_pX& tmp) const;  



  // expects zp context to be set externally
  void iFFT(NTL::zz_pX &x, const NTL::vec_long& y) const; // x = FFT^{-1}(y)

  // returns thread-local scratch space
  // DIRT: this zz_pX is used for several zz_p moduli,
  // which is not officially sanctioned by NTL, but should be OK.
  static NTL::zz_pX& getScratch_zz_pX();

  static NTL::Vec<long>& getScratch_vec_long();

  // returns thread-local scratch space
  // DIRT: this use a couple of internal, undocumented
  // NTL interfaces
  static NTL::fftRep& getScratch_fftRep(long k);
};


#endif // ifndef HELIB_CMODULUS_H
