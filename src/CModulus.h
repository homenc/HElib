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
 * Two classes are defined here, Cmodulus for a small moduli (long) and
 * CModulus for a large ones (ZZ). These classes are otherwise identical
 * hence they are implemented using a class template.
 **/
#include "PAlgebra.h"
#include "bluestein.h"
#include "cloned_ptr.h"
//NTL_CLIENT

/**
 * @class CMOD_zz_p
 * @brief typedefs for smallint Cmodulus
 **/
class CMOD_zz_p {
public:
  typedef long zz;
  typedef zz_p zp;
  typedef zz_pX zpx;
  typedef vec_long zzv;
  typedef fftRep fftrep;
  typedef zz_pContext zpContext;
  typedef zz_pBak zpBak;
  typedef zz_pXModulus zpxModulus;

};

/**
* @class CMOD_ZZ_p
* @brief typedefs for bigint CModulus
**/
class CMOD_ZZ_p {
public:
  typedef ZZ zz;
  typedef ZZ_p zp;
  typedef ZZ_pX zpx;
  typedef vec_ZZ zzv;
  typedef FFTRep fftrep;
  typedef ZZ_pContext zpContext ;
  typedef ZZ_pBak zpBak;
  typedef ZZ_pXModulus zpxModulus;
};

#define INJECT_TYPE(type,subtype) typedef typename type::subtype subtype


/**
* @class Cmod
* @brief template class for both bigint and smallint implementations
*
* This is a wrapper around the bluesteinFFT routines, for one modulus q.
* Two classes are defined here, Cmodulus for a small moduli (long) and
* CModulus for a large ones (ZZ). These classes are otherwise identical
* hence they are implemented using a class template.
*
* On initialization, it initizlies NTL's zz_pContext/ZZ_pContext for this q
* and computes a 2m-th root of unity r mod q and also r^{-1} mod q.
* Thereafter this class provides FFT and iFFT routines that converts between
* time & frequency domains. Some tables are computed the first time that
* each dierctions is called, which are then used in subsequent computations.
* 
* The "time domain" polynomials are represented as ZZX, whic are reduced
* modulo Phi_m(X). The "frequency domain" are jusr vectors of integers
* (vec_long or vec_ZZ), that store only the evaluation in primitive m-th
* roots of unity.
**/
template <class type>
class Cmod {
  INJECT_TYPE(type,zz);
  INJECT_TYPE(type,zp);
  INJECT_TYPE(type,zpx);
  INJECT_TYPE(type,zzv);
  INJECT_TYPE(type,fftrep);
  INJECT_TYPE(type,zpContext);
  INJECT_TYPE(type,zpBak);
  INJECT_TYPE(type,zpxModulus);

  zz          q;       // the modulus
  zpContext   context; // NTL's tables for this modulus

  const PAlgebra* zMStar;  // points to the Zm* structure, m is FFT size

  zz          m_inv;   // m^{-1} mod q

  zz          root;    // 2m-th root of unity modulo q
  zz          rInv;    // root^{-1} mod q

  zpx*        powers;  // tables for forward FFT
  mutable Vec<mulmod_precon_t> powers_aux;
  fftrep*     Rb;
  mutable fftrep_aux Rb_aux;
  fftrep*     Ra;

  zpx*        ipowers; // tables for backward FFT
  mutable Vec<mulmod_precon_t> ipowers_aux;
  fftrep*     iRb;
  mutable fftrep_aux iRb_aux;

  zpxModulus* phimx; // PhimX modulo q, for faster division w/ remainder
  zpx*        scratch; // temporary space, to satisfy NTL's rules

  // Allocate memory and compute roots
  void privateInit(const PAlgebra&, const zz& rt);

  void freeSpace() 
  {
    if (powers!=NULL)  { delete powers;  powers=NULL; }
    if (Rb!=NULL)      { delete Rb;      Rb=NULL; }
    if (Ra!=NULL)      { delete Ra;      Ra=NULL; }
    if (ipowers!=NULL) { delete ipowers; ipowers=NULL; }
    if (iRb!=NULL)     { delete iRb;     iRb = NULL; }
    if (phimx!=NULL)   { delete phimx;     phimx = NULL; }
    if (scratch!=NULL)   { delete scratch; scratch = NULL; }
  }

 public:

  // Destructor and constructors

  ~Cmod() { freeSpace(); } // destructor

  // Default constructor
  Cmod(): zMStar(NULL), powers(NULL), Rb(NULL), Ra(NULL), ipowers(NULL), iRb(NULL), 
          phimx(NULL), scratch(NULL) {}

  Cmod(const Cmod &other):
    zMStar(NULL),powers(NULL),Rb(NULL),Ra(NULL),ipowers(NULL),iRb(NULL),phimx(NULL),scratch(NULL)
  { *this = other; }

  // Specify m and q, and optionally also the root
  Cmod(const PAlgebra &zms, const zz &qq, const zz &rt);

  // Copy operator
  Cmod& operator=(const Cmod &other);

  // utility methods

  const PAlgebra &getZMStar() const { return *zMStar; }
  unsigned long getM() const    { return zMStar->getM(); }
  unsigned long getPhiM() const { return zMStar->getPhiM(); }
  const zz& getQ() const          { return q; }
  const zz& getRoot() const       { return root; }
  const zpxModulus& getPhimX() const  { return *phimx; }
  zpx& getScratch() const { return *scratch; }

  //! @brief Restore NTL's current modulus
  void restoreModulus() const {context.restore();}

  // FFT routines

  // sets zp context internally
  void FFT(zzv &y, const ZZX& x) const;  // y = FFT(x)

  // expects zp context to be set externally
  void iFFT(zpx &x, const zzv& y) const; // x = FFT^{-1}(y)
};

typedef Cmod<CMOD_zz_p> Cmodulus;
typedef Cmod<CMOD_ZZ_p> CModulus;

#endif // ifdef _CModulus_H_
