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
/* CModulus.cpp - supports forward and backward length-m FFT transformations
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
 */
#include <cassert>

#include "NumbTh.h"
#include "CModulus.h"
#include "timing.h"

// Some simple functions that should have been provided by NTL but are not
inline bool IsZero(long i) { return (i==0); }
inline void conv(NTL::vec_zz_p& to, NTL::vec_long& from)
{
  to.SetLength(from.length());
  for (long i=0; i<from.length(); i++) conv(to[i], from[i]);
  // It is assumed that the NTL modulus is already set
}
inline void conv(NTL::vec_long& to, NTL::vec_zz_p& from)
{
  to.SetLength(from.length());
  for (long i=0; i<from.length(); i++) to[i]=rep(from[i]);
  // It is assumed that the NTL modulus is already set
}

inline void SetCoeff(NTL::ZZ_pX& poly, long idx, const NTL::ZZ& val)
{ SetCoeff(poly, idx, to_ZZ_p(val)); }

// It is assumed that m,q,context, and root are already set. If root is set
// to zero, it will be computed by the compRoots() method. Then rInv is
// computed as the inverse of root.

//void Cmod<zz,zp,zpx,zzv,fftrep,zpContext,zpBak>::privateInit(const PAlgebra& zms, const zz& rt)


ZZ_pContext BuildContext(const ZZ& p, long maxroot)
  { return ZZ_pContext(p); }

zz_pContext BuildContext(long p, long maxroot)
  { return zz_pContext(p, maxroot); }


// Constructor: it is assumed that zms is already set with m>1
template <class type> Cmod<type>::
Cmod(const PAlgebra &zms, const zz &qq, const zz &rt)
{
  assert(zms.getM()>1);

  zMStar = &zms;
  q = qq;
  root = rt;

  zz mm;
  mm = zms.getM();
  m_inv = InvMod(mm, q);

  zz_pBak bak; bak.save(); // backup the current modulus
  context = BuildContext(qq, NextPowerOfTwo(zms.getM()) + 1);
  context.restore();       // set NTL's current modulus to q

  if (IsZero(root)) { // Find a 2m-th root of unity modulo q, if not given
    zp rtp;
    long e = 2*zms.getM();
    FindPrimitiveRoot(rtp,e); // NTL routine, relative to current modulus
    if (IsZero(rtp)) // sanity check
      Error("Cmod::compRoots(): no 2m'th roots of unity mod q");
    root = rep(rtp);
  }
  rInv = InvMod(root,q); // set rInv = root^{-1} mod q

  // Allocate memory (relative to current modulus that was defined above).
  // These objects will be initialized when anyone calls FFT/iFFT.

  zpx phimx_poly;
  conv(phimx_poly, zms.getPhimX());

  powers  = new zpx();
  Rb      = new fftrep();
  Ra      = new fftrep();
  ipowers = new zpx();
  iRb     = new fftrep();
  phimx   = new zpxModulus(phimx_poly);
  scratch = new zpx();
}

template <class type>
Cmod<type>& Cmod<type>::operator=(const Cmod &other)
{
  if (this == &other) return *this;

  zMStar  =  other.zMStar; // Yes, really copy this pointer
  q       = other.q;
  m_inv   = other.m_inv;

  context = other.context;
  zz_pBak bak; bak.save(); // backup the current modulus
  context.restore();       // Set NTL's current modulus to q

  root = other.root;
  rInv = other.rInv;

  powers_aux = other.powers_aux;
  ipowers_aux = other.ipowers_aux;
  Rb_aux = other.Rb_aux;
  iRb_aux = other.iRb_aux;

  // copy data, not pointers in these fields
  freeSpace(); // just in case
  if (other.powers)  powers  = new zpx(*(other.powers));
  if (other.Rb)      Rb      = new fftrep(*(other.Rb));
  if (other.Ra)      Ra      = new fftrep(*(other.Ra));
  if (other.ipowers) ipowers = new zpx(*(other.ipowers));
  if (other.iRb)     iRb     = new fftrep(*(other.iRb));
  if (other.phimx)   phimx   = new zpxModulus(*(other.phimx));
  if (other.scratch) scratch   = new zpx(*(other.scratch));

  return *this;
}

template <class type>
void Cmod<type>::FFT(zzv &y, const ZZX& x) const
{
  FHE_TIMER_START;
  zpBak bak; bak.save();
  context.restore();
  zp rt;
  zpx& tmp = getScratch();

  conv(tmp,x);      // convert input to zpx format
  conv(rt, root);  // convert root to zp format

  BluesteinFFT(tmp, getM(), rt, *powers, powers_aux, *Rb, Rb_aux, *Ra); // call the FFT routine

  // copy the result to the output vector y, keeping only the
  // entries corresponding to primitive roots of unity
  y.SetLength(zMStar->getPhiM());
  long i,j;
  long m = getM();
  for (i=j=0; i<m; i++)
    if (zMStar->inZmStar(i)) y[j++] = rep(coeff(tmp,i));
  FHE_TIMER_STOP;
}

template <class type>
void Cmod<type>::iFFT(zpx &x, const zzv& y)const
{
  FHE_TIMER_START;
  zpBak bak; bak.save();
  context.restore();
  zp rt;

  long m = getM();

  // convert input to zpx format, initializing only the coeffs i s.t. (i,m)=1
  x.SetMaxLength(m);
  long i,j;
  for (i=j=0; i<m; i++)
    if (zMStar->inZmStar(i)) SetCoeff(x, i, y[j++]);
  x.normalize();
  conv(rt, rInv);  // convert rInv to zp format

  BluesteinFFT(x, m, rt, *ipowers, ipowers_aux, *iRb, iRb_aux, *Ra); // call the FFT routine

  // reduce the result mod (Phi_m(X),q) and copy to the output polynomial x
FHE_NTIMER_START("iFFT:division")
  rem(x, x, *phimx); // out %= (Phi_m(X),q)
FHE_NTIMER_STOP("iFFT:division")

  // normalize
  zp mm_inv;
  conv(mm_inv, m_inv);
  x *= mm_inv; 


  FHE_TIMER_STOP;
}

// instantiating the template classes
template class Cmod<CMOD_zz_p>; // small q
template class Cmod<CMOD_ZZ_p>; // large q

