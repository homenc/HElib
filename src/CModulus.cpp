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
/* CModulus.cpp - supports forward and backward length-m FFT transformations
 *
 * This is a wrapper around the bluesteinFFT routines, for one modulus q.
 *
 * On initialization, it initizlies NTL's zz_pContext for this q
 * and computes a 2m-th root of unity r mod q and also r^{-1} mod q.
 * Thereafter this class provides FFT and iFFT routines that converts between
 * time & frequency domains. Some tables are computed the first time that
 * each directions is called, which are then used in subsequent computations.
 * 
 * The "time domain" polynomials are represented as ZZX, which are reduced
 * modulo Phi_m(X). The "frequency domain" are jusr vectors of integers
 * (vec_long), that store only the evaluation in primitive m-th
 * roots of unity.
 */
#include "CModulus.h"
#include "timing.h"

// It is assumed that m,q,context, and root are already set. If root is set
// to zero, it will be computed by the compRoots() method. Then rInv is
// computed as the inverse of root.


zz_pContext BuildContext(long p, long maxroot) {
   if (maxroot <= CalcMaxRoot(p))
      return zz_pContext(INIT_USER_FFT, p);
   else
      return zz_pContext(p, maxroot);
}


// Constructor: it is assumed that zms is already set with m>1
// If q == 0, then the current context is used
Cmodulus::Cmodulus(const PAlgebra &zms, long qq, long rt)
{
  assert(zms.getM()>1);
  bool explicitModulus = true;

  if (qq == 0) {
    q = zz_p::modulus();
    explicitModulus = false;
  }
  else
    q = qq;

  qinv = PrepMulMod(q);

  zMStar = &zms;
  root = rt;

  long mm;
  mm = zms.getM();
  m_inv = InvMod(mm, q);

  zz_pBak bak; 

  if (zms.getPow2()) {
    // special case when m is a power of 2

    assert( explicitModulus );
    bak.save();

    RandomState state;  SetSeed(conv<ZZ>("84547180875373941534287406458029"));
    // DIRT: this ensures the roots are deterministically generated
    //    inside the zz_pContext constructor
    context = zz_pContext(INIT_USER_FFT, q);
    state.restore();

    context.restore();

    powers.set_ptr(new zz_pX);
    ipowers.set_ptr(new zz_pX);


    long k = zms.getPow2();
    long phim = 1L << (k-1); 

    assert(k <= zz_pInfo->MaxRoot); 
    // rootTables get initialized 0..zz_pInfo->Maxroot

#ifdef FHE_OPENCL
    altFFTInfo = MakeSmart<AltFFTPrimeInfo>();
    InitAltFFTPrimeInfo(*altFFTInfo, *zz_pInfo->p_info, k-1);
#endif

    long w0 = zz_pInfo->p_info->RootTable[0][k];
    long w1 = zz_pInfo->p_info->RootTable[1][k];

    powers->rep.SetLength(phim);
    powers_aux.SetLength(phim);
    for (long i = 0, w = 1; i < phim; i++) {
      powers->rep[i] = w;
      powers_aux[i] = PrepMulModPrecon(w, q);
      w = MulMod(w, w0, q);
    }

    ipowers->rep.SetLength(phim);
    ipowers_aux.SetLength(phim);
    for (long i = 0, w = 1; i < phim; i++) {
      ipowers->rep[i] = w;
      ipowers_aux[i] = PrepMulModPrecon(w, q);
      w = MulMod(w, w1, q);
    }

  
    return;
  }

  if (explicitModulus) {
    bak.save(); // backup the current modulus
    context = BuildContext(q, NextPowerOfTwo(zms.getM()) + 1);
    context.restore();       // set NTL's current modulus to q
  }
  else
    context.save();

  if (root==0) { // Find a 2m-th root of unity modulo q, if not given
    zz_p rtp;
    long e = 2*zms.getM();
    FindPrimitiveRoot(rtp,e); // NTL routine, relative to current modulus
    if (rtp==0) // sanity check
      Error("Cmod::compRoots(): no 2m'th roots of unity mod q");
    root = rep(rtp);
  }
  rInv = InvMod(root,q); // set rInv = root^{-1} mod q

  // Allocate memory (relative to current modulus that was defined above).
  // These objects will be initialized when anyone calls FFT/iFFT.

  zz_pX phimx_poly;
  conv(phimx_poly, zms.getPhimX());

  powers.set_ptr(new zz_pX);
  Rb.set_ptr(new fftRep);
  ipowers.set_ptr(new zz_pX);
  iRb.set_ptr(new fftRep);
  phimx.set_ptr(new zz_pXModulus1(zms.getM(), phimx_poly));

  BluesteinInit(mm, conv<zz_p>(root), *powers, powers_aux, *Rb);
  BluesteinInit(mm, conv<zz_p>(rInv), *ipowers, ipowers_aux, *iRb);
}

Cmodulus& Cmodulus::operator=(const Cmodulus &other)
{
  if (this == &other) return *this;

  zMStar  =  other.zMStar; // Yes, really copy this pointer
  q       = other.q;
  qinv    = other.qinv;
  m_inv   = other.m_inv;

  context = other.context;
  zz_pBak bak; bak.save(); // backup the current modulus
  context.restore();       // Set NTL's current modulus to q

  // NOTE: newer versions of NTL allow fftRep's and zz_pXModulus's to be copied
  // "out of context" (versions after 7.0.*). However, those copies
  // are not intended to allow copies out of one context into another,
  // so we still need to use copied_ptr's (but not context restoration).
  // All of this is fairly academic, as I don't think we really
  // copy FHEcontexts around anywhere. Also, it would be cleaner
  // to make the vector in FHEcontext be a vector of copied_ptr<Cmodulus>

  root = other.root;
  rInv = other.rInv;

  powers_aux = other.powers_aux;
  ipowers_aux = other.ipowers_aux;

  // copy data, not pointers in these fields
  powers = other.powers;
  Rb = other.Rb;
  ipowers = other.ipowers;
  iRb = other.iRb;
  phimx = other.phimx;

#ifdef FHE_OPENCL
  altFFTInfo = other.altFFTInfo;
#endif



  return *this;
}


void Cmodulus::FFT_aux(vec_long &y, zz_pX& tmp) const
{

  if (zMStar->getPow2()) {
    // special case when m is a power of 2

    long k = zMStar->getPow2();
    long phim = (1L << (k-1));
    long dx = deg(tmp);
    long p = zz_p::modulus();

    const zz_p *powers_p = (*powers).rep.elts();
    const mulmod_precon_t *powers_aux_p = powers_aux.elts();

    y.SetLength(phim);
    long *yp = y.elts();

    zz_p *tmp_p = tmp.rep.elts();

    for (long i = 0; i <= dx; i++)
      yp[i] = MulModPrecon(rep(tmp_p[i]), rep(powers_p[i]), p, powers_aux_p[i]);
    for (long i = dx+1; i < phim; i++)
      yp[i] = 0;

#ifdef FHE_OPENCL
    AltFFTFwd(yp, yp, k-1, *altFFTInfo);
#else
    FFTFwd(yp, yp, k-1, *zz_pInfo->p_info);
#endif

    return;
  }
    

  zz_p rt;
  conv(rt, root);  // convert root to zp format

  BluesteinFFT(tmp, getM(), rt, *powers, powers_aux, *Rb); // call the FFT routine

  // copy the result to the output vector y, keeping only the
  // entries corresponding to primitive roots of unity
  y.SetLength(zMStar->getPhiM());
  long i,j;
  long m = getM();
  for (i=j=0; i<m; i++)
    if (zMStar->inZmStar(i)) y[j++] = rep(coeff(tmp,i));
}

void Cmodulus::FFT(vec_long &y, const ZZX& x) const
{
  FHE_TIMER_START;
  zz_pBak bak; bak.save();
  context.restore();

  zz_pX& tmp = Cmodulus::getScratch_zz_pX();
  { FHE_NTIMER_START(FFT_remainder);
    convert(tmp,x);      // convert input to zpx format
  }

  FFT_aux(y, tmp);
};

void Cmodulus::FFT(vec_long &y, const zzX& x) const
{
  FHE_TIMER_START;
  zz_pBak bak; bak.save();
  context.restore();

  zz_pX& tmp = Cmodulus::getScratch_zz_pX();
  { FHE_NTIMER_START(FFT_remainder);
    convert(tmp,x);      // convert input to zpx format
  }

  FFT_aux(y, tmp);
};



void Cmodulus::iFFT(zz_pX &x, const vec_long& y)const
{
  FHE_TIMER_START;
  zz_pBak bak; bak.save();
  context.restore();

  if (zMStar->getPow2()) {
    // special case when m is a power of 2

    long k = zMStar->getPow2();
    long phim = (1L << (k-1));
    long p = zz_p::modulus();

    const zz_p *ipowers_p = (*ipowers).rep.elts();
    const mulmod_precon_t *ipowers_aux_p = ipowers_aux.elts();

    const long *yp = y.elts();

    vec_long& tmp = Cmodulus::getScratch_vec_long();
    tmp.SetLength(phim);
    long *tmp_p = tmp.elts();

#ifdef FHE_OPENCL
    AltFFTRev1(tmp_p, yp, k-1, *altFFTInfo);
#else
    FFTRev1(tmp_p, yp, k-1, *zz_pInfo->p_info);
#endif

    x.rep.SetLength(phim);
    zz_p *xp = x.rep.elts();

    for (long i = 0; i < phim; i++)
      xp[i].LoopHole() = MulModPrecon(tmp_p[i], rep(ipowers_p[i]), p, ipowers_aux_p[i]);


    x.normalize();

    return;
  }



  zz_p rt;
  long m = getM();

  // convert input to zpx format, initializing only the coeffs i s.t. (i,m)=1
  x.rep.SetLength(m);
  long i,j;
  for (i=j=0; i<m; i++)
    if (zMStar->inZmStar(i)) x.rep[i].LoopHole() = y[j++]; // DIRT: y[j] already reduced
  x.normalize();
  conv(rt, rInv);  // convert rInv to zp format

  BluesteinFFT(x, m, rt, *ipowers, ipowers_aux, *iRb); // call the FFT routine

  // reduce the result mod (Phi_m(X),q) and copy to the output polynomial x
  { FHE_NTIMER_START(iFFT_division);
    rem(x, x, *phimx); // out %= (Phi_m(X),q)
  }

  // normalize
  zz_p mm_inv;
  conv(mm_inv, m_inv);
  x *= mm_inv; 
}


zz_pX& Cmodulus::getScratch_zz_pX() 
{
   NTL_THREAD_LOCAL static zz_pX scratch;
   return scratch;
}

Vec<long>& Cmodulus::getScratch_vec_long()
{
   NTL_THREAD_LOCAL static Vec<long> scratch;
   return scratch;
}


fftRep& Cmodulus::getScratch_fftRep(long k)
{
  NTL_THREAD_LOCAL static fftRep rep;
  NTL_THREAD_LOCAL static long MaxK[4] = {-1, -1, -1, -1};

  long NumPrimes = zz_pInfo->NumPrimes;

  for (long i = 0; i < NumPrimes; i++) {
    if (k > MaxK[i]) {
      rep.tbl[i].SetLength(1L << k);
      MaxK[i] = k;
    }
  }

  rep.NumPrimes = NumPrimes;
  rep.k = rep.MaxK = k;

  return rep;
}




