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
/* blockMatmul.cpp - Data-movement operations on arrays of slots
 */
#include <algorithm>
#include <NTL/BasicThreadPool.h>
#include "matmul.h"

// A useful helper function to get information from cache
bool
getDataFromCache(CachedConstants& cache, long i,
		 CachedConstants::CacheTag tag, const FHEcontext& context,
		 NTL::ZZX*& zzxPtr, DoubleCRT*& dcrtPtr)
{
  if (cache.isZero(i)) return true; // zero constant

  if (cache.isDCRT(i)) dcrtPtr = cache.getDCRT(i);
  else if (cache.isZZX(i)) {
    zzxPtr = cache.getZZX(i);
    if (tag == CachedConstants::tagDCRT) { // upgrade cache to DoubleCRT
      // DIRT: this "upgrade" logic may not be thread safe
      dcrtPtr = new DoubleCRT(*zzxPtr, context);
      cache.setAt(i,dcrtPtr);
      zzxPtr = NULL;
    }
  }
  else throw std::logic_error("cached constant is NULL");
  return false;
}

// helper routines

void CachedMatrixConvert(CachedDCRTPtxtMatrix& v, 
			 const CachedPtxtMatrix& w, const FHEcontext& context)
{
  long n = w.length();
  v.SetLength(n);
  for (long i = 0; i < n; i++)
    if (w[i]) v[i] = DCRTptr(new DoubleCRT(*w[i], context));
    // DoubleCRT defined relative to all primes, even the "special" ones
}

// Utility function to ocnvert cache formats
void CachedBlockMatrixConvert(CachedDCRTPtxtBlockMatrix& v, 
     const CachedPtxtBlockMatrix& w, const FHEcontext& context)
{
  long m = w.NumRows();
  long n = w.NumCols();

  v.SetDims(m, n);
  for (long i = 0; i < m; i++)
    for (long j = 0; j < n; j++)
    if (w[i][j]) v[i][j] = DCRTptr(new DoubleCRT(*w[i][j], context));
    // DoubleCRT defined relative to all primes, even the "special" ones
    // NOTE: Fixes potential bug in original version
}


template<class type>
class blockMat_mul_impl {
public:
  PA_INJECT(type)

  // This code has a complexity of N+d (instead of N*d) where N is the
  // number of nonzero diagonal blocks. However, it requires space for
  // d extra ciphertexts
  static void apply(const EncryptedArrayDerived<type>& ea, Ctxt& ctxt, 
    const PlaintextBlockMatrixBaseInterface& mat) 
  {
    assert(&ea == &mat.getEA().getDerived(type()));
    assert(&ea.getContext() == &ctxt.getContext());

    const PAlgebra& zMStar = ea.getContext().zMStar;
    long p = zMStar.getP(); 
    long m = zMStar.getM();
    const RXModulus& F = ea.getTab().getPhimXMod();

    RBak bak; bak.save(); ea.getTab().restoreContext();

    const PlaintextBlockMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

    ctxt.cleanUp(); // not sure, but this may be a good idea

    long nslots = ea.size();
    long d = ea.getDegree();

    Vec< shared_ptr<Ctxt> > acc;
    acc.SetLength(d);
    for (long k = 0; k < d; k++)
      acc[k] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

    mat_R entry;
    entry.SetDims(d, d);

    vector<RX> entry1;
    entry1.resize(d);
    
    vector< vector<RX> > diag;
    diag.resize(nslots);
    for (long j = 0; j < nslots; j++) diag[j].resize(d);

    for (long i = 0; i < nslots; i++) { // process diagonal i
      bool zDiag = true;
      long nzLast = -1;

      for (long j = 0; j < nslots; j++) {
        bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j);
        assert(zEntry || (entry.NumRows() == d && entry.NumCols() == d));
        // get(...) returns true if the entry is empty, false otherwise

        if (!zEntry && IsZero(entry)) zEntry=true; // zero is an empty entry too

        if (!zEntry) {    // non-empty entry
          zDiag = false;  // mark diagonal as non-empty

          // clear entries between last nonzero entry and this one
          for (long jj = nzLast+1; jj < j; jj++) {
            for (long k = 0; k < d; k++)
              clear(diag[jj][k]);
          }
          nzLast = j;

          // recode entry as a vector of polynomials
          for (long k = 0; k < d; k++) conv(entry1[k], entry[k]);

          // compute the lin poly coeffs
          ea.buildLinPolyCoeffs(diag[j], entry1);
        }
      }
      if (zDiag) continue; // zero diagonal, continue

      // clear trailing zero entries    
      for (long jj = nzLast+1; jj < nslots; jj++) {
        for (long k = 0; k < d; k++)
          clear(diag[jj][k]);
      }

      // now diag[j] contains the lin poly coeffs

      Ctxt shCtxt = ctxt;
      ea.rotate(shCtxt, i); 
      shCtxt.cleanUp();

      RX cpoly1, cpoly2;
      ZZX cpoly;

      // apply the linearlized polynomial
      for (long k = 0; k < d; k++) {

        // compute the constant
        bool zConst = true;
        vector<RX> cvec;
        cvec.resize(nslots);
        for (long j = 0; j < nslots; j++) {
          cvec[j] = diag[j][k];
          if (!IsZero(cvec[j])) zConst = false;
        }

        if (zConst) continue;

        ea.encode(cpoly, cvec);
        conv(cpoly1, cpoly);

        // apply inverse automorphism to constant
        plaintextAutomorph(cpoly2, cpoly1, PowerMod(mcMod(p, m), mcMod(-k, d), m), zMStar.getM(), F);
        conv(cpoly, cpoly2);
        Ctxt shCtxt1 = shCtxt;
        shCtxt1.multByConstant(cpoly);
        *acc[k] += shCtxt1;
      }
    }

    Ctxt res(ZeroCtxtLike, ctxt);

    for (long k = 0; k < d; k++) {
      acc[k]->frobeniusAutomorph(k);
      res += *acc[k];
    }

    ctxt = res;
  }
};


void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const PlaintextBlockMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  ea.dispatch<blockMat_mul_impl>(Fwd(ctxt), mat);
}



// A class that implements the basic caching functionality for (sparse)
// matrix-vector products
template<class type>
class compBlockMat_impl{
public:
  PA_INJECT(type)

  // This "apply" function only computes the cached block matrix cmat
  static void apply(const EncryptedArrayDerived<type>& ea, 
    CachedPtxtBlockMatrix& cmat, const PlaintextBlockMatrixBaseInterface& mat) 
  {
    assert(&ea == &mat.getEA().getDerived(type()));
    const PAlgebra& zMStar = ea.getContext().zMStar;
    long p = zMStar.getP(); 
    long m = zMStar.getM();
    const RXModulus& F = ea.getTab().getPhimXMod();

    RBak bak; bak.save(); ea.getTab().restoreContext();

    // Get the derived type
    const PlaintextBlockMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

    long nslots = ea.size();
    long d = ea.getDegree();

    mat_R entry;
    entry.SetDims(d, d);
    vector<RX> entry1;
    entry1.resize(d);
    
    vector< vector<RX> > diag;
    diag.resize(nslots);
    for (long j = 0; j < nslots; j++) diag[j].resize(d);
    cmat.SetDims(nslots, d);

    for (long i = 0; i < nslots; i++) { // process diagonal i
      bool zDiag = true;
      long nzLast = -1;

      for (long j = 0; j < nslots; j++) {
        bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j);
        assert(zEntry || (entry.NumRows() == d && entry.NumCols() == d));
        // get(...) returns true if the entry is empty, false otherwise

        if (!zEntry && IsZero(entry)) zEntry=true; // zero is an empty entry too

        if (!zEntry) {    // non-empty entry
          zDiag = false;  // mark diagonal as non-empty

          // clear entries between last nonzero entry and this one
          for (long jj = nzLast+1; jj < j; jj++) {
            for (long k = 0; k < d; k++)
              clear(diag[jj][k]);
          }
          nzLast = j;

          // recode entry as a vector of polynomials
          for (long k = 0; k < d; k++) conv(entry1[k], entry[k]);

          // compute the lin poly coeffs
          ea.buildLinPolyCoeffs(diag[j], entry1);
        }
      }
      if (zDiag) continue; // zero diagonal, continue

      // clear trailing zero entries    
      for (long jj = nzLast+1; jj < nslots; jj++) {
        for (long k = 0; k < d; k++)
          clear(diag[jj][k]);
      }
      // now diag[j] contains the lin poly coeffs

      RX cpoly1, cpoly2;
      ZZX cpoly;

      // apply the linearlized polynomial
      for (long k = 0; k < d; k++) {

        // compute the constant
        bool zConst = true;
        vector<RX> cvec;
        cvec.resize(nslots);
        for (long j = 0; j < nslots; j++) {
          cvec[j] = diag[j][k];
          if (!IsZero(cvec[j])) zConst = false;
        }
        if (zConst) continue;

        ea.encode(cpoly, cvec);
        conv(cpoly1, cpoly);

        // apply inverse automorphism to constant
        plaintextAutomorph(cpoly2, cpoly1, PowerMod(mcMod(p, m), mcMod(-k, d), m), zMStar.getM(), F);
        conv(cpoly, cpoly2);
        cmat[i][k] = ZZXptr(new ZZX(cpoly));
      }
    }
  }
};

void compMat(const EncryptedArray& ea, CachedPtxtBlockMatrix& cmat, 
	     const PlaintextBlockMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  ea.dispatch<compBlockMat_impl>(Fwd(cmat), mat);
}

void compMat(const EncryptedArray& ea, CachedDCRTPtxtBlockMatrix& cmat, 
	     const PlaintextBlockMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  CachedPtxtBlockMatrix zzxMat;
  compMat(ea, zzxMat, mat);
  CachedBlockMatrixConvert(cmat, zzxMat, ea.getContext());
}





template<class CachedMatrix>
void blockMat_mul_tmpl(const EncryptedArray& ea, Ctxt& ctxt, 
		       const CachedMatrix& cmat)
{
  FHE_TIMER_START;
  ctxt.cleanUp(); // not sure, but this may be a good idea

  long nslots = ea.size();
  long d = ea.getDegree();

  Vec< shared_ptr<Ctxt> > acc;
  acc.SetLength(d);
  for (long k = 0; k < d; k++)
    acc[k] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

  for (long i = 0; i < nslots; i++) { // process diagonal i
    // apply the linearlized polynomial
    Ctxt shCtxt = ctxt;
    ea.rotate(shCtxt, i); 
    shCtxt.cleanUp();

    for (long k = 0; k < d; k++) {
      if (!cmat[i][k]) continue; // a zero constant
      Ctxt shCtxt1 = shCtxt;
      shCtxt1.multByConstant(*cmat[i][k]);
      *acc[k] += shCtxt1;
    }
  }

  Ctxt res(ZeroCtxtLike, ctxt);
  for (long k = 0; k < d; k++) {
    acc[k]->frobeniusAutomorph(k);
    res += *acc[k];
  }
  ctxt = res;
}
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedPtxtBlockMatrix& cmat)
{
  blockMat_mul_tmpl(ea, ctxt, cmat);
}
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedDCRTPtxtBlockMatrix& cmat)
{
  blockMat_mul_tmpl(ea, ctxt, cmat);
}



// Applying matmul to plaintext, useful for debugging
template<class type>
class blockMat_mul_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, 
    const PlaintextBlockMatrixBaseInterface& mat)
  {
    PA_BOILER

    const PlaintextBlockMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

    vector<RX> res;
    res.resize(n);
    for (long j = 0; j < n; j++) {
      vec_R acc, tmp, tmp1;
      mat_R val;

      acc.SetLength(d);
      for (long i = 0; i < n; i++) {
         if (!mat1.get(val, i, j)) {
            VectorCopy(tmp1, data[i], d);
            mul(tmp, tmp1, val);
            add(acc, acc, tmp);
         }
      }
      conv(res[j], acc);
    }

    data = res;
  }
}; 
void mat_mul(const EncryptedArray& ea, NewPlaintextArray& pa, 
	     const PlaintextBlockMatrixBaseInterface& mat)
{
  ea.dispatch<blockMat_mul_pa_impl>(Fwd(pa), mat); 
}
