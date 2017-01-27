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
/* matmul.cpp - Data-movement operations on arrays of slots
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


template<class type> class mat_mul_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, 
    Ctxt& ctxt, const PlaintextMatrixBaseInterface& mat) 
  {
    assert(&ea == &mat.getEA().getDerived(type()));
    assert(&ea.getContext() == &ctxt.getContext());

    RBak bak; bak.save(); ea.getTab().restoreContext();

    // Get the derived type
    const PlaintextMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

    ctxt.cleanUp(); // not sure, but this may be a good idea

    Ctxt res(ctxt.getPubKey(), ctxt.getPtxtSpace()); // fresh encryption of zero

    long nslots = ea.size();
    long d = ea.getDegree();

    RX entry;
    vector<RX> diag;
    diag.resize(nslots);

    // Process the diagonals one at a time
    for (long i = 0; i < nslots; i++) {  // process diagonal i
      bool zDiag = true; // is this a zero diagonal?
      long nzLast = -1;  // index of last non-zero entry on this diagonal

      // Compute constants for each entry on this diagonal
      for (long j = 0; j < nslots; j++) { // process entry j
        bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j); // callback
        assert(zEntry || deg(entry) < d);

        if (!zEntry && IsZero(entry)) zEntry = true; // check for zero

        if (!zEntry) { // non-zero diagonal entry

          zDiag = false; // diagonal is non-zero

          // clear entries between last nonzero entry and this one
          for (long jj = nzLast+1; jj < j; jj++) clear(diag[jj]);
          nzLast = j;

          diag[j] = entry;
        }
      }
      
      if (zDiag) continue; // zero diagonal, continue

      // clear trailing zero entries
      for (long jj = nzLast+1; jj < nslots; jj++) clear(diag[jj]);

      // Now we have the constants for all the diagonal entries, encode the
      // diagonal as a single polynomial with these constants in the slots
      ZZX cpoly;
      ea.encode(cpoly, diag);

      // rotate by i, multiply by the polynomial, then add to the result
      Ctxt shCtxt = ctxt;
      ea.rotate(shCtxt, i); // rotate by i
      shCtxt.multByConstant(cpoly);
      res += shCtxt;
    }
    ctxt = res;
  }
};
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
	     const PlaintextMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  ea.dispatch<mat_mul_impl>(Fwd(ctxt), mat);
}




template<class type> class mat_mul_dense_impl {
public:
  PA_INJECT(type)

  // A recursive matrix-by-vector multiply, used by the dense matrix code.
  // This routine is optimized to use only the rotate1D routine rather
  // than the more expensive linear-array rotations.
  static void rec_mul(const EncryptedArrayDerived<type>& ea,
    long dim, Ctxt& res, const Ctxt& pdata, const vector<long>& idx,
    const PlaintextMatrixInterface<type>& mat,
    const vector<long>& dimx) 
  {
    long ndims = ea.dimension();
    long nslots = ea.size();

    if (dim >= ndims) { // Last dimension (recursion edge condition)
      vector<RX> pmat;  // the plaintext diagonal
      pmat.resize(nslots);
      bool zDiag = true; // is this a zero diagonal
      for (long j = 0; j < nslots; j++) {
        long i = idx[j];
        RX val;
        if (mat.get(val, i, j)) // returns true if the entry is zero
          clear(pmat[j]);
        else {           // not a zero entry
          pmat[j] = val;
          zDiag = false; // not a zero diagonal
        }
      }
      if (zDiag) return; // zero diagonal, nothing to do

      // Now we have the constants for all the diagonal entries, encode the
      // diagonal as a single polynomial with these constants in the slots
      ZZX epmat;
      ea.encode(epmat, pmat);

      // multiply by the polynomial, then add to the result
      Ctxt tmp = pdata;
      tmp.multByConstant(epmat);
      res += tmp;
    }
    else { // not the last dimension, make a recursive call
      long sdim = ea.sizeOfDimension(dimx[dim]);

      // compute "in spirit" sum_i (pdata >> i) * i'th-diagonal, but
      // adjust the indexes so that we only need to rotate the cipehrtext
      // along the different dimensions separately
      for (long offset = 0; offset < sdim; offset++) {
        Ctxt pdata1 = pdata;
        vector<long> idx1;
        ea.rotate1D(pdata1, dimx[dim], offset);
        ea.EncryptedArrayBase::rotate1D(idx1, idx, dimx[dim], offset);
        // indexes adjusted, make the recursive call
        rec_mul(ea, dim+1, res, pdata1, idx1, mat, dimx);
      }
    }
  }

  // helper class to sort dimensions, so that
  //    - bad dimensions come before good dimensions (primary sort key)
  //    - small dimensions come before large dimesnions (secondary sort key)
  // this is a good order to process the dimensions in the recursive mat_mul_dense
  // routine: it ensures that the work done at the leaves of the recursion is
  // minimized, and that the work done at the non-leaves is dominated by the
  // work done at the leaves.

  struct MatMulDimComp {
    const EncryptedArrayDerived<type> *ea;
    MatMulDimComp(const EncryptedArrayDerived<type> *_ea) : ea(_ea) {}

    bool operator()(long i, long j) { 
      return (!ea->nativeDimension(i) && ea->nativeDimension(j)) ||
             (  (ea->nativeDimension(i) == ea->nativeDimension(j)) &&
                (ea->sizeOfDimension(i) < ea->sizeOfDimension(j))  );
    }
  };

  // Multiply a ciphertext vector by a plaintext dense matrix
  static void apply(const EncryptedArrayDerived<type>& ea, Ctxt& ctxt, 
    const PlaintextMatrixBaseInterface& mat) 
  {
    assert(&ea == &mat.getEA().getDerived(type()));
    assert(&ea.getContext() == &ctxt.getContext());

    RBak bak; bak.save(); ea.getTab().restoreContext();

    // Get the derived type
    const PlaintextMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

    ctxt.cleanUp(); // not sure, but this may be a good idea

    Ctxt res(ctxt.getPubKey(), ctxt.getPtxtSpace()); // fresh encryption of zero
    vector<long> idx;
    idx.resize(ea.size());
    for (long i = 0; i < ea.size(); i++)
       idx[i] = i;

    vector<long> dimx;
    dimx.resize(ea.dimension());
    for (long i = 0; i < ea.dimension(); i++)
      dimx[i] = i;

    sort(dimx.begin(), dimx.end(), MatMulDimComp(&ea));
    // sort the dimenesions so that bad ones come before good,
    // and then small ones come before large

    // call the recursive procedure to do the actual work
    rec_mul(ea, 0, res, ctxt, idx, mat1, dimx);

    ctxt = res; // copy the result back to ctxt
  }
};
void mat_mul_dense(const EncryptedArray& ea, Ctxt& ctxt, 
		   const PlaintextMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  ea.dispatch<mat_mul_dense_impl>(Fwd(ctxt), mat);
}



// A class that implements the basic caching functionality for (sparse)
// matrix-vector products
template<class type> class compMat_impl {
public:
  PA_INJECT(type)

  // This "apply" function only computes the cached matrix cmat
  static void apply(const EncryptedArrayDerived<type>& ea, 
    CachedPtxtMatrix& cmat, const PlaintextMatrixBaseInterface& mat) 
  {
    assert(&ea == &mat.getEA().getDerived(type()));

    RBak bak; bak.save(); ea.getTab().restoreContext();

    // Get the derived type
    const PlaintextMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

    long nslots = ea.size();
    long d = ea.getDegree();
    RX entry;
    vector<RX> diag;
    diag.resize(nslots);
    cmat.SetLength(nslots);

    // Process the diagonals one at a time
    for (long i = 0; i < nslots; i++) {  // process diagonal i
      bool zDiag = true; // is this a zero diagonal?
      long nzLast = -1;  // index of last non-zero entry on this diagonal

      // Compute constants for each entry on this diagonal
      for (long j = 0; j < nslots; j++) { // process entry j
        bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j); // callback
        assert(zEntry || deg(entry) < d);

        if (!zEntry && IsZero(entry)) zEntry = true; // check for zero

        if (!zEntry) { // non-zero diagonal entry

          zDiag = false; // diagonal is non-zero

          // clear entries between last nonzero entry and this one
          for (long jj = nzLast+1; jj < j; jj++) clear(diag[jj]);
          nzLast = j;

          diag[j] = entry;
        }
      }
      
      if (zDiag) continue; // zero diagonal, continue

      // clear trailing zero entries
      for (long jj = nzLast+1; jj < nslots; jj++) clear(diag[jj]);

      // Now we have the constants for all the diagonal entries, encode the
      // diagonal as a single polynomial with these constants in the slots
      ZZX cpoly;
      ea.encode(cpoly, diag);
      cmat[i] = ZZXptr(new ZZX(cpoly));
    }
  }
};


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

void compMat(const EncryptedArray& ea, CachedPtxtMatrix& cmat, 
	     const PlaintextMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  ea.dispatch<compMat_impl>(Fwd(cmat), mat);
}

void compMat(const EncryptedArray& ea, CachedDCRTPtxtMatrix& cmat, 
	     const PlaintextMatrixBaseInterface& mat)
{
  FHE_TIMER_START;
  CachedPtxtMatrix zzxMat;
  compMat(ea, zzxMat, mat);
  CachedMatrixConvert(cmat, zzxMat, ea.getContext());
}



template<class CachedMatrix>
void mat_mul_tmpl(const EncryptedArray& ea, Ctxt& ctxt,
		  const CachedMatrix& cmat)
{
  FHE_TIMER_START;
  ctxt.cleanUp(); // not sure, but this may be a good idea
  Ctxt res(ctxt.getPubKey(), ctxt.getPtxtSpace()); // fresh encryption of zero

  // Process the diagonals one at a time
  long nslots = ea.size();
  for (long i = 0; i < nslots; i++) {  // process diagonal i
    if (!cmat[i]) continue; // a zero diagonal

    // rotate by i, multiply, and add to the result
    Ctxt shCtxt = ctxt;
    ea.rotate(shCtxt, i); // rotate by i
    shCtxt.multByConstant(*cmat[i]);
    res += shCtxt;
  }
  ctxt = res;
}
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedPtxtMatrix& cmat)
{
  mat_mul_tmpl(ea, ctxt, cmat);
}
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedDCRTPtxtMatrix& cmat)
{
  mat_mul_tmpl(ea, ctxt, cmat);
}


// Applying matmul to plaintext, useful for debugging
template<class type>
class mat_mul_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, 
    const PlaintextMatrixBaseInterface& mat)
  {
    PA_BOILER

    const PlaintextMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

    vector<RX> res;
    res.resize(n);
    for (long j = 0; j < n; j++) {
      RX acc, val, tmp; 
      acc = 0;
      for (long i = 0; i < n; i++) {
        if (!mat1.get(val, i, j)) {
          NTL::mul(tmp, data[i], val);
          NTL::add(acc, acc, tmp);
        }
      }
      rem(acc, acc, G);
      res[j] = acc;
    }

    data = res;
  }
}; 
void mat_mul(const EncryptedArray& ea, NewPlaintextArray& pa, 
	     const PlaintextMatrixBaseInterface& mat)
{
  ea.dispatch<mat_mul_pa_impl>(Fwd(pa), mat); 
}
