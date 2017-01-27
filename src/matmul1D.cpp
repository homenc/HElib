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
/* matmul1D.cpp - Data-movement operations on arrays of slots
 */
#include <algorithm>
#include <NTL/BasicThreadPool.h>
#include "matmul.h"

// Useful utility functions from matmul.cpp

// get data from cache, returns true if cache[i] is zero
extern bool
getDataFromCache(CachedConstants& cache, long i,
		 CachedConstants::CacheTag tag, const FHEcontext& context,
		 NTL::ZZX*& zzxPtr, DoubleCRT*& dcrtPtr);

// convert from ZZX to DoubleCRT cache
extern void CachedMatrixConvert(CachedDCRTPtxtMatrix& v, 
			 const CachedPtxtMatrix& w, const FHEcontext& context);

template<class type, class RX>
static bool processDiagonal(vector<RX>& diag, long dim, long i,
  vector<RX>& tmpDiag, const PlaintextMatrixInterface<type>& mat,
  const PAlgebra& zMStar, long extDeg, bool special)
{
  long D = tmpDiag.size();
  long nslots = diag.size();
  bool zDiag = true; // is this a zero diagonal?
  long nzLast = -1;  // index of last non-zero entry
  RX entry;

  // Process the entries in this diagonal one at a time
  for (long j = 0; j < D; j++) { // process entry j
    bool zEntry = mat.get(entry, mcMod(j-i, D), j); // entry [i,j-i mod D]
    assert(zEntry || deg(entry) < extDeg);
    // get(...) returns true if the entry is empty, false otherwise

    if (!zEntry && IsZero(entry)) zEntry = true; // zero is an empty entry too

    if (!zEntry) {   // not a zero entry
      zDiag = false; // mark diagonal as non-empty

      // clear entries between last nonzero entry and this one
      for (long jj = nzLast+1; jj < j; jj++) clear(tmpDiag[jj]);
      nzLast = j;

      tmpDiag[j] = entry;
    }
  }    
  if (zDiag) return true; // zero diagonal, nothing to do

  // clear trailing zero entries
  for (long jj = nzLast+1; jj < D; jj++) clear(tmpDiag[jj]);

  if (special) diag.assign(nslots, tmpDiag[0]); // order-1 dimension
  else for (long j = 0; j < nslots; j++)
    diag[j] = tmpDiag[ zMStar.coordinate(dim,j) ];
    // rearrange the indexes based on the current dimension

  return false; // a nonzero diagonal
}

// A class that implements the basic (sparse) 1D matrix-vector functions
template<class type> class mat_mul1D_impl{
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, Ctxt& ctxt,
    const PlaintextMatrixBaseInterface& mat, long dim) 
  {
    const PAlgebra& zMStar = ea.getContext().zMStar;

    assert(&ea == &mat.getEA().getDerived(type()));
    assert(&ea.getContext() == &ctxt.getContext());
    assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

    // special case for the extra dimension
    bool special = (dim == LONG(zMStar.numOfGens()));
    long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator

    RBak bak; bak.save(); ea.getTab().restoreContext(); // backup the NTL modulus

    // Get the derived type
    const PlaintextMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

    ctxt.cleanUp();  // not sure, but this may be a good idea
    Ctxt res(ZeroCtxtLike, ctxt); // fresh encryption of zero

    long nslots = ea.size();
    long d = ea.getDegree();
    vector<RX> diag, diag1;
    diag.resize(D);
    diag1.resize(nslots);
    ZZX cpoly;

    // Process the diagonals one at a time
    for (long i = 0; i < D; i++) { // process diagonal i
      if (processDiagonal(diag1, dim, i, diag, mat1, zMStar, d, special))
        continue; // zero diagonal

      // encode as a polynomial, then multiply and add
      ea.encode(cpoly, diag1);
      Ctxt shCtxt = ctxt;
      if (i != 0) ea.rotate1D(shCtxt, dim, i);   
      shCtxt.multByConstant(cpoly);
      res += shCtxt;
    }

    ctxt = res;
  }
};

void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
	       const PlaintextMatrixBaseInterface& mat, long dim) 
{
  FHE_TIMER_START;
  ea.dispatch<mat_mul1D_impl>(Fwd(ctxt), mat, dim);
}



// A class that implements the basic caching functionality for (sparse)
// 1D matrix-vector products
template<class type> class compMat1D_impl {
public:
  PA_INJECT(type)

  // This "apply" function only computes the cached 1D matrix cmat
  static void apply(const EncryptedArrayDerived<type>& ea, 
    CachedPtxtMatrix& cmat, const PlaintextMatrixBaseInterface& mat, long dim) 
  {
    const PAlgebra& zMStar = ea.getContext().zMStar;

    assert(&ea == &mat.getEA().getDerived(type()));
    assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

    // special case for the extra dimension
    bool special = (dim == LONG(zMStar.numOfGens()));
    long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator

    RBak bak; bak.save(); ea.getTab().restoreContext(); // backup the NTL modulus

    // Get the derived type
    const PlaintextMatrixInterface<type>& mat1 = 
      dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

    long nslots = ea.size();
    long d = ea.getDegree();
    vector<RX> diag, diag1;
    diag.resize(D);
    diag1.resize(nslots);
    cmat.SetLength(D);
    ZZX cpoly;

    // Process the diagonals one at a time
    for (long i = 0; i < D; i++) { // process diagonal i
      if (processDiagonal(diag1, dim, i, diag, mat1, zMStar, d, special))
        continue; // zero diagonal

      // encode as a polynomial, then multiply and add
      ea.encode(cpoly, diag1);
      cmat[i] = ZZXptr(new ZZX(cpoly));
    }
  }
};


void compMat1D(const EncryptedArray& ea, CachedPtxtMatrix& cmat, 
  const PlaintextMatrixBaseInterface& mat, long dim) 
{
  FHE_TIMER_START;
  ea.dispatch<compMat1D_impl>(Fwd(cmat), mat, dim);
}

void compMat1D(const EncryptedArray& ea, CachedDCRTPtxtMatrix& cmat, 
  const PlaintextMatrixBaseInterface& mat, long dim)
{
  FHE_TIMER_START;
  CachedPtxtMatrix zzxMat;
  compMat1D(ea, zzxMat, mat, dim);
  CachedMatrixConvert(cmat, zzxMat, ea.getContext());
}




#ifdef FHE_BOOT_THREADS

template<class CachedMatrix>
static void mat_mul1D_tmpl(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedMatrix& cmat, long dim)
{
  FHE_TIMER_START;
  const FHEcontext& context = ctxt.getContext();
  const PAlgebra& zMStar = context.zMStar;
  assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator

  ctxt.cleanUp();  // not sure, but this may be a good idea
  Ctxt res(ZeroCtxtLike, ctxt); // fresh encryption of zero

  Vec< shared_ptr<Ctxt> > tvec;
  tvec.SetLength(D);
  for (long i = 0; i < D; i++)
    tvec[i] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

  NTL_EXEC_RANGE(D, first, last)
      for (long i = first; i < last; i++) { // process diagonal i
        if (!cmat[i]) continue;      // zero diagonal
    
        // rotate and multiply
        (*tvec[i]) = ctxt;
        if (i != 0) ea.rotate1D(*tvec[i], dim, i);   
        tvec[i]->multByConstant(*cmat[i]);
      }
  NTL_EXEC_RANGE_END

  for (long i = 0; i < D; i++)
    res += *tvec[i];

  ctxt = res;
}

#else

template<class CachedMatrix>
static void mat_mul1D_tmpl(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedMatrix& cmat, long dim)
{
  FHE_TIMER_START;
  const FHEcontext& context = ctxt.getContext();
  const PAlgebra& zMStar = context.zMStar;
  assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator

  ctxt.cleanUp();  // not sure, but this may be a good idea
  Ctxt res(ZeroCtxtLike, ctxt); // fresh encryption of zero

  // Process the diagonals one at a time
  for (long i = 0; i < D; i++) { // process diagonal i
    if (!cmat[i]) continue;      // zero diagonal

    // rotate, multiply and add
    Ctxt shCtxt = ctxt;
    if (i != 0) ea.rotate1D(shCtxt, dim, i);   
    shCtxt.multByConstant(*cmat[i]);
    res += shCtxt;
  }
  ctxt = res;
}
#endif

void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedPtxtMatrix& cmat, long dim)
{ mat_mul1D_tmpl(ea, ctxt, cmat, dim); }

void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedDCRTPtxtMatrix& cmat, long dim)
{ mat_mul1D_tmpl(ea, ctxt, cmat, dim); }



/*********************************************************************/
/* mat_multi1D: Similar to mat_mul1D but different submatrices have
 * different transformations. 
 */

// This implementation uses one set of procedures to handle all of the
// caching options (none/ZZX/DCRT), at some point we should migrate
// all the stuff above to the same format.


// Process a single diagonal with index idx along dimenssion dim,
// making calls to the get(i,j,k) method and storing the result in the
// diag vector. If all the calls to get(i,j,k) return empty entries
// then return true and do not touch the diag vector. Otherwise, upon
// return diag[0...nSlots-1] contain the corresponding diagonal
// entries in all the dimension-dim matrices handled by mats
template<class type, class RX> static bool
processDiagonal(vector<RX>& diag, long dim, long idx,
                PlaintextMultiMatrixInterface<type>& mats,
                const PAlgebra& zMStar, long extDeg, bool special)
{
  long nSlots = diag.size();
  long nParts = mats.size();
  long d = special? 1 : zMStar.OrderOf(dim); // size of dimenssion dim

  bool zDiag = true; // is this a zero diagonal?
  long nzLast = -1;  // index of last non-zero entry
  RX entry;

  // Process the entries in this diagonal one at a time
  long blockIdx, innerIdx;
  for (long j=0; j<nSlots; j++) {
    if (special) {
      blockIdx=j; innerIdx = 0;
    } else {
      std::pair<long,long> idxes = zMStar.breakIndexByDim(j, dim);
      blockIdx = idxes.first;  // index into mats,
      innerIdx = idxes.second; // index along diemnssion dim
    }
    // process entry j
    bool zEntry = mats.get(entry, mcMod(innerIdx-idx, d), innerIdx, blockIdx);
    // entry [i,j-i mod d] in the block corresponding to blockIdx
    // get(...) returns true if the entry is empty, false otherwise

    // If non-zero, make sure the degree is not too large
    assert(zEntry || deg(entry) < extDeg);

    if (!zEntry && IsZero(entry)) zEntry = true; // zero is an empty entry too

    if (!zEntry) {   // not a zero entry
      zDiag = false; // mark diagonal as non-empty

      // clear entries between last nonzero entry and this one
      for (long jj = nzLast+1; jj < j; jj++) clear(diag[jj]);
      nzLast = j;
      diag[j] = entry;
    }
  }    
  if (zDiag) return true; // zero diagonal, nothing to do

  // clear trailing zero entries
  for (long jj = nzLast+1; jj < nSlots; jj++) clear(diag[jj]);

  return false; // a nonzero diagonal
}



// A class that implements the basic caching functionality for (sparse)
// 1D matrix-vector products
template<class type>
class mat_multi1D_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type> &ea,Ctxt &ctxt,long dim,
                    const PlaintextMatrixBaseInterface& mats,
		      CachedConstants::CacheTag tag)
  {
    const PAlgebra& zMStar = ea.getContext().zMStar;

    assert(&ea == &(mats.getEA().getDerived(type())));
    assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

    // special case for the extra dimension
    bool special = (dim == LONG(zMStar.numOfGens()));
    long D = special ? 1 : ea.sizeOfDimension(dim); // order of current dim
    long d = ea.getDegree();

    RBak bak; bak.save(); ea.getTab().restoreContext(); // backup NTL modulus

    // Get the derived type (DIRT: need to get rid of const modifier)
    PlaintextMultiMatrixInterface<type>& mats1 =
      (PlaintextMultiMatrixInterface<type>&)
      dynamic_cast< const PlaintextMultiMatrixInterface<type>& >( mats );

    CachedConstants& cache = mats1.getCache();
    std::vector<RX> diag(ea.size()); // scratch space
    ZZX cpoly;

    // It is assumed that a cache of the right size is ready to use
    bool cacheAvailable = (cache.size() == D);
    if (!cacheAvailable) {
      cache.clear(); // ensure that there is no stale cache
      if (tag == CachedConstants::tagZZX || tag == CachedConstants::tagDCRT)
	cache.resize(D);
    }

    // Process the diagonals one at a time. The code below will use only
    // rotate-by-one operations if this is a good dimenssion and the
    // matrix is dense, hence it may require fewer key-switching matrices

    long lastShift = 0;
    Ctxt shCtxt = ctxt;                  // temporary shifted ciphertext
    Ctxt acc = Ctxt(ZeroCtxtLike, ctxt); // initialize to zero

    for (long i = 0; i < D; i++) {       // process diagonal i
      bool zero = false;
      ZZX* zzxPtr = NULL;
      DoubleCRT *dcrtPtr = NULL;

      // If cached constant is available, use it
      if (cacheAvailable && !cache.isEmpty(i)) {
	zero = getDataFromCache(cache, i, tag,
				ctxt.getContext(), zzxPtr, dcrtPtr);
      }
      else { // no cached constant, need to compute it
        zero = processDiagonal(diag, dim, i, mats1, zMStar, d, special);
	if (!zero) {
          ea.encode(cpoly, diag); // encode as ZZX
	  if (tag == CachedConstants::tagDCRT)      // allocate a new DoubleCRT
	    dcrtPtr = new DoubleCRT(cpoly, ctxt.getContext());
	  else if (tag != CachedConstants::tagEmpty)// allocate a new ZZX
	    zzxPtr = new NTL::ZZX(cpoly);
	  else                                      // just use temporary ZZX
	    zzxPtr = &cpoly;
	}
      }

      // Depending on zero, zzxPtr, dcrtPtr, update the accumulated sum
      if (!zero) {
	if (i > 0) { // rotate the ciphertext
	  if (ea.nativeDimension(dim)) { // rotate the previous version
	    ea.rotate1D(ctxt, dim, i-lastShift);
	    shCtxt = ctxt;
	  } else {                       // rotate the original ciphertext
	    shCtxt = ctxt;
	    ea.rotate1D(shCtxt, dim, i);
	  }
	  lastShift = i;
	}
	if (dcrtPtr != NULL)
	  shCtxt.multByConstant(*dcrtPtr);
	else if (zzxPtr != NULL)
	  shCtxt.multByConstant(*zzxPtr);
	acc += shCtxt;
      }
      // The implementation above incurs an extra mult-by-constant due
      // to the masks in rotate1D when applied in a "bad dimension".
      // These masks can be folded into the constants here, but then
      // we would need to store two constants for each (e,f) rather
      // than one, namely const*mask and const*(1-mask).
      // We should implement that optimization at some point.

      // update the cache if needed
      if (tag != CachedConstants::tagEmpty
	                          && (!cacheAvailable || cache.isEmpty(i))) {
	if (zero)                 cache.setZero(i);
	else if (dcrtPtr != NULL) cache.setAt(i,dcrtPtr);
	else if (tag == CachedConstants::tagZZX && zzxPtr != NULL)
	  cache.setAt(i,zzxPtr);
      }
    }// end of this diagonal

    ctxt = acc; // return the result in ctxt

  }  // end of apply(...)
};   // end of class mat_multi1D_impl


//! @brief Multiply ctx by plaintext matrix
void mat_multi1D(Ctxt& ctxt, const EncryptedArray& ea, long dim,
                 const PlaintextMatrixBaseInterface& mats,
                 CachedConstants::CacheTag tag)
{
  FHE_TIMER_START;
  ea.dispatch<mat_multi1D_impl>(Fwd(ctxt), dim, mats, tag);
}
