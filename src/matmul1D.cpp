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

// A class that implements the basic (sparse) 1D matrix-vector functions
template<class type> class matmul1D_impl {
  PA_INJECT(type)
  const MatrixCacheType buildCache;
  std::unique_ptr<CachedzzxMatrix> zCache;
  std::unique_ptr<CachedDCRTMatrix> dCache;

  MatMul<type>& mat;
  const EncryptedArrayDerived<type>& ea;

public:
  matmul1D_impl(MatMulBase& _mat, MatrixCacheType tag)
    : buildCache(tag), mat(dynamic_cast< MatMul<type>& >(_mat)),
      ea(_mat.getEA().getDerived(type()))
  {
    if (buildCache==cachezzX)
      zCache.reset(new CachedzzxMatrix(NTL::INIT_SIZE,ea.size()));
    else if (buildCache==cacheDCRT)
      dCache.reset(new CachedDCRTMatrix(NTL::INIT_SIZE,ea.size()));
  }

  // Get the i'th diagonal along dimension dim, encoded as a
  // single constant. All blocks use the same transofmration.
  // Returns true if this is a zero diagonal, false otherwise
  bool processDiagonal1(zzX& cPoly, long dim, long i, long D)
  {  
    vector<RX> tmpDiag(D);
    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry
    RX entry;

    // Process the entries in this diagonal one at a time
    for (long j = 0; j < D; j++) { // process entry j
      bool zEntry = mat.get(entry, mcMod(j-i, D), j); // entry [j-i mod D, j]
      assert(zEntry || deg(entry) < ea.getDegree());
      // get(...) returns true if the entry is empty, false otherwise

      if (!zEntry && IsZero(entry)) zEntry = true;// zero is an empty entry too

      if (!zEntry) {   // not a zero entry
        zDiag = false; // mark diagonal as non-empty

        // clear entries between last nonzero entry and this one
        for (long jj = nzLast+1; jj < j; jj++) clear(tmpDiag[jj]);
        tmpDiag[j] = entry;
        nzLast = j;
      }
    }    
    if (zDiag) return true; // zero diagonal, nothing to do

    // clear trailing zero entries
    for (long jj = nzLast+1; jj < D; jj++) clear(tmpDiag[jj]);
    
    vector<RX> diag(ea.size());
    if (D==1) diag.assign(ea.size(), tmpDiag[0]); // dimension of size one
    else for (long j = 0; j < ea.size(); j++)
           diag[j] = tmpDiag[ ea.coordinate(dim,j) ];
           // rearrange the indexes based on the current dimension

    ea.encode(cPoly, diag);
    return false; // a nonzero diagonal
  }

  // Get the i'th diagonal along dimension dim, encoded as a
  // single constant. Different blocks use different transofmrations.
  // Returns true if this is a zero diagonal, false otherwise
  bool processDiagonal2(zzX& poly, long dim, long idx, long D)
  {
    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry
    RX entry;

    // Process the entries in this diagonal one at a time
    long blockIdx, innerIdx;
    vector<RX> diag(ea.size());
    for (long j=0; j<ea.size(); j++) {
      if (D==1) {
	blockIdx=j; innerIdx = 0;
      } else {
	std::tie(blockIdx, innerIdx) // std::pair<long,long> idxes
	  = ea.getContext().zMStar.breakIndexByDim(j, dim);
	//	blockIdx = idxes.first;  // which transformation
	//	innerIdx = idxes.second; // index along diemnssion dim
      }
      // process entry j
      bool zEntry=mat.multiGet(entry,mcMod(innerIdx-idx,D),innerIdx,blockIdx);
      // entry [i,j-i mod D] in the block corresponding to blockIdx
      // multiGet(...) returns true if the entry is empty, false otherwise

      // If non-zero, make sure the degree is not too large
      assert(zEntry || deg(entry) < ea.getDegree());

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
    for (long jj = nzLast+1; jj < ea.size(); jj++) clear(diag[jj]);
    ea.encode(poly, diag);
    return false; // a nonzero diagonal
  }

  void multiply(Ctxt* ctxt, long dim, bool oneTransform) 
  {
    assert(dim >= 0 && dim <= ea.dimension());
    RBak bak; bak.save(); ea.getTab().restoreContext(); // backup NTL modulus

    std::unique_ptr<Ctxt> res, shCtxt;
    if (ctxt!=nullptr) { // we need to do an actual multiplication
      ctxt->cleanUp(); // not sure, but this may be a good idea
      res.reset(new Ctxt(ZeroCtxtLike, *ctxt));
      shCtxt.reset(new Ctxt(*ctxt));
    }

    // Check if we have the relevant constant in cache
    CachedzzxMatrix* zcp;
    CachedDCRTMatrix* dcp;
    mat.getCache(&zcp, &dcp);

    // Process the diagonals one at a time
    zzX cpoly;
    long lastRotate = 0;
    long D = (dim==ea.dimension())? 1 : ea.sizeOfDimension(dim);
    for (long i = 0; i < D; i++) { // process diagonal i
      zzX* zxPtr=nullptr;
      DoubleCRT* dxPtr=nullptr;

      if (dcp != nullptr)         // DoubleCRT cache exists
	dxPtr = (*dcp)[i].get();
      else if (zcp != nullptr)    // zzx cache exists but no DoubleCRT
	zxPtr = (*zcp)[i].get();
      else { // no cache, compute const
	bool zero = oneTransform? processDiagonal1(cpoly, dim, i, D)
	                        : processDiagonal2(cpoly, dim, i, D);
        if (!zero) zxPtr = &cpoly; // if it is not a zero value, point to it
      }

      // if zero diagonal, nothing to do for this iteration
      if (zxPtr==nullptr && dxPtr==nullptr)
        continue;

      // Non-zero diagonal, store it in cache and/or multiply/add it

      if (ctxt!=nullptr && res!=nullptr) {
        // rotate by i, multiply by the polynomial, then add to the result
        if (i>0) {
          if (ea.nativeDimension(dim)) {
            ea.rotate1D(*ctxt, dim, i-lastRotate);
            *shCtxt = *ctxt;
          } else {
            *shCtxt = *ctxt;
            ea.rotate1D(*shCtxt, dim, i); // rotate by i
          }
          lastRotate = i;
	} // if i==0 we already have *shCtxt == *ctxt

        if (dxPtr!=nullptr) shCtxt->multByConstant(*dxPtr);
	else                shCtxt->multByConstant(*zxPtr);
	*res += *shCtxt;
      }
      if (buildCache==cachezzX) {
        (*zCache)[i].reset(new zzX(*zxPtr));
      }
      else if (buildCache==cacheDCRT) {
        (*dCache)[i].reset(new DoubleCRT(*zxPtr, ea.getContext()));
      }
    } // end of loop over diagonals

    if (ctxt!=nullptr && res!=nullptr) // copy result back to ctxt
      *ctxt = *res;

    // "install" the cache (if needed)
    if (buildCache == cachezzX)
      mat.installzzxcache(zCache);
    else if (buildCache == cacheDCRT)
      mat.installDCRTcache(dCache);
  } // end of multiply(...)
};

// Wrapper functions around the implemenmtation class
static void matmul1d(Ctxt* ctxt, MatMulBase& mat, long dim,
		     bool oneTransform, MatrixCacheType buildCache)
{
  MatMulLock locking(mat, buildCache);

  // If locking.getType()!=cacheEmpty then we really do need to
  // build the cache, and we also have the lock for it.

  if (locking.getType()==cacheEmpty && ctxt==nullptr) //  nothing to do
    return;

  switch (mat.getEA().getTag()) {
    case PA_GF2_tag: {
      matmul1D_impl<PA_GF2> M(mat, locking.getType());
      M.multiply(ctxt, dim, oneTransform);
      break;
    }
    case PA_zz_p_tag: {
      matmul1D_impl<PA_zz_p> M(mat, locking.getType());
      M.multiply(ctxt, dim, oneTransform);
      break;
    }
    default:
      throw std::logic_error("matmul1d: neither PA_GF2 nor PA_zz_p");
  }
}
void buildCache4MatMul1D(MatMulBase& mat, long dim, MatrixCacheType buildCache)
{ matmul1d(nullptr, mat, dim, true, buildCache); }

void matMul1D(Ctxt& ctxt, MatMulBase& mat,long dim, MatrixCacheType buildCache)
{ matmul1d(&ctxt, mat, dim, true, buildCache); }

void buildCache4MatMulti1D(MatMulBase& mat,long dim,MatrixCacheType buildCache)
{ matmul1d(nullptr, mat, dim, false, buildCache); }

void matMulti1D(Ctxt& ctxt,MatMulBase& mat,long dim,MatrixCacheType buildCache)
{ matmul1d(&ctxt, mat, dim, false, buildCache); }


/********************************************************************
 ********************************************************************/
// Applying matmul to plaintext, useful for debugging

template<class type>
class matmul1D_pa_impl {
public:
  PA_INJECT(type)

  static void multiply(NewPlaintextArray& pa, MatMul<type>& mat,
		       long dim, bool oneTrans)
  {
    const EncryptedArrayDerived<type>& ea = mat.getEA().getDerived(type());
    RBak bak; bak.save(); ea.getTab().restoreContext();

    long n = ea.size();
    long D = ea.sizeOfDimension(dim);

    vector< vector<RX> > data1(n/D);
    for (long k = 0; k < n/D; k++)
      data1[k].resize(D);

    // copy the data into a vector of 1D vectors
    vector<RX>& data = pa.getData<type>();
    for (long i = 0; i < n; i++) {
      long k,j;
      std::tie(k,j) = ea.getContext().zMStar.breakIndexByDim(i, dim);
      data1[k][j] = data[i];       // k= along dim, j = the rest of i
    }

    // multiply each one of the vectors by the same matrix
    for (long k = 0; k < n/D; k++) {
      for (long j = 0; j < D; j++) { // simple matrix-vector multiplication
	std::pair<long,long> p(k,j);
	long idx = ea.getContext().zMStar.assembleIndexByDim(p, dim);

	RX acc, val, tmp; 
	acc = 0;
        for (long i = 0; i < D; i++) {
          bool zero = oneTrans? mat.get(val, i, j) : mat.multiGet(val,i,j,k);
          if (!zero) {
            NTL::mul(tmp, data1[k][i], val);
            NTL::add(acc, acc, tmp);
          }
        }
        rem(data[idx], acc, ea.getG()); // store the result in the data array
      }
    }
  }
};
// A wrapper around the implementation class
static void matmul1d(NewPlaintextArray& pa, MatMulBase& mat,
		     long dim, bool oneTrans)
{
  switch (mat.getEA().getTag()) {
    case PA_GF2_tag: {
      matmul1D_pa_impl<PA_GF2>::multiply(pa,
                      dynamic_cast< MatMul<PA_GF2>& >(mat), dim, oneTrans);
      return;
    }
    case PA_zz_p_tag: {
      matmul1D_pa_impl<PA_zz_p>::multiply(pa,
                      dynamic_cast<MatMul<PA_zz_p>&>(mat), dim, oneTrans);
      return;
    }
    default:
      throw std::logic_error("matMul1D: neither PA_GF2 nor PA_zz_p");
  }
}

void matMul1D(NewPlaintextArray& pa, MatMulBase& mat, long dim)
{
  matmul1d(pa, mat, dim, true);
}
void matMulti1D(NewPlaintextArray& pa, MatMulBase& mat, long dim)
{
  matmul1d(pa, mat, dim, false);
}
