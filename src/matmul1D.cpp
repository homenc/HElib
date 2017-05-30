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
/* matmul1D.cpp - Data-movement operations on arrays of slots
 */
#include <algorithm>
#include <NTL/BasicThreadPool.h>
#include "EncryptedArray.h"
#include "matmul.h"

struct PtxtPtr {
  enum Type { ZERO=0, ZZX=1, DCRT=2 };
  zzX*       zp;
  DoubleCRT* dp;

  PtxtPtr(): zp(nullptr), dp(nullptr) {}
};

// A class that implements the basic (sparse) 1D matrix-vector functions
template<class type> class matmul1D_impl {
  PA_INJECT(type)
  const MatrixCacheType buildCache;
  std::unique_ptr<CachedzzxMatrix> zCache;
  std::unique_ptr<CachedDCRTMatrix> dCache;

  MatMul<type>& mat;
  const EncryptedArrayDerived<type>& ea;

public:
  matmul1D_impl(MatMulBase& _mat, long dim, MatrixCacheType tag)
    : buildCache(tag), mat(dynamic_cast< MatMul<type>& >(_mat)),
      ea(_mat.getEA().getDerived(type()))
  {
    if (buildCache==cachezzX) {
      zCache.reset(new CachedzzxMatrix(NTL::INIT_SIZE,ea.sizeOfDimension(dim)));
    } else if (buildCache==cacheDCRT) {
      dCache.reset(new CachedDCRTMatrix(NTL::INIT_SIZE,ea.sizeOfDimension(dim)));
    }
  }

  // Get the i'th diagonal along dimension dim, encoded as a
  // single constant. All blocks use the same transofmration.
  // Returns true if this is a zero diagonal, false otherwise
  bool processDiagonal1(zzX& cPoly, long dim, long i, long D, long rotAmt)
  {
    vector<RX> tmpDiag(D);
    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry
    RX entry;

    // Process the entries in this diagonal one at a time
    for (long j = 0; j < D; j++) { // process entry j
      long rotJ = (j+rotAmt) % D;  // need to rotate constant by rotAmt
      bool zEntry = mat.get(entry, mcMod(rotJ-i, D), rotJ); // entry [j-i mod D, j]
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
  bool processDiagonal2(zzX& poly, long dim, long idx, long D, long rotAmt)
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
        innerIdx = (innerIdx+rotAmt) % D;  // need to rotate constant by rotAmt
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

  // Returns in vec all the constants needed for a single giant-step.
  // The return value is ZERO(=0) if all these constants are zero.
  // Otherwise it is ZZX(=1) if using zzX or DCRT(=2) if using DoubleCRT
  PtxtPtr::Type getConsts(std::vector<PtxtPtr>& vec, std::vector<zzX>& polys,
                          CachedzzxMatrix* zcp, CachedDCRTMatrix* dcp,
                          long dim, long jg, long D, bool oneTrans)
  {
    PtxtPtr::Type typ = PtxtPtr::ZERO;
    long g = polys.size(); // how many constants do we need
    vec.assign(g,PtxtPtr());

    if (dcp!=nullptr) {        // we have DCRT cache
      for (long i=0, idx=jg; i<g && idx<D; i++, idx++) {
        if ((vec[i].dp = (*dcp)[idx].get()) != nullptr)
          typ = PtxtPtr::DCRT;
      }
    } else if (zcp!=nullptr) { // we have ZZX cache
      for (long i=0, idx=jg; i<g && idx<D; i++, idx++) {
        if ((vec[i].zp = (*zcp)[idx].get()) != nullptr)
          typ = PtxtPtr::ZZX;
      }
    } else {                   // no cache, compute consts
      for (long i=0, idx=jg; i<g && idx<D; i++, idx++) {
        bool zero = oneTrans? processDiagonal1(polys[i], dim, idx, D, i)
                            : processDiagonal2(polys[i], dim, idx, D, i);
        if (!zero) {
          vec[i].zp = &(polys[i]); // if not a zero, point to it
          typ = PtxtPtr::ZZX;
        }
      }
    }
    return typ;
  }

  void multiply(Ctxt* ctxt, long dim, bool oneTransform) 
  {
    assert(dim >= 0 && dim <= ea.dimension());
    RBak bak; bak.save(); ea.getTab().restoreContext(); // backup NTL modulus

    long D = (dim==ea.dimension())? 1 : ea.sizeOfDimension(dim);
    long g = mat.getGstep();
    if (g<1 || g>=D) g=1; // sanity check
    long dDivg = divc(D,g);

    // Process the diagonals in baby-step/giant-step ordering.
    //   sum_{i=0}^{d-1} const_i rot^i(X)
    //   = \sum_{i=0}^{g-1} \sum_{j=0}^{d/g -1} const_{i+g*j} rot^{i+g*j}(X)
    //   = \sum_{i=0}^{g-1} rot^i(sum_j rot^{-i}(const_{i+g*j}) rot^{g*j}(X))
    //
    // so for i=0..g-1 we let
    //    Y_i = sum_j rot^{-i}(const_{i+g*j}) rot^{g*j}(X)
    // then compute \sum_{i=0}^{g-1} rot^i(Y_i).
    //
    // Computing the Y_i's, we initialize an accumulator for each Y_i,
    // then compute the rotations X_j = rot^{g*j}(X), j=0,...,d/g-1.
    // Each X_j is multiplied by all the constants rot^{-i}(const_{i+g*j}),
    // i=0,...,g-1, and the i'th product is added to the accumulator for
    // the corresponding Y_i.

    long lastRotate = 0;
    std::vector<Ctxt> acc; // accumulators
    if (ctxt!=nullptr) {   // we need to do an actual multiplication
      ctxt->cleanUp();     // not sure, but this may be a good idea
      Ctxt tmp(ZeroCtxtLike, *ctxt);
      acc.resize(g, tmp);
    }

    // Check if we have the relevant constant in cache
    CachedzzxMatrix* zcp;
    CachedDCRTMatrix* dcp;
    mat.getCache(&zcp, &dcp);

    // Process the diagonals in giant-step/baby-step order
    std::vector<zzX> cpolys(g); // scratch space for encoding consts
    std::vector<PtxtPtr> ptrs;  // pointers to these constants
    for (long j = 0; j < dDivg; j++) { // giant steps
      long jg = j*g;            // beginning index of this giant step

      // get all the constants rot^{-i}(const_{i+g*j}) for this step
      PtxtPtr::Type ty = getConsts(ptrs, cpolys, zcp, dcp,
                                   dim, jg, D, oneTransform);

      if (ty==PtxtPtr::ZERO) continue; // all consts are zero

      // Store constants in cache and/or multiply/add them

      if (ctxt!=nullptr) {  // rotate by jg, multiply & add
        Ctxt shCtxt(*ctxt); // temporary to hold current rot^{g*j}(X)

        if (j>0) {
          shCtxt = *ctxt;
          ea.rotate1D(shCtxt, dim, jg); // rotate by i
	} // if j==0 we already have *shCtxt == *ctxt

        if (ty==PtxtPtr::DCRT) for (long i=0; i<min(g,D-jg); i++) {
            if (ptrs[i].dp!=nullptr) {
              Ctxt tmp(shCtxt);
              tmp.multByConstant(*(ptrs[i].dp));
              acc[i] += tmp;
            }
          }
        else /*ty==PtxtPtr::ZZX*/ for (long i=0; i<min(g,D-jg); i++) {
            if (ptrs[i].zp!=nullptr) {
              Ctxt tmp(shCtxt);
              tmp.multByConstant(*(ptrs[i].zp));
              acc[i] += tmp;
            }
          }
      }

      if (buildCache==cachezzX) for (long i=0; i<min(g,D-jg); i++) {
          if (ptrs[i].zp!=nullptr)
            (*zCache)[i+jg].reset(new zzX(*(ptrs[i].zp)));
        }
      else if (buildCache==cacheDCRT) for (long i=0; i<min(g,D-jg); i++) {
          if (ptrs[i].zp!=nullptr)
            (*dCache)[i+jg].reset(new DoubleCRT(*(ptrs[i].zp),ea.getContext()));
        }
      // end of giant-step loop
    }
    // Compute the result as \sum_{i=0}^{g-1} rho^i(Y_i)
    if (ctxt!=nullptr) {
      *ctxt = acc[0];
      for (long i = 1; i < g; i++) {
        ea.rotate1D(acc[i], dim, i);
        *ctxt += acc[i];
      }
    }
    // "install" the cache if needed
    if (buildCache == cachezzX)
      mat.installzzxcache(zCache);
    else if (buildCache == cacheDCRT)
      mat.installDCRTcache(dCache);
  } // end of multiply(...)
};

// Wrapper functions around the implemenmtation class
static void matmul1d(Ctxt* ctxt, MatMulBase& mat, long dim, bool oneTransform,
                     MatrixCacheType buildCache)
{
  MatMulLock locking(mat, buildCache);

  // If locking.getType()!=cacheEmpty then we really do need to
  // build the cache, and we also have the lock for it.

  if (locking.getType()==cacheEmpty && ctxt==nullptr) //  nothing to do
    return;

  switch (mat.getEA().getTag()) {
    case PA_GF2_tag: {
      matmul1D_impl<PA_GF2> M(mat, dim, locking.getType());
      M.multiply(ctxt, dim, oneTransform);
      break;
    }
    case PA_zz_p_tag: {
      matmul1D_impl<PA_zz_p> M(mat, dim, locking.getType());
      M.multiply(ctxt, dim, oneTransform);
      break;
    }
    default:
      throw std::logic_error("matmul1d: neither PA_GF2 nor PA_zz_p");
  }
}
void buildCache4MatMul1D(MatMulBase& mat, long dim,
                         MatrixCacheType buildCache)
{ matmul1d(nullptr, mat, dim, true, buildCache); }

void matMul1D(Ctxt& ctxt, MatMulBase& mat,long dim,
              MatrixCacheType buildCache)
{ matmul1d(&ctxt, mat, dim, true, buildCache); }

void buildCache4MatMulti1D(MatMulBase& mat,long dim,
                           MatrixCacheType buildCache)
{ matmul1d(nullptr, mat, dim, false, buildCache); }

void matMulti1D(Ctxt& ctxt,MatMulBase& mat,long dim,
                MatrixCacheType buildCache)
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
