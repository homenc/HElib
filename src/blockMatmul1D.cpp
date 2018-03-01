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
/* blockMatmul1D.cpp - Data-movement operations on arrays of slots
 */
#include <algorithm>
#include <NTL/BasicThreadPool.h>
#include "matmul.h"

// A class that implements the basic (sparse) 1D matrix-vector functions
template<class type> class blockMatmul1D_impl {
  PA_INJECT(type)
  const MatrixCacheType buildCache;
  std::unique_ptr<CachedzzxMatrix> zCache;
  std::unique_ptr<CachedDCRTMatrix> dCache;

  BlockMatMul<type>& mat;
  const EncryptedArrayDerived<type>& ea;

public:
  blockMatmul1D_impl(MatMulBase& _mat, MatrixCacheType tag, long dim)
    : buildCache(tag), mat(dynamic_cast< BlockMatMul<type>& >(_mat)),
      ea(_mat.getEA().getDerived(type()))
  {
    long n = ea.sizeOfDimension(dim) * ea.getDegree();
    if (buildCache==cachezzX)
      zCache.reset(new CachedzzxMatrix(NTL::INIT_SIZE, n));
    else if (buildCache==cacheDCRT)
      dCache.reset(new CachedDCRTMatrix(NTL::INIT_SIZE, n));
  }

  // Extract one "column" from a matrix that was built with buildLinPolyCoeffs
  bool shiftedColumnInDiag(zzX& zpoly, long f,
                           const std::vector< std::vector<RX> >& diag)
  {
    // extract "column" with index f and store it in cvev
    bool zero = true;
    vector<RX> cvec(ea.size());
    for (long j = 0; j < ea.size(); j++) {
      cvec[j] = diag[j][f];
      if (!NTL::IsZero(cvec[j])) zero = false;
    }
    if (zero) { // all are zeros
      zpoly.kill();
      return true;
    }
    ea.encode(zpoly, cvec);

    if (f>0) {
      long p = ea.getContext().zMStar.getP();
      long m = ea.getContext().zMStar.getM();
      long d = ea.getDegree();
      long exp = PowerMod(mcMod(p, m), d-f, m); // apply inverse automorphism
      const auto& F = ea.getTab().getPhimXMod();
      RX rpoly1, rpoly2;

      convert(rpoly1, zpoly);
      plaintextAutomorph(rpoly2,rpoly1, exp, m, F);
      convert(zpoly, rpoly2);
    }
    return false;
  }

  // Process a single block diagonal with index idx along dimenssion dim,
  // making calls to the get(i,j) method and storing the result in
  // the diag matrix (# of vectors = extenssion-degree).
  bool processDiagonal1(std::vector< std::vector<RX> >& diag,
			long dim, long i, long D, long d)
  {
    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry

    mat_R entry(INIT_SIZE, d, d);
    std::vector<RX> entry1(d);
    std::vector< std::vector<RX> > tmpDiag(D);

    // Process the entries in this diagonal one at a time
    for (long j = 0; j < D; j++) { // process entry j
      bool zEntry = mat.get(entry, mcMod(j-i, D), j); // entry [j-i mod D, j]
      // get(...) returns true if the entry is empty, false otherwise

      if (!zEntry && IsZero(entry)) zEntry = true;// zero is an empty entry too
      assert(zEntry || (entry.NumRows() == d && entry.NumCols() == d));

      if (!zEntry) {   // not a zero entry
        zDiag = false; // mark diagonal as non-empty

	for (long jj = nzLast+1; jj < j; jj++) {// clear from last nonzero entry
          tmpDiag[jj].assign(d, RX());
        }
        nzLast = j; // current entry is the last nonzero one

        // recode entry as a vector of polynomials
        for (long k = 0; k < d; k++) conv(entry1[k], entry[k]);

        // compute the linearlized polynomial coefficients
	ea.buildLinPolyCoeffs(tmpDiag[j], entry1);
      }
    }
    if (zDiag) return true; // zero diagonal, nothing to do

    // clear trailing zero entries
    for (long jj = nzLast+1; jj < D; jj++) {
      tmpDiag[jj].assign(d, RX());
    }

    if (D==1) diag.assign(ea.size(), tmpDiag[0]); // dimension of size one
    else for (long j = 0; j < ea.size(); j++)
           diag[j] = tmpDiag[ ea.coordinate(dim,j) ];
           // rearrange the indexes based on the current dimension

    return false; // a nonzero diagonal
  }

  // Process a single block diagonal with index idx along dimenssion dim,
  // making calls to the multiGet(i,j,k) method and storing the result in
  // the diag matrix (# of vectors = extenssion-degree).
  bool processDiagonal2(std::vector< std::vector<RX> >& diag,
			long dim, long idx, long D, long d)
  {
    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry

    mat_R entry(INIT_SIZE, d, d);
    std::vector<RX> entry1(d);

    // Get the slots in this diagonal one at a time
    long blockIdx, rowIdx, colIdx;
    for (long j = 0; j < ea.size(); j++) { // process entry j
      if (dim == ea.dimension()) { // "special" last dimenssion of size 1
	rowIdx = colIdx = 0; blockIdx=j;
      } else {
        std::tie(blockIdx, colIdx)
	  = ea.getContext().zMStar.breakIndexByDim(j, dim);
	rowIdx = mcMod(colIdx-idx,D);
      }
      bool zEntry = mat.multiGet(entry,rowIdx,colIdx,blockIdx);
      // entry [i,j-i mod D] in the block corresponding to blockIdx
      // multiGet(...) returns true if the entry is empty, false otherwise

      if (!zEntry && IsZero(entry)) zEntry=true; // zero is an empty entry too
      assert(zEntry ||
             (entry.NumRows() == d && entry.NumCols() == d));

      if (!zEntry) {    // non-empty entry
	zDiag = false;  // mark diagonal as non-empty

	for (long jj = nzLast+1; jj < j; jj++) {// clear from last nonzero entry
	  for (long k = 0; k < d; k++)
	    clear(diag[jj][k]);
        }
	nzLast = j; // current entry is the last nonzero one

	// recode entry as a vector of polynomials
	for (long k = 0; k < d; k++) conv(entry1[k], entry[k]);

        // compute the linearlized polynomial coefficients
	ea.buildLinPolyCoeffs(diag[j], entry1);
      }
    }
    if (zDiag) return true; // zero diagonal, nothing to do

    // clear trailing zero entries
    for (long jj = nzLast+1; jj < ea.size(); jj++)
      for (long k = 0; k < d; k++)
	clear(diag[jj][k]);

    return false; // a nonzero diagonal
  }

  void multiply(Ctxt* ctxt, long dim, bool oneTransform) 
  {
    assert(dim >= 0 && dim <= ea.dimension());
    RBak bak; bak.save(); ea.getTab().restoreContext(); // backup NTL modulus

    long nslots = ea.size();
    long d = ea.getDegree();
    long D = (dim == ea.dimension()) ? 1 : ea.sizeOfDimension(dim);

    const PAlgebra& zMStar = ea.getContext().zMStar;
    long p = zMStar.getP(); 
    long m = zMStar.getM();

    std::vector< std::vector<RX> > diag(nslots); // scratch space
    for (long j = 0; j < nslots; j++) diag[j].resize(d);

    std::vector<Ctxt> acc;
    std::unique_ptr<Ctxt> shCtxt;
    if (ctxt!=nullptr) { // we need to do an actual multiplication
      ctxt->cleanUp(); // not sure, but this may be a good idea
      acc.assign(d, Ctxt(ZeroCtxtLike, *ctxt));
      shCtxt.reset(new Ctxt(*ctxt));
    }

    // Check if we have the relevant constant in cache
    CachedzzxMatrix* zcp;
    CachedDCRTMatrix* dcp;
    mat.getCache(&zcp, &dcp);

    // Process the diagonals one at a time
    for (long e = 0; e < D; e++) { // process diagonal e
      bool zeroDiag = true;
      // For each diagonal e, we update the d accumulators y_0,..,y_{d-1}
      // with y_f += \sigma^{-f}(\lambda_{e,f}) * \rho^e(x)

      std::vector<zzX> zpoly(d, zzX());
      std::vector<zzX*> zzxPtr(d,NULL);
      std::vector<DoubleCRT*> dcrtPtr(d,NULL);

      if (dcp!=nullptr) for (long f=0; f<d; f++) { // DoubleCRT cache exists
          dcrtPtr[f] = (*dcp)[d*e +f].get();
          if (dcrtPtr[f]!=nullptr) zeroDiag = false;
      }
      else if (zcp!=nullptr) for (long f=0; f<d; f++) {// zzX but no DoubleCRT
          zzxPtr[f] = (*zcp)[d*e +f].get();
          if (zzxPtr[f]!=nullptr) zeroDiag = false;
      } else {
       zeroDiag = oneTransform? this->processDiagonal1(diag, dim, e, D, d)
                               : this->processDiagonal2(diag, dim, e, D, d);

        // extract the "columns" from diag and encode them in zpoly
        if (!zeroDiag) for (long f=0; f<d; f++) {
          if (!shiftedColumnInDiag(zpoly[f], f, diag))// returns true on zero
            zzxPtr[f] = &(zpoly[f]);
        }
      }
      if (zeroDiag) continue; // nothing to do for this diagonal
      // done preparing all the zzxPtr, dcrtPtr variables

      // Rotate the ciphertext to position corresponding to diagonal e.
      // The code below uses only rotate-by-one operations if this is
      // a good dimenssion and the matrix is dense, hence it may require
      // fewer key-switching matrices

      if (ctxt!=nullptr) {
	if (e > 0) { // rotate the ciphertext
          *shCtxt = *ctxt;
          ea.rotate1D(*shCtxt, dim, e);
	} // if (e>0)
	// The implementation above incurs an extra mult-by-constant due
	// to the masks in rotate1D when applied in a "bad dimension".
	// These masks can be folded into the constants here, but then
	// we would need two constants for each (e,f) rather than one,
	// namely const*mask and const*(1-mask).
	// We should implement that optimization at some point.

	// Depending on zzxPtr, dcrtPtr, update the accumulated sums
	for (long f=0; f<d; f++) if (dcrtPtr[f]!=NULL || zzxPtr[f]!=NULL) {
            Ctxt tmp1(*shCtxt);
            if (dcrtPtr[f] != NULL) tmp1.multByConstant(*(dcrtPtr[f]));
            else                    tmp1.multByConstant(*(zzxPtr[f]));
            acc[f] += tmp1;
          }
      } // if (ctxt!=nullptr)
 
      // allocate constants and store in the cache, if needed
     if (buildCache==cachezzX) for (long f=0; f<d; f++) {
          (*zCache)[d*e +f].reset( new zzX(*(zzxPtr[f])) );
      }
      else if (buildCache==cacheDCRT) for (long f=0; f<d; f++) {
          (*dCache)[d*e +f].reset(new DoubleCRT(*(zzxPtr[f]),ea.getContext()));
      }
    } // end of e'th diagonal

    // Finally, compute the result as \sum_{f=0}^{d-1} \sigma_f(y_f)
    if (ctxt!=nullptr) {
      *ctxt = acc[0];
      for (long f = 1; f < d; f++) {
	acc[f].frobeniusAutomorph(f);
	*ctxt += acc[f];
      }
    }
    // "install" the cache (if needed)
    if (buildCache == cachezzX)
      mat.installzzxcache(zCache);
    else if (buildCache == cacheDCRT)
      mat.installDCRTcache(dCache);
  } // end of multiply(...)
};

// Wrapper functions around the implemenmtation class
static void blockMatmul1d(Ctxt* ctxt, MatMulBase& mat, long dim,
		          MatrixCacheType buildCache, bool oneTransform)
{
  MatMulLock locking(mat, buildCache);

  // If locking.getType()!=cacheEmpty then we really do need to
  // build the cache, and we also have the lock for it.

  if (locking.getType()==cacheEmpty && ctxt==nullptr) //  nothing to do
    return;

  switch (mat.getEA().getTag()) {
    case PA_GF2_tag: {
      blockMatmul1D_impl<PA_GF2> M(mat, locking.getType(), dim);
      M.multiply(ctxt, dim, oneTransform);
      break;
    }
    case PA_zz_p_tag: {
      blockMatmul1D_impl<PA_zz_p> M(mat, locking.getType(), dim);
      M.multiply(ctxt, dim, oneTransform);
      break;
    }
    default:
      throw std::logic_error("matmul1d: neither PA_GF2 nor PA_zz_p");
  }
}
void buildCache4BlockMatMul1D(MatMulBase& mat,
			      long dim, MatrixCacheType buildCache)
{ blockMatmul1d(nullptr, mat, dim, buildCache, true); }

void blockMatMul1D(Ctxt& ctxt, MatMulBase& mat, long dim,
               MatrixCacheType buildCache)
{ blockMatmul1d(&ctxt, mat, dim, buildCache, true); }

void buildCache4BlockMatMulti1D(MatMulBase& mat,
				long dim, MatrixCacheType buildCache)
{ blockMatmul1d(nullptr, mat, dim, buildCache, false); }

void blockMatMulti1D(Ctxt& ctxt, MatMulBase& mat, long dim,
		     MatrixCacheType buildCache)
{ blockMatmul1d(&ctxt, mat, dim, buildCache, false); }



// Versions for plaintext rather than ciphertext, useful for debugging
template<class type> class blockMatmul1D_pa_impl {
public:
  PA_INJECT(type)

  static void multiply(NewPlaintextArray& pa, BlockMatMul<type>& mat,
		       long dim, bool oneTrans=false)
  {
    const EncryptedArrayDerived<type>& ea = mat.getEA().getDerived(type());
    const PAlgebra& zMStar = ea.getContext().zMStar;
    RBak bak; bak.save(); ea.getTab().restoreContext();

    long n = ea.size();
    long D = ea.sizeOfDimension(dim);
    long d = ea.getDegree();

    vector< vector<RX> > data1(n/D);
    for (long k = 0; k < n/D; k++)
      data1[k].resize(D);

    // copy the data into a vector of 1D vectors
    vector<RX>& data = pa.getData<type>();
    for (long i = 0; i < n; i++) {
      long k,j;
      std::tie(k,j) = zMStar.breakIndexByDim(i, dim);
      data1[k][j] = data[i];       // k= along dim, j = the rest of i
    }
    for (long k = 0; k < n/D; k++) { // multiply each vector by a matrix
      for (long j = 0; j < D; j++) { // matrix-vector multiplication
	vec_R acc, tmp, tmp1;
	mat_R val;
	acc.SetLength(d);
	for (long i = 0; i < D; i++) {
          bool zero = oneTrans? mat.get(val, i, j)
                              : mat.multiGet(val, i, j, k);
	  if (!zero) { // if non-zero, multiply and add
            VectorCopy(tmp1, data1[k][i], d);
            mul(tmp, tmp1, val);
            add(acc, acc, tmp);
	  }
	}
	long idx = zMStar.assembleIndexByDim(make_pair(k,j), dim);
        conv(data[idx], acc);
      }
    }
  }
}; 
void blockMatMul1D(NewPlaintextArray& pa, MatMulBase& mat, long dim)
{
  switch (mat.getEA().getTag()) {
    case PA_GF2_tag: {
      BlockMatMul<PA_GF2>& mat1= dynamic_cast< BlockMatMul<PA_GF2>& >(mat);
      blockMatmul1D_pa_impl<PA_GF2>::multiply(pa, mat1, dim, true);
      break;
    }
    case PA_zz_p_tag: {
      BlockMatMul<PA_zz_p>& mat1= dynamic_cast< BlockMatMul<PA_zz_p>& >(mat);
      blockMatmul1D_pa_impl<PA_zz_p>::multiply(pa, mat1, dim, true);
      break;
    }
    default:
      throw std::logic_error("blockMatMul1D: neither PA_GF2 nor PA_zz_p");
  }
}
void blockMatMulti1D(NewPlaintextArray& pa, MatMulBase& mat, long dim)
{
  switch (mat.getEA().getTag()) {
    case PA_GF2_tag: {
      BlockMatMul<PA_GF2>& mat1= dynamic_cast< BlockMatMul<PA_GF2>& >(mat);
      blockMatmul1D_pa_impl<PA_GF2>::multiply(pa, mat1, dim, false);
      break;
    }
    case PA_zz_p_tag: {
      BlockMatMul<PA_zz_p>& mat1= dynamic_cast< BlockMatMul<PA_zz_p>& >(mat);
      blockMatmul1D_pa_impl<PA_zz_p>::multiply(pa, mat1, dim, false);
      break;
    }
    default:
      throw std::logic_error("blockMatMulti1D: neither PA_GF2 nor PA_zz_p");
  }
}
