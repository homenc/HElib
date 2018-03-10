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
/* blockMatmul.cpp - Data-movement operations on arrays of slots
 */
#include <algorithm>
#include <NTL/BasicThreadPool.h>
#include "matmul.h"

// A class that implements (sparse) linear transformation over the base ring
template<class type> class blockMatmul_impl {
  PA_INJECT(type)
  const MatrixCacheType buildCache;
  std::unique_ptr<CachedzzxMatrix> zCache;
  std::unique_ptr<CachedDCRTMatrix> dCache;

  BlockMatMul<type>& mat;
  const EncryptedArrayDerived<type>& ea;

public:
  blockMatmul_impl(MatMulBase& _mat, MatrixCacheType tag)
    : buildCache(tag), mat(dynamic_cast< BlockMatMul<type>& >(_mat)),
      ea(_mat.getEA().getDerived(type()))
  {
    long n = ea.size() * ea.getDegree();
    if (buildCache==cachezzX)
      zCache.reset(new CachedzzxMatrix(NTL::INIT_SIZE, n));
    else if (buildCache==cacheDCRT)
      dCache.reset(new CachedDCRTMatrix(NTL::INIT_SIZE, n));
  }

  // Extract one "column" from a matrix that was built with buildLinPolyCoeffs
  bool shiftedColumnInDiag(zzX& zpoly, long k,
                           const std::vector< std::vector<RX> >& diag)
  {
    // extract "column" with index k and store it in cvev
    bool zero = true;
    vector<RX> cvec(ea.size());
    for (long j = 0; j < ea.size(); j++) {
      cvec[j] = diag[j][k];
      if (!NTL::IsZero(cvec[j])) zero = false;
    }
    if (zero) { // all are zeros
      convert(zpoly, RX::zero());
      return true;
    }
    ea.encode(zpoly, cvec);

    if (k>0) {
      long p = ea.getContext().zMStar.getP();
      long m = ea.getContext().zMStar.getM();
      long d = ea.getDegree();
      long exp = PowerMod(mcMod(p, m), d-k, m); // apply inverse automorphism
      const auto& F = ea.getTab().getPhimXMod();
      RX rpoly1, rpoly2;

      convert(rpoly1, zpoly);
      plaintextAutomorph(rpoly2,rpoly1, exp, m, F);
      convert(zpoly, rpoly2);
    }
    return false;
  }

  // Process a single block diagonal with index idx, making calls
  // to the get(i,j) method and storing the result in the diag argument
  bool processDiagonal(std::vector< std::vector<RX> >& diag,
                       long idx, long nslots, long d)
  {
    bool zDiag = true;
    long nzLast = -1;
    std::vector<RX> entry1(d);
    mat_R entry(INIT_SIZE, d, d);
    for (long j = 0; j < nslots; j++) {
      bool zEntry = mat.get(entry, mcMod(j-idx, nslots), j);
      if (!zEntry && IsZero(entry)) zEntry=true; // zero is an empty entry too
      // get(...) returns true if the entry is empty, false otherwise

      assert(zEntry || (entry.NumRows() == d && entry.NumCols() == d));

      if (!zEntry) {    // non-empty entry
        zDiag = false;  // mark diagonal as non-empty

        // clear entries between last nonzero entry and this one
        for (long jj = nzLast+1; jj < j; jj++)
          for (long k = 0; k < d; k++)
            clear(diag[jj][k]);
	nzLast = j;

	// recode entry as a vector of polynomials
	for (long k = 0; k < d; k++) conv(entry1[k], entry[k]);

	ea.buildLinPolyCoeffs(diag[j], entry1);
      } // now diag[j] contains the lin poly coeffs
    }

    if (!zDiag) // clear trailing zero entries    
      for (long jj = nzLast+1; jj < nslots; jj++)
	for (long k = 0; k < d; k++)
	  clear(diag[jj][k]);

    return zDiag;
  }

  // This code has a complexity of N+d (instead of N*d) where N is the
  // number of nonzero diagonal blocks. However, it requires space for
  // d extra ciphertexts
  void multiply(Ctxt* ctxt) 
  {
    RBak bak; bak.save(); ea.getTab().restoreContext();

    const RXModulus& F = ea.getTab().getPhimXMod();
    const PAlgebra& zMStar = ea.getContext().zMStar;
    long p = zMStar.getP(); 
    long m = zMStar.getM();
    long nslots = ea.size();
    long d = ea.getDegree();

    std::vector< std::vector<RX> > diag(nslots);
    for (long j = 0; j < nslots; j++) diag[j].resize(d);

    std::vector<Ctxt> acc;
    std::unique_ptr<Ctxt> shCtxt;
    if (ctxt!=nullptr) { // we need to do an actual multiplication
      ctxt->cleanUp();   // not sure, but this may be a good idea
      acc.assign(d, Ctxt(ZeroCtxtLike, *ctxt));
      shCtxt.reset(new Ctxt(*ctxt));
    }

    // Check if we have the relevant constant in cache
    CachedzzxMatrix* zcp;
    CachedDCRTMatrix* dcp;
    mat.getCache(&zcp, &dcp);

    long lastShift = 0;
    bool sequential = (ea.dimension()==1) && ea.nativeDimension(0);
    // if just a single native dimension, then rotate adds only little noise

    for (long i = 0; i < nslots; i++) { // Process the diagonals one at a time
      bool zDiag = true;
      // For each diagonal i, we update the d accumulators y_0,..,y_{d-1}
      // with y_k += \sigma^{-k}(\lambda_{i,k}) * \rho^i(x)

      std::vector<zzX> zpoly(d);
      std::vector<zzX*> zzxPtr(d,nullptr);
      std::vector<DoubleCRT*> dcrtPtr(d,nullptr);

      if (dcp!=nullptr) for (long k=0; k<d; k++) { // DoubleCRT cache exists
          dcrtPtr[k] = (*dcp)[d*i +k].get();
          if (dcrtPtr[k]!=nullptr) zDiag = false;
      }
      else if (zcp!=nullptr) for (long k=0; k<d; k++) {// zzX but no DoubleCRT
          zzxPtr[k] = (*zcp)[d*i +k].get();
          if (zzxPtr[k]!=nullptr) zDiag = false;
      } else {
        zDiag = processDiagonal(diag, i, nslots, d);

        // extract the "columns" from diag and encode them in zpoly
        if (!zDiag) for (long k=0; k<d; k++) {
            if (!shiftedColumnInDiag(zpoly[k],k,diag))// compute the constant
            zzxPtr[k] = &(zpoly[k]);                  // returns true on zero
        }
      }
      if (zDiag) continue; // nothing to do for this diagonal
      // done preparing all the zzxPtr, dcrtPtr variables

      if (ctxt!=nullptr) {
        if (i > 0) { // rotate the ciphertext
          if (sequential)  // rotate the previous version
            ea.rotate(*shCtxt, i-lastShift);
          else {           // rotate the original ciphertext
            *shCtxt = *ctxt;
            ea.rotate(*shCtxt, i);
	  }
	  lastShift = i;
	} // if (i>0)

	// Depending on zzxPtr, dcrtPtr, update the accumulated sums
	for (long k=0; k<d; k++) if (dcrtPtr[k]!=NULL || zzxPtr[k]!=NULL) {
            Ctxt tmp1(*shCtxt);
            if (dcrtPtr[k] != NULL) tmp1.multByConstant(*(dcrtPtr[k]));
            else                    tmp1.multByConstant(*(zzxPtr[k]));
            acc[k] += tmp1;
          }
      } // if (ctxt!=nullptr)

      // allocate constants and store in the cache, if needed
      if (buildCache==cachezzX) for (long k=0; k<d; k++) {
          (*zCache)[d*i +k].reset( new zzX(*(zzxPtr[k])) );
      }
      else if (buildCache==cacheDCRT) for (long k=0; k<d; k++) {
          (*dCache)[d*i +k].reset(new DoubleCRT(*(zzxPtr[k]),ea.getContext()));
      }
    } // end of i'th diagonal

    // Finally, compute the result as \sum_{k=0}^{d-1} \sigma_k(y_k)
    if (ctxt!=nullptr) {
      *ctxt = acc[0];
      for (long k = 1; k < d; k++) {
	acc[k].frobeniusAutomorph(k);
	*ctxt += acc[k];
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
static void block_matmul(Ctxt* ctxt,MatMulBase& mat,MatrixCacheType buildCache)
{
  MatMulLock locking(mat, buildCache);

  // If locking.getType()!=cacheEmpty then we really do need to
  // build the cache, and we also have the lock for it.

  if (locking.getType()==cacheEmpty && ctxt==nullptr) //  nothing to do
    return;

  switch (mat.getEA().getTag()) {
    case PA_GF2_tag: {
      blockMatmul_impl<PA_GF2> M(mat, locking.getType());
      M.multiply(ctxt);
      break;
    }
    case PA_zz_p_tag: {
      blockMatmul_impl<PA_zz_p> M(mat, locking.getType());
      M.multiply(ctxt);
      break;
    }
    default:
      throw std::logic_error("block_matmul: neither PA_GF2 nor PA_zz_p");
  }
}
void buildCache4BlockMatMul(MatMulBase& mat, MatrixCacheType buildCache)
{ block_matmul(nullptr, mat, buildCache); }

void blockMatMul(Ctxt& ctxt, MatMulBase& mat, MatrixCacheType buildCache)
{ block_matmul(&ctxt, mat, buildCache); }



// Applying matmul to plaintext, useful for debugging
template<class type>
class blockMatmul_pa_impl {
public:
  PA_INJECT(type)

  static void multiply(NewPlaintextArray& pa, BlockMatMul<type>& mat)
  {
    const EncryptedArrayDerived<type>& ea = mat.getEA().getDerived(type());
    const PAlgebra& zMStar = ea.getContext().zMStar;
    RBak bak; bak.save(); ea.getTab().restoreContext();

    long n = ea.size();
    long d = ea.getDegree();
    vector<RX>& data = pa.getData<type>();

    vector<RX> res(n);
    for (long j = 0; j < n; j++) {
      vec_R acc, tmp, tmp1;
      mat_R val;

      acc.SetLength(d);
      for (long i = 0; i < n; i++) {
         if (!mat.get(val, i, j)) {
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
void blockMatMul(NewPlaintextArray& pa, MatMulBase& mat)
{
  switch (mat.getEA().getTag()) {
    case PA_GF2_tag: {
      BlockMatMul<PA_GF2>& mat1= dynamic_cast< BlockMatMul<PA_GF2>& >(mat);
      blockMatmul_pa_impl<PA_GF2>::multiply(pa, mat1);
      break;
    }
    case PA_zz_p_tag: {
      BlockMatMul<PA_zz_p>& mat1= dynamic_cast< BlockMatMul<PA_zz_p>& >(mat);
      blockMatmul_pa_impl<PA_zz_p>::multiply(pa, mat1);
      break;
    }
    default:
      throw std::logic_error("blockMatMul: neither PA_GF2 nor PA_zz_p");
  }
}
