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
/* matmul.cpp - Data-movement operations on arrays of slots
 */
#include <algorithm>
#include <NTL/BasicThreadPool.h>
#include "matmul.h"

#ifdef DEBUG
void printCache(const CachedzzxMatrix& cache);
void printCache(const CachedDCRTMatrix& cache);
#endif


/********************************************************************
 ************* Helper routines to handle caches *********************/

bool MatMulBase::lockCache(MatrixCacheType ty)
{
  // Check if we need to build cache (1st time)
  if (ty==cacheEmpty || hasDCRTcache()
      || (ty==cachezzX && haszzxcache())) // no need to build
    return false;

  cachelock.lock();  // Take the lock

  // Check again if we need to build the cache
  if (ty==cacheEmpty || hasDCRTcache()
      || (ty==cachezzX && haszzxcache())) { // no need to build
    cachelock.unlock();
    return false; // no need to build
  }

  // We have the lock, can upgrade zzx to dcrt if needed
  if (ty==cacheDCRT && haszzxcache()) {
    upgradeCache();     // upgrade zzx to DCRT
    cachelock.unlock(); // release the lock
    return false; // already built
  }
  return true; // need to build, and we have the lock
}

// This function assumes that we have the lock
void MatMulBase::upgradeCache()
{
  std::unique_ptr<CachedDCRTMatrix> dCache(new CachedDCRTMatrix());
  dCache->SetLength(zzxCache->length());
  for (long i=0; i<zzxCache->length(); i++) {
    if ((*zzxCache)[i] != nullptr)
      (*dCache)[i].reset(new DoubleCRT(*(*zzxCache)[i], ea.getContext()));
  }
  dcrtCache.swap(dCache);
}
/********************************************************************/


/********************************************************************
 * An implementation class for dense matmul.
 *
 * Using such a class (rather than plain functions) is convenient
 * since you only need to PA_INJECT once (and the injected names do
 * not pollute the external namespace).
 * We also use it here to store some pieces of data between recursive
 * calls, as well as pointers to temporary caches as we build them.
 *********************************************************************/
template<class type> class matmul_impl {
  PA_INJECT(type)
  const MatrixCacheType buildCache;
  std::unique_ptr<CachedzzxMatrix> zCache;
  std::unique_ptr<CachedDCRTMatrix> dCache;

  Ctxt* res;
  MatMul<type>& mat;
  std::vector<long> dims;
  const EncryptedArrayDerived<type>& ea;

  // helper class to sort dimensions, so that
  //    - bad dimensions come before good dimensions (primary sort key)
  //    - small dimensions come before large dimesnions (secondary sort key)
  // this is a good order to process the dimensions in the recursive
  // mat_mul_dense routine: it ensures that the work done at the leaves
  // of the recursion is minimized, and that the work done at the non-leaves
  // is dominated by the work done at the leaves.
  struct MatMulDimComp {
    const EncryptedArrayDerived<type> *ea;
    MatMulDimComp(const EncryptedArrayDerived<type> *_ea) : ea(_ea) {}

    bool operator()(long i, long j) { 
      return (!ea->nativeDimension(i) && ea->nativeDimension(j)) ||
             (  (ea->nativeDimension(i) == ea->nativeDimension(j)) &&
                (ea->sizeOfDimension(i) < ea->sizeOfDimension(j))  );
    }
  };

public:
  matmul_impl(Ctxt* c, MatMulBase& _mat, MatrixCacheType tag)
    : buildCache(tag), res(c),
      mat(dynamic_cast< MatMul<type>& >(_mat)),
      ea(_mat.getEA().getDerived(type()))
  {
    if (buildCache==cachezzX)
      zCache.reset(new CachedzzxMatrix(NTL::INIT_SIZE,ea.size()));
    else if (buildCache==cacheDCRT)
      dCache.reset(new CachedDCRTMatrix(NTL::INIT_SIZE,ea.size()));

    // dims stores the order of dimensions
    dims.resize(ea.dimension());
    for (long i = 0; i < ea.dimension(); i++)
      dims[i] = i;
    sort(dims.begin(), dims.end(), MatMulDimComp(&ea));
    // sort the dimenesions so that bad ones come before good,
    // and then small ones come before large
  }

  // Get a diagonal encoded as a single constant
  bool processDiagonal(zzX& epmat, const vector<long>& idxes)
  {
    vector<RX> pmat;  // the plaintext diagonal
    pmat.resize(ea.size());
    bool zDiag = true; // is this a zero diagonal
    for (long j = 0; j < ea.size(); j++) {
      long i = idxes[j];
      RX val;
      if (mat.get(val, i, j)) // returns true if the entry is zero
	clear(pmat[j]);
      else {           // not a zero entry
	pmat[j] = val;
	zDiag = false; // not a zero diagonal
      }
    }
    // Now we have the constants for all the diagonal entries, encode the
    // diagonal as a single polynomial with these constants in the slots
    if (!zDiag) {
      ea.encode(epmat, pmat);
    }
    return zDiag;
  }


  // A recursive matrix-by-vector multiply, used by the dense matrix code.
  // This routine is optimized to use only the rotate1D routine rather
  // than the more expensive linear-array rotations.
  long rec_mul(const Ctxt* pdata, long dim, long idx,
               const vector<long>& idxes)
  {
    if (dim >= ea.dimension()) { // Last dimension (recursion edge condition)
      zzX pt;
      zzX* zxPtr=nullptr;
      DoubleCRT* dxPtr=nullptr;

      // Check if we have the relevant constant in cache
      CachedzzxMatrix* zcp;
      CachedDCRTMatrix* dcp;
      mat.getCache(&zcp, &dcp);
      if (dcp != nullptr)         // DoubleCRT cache exists
	dxPtr = (*dcp)[idx].get();
      else if (zcp != nullptr)    // zzx cache exists but no DoubleCRT
	zxPtr = (*zcp)[idx].get();
      else if (!processDiagonal(pt, idxes)) // no cache, compute const
	zxPtr = &pt; // if it is not a zero value, point to it

      // if constant is zero, return without doing anything
      if (zxPtr==nullptr && dxPtr==nullptr)
	return idx+1;

      // Constant is non-zero, store it in cache and/or multiply/add it

      if (pdata!=nullptr && res!=nullptr) {
        Ctxt tmp = *pdata;
        if (dxPtr!=nullptr) tmp.multByConstant(*dxPtr); // mult by DCRT
        else                tmp.multByConstant(*zxPtr); // mult by zzx
        *res += tmp;
      }

      if (buildCache==cachezzX) {
        (*zCache)[idx].reset(new zzX(*zxPtr));
      }
      else if (buildCache==cacheDCRT) {
        (*dCache)[idx].reset(new DoubleCRT(*zxPtr, ea.getContext()));
      }
      return idx+1;      
    }

    // not the last dimension, make a recursive call
    long sdim = ea.sizeOfDimension(dims[dim]);

    // compute "in spirit" sum_i (pdata >> i) * i'th-diagonal, but
    // adjust the indexes so that we only need to rotate the cipehrtext
    // along the different dimensions separately
    for (long offset = 0; offset < sdim; offset++) {
      vector<long> idxes1;
      ea.EncryptedArrayBase::rotate1D(idxes1, idxes, dims[dim], offset);
      if (pdata!=nullptr && res!=nullptr) {
	Ctxt pdata1 = *pdata;
	ea.rotate1D(pdata1, dims[dim], offset);
	// indexes adjusted, make the recursive call
	idx = rec_mul(&pdata1, dim+1, idx, idxes1);
      }
      else // don't bother with the ciphertext
	idx = rec_mul(pdata, dim+1, idx, idxes1);
    }

    return idx;
  }

  // Multiply a ciphertext vector by a plaintext dense matrix
  // and/or build a cache with the multiplication constants
  void multilpy(Ctxt* ctxt) 
  {
    RBak bak; bak.save(); ea.getTab().restoreContext();
    // idxes describes a genealized diagonal, {(i,idx[i])}_i
    // initially just the identity, idx[i]==i
    vector<long> idxes(ea.size());
    for (long i = 0; i < ea.size(); i++) idxes[i] = i;

    // call the recursive procedure to do the actual work
    rec_mul(ctxt, 0, 0, idxes);

    if (ctxt!=nullptr && res!=nullptr)
      *ctxt = *res; // copy the result back to ctxt

    // "install" the cache (if needed)
    if (buildCache == cachezzX)
      mat.installzzxcache(zCache);
    else if (buildCache == cacheDCRT)
      mat.installDCRTcache(dCache);
  }
};

// Wrapper functions around the implemenmtation class
static void mat_mul(Ctxt* ctxt, MatMulBase& mat, MatrixCacheType buildCache)
{
  MatMulLock locking(mat, buildCache);

  // If locking.getType()!=cacheEmpty then we really do need to
  // build the cache, and we also have the lock for it.

  if (locking.getType() == cacheEmpty && ctxt==nullptr) //  nothing to do
    return;

  std::unique_ptr<Ctxt> res;
  if (ctxt!=nullptr) { // we need to do an actual multiplication
    ctxt->cleanUp(); // not sure, but this may be a good idea
    res.reset(new Ctxt(ZeroCtxtLike, *ctxt));
  }

  switch (mat.getEA().getTag()) {
    case PA_GF2_tag: {
      matmul_impl<PA_GF2> M(res.get(), mat, locking.getType());
      M.multilpy(ctxt);
      break;
    }
    case PA_zz_p_tag: {
      matmul_impl<PA_zz_p> M(res.get(), mat, locking.getType());
      M.multilpy(ctxt);
      break;
    }
    default:
      throw std::logic_error("mat_mul: neither PA_GF2 nor PA_zz_p");
  }
}
void buildCache4MatMul(MatMulBase& mat, MatrixCacheType buildCache)
{ mat_mul(nullptr, mat, buildCache); }

void matMul(Ctxt& ctxt, MatMulBase& mat, MatrixCacheType buildCache)
{ mat_mul(&ctxt, mat, buildCache); }



/********************************************************************
 * An implementation class for sparse diagonal matmul.
 *
 * Using such a class (rather than plain functions) is convenient
 * since you only need to PA_INJECT once (and the injected names do
 * not pollute the external namespace). We also use it here to store
 * pointers to temporary caches as we build them.
 ********************************************************************/
template<class type> class matmul_sparse_impl {
  PA_INJECT(type)
  const MatrixCacheType buildCache;
  std::unique_ptr<CachedzzxMatrix> zCache;
  std::unique_ptr<CachedDCRTMatrix> dCache;

  MatMul<type>& mat;
  const EncryptedArrayDerived<type>& ea;
public:
  matmul_sparse_impl(MatMulBase& _mat, MatrixCacheType tag)
    : buildCache(tag), mat(dynamic_cast< MatMul<type>& >(_mat)),
      ea(_mat.getEA().getDerived(type()))
  {
    if (buildCache==cachezzX)
      zCache.reset(new CachedzzxMatrix(NTL::INIT_SIZE,ea.size()));
    else if (buildCache==cacheDCRT)
      dCache.reset(new CachedDCRTMatrix(NTL::INIT_SIZE,ea.size()));
  }

  // Get a diagonal encoded as a single constant
  bool processDiagonal(zzX& cPoly, long diagIdx, long nslots)
  {
    bool zDiag = true; // is this a zero diagonal
    vector<RX> diag(ea.size()); // the plaintext diagonal

    for (long j = 0; j < nslots; j++) { // process entry j
      long ii = mcMod(j-diagIdx, nslots);    // j-diagIdx mod nslots
      bool zEntry = mat.get(diag[j], ii, j); // callback
      assert(zEntry || deg(diag[j]) < ea.getDegree());

      if (zEntry) clear(diag[j]);
      else if (!IsZero(diag[j]))
        zDiag = false; // diagonal is non-zero
    }
    // Now we have the constants for all the diagonal entries, encode the
    // diagonal as a single polynomial with these constants in the slots
    if (!zDiag) {
      ea.encode(cPoly, diag);
    }
    return zDiag;
  }

  void multiply(Ctxt* ctxt) 
  {
    RBak bak; bak.save(); ea.getTab().restoreContext();

    long nslots = ea.size();
    bool sequential = (ea.dimension()==1) && ea.nativeDimension(0);
    // if just a single native dimension, then rotate adds only little noise

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
    long lastRotate = 0;
    for (long i = 0; i < nslots; i++) {  // process diagonal i
      zzX cpoly;
      zzX* zxPtr=nullptr;
      DoubleCRT* dxPtr=nullptr;

      if (dcp != nullptr)         // DoubleCRT cache exists
	dxPtr = (*dcp)[i].get();
      else if (zcp != nullptr)    // zzx cache exists but no DoubleCRT
	zxPtr = (*zcp)[i].get();
      else { // no cache, compute const
        if (!processDiagonal(cpoly,i,nslots)) { // returns true if zero
          zxPtr = &cpoly;   // if it is not a zero value, point to it
	}
      }

      // if zero diagonal, nothing to do for this iteration
      if (zxPtr==nullptr && dxPtr==nullptr)
        continue;

      // Non-zero diagonal, store it in cache and/or multiply/add it

      if (ctxt!=nullptr && res!=nullptr) {
        // rotate by i, multiply by the polynomial, then add to the result
        if (i>0) {
          if (sequential) {
            ea.rotate(*ctxt, i-lastRotate);
            *shCtxt = *ctxt;
          } else {
            *shCtxt = *ctxt;
            ea.rotate(*shCtxt, i); // rotate by i
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
    }

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
static void
mat_mul_sparse(Ctxt* ctxt, MatMulBase& mat, MatrixCacheType buildCache)
{
  MatMulLock locking(mat, buildCache);

  // If locking.getType()!=cacheEmpty then we really do need to
  // build the cache, and we also have the lock for it.

  if (locking.getType() == cacheEmpty && ctxt==nullptr) //  nothing to do
    return;

  switch (mat.getEA().getTag()) {
    case PA_GF2_tag: {
      matmul_sparse_impl<PA_GF2> M(mat, locking.getType());
      M.multiply(ctxt);
      break;
    }
    case PA_zz_p_tag: {
      matmul_sparse_impl<PA_zz_p> M(mat, locking.getType());
      M.multiply(ctxt);
      break;
    }
    default:
      throw std::logic_error("mat_mul_sparse: neither PA_GF2 nor PA_zz_p");
  }
}
// Same as matMul but optimized for matrices with few non-zero diagonals
void matMul_sparse(Ctxt& ctxt, MatMulBase& mat,
                   MatrixCacheType buildCache)
{ mat_mul_sparse(&ctxt, mat, buildCache); }

// Build a cache without performing multiplication
void buildCache4MatMul_sparse(MatMulBase& mat, MatrixCacheType buildCache)
{ mat_mul_sparse(nullptr, mat, buildCache); }


/********************************************************************
 ********************************************************************/
// Applying matmul to plaintext, useful for debugging
template<class type> class matmul_pa_impl {
public:
  PA_INJECT(type)

  static void matmul(NewPlaintextArray& pa, MatMul<type>& mat)
  {
    const EncryptedArrayDerived<type>& ea = mat.getEA().getDerived(type());
    PA_BOILER

    vector<RX> res;
    res.resize(n);
    for (long j = 0; j < n; j++) {
      RX acc, val, tmp; 
      acc = 0;
      for (long i = 0; i < n; i++) {
        if (!mat.get(val, i, j)) {
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
// A wrapper around the implementation class
void matMul(NewPlaintextArray& pa, MatMulBase& mat)
{
  switch (mat.getEA().getTag()) {
    case PA_GF2_tag: {
      matmul_pa_impl<PA_GF2>::matmul(pa, dynamic_cast< MatMul<PA_GF2>& >(mat));
      return;
    }
    case PA_zz_p_tag: {
      matmul_pa_impl<PA_zz_p>::matmul(pa,dynamic_cast<MatMul<PA_zz_p>&>(mat));
      return;
    }
    default:
      throw std::logic_error("mat_mul: neither PA_GF2 nor PA_zz_p");
  }
}


#ifdef DEBUG
void printCache(const CachedzzxMatrix& cache)
{
  std::cerr << " zzxCache=[";
  for (long i=0; i<cache.length(); i++) {
    if (cache[i]==nullptr)
      std::cerr << "null ";
    else {
      std::cerr << (*(cache[i])) << " ";
    }
  }
  std:cerr << "]\n";
}

void printCache(const CachedDCRTMatrix& cache)
{
  std::cerr << "dcrtCache=[";
  for (long i=0; i<cache.length(); i++) {
    if (cache[i]==nullptr)
      std::cerr << "null ";
    else {
      ZZX poly;
      cache[i]->toPoly(poly);
      std::cerr << poly << " ";
    }
  }
  std:cerr << "]\n";
}
#endif
