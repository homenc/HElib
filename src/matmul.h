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
#ifndef FHE_matrix_H_
#define FHE_matrix_H_
/**
 * @file matmul.h
 * @brief some matrix / linear algenra stuff
 */
#include <tuple>
#include "EncryptedArray.h"
#include "CachedConstants.h"

/********************************************************************/
/****************** Linear transformation classes *******************/

//! @class PlaintextBlockMatrixBaseInterface
//! @brief An abstract interface for linear transformations.
//!
//! A block matrix implements linear transformation over the base field/ring
//! (e.g., Z_2, Z_3, Z_{2^8}, etc.) Any class implementing this interface
//! should be linked to a specific EncryptedArray object, a reference to which
//! is returned by the getEA() method -- this method will generally be invoked
//! by an EncryptedArray object to verify consistent use.

class PlaintextBlockMatrixBaseInterface {
  CachedConstants cache;
public:
  virtual const EncryptedArray& getEA() const = 0;
  virtual ~PlaintextBlockMatrixBaseInterface() {}

  CachedConstants& getCache() { return cache; }
};


//! @class PlaintextBlockMatrixInterface
//! @brief A somewhat less abstract interface for linear transformations.
//! 
//! A block matrix implements linear transformation over the base field/ring
//! (e.g., Z_2, Z_3, Z_{2^8}, etc.) An object of this type specifies a
//! particular linear transformation, represented as a matrix of small
//! matrices. The method get(out,i,j) copies the small matrix ar row i,
//! column j into the variable out. The type of out is mat_R (so either
//! mat_GF2 or mat_zz_p). The dimenssion of our is d-by-d, where d is
//! the extension degree of the slots (per the EncryptedArray object).
//! A return value of true means that the entry is zero, and out is not
//! touched.

template<class type> 
class  PlaintextBlockMatrixInterface : public PlaintextBlockMatrixBaseInterface {
public:
  PA_INJECT(type)

  virtual bool get(mat_R& out, long i, long j) const = 0;
};



/**************** End linear transformation classes *****************/
/********************************************************************/


typedef NTL::Vec<ZZXptr> CachedPtxtMatrix;
typedef NTL::Mat<ZZXptr> CachedPtxtBlockMatrix;
typedef NTL::Vec<DCRTptr> CachedDCRTPtxtMatrix;
typedef NTL::Mat<DCRTptr> CachedDCRTPtxtBlockMatrix;
// ZZXptr is std::shared_ptr<NTL::ZZX>
// DCRTptr is std::shared_ptr<DoubleCRT>


///@{
//! @name Matrix multiplication routines

/*
//! @brief Multiply ctx by plaintext matrix. Ctxt is treated as
//! a row matrix v, and replaced by an encryption of v * mat.
//! Optimized for sparse diagonals
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const PlaintextMatrixBaseInterface& mat);
void compMat(const EncryptedArray& ea, CachedPtxtMatrix& cmat, 
  const PlaintextMatrixBaseInterface& mat);
void compMat(const EncryptedArray& ea, CachedDCRTPtxtMatrix& cmat, 
  const PlaintextMatrixBaseInterface& mat);
// The latter two functions used a "cached" version of the matrix,
// which was built from the PlaintextMatrixBaseInterface object,
// where all the constants are already encoded as ZZX or DoubleCRT.
*/


//! @brief Multiply ctx by plaintext block matrix (over the base field/ring).
//! Ctxt is treated as a row matrix v, and replaced by an encryption of v*mat.
//! Optimized for sparse diagonals
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const PlaintextBlockMatrixBaseInterface& mat);
void compMat(const EncryptedArray& ea, CachedPtxtBlockMatrix& cmat, 
  const PlaintextBlockMatrixBaseInterface& mat);
void compMat(const EncryptedArray& ea, CachedDCRTPtxtBlockMatrix& cmat, 
  const PlaintextBlockMatrixBaseInterface& mat);
// The latter two functions used a "cached" version of the matrix,
// which was built from the PlaintextMatrixBaseInterface object,
// where all the constants are already encoded as ZZX or DoubleCRT.

///@}


//! @brief Functions for various flavors of cached matrix multiplication
//!
//! To save time, the constants that need to be computed during a matrix-vector
//! multiply can be precomputed and stored.  One can store these either
//! as polynomials (ZZX's) using a CachedPtxtMatrix, or as DoubleCRT's
//! using a CachedDCRTMatrix. These caches are computed using EncryptedArray
//! methods.


//! @brief Functions corresponding to EncryptedArray::mat_mul
//! Caches are computed using EncryptedArray::compMat
//! These functions are designed to work with sparse matrices.
//! The first two work with matrices over Z_{p^r}[X]/(G).
//! The second two work with "block" matrices whose entries are themselves
//! matrices over Z_{p^r}
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedPtxtBlockMatrix& cmat);
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedDCRTPtxtBlockMatrix& cmat);

///@}

// Versions of the matrix-vector functions that work on plaintext
// rather than cipehrtext, useful for debugging.
void mat_mul(const EncryptedArray& ea, NewPlaintextArray& pa, 
  const PlaintextBlockMatrixBaseInterface& mat);

/********************************************************************/
/********************************************************************/
/********************************************************************/
/********************************************************************/
#include <cstddef>
#include <mutex>

typedef std::shared_ptr< zzX > zzxptr; // zzX=Vec<long> defined in NumbTh.h
typedef NTL::Vec<zzxptr> CachedzzxMatrix;

enum MatrixCacheType : int { cacheEmpty=0, cachezzX=1, cacheDCRT=2 };

//! @class MatMulBase
//! @brief An abstract interface for linear transformations.
//!
//! A matrix implements linear transformation over an extension
//! field/ring, e.g., GF(2^d) or Z_{2^8}[X]/G(X) for irreducible G.
class MatMulBase {
  const EncryptedArray& ea;
  std::unique_ptr<CachedzzxMatrix> zzxCache;
  std::unique_ptr<CachedDCRTPtxtMatrix> dcrtCache;
  std::mutex cachelock;
public:
  MatMulBase(const EncryptedArray& _ea): ea(_ea) {}
  virtual ~MatMulBase() {}

  const EncryptedArray& getEA() const { return ea; }

  bool haszzxcache() const { return (bool)zzxCache; }   // check if not null
  bool hasDCRTcache() const { return (bool)dcrtCache; } // check if not null
  void getCache(CachedzzxMatrix** zcp, CachedDCRTPtxtMatrix** dcp) const
  { *zcp = zzxCache.get(); *dcp = dcrtCache.get(); }

  bool lockCache(MatrixCacheType ty);
  void releaseCache() { cachelock.unlock(); }

  void upgradeCache(); // build DCRT cache from zzx cache

  void installzzxcache(std::unique_ptr<CachedzzxMatrix>& zc)
  { zzxCache.swap(zc); }
  void installDCRTcache(std::unique_ptr<CachedDCRTPtxtMatrix>& dc)
  { dcrtCache.swap(dc); }
};

//! @class MatMul
//! @brief Linear transformation interfaces, specialized to GF2/zz_p
//!
//! An implementation that specifies particular transformation(s) is
//! derived from MatMul<PA_GF2> or MatMul<PA_zz_p>. It should implement
//! either get(i,j) or multiGet(i,j,k). They return the element mat[i,j])
//! (in the k'th trasformation in the latter case), which is in some
//! field/ring (e.g., GF(p^d)), and is returned as a GF2X or a zz_pX.
template<class type>
class MatMul : public MatMulBase { // type is PA_GF2 or PA_zz_p
public:
  PA_INJECT(type)
  MatMul(const EncryptedArray& _ea): MatMulBase(_ea) {}

  // Should return true when the entry is a zero. An application must
  // implement (at least) one of these get functions, calling the base
  // methods below will throw an exception.
  virtual bool get(RX& out, long i, long j) const {
    throw std::logic_error("MatMul: get() not implemented in base class");
    return true;
  }
  virtual bool multiGet(RX& out, long i, long j, long k) const {
    throw std::logic_error("MatMul: multiGet() not implemented in base class");
    return true;
  }
};

//! @brief Multiply ctx by plaintext matrix. Ctxt is treated as
//! a row matrix v, and replaced by an encryption of v * mat.
//! If buildCache != cacheEmpty and the cache is not available,
//! then it will be built (However, a zzx cahce is never built
//! if the dcrt cache exists).
void matMul(Ctxt& ctxt, MatMulBase& mat,
            MatrixCacheType buildCache=cacheEmpty);

//! Build a cache without performing multiplication
void buildCache4MatMul(MatMulBase& mat, MatrixCacheType buildCache);


//! Same as mat_mul but optimized for matrices with few non-zero diagonals
void matMul_sparse(Ctxt& ctxt, MatMulBase& mat,
                   MatrixCacheType buildCache=cacheEmpty);

//! Build a cache without performing multiplication
void buildCache4MatMul_sparse(MatMulBase& mat, MatrixCacheType buildCache);

//FIXME: With the interfaces above, an application can call buildCache4MatMul
// and then use the cache with matMul_sparse (or vise versa), and currently
// there is no run-time check to detect that we have the wrong cache.


///@{
//! @name 1D Matrix multiplication routines
//! A single ciphertext holds many vectors, all of length equal to the
//! the size of the relevant dimenssion. Each vector is multiplied by
//! a potentially different matrix, all products done in SIMD.


//! @brief Multiply ctx by plaintext matrix. Ctxt is treated as a row
//! matrix v, and replaced by an encryption of v * mat, where mat is a
//! D x D matrix (where D is the order of generator dim).  We allow
//! dim to be one greater than the number of generators in zMStar, as
//! if there were an implicit generator of order 1, this is convenient
//! in some applications.
//! If buildCache != cacheEmpty and the cache is not available, then
//! it will be built (However, a zzx cahce is never built if the dcrt
//! cache exists).

void matMul1D(Ctxt& ctxt, MatMulBase& mat, long dim,
               MatrixCacheType buildCache=cacheEmpty);
//! Build a cache without performing multiplication
void buildCache4MatMul1D(MatMulBase& mat, long dim,MatrixCacheType buildCache);

void matMulti1D(Ctxt& ctxt, MatMulBase& mat, long dim,
               MatrixCacheType buildCache=cacheEmpty);
//! Build a cache without performing multiplication
void buildCache4MatMulti1D(MatMulBase& mat,long dim,MatrixCacheType buildCache);

// Versions for plaintext rather than cipehrtext, useful for debugging
void matMul(NewPlaintextArray& pa, MatMulBase& mat);
void matMul1D(NewPlaintextArray& pa, MatMulBase& mat, long dim);
void matMulti1D(NewPlaintextArray& pa, MatMulBase& mat, long dim);
///@}


//! @class BlockMatMul
//! @brief Linear transformation interfaces, specialized to GF2/zz_p
//!
//! An implementation that specifies particular transformation(s) is
//! derived from BlockMatMul<PA_GF2> or BlockMatMul<PA_zz_p>. It should
//! implement either get(i,j) or multiGet(i,j,k). They return the element
//! mat[i,j]) (in the k'th trasformation in the latter case), which is a
//! small matrix over the base ring (e.g., Z_p), and is returned as a
//! NTL::Mat<GF2> or NTL::Mat<zz_p>.
template<class type>
class BlockMatMul : public MatMulBase { // type is PA_GF2 or PA_zz_p
public:
  PA_INJECT(type)
  BlockMatMul(const EncryptedArray& _ea): MatMulBase(_ea) {}

  // Should return true when the entry is a zero. An application must
  // implement (at least) one of these get functions, calling the base
  // methods below will throw an exception.
  virtual bool get(mat_R& out, long i, long j) const {
    throw std::logic_error("MatMul: get() not implemented in base class");
    return true;
  }
  virtual bool multiGet(mat_R& out, long i, long j, long k) const {
    throw std::logic_error("MatMul: multiGet() not implemented in base class");
    return true;
  }
};

//! @brief Multiply ctx by plaintext matrix. Ctxt is treated as
//! a row matrix v, and replaced by an encryption of v * mat.
//! If buildCache != cacheEmpty and the cache is not available,
//! then it will be built (However, a zzx cahce is never built
//! if the dcrt cache exists).
void blockMatMul(Ctxt& ctxt, MatMulBase& mat,
		 MatrixCacheType buildCache=cacheEmpty);

//! Build a cache without performing multiplication
void buildCache4BlockMatMul(MatMulBase& mat, MatrixCacheType buildCache);


///@{
//! @name 1D block Matrix multiplication routines
//! A single ciphertext holds many vectors, all of length equal to the
//! the size of the relevant dimenssion (times slot size). Each vector is
//! multiplied by a potentially different matrix, all products done in SIMD.


//! @brief Multiply ctx by plaintext matrix. Ctxt is treated as a row
//! matrix v, and replaced by an encryption of v * mat, where mat is a
//! Dd x Dd matrix (where D is the order of generator dim, and d is the
//! size of slots). We allow dim to be one greater than the number of
//! generators in zMStar, as if there were an implicit generator of order 1,
//! this is convenient in some applications.
//! If buildCache != cacheEmpty and the cache is not available, then
//! it will be built (However, a zzx cahce is never built if the dcrt
//! cache exists).

void blockMatMul1D(Ctxt& ctxt, MatMulBase& mat, long dim,
               MatrixCacheType buildCache=cacheEmpty);
//! Build a cache without performing multiplication
void buildCache4BlockMatMul1D(MatMulBase& mat,
			      long dim, MatrixCacheType buildCache);

void blockMatMulti1D(Ctxt& ctxt, MatMulBase& mat, long dim,
		     MatrixCacheType buildCache=cacheEmpty);
//! Build a cache without performing multiplication
void buildCache4BlockMatMulti1D(MatMulBase& mat,
				long dim, MatrixCacheType buildCache);

// Versions for plaintext rather than cipehrtext, useful for debugging
void blockMatMul(NewPlaintextArray& pa, MatMulBase& mat);
void blockMatMul1D(NewPlaintextArray& pa, MatMulBase& mat, long dim);
void blockMatMulti1D(NewPlaintextArray& pa, MatMulBase& mat, long dim);
///@}


/*********************************************************************
 * MatMulLock: A helper class that handles the lock in MatMulBase.
 * It is used in a function as follows:
 *
 *   void some_Function(...)
 *   {
 *       MatMulLock locking(mat, cachetype);
 *       if (locking.getType()!=cacheEmpty) {
 *           // we need to build a cache and we have the lock
 *           ...
 *       }
 *   } // Lock is released by MatMulLock destructor
 *
 ********************************************************************/
class MatMulLock {
  MatMulBase& mat;
  MatrixCacheType ctype;
public:
  MatMulLock(MatMulBase& _mat, MatrixCacheType _ctype): mat(_mat),ctype(_ctype)
  {
    if (ctype != cacheEmpty) { // build a cache if it is not there already
      if (!mat.lockCache(ctype))
	ctype = cacheEmpty; // no need to build
    }
  }
  MatrixCacheType getType() const { return ctype; }
  ~MatMulLock() { if (ctype != cacheEmpty) mat.releaseCache(); }
};

#endif /* ifdef FHE_matrix_H_ */
