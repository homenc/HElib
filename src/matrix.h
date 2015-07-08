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
 * @file matrix.h
 * @brief some matrix / linear algenra stuff
 */

#include "EncryptedArray.h"


/********************************************************************/
/****************** Linear transformation classes *******************/

//! @class PlaintextMatrixBaseInterface
//! @brief An abstract interface for linear transformations.
//!
//! A matrix implements linear transformation over an extension field/ring
//! (e.g., GF(2^d) or Z_{2^8}[X]/G(X) for irreducible G). Any class
//! implementing this interface should be linked to a specific EncryptedArray
//! object, a reference to which is returned by the getEA() method -- this
//! method will generally be invoked by an EncryptedArray object to verify
//! consistent use.

class PlaintextMatrixBaseInterface {
public:
  virtual const EncryptedArray& getEA() const = 0;

  virtual ~PlaintextMatrixBaseInterface() {}
};


//! @class PlaintextMatrixInterface
//! @brief A somewhat less abstract interface for linear transformations.
//! 
//! A matrix implements linear transformation over an extension field/ring
//! (e.g., GF(2^d) or Z_{2^8}[X]/G(X) for irreducible G). The method
//! get(out, i, j) copies the element at row i column j of a matrix into
//! the variable out. The type of out is RX, which is GF2X if type is PA_GF2,
//! and zz_pX if type is PA_zz_p. A return value of true means that the
//! entry is zero, and out is not touched.
template<class type> 
class  PlaintextMatrixInterface : public PlaintextMatrixBaseInterface {
public:
  PA_INJECT(type)

  virtual bool get(RX& out, long i, long j) const = 0;
};



//! @class PlaintextBlockMatrixBaseInterface
//! @brief An abstract interface for linear transformations.
//!
//! A block matrix implements linear transformation over the base field/ring
//! (e.g., Z_2, Z_3, Z_{2^8}, etc.) Any class implementing this interface
//! should be linked to a specific EncryptedArray object, a reference to which
//! is returned by the getEA() method -- this method will generally be invoked
//! by an EncryptedArray object to verify consistent use.

class PlaintextBlockMatrixBaseInterface {
public:
  virtual const EncryptedArray& getEA() const = 0;

  virtual ~PlaintextBlockMatrixBaseInterface() {}
};


//! @class PlaintextBlockMatrixInterface
//! @brief A somewhat less abstract interface for linear transformations.
//! 
//! A block matrix implements linear transformation over the base field/ring
//! (e.g., Z_2, Z_3, Z_{2^8}, etc.) The method get(out, i, j) copies the
//! element at row i column j of a matrix into the variable out. The type
//! of out is mat_R (so either mar_GF2 or mat_zz_p.  A return value of true
//! means that the entry is zero, and out is not touched.

template<class type> 
class  PlaintextBlockMatrixInterface : public PlaintextBlockMatrixBaseInterface {
public:
  PA_INJECT(type)

  virtual bool get(mat_R& out, long i, long j) const = 0;
};


/**************** End linear transformation classes *****************/
/********************************************************************/








typedef Vec<ZZXptr> CachedPtxtMatrix;
typedef Mat<ZZXptr> CachedPtxtBlockMatrix;
typedef Vec<DCRTptr> CachedDCRTPtxtMatrix;
typedef Mat<DCRTptr> CachedDCRTPtxtBlockMatrix;


///@{
//! @name Matrix multiplication routines


//! @brief Multiply ctx by plaintext matrix. Ctxt is treated as
//! a row matrix v, and replaced by an encryption of v * mat.
//! Optimized for dense matrices
void mat_mul_dense(const EncryptedArray& ea, Ctxt& ctxt, 
  const PlaintextMatrixBaseInterface& mat);



//! @brief Multiply ctx by plaintext matrix. Ctxt is treated as
//! a row matrix v, and replaced by an encryption of v * mat.
//! Optimized for sparse diagonals
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const PlaintextMatrixBaseInterface& mat);
void compMat(const EncryptedArray& ea, CachedPtxtMatrix& cmat, 
  const PlaintextMatrixBaseInterface& mat);
void compMat(const EncryptedArray& ea, CachedDCRTPtxtMatrix& cmat, 
  const PlaintextMatrixBaseInterface& mat);





//! @brief Multiply ctx by plaintext block matrix (over the base field/ring).
//! Ctxt is treated as a row matrix v, and replaced by an encryption of v*mat.
//! Optimized for sparse diagonals
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const PlaintextBlockMatrixBaseInterface& mat);
void compMat(const EncryptedArray& ea, CachedPtxtBlockMatrix& cmat, 
  const PlaintextBlockMatrixBaseInterface& mat);
void compMat(const EncryptedArray& ea, CachedDCRTPtxtBlockMatrix& cmat, 
  const PlaintextBlockMatrixBaseInterface& mat);



//! @brief Multiply ctx by plaintext matrix.
//! Ctxt is treated as a row matrix v, and replaced by en encryption of
//! v * mat' where mat' is the block-diagonal matrix defined by mat in
//! dimension dim. Here, mat should represent a D x D matrix, where D is
//! the order of generator dim.
//! We also allow dim to be one greater than the number of generators in
//! zMStar, as if there were an implicit generator of order 1, this is
//! convenient in some applications.
void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
  const PlaintextMatrixBaseInterface& mat, long dim); 
void compMat1D(const EncryptedArray& ea, CachedPtxtMatrix& cmat,
  const PlaintextMatrixBaseInterface& mat, long dim); 
void compMat1D(const EncryptedArray& ea, CachedDCRTPtxtMatrix& cmat, 
  const PlaintextMatrixBaseInterface& mat, long dim);



void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
  const PlaintextBlockMatrixBaseInterface& mat, long dim);
void compMat1D(const EncryptedArray& ea, CachedPtxtBlockMatrix& cmat,
  const PlaintextBlockMatrixBaseInterface& mat, long dim); 
void compMat1D(const EncryptedArray& ea, CachedDCRTPtxtBlockMatrix& cmat, 
  const PlaintextBlockMatrixBaseInterface& mat, long dim);







//! @brief Free functions for various flavors of cached matrix multiplication
//!
//! To save time, the constants that need to be computed during a matrix-vector
//! multiply can be precomputed and stored.  One can store these either
//! as polynomials (ZZX's) using a CachedPtxtMatrix, or as DoubleCRT's
//! using a CachedDCRTMatrix.  These caches are computed using EncryptedArray
//! methods.


//! @brief Functions corresponding to EncryptedArray::mat_mul
//! Caches are computed using EncryptedArray::compMat
//! These functions are designed to work with sparse matrices.
//! The first two work with matrices over Z_{p^r}[X]/(G).
//! The second two work with "block" matrices whose entries are themselves
//! matrices over Z_{p^r}
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedPtxtMatrix& cmat);
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedDCRTPtxtMatrix& cmat);

void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedPtxtBlockMatrix& cmat);
void mat_mul(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedDCRTPtxtBlockMatrix& cmat);

//! @brief Functions corresponding to EncryptedArray::mat_mul1D
//! Caches are computed using EncryptedArray::compMat1D
//! These functions apply a matrix concurrently to all the hypercolumns
//! in a single dimension.
//! The first two work with matrices over Z_{p^r}[X]/(G).
//! The second two work with "block" matrices whose entries are themselves
//! matrices over Z_{p^r}
void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedPtxtMatrix& cmat, long dim);

void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedDCRTPtxtMatrix& cmat, long dim);

void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedPtxtBlockMatrix& cmat, long dim);

void mat_mul1D(const EncryptedArray& ea, Ctxt& ctxt, 
  const CachedDCRTPtxtBlockMatrix& cmat, long dim);




///@}



void mat_mul(const EncryptedArray& ea, NewPlaintextArray& pa, 
  const PlaintextMatrixBaseInterface& mat);
void mat_mul(const EncryptedArray& ea, NewPlaintextArray& pa, 
  const PlaintextBlockMatrixBaseInterface& mat);




#endif /* ifdef FHE_matrix_H_ */

