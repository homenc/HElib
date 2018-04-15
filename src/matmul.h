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
#ifndef _matmul_H
#define _matmul_H

#include "EncryptedArray.h"


class MatMulFullExec;

// Abstract base class for representing a linear transformation on a full vector.
class MatMulFull {
public:
  virtual ~MatMulFull() {}
  virtual const EncryptedArray& getEA() const = 0; 
  typedef MatMulFullExec ExecType;
  
};
  
// Concrete derived class that defines the matrix entries.
template<class type>
class MatMulFull_derived : public MatMulFull { 
public:
  PA_INJECT(type)

  // Get (i, j) entry of matrix.
  // Should return true when the entry is a zero. 
  virtual bool get(RX& out, long i, long j) const = 0;
};

//====================================

class BlockMatMulFullExec;

// Abstract base class for representing a block linear transformation on a full vector.
class BlockMatMulFull {
public:
  virtual ~BlockMatMulFull() {}
  virtual const EncryptedArray& getEA() const = 0;
  typedef BlockMatMulFullExec ExecType;

};
  
// Concrete derived class that defines the matrix entries.
template<class type>
class BlockMatMulFull_derived : public BlockMatMulFull { 
public:
  PA_INJECT(type)

  // Get (i, j) entry of matrix.
  // Each entry is a d x d matrix over the base ring.
  // Should return true when the entry is a zero. 
  virtual bool get(mat_R& out, long i, long j) const = 0;
};

//====================================

class MatMul1DExec;

// Abstract base class for representing a 1D linear transformation.
class MatMul1D {
public:
  virtual ~MatMul1D() {}
  virtual const EncryptedArray& getEA() const = 0;
  virtual long getDim() const = 0;
  typedef MatMul1DExec ExecType;
};

// An intermediate class that is mainly intended for internal use.
template<class type>
class MatMul1D_partial : public MatMul1D {
public:
  PA_INJECT(type)

  // Get the i'th diagonal, encoded as a single constant.
  // MatMul1D_derived (below) supplies a default implementation,
  // which can be overriden in special circumstances.
  virtual void 
  processDiagonal(RX& poly, long i,
                  const EncryptedArrayDerived<type>& ea) const = 0;

};

// Concrete derived class that defines the matrix entries.
template<class type>
class MatMul1D_derived : public MatMul1D_partial<type> { 
public:
  PA_INJECT(type)

  // Should return true if their are multiple (different) transforms
  // among the various components.
  virtual bool multipleTransforms() const = 0;

  // Get coordinate (i, j) of the kth component.
  // Should return true when the entry is a zero. 
  virtual bool get(RX& out, long i, long j, long k) const = 0;

  void 
  processDiagonal(RX& poly, long i,
                  const EncryptedArrayDerived<type>& ea) const override;
};

//====================================

class BlockMatMul1DExec;

// Abstract base class for representing a block 1D linear transformation.
class BlockMatMul1D {
public:
  virtual ~BlockMatMul1D() {}
  virtual const EncryptedArray& getEA() const = 0;
  virtual long getDim() const = 0;
  typedef BlockMatMul1DExec ExecType;
};


// An intermediate class that is mainly intended for internal use.
template<class type>
class BlockMatMul1D_partial : public BlockMatMul1D {
public:
  PA_INJECT(type)

  // Get the i'th diagonal, encoded as a vector of d constants,
  // where d is the order of p.
  // BlockMatMul1D_derived (below) supplies a default implementation,
  // which can be overriden in special circumstances.
  virtual bool
  processDiagonal(vector<RX>& poly, long i,
                  const EncryptedArrayDerived<type>& ea) const = 0;

};

// Concrete derived class that defines the matrix entries.
template<class type>
class BlockMatMul1D_derived : public BlockMatMul1D_partial<type> { 
public:
  PA_INJECT(type)

  // Should return true if their are multiple (different) transforms
  // among the various components.
  virtual bool multipleTransforms() const = 0;

  // Get coordinate (i, j) of the kth component.
  // Each entry is a d x d matrix over the base ring.
  // Should return true when the entry is a zero. 
  virtual bool get(mat_R& out, long i, long j, long k) const = 0;

  bool
  processDiagonal(vector<RX>& poly, long i,
                  const EncryptedArrayDerived<type>& ea) const override;
};

//====================================


struct ConstMultiplier; 
// Defined in matmul.cpp.
// Holds a constant by which a ciphertext can be multiplied.
// Internally, it is represented as either zzX or a DoubleCRT.
// The former occupies less space, but the latter makes for
// much faster multiplication.

struct ConstMultiplierCache {
  std::vector<std::shared_ptr<ConstMultiplier>> multiplier;

  // Upgrade zzX constants to DoubleCRT constants.
  void upgrade(const FHEcontext& context);
};

//====================================


// Abstract base case for multiplying an encrypted vector by a plaintext matrix.
class MatMulExecBase {
public:
  virtual ~MatMulExecBase() { }

  virtual const EncryptedArray& getEA() const = 0;

  // Upgrade zzX constants to DoubleCRT constants.
  virtual void upgrade() = 0;

  // If ctxt enctrypts a row vector v, then this replaces ctxt
  // by an encryption of the row vector v*mat, where mat is 
  // a matrix provided to the constructor of one of the
  // concrete subclasses MatMul1DExec, BlockMatMul1DExec,
  // MatMulFullExec, BlockMatMulFullExec, defined below.
  virtual void mul(Ctxt& ctxt) const = 0;
};

//====================================

// Class used to multiply an encrypted row vector by a 1D linear transformation.
class MatMul1DExec : public MatMulExecBase {
public:

  const EncryptedArray& ea;

  long dim;
  long D;
  bool native;
  bool minimal;
  long g;

  ConstMultiplierCache cache;
  ConstMultiplierCache cache1; // only for non-native dimension


  // The constructor encodes all the constants for a given
  // matrix in zzX format.
  // The mat argument defines the entries of the matrix.
  // Use the upgrade method (below) to convert to DoubleCRT format.
  // If the minimal flag is set to true, a strategy that relies
  // on a minimal number of key switching matrices will be used;
  // this is intended for use in conjunction with the 
  // addMinimal{1D,Frb}Matrices routines decalred in FHE.h.
  // If the minimal flag is false, it is best to use the
  // addSome{1D,Frb}Matrices routines declared in FHE.h.
  explicit
  MatMul1DExec(const MatMul1D& mat, bool minimal=false);

  // Replaces an encryption of row vector v by encryption of v*mat
  void mul(Ctxt& ctxt) const override;

  // Upgrades encoded constants from zzX to DoubleCRT.
  void upgrade() override { 
    cache.upgrade(ea.getContext()); 
    cache1.upgrade(ea.getContext()); 
  }

  const EncryptedArray& getEA() const override { return ea; }
};

//====================================

// Class used to multiply an encrypted row vector by a block 1D linear transformation.
class BlockMatMul1DExec : public MatMulExecBase {
public:

  const EncryptedArray& ea;

  long dim;
  long D;
  long d;
  bool native;
  long strategy;

  ConstMultiplierCache cache;
  ConstMultiplierCache cache1; // only for non-native dimension


  // The constructor encodes all the constants for a given
  // matrix in zzX format.
  // The mat argument defines the entries of the matrix.
  // Use the upgrade method (below) to convert to DoubleCRT format.
  // If the minimal flag is set to true, a strategy that relies
  // on a minimal number of key switching matrices will be used;
  // this is intended for use in conjunction with the 
  // addMinimal{1D,Frb}Matrices routines decalred in FHE.h.
  // If the minimal flag is false, it is best to use the
  // addSome{1D,Frb}Matrices routines declared in FHE.h.
  explicit
  BlockMatMul1DExec(const BlockMatMul1D& mat, bool minimal=false);

  // Replaces an encryption of row vector v by encryption of v*mat
  void mul(Ctxt& ctxt) const override;

  // Upgrades encoded constants from zzX to DoubleCRT.
  void upgrade() override { 
    cache.upgrade(ea.getContext()); 
    cache1.upgrade(ea.getContext()); 
  }

  const EncryptedArray& getEA() const override { return ea; }
};

//====================================

// Class used to multiply an encrypted row vector by a full linear transformation.
class MatMulFullExec : public MatMulExecBase {
public:

  const EncryptedArray& ea;
  bool minimal;
  std::vector<long> dims;
  std::vector<MatMul1DExec> transforms;

  // The constructor encodes all the constants for a given
  // matrix in zzX format.
  // The mat argument defines the entries of the matrix.
  // Use the upgrade method (below) to convert to DoubleCRT format.
  // If the minimal flag is set to true, a strategy that relies
  // on a minimal number of key switching matrices will be used;
  // this is intended for use in conjunction with the 
  // addMinimal{1D,Frb}Matrices routines decalred in FHE.h.
  // If the minimal flag is false, it is best to use the
  // addSome{1D,Frb}Matrices routines declared in FHE.h.
  explicit
  MatMulFullExec(const MatMulFull& mat, bool minimal=false);

  // Replaces an encryption of row vector v by encryption of v*mat
  void mul(Ctxt& ctxt) const override;

  // Upgrades encoded constants from zzX to DoubleCRT.
  void upgrade() override { 
    for (auto& t: transforms) t.upgrade();
  }

  const EncryptedArray& getEA() const override { return ea; }

  // This really should be private.
  long rec_mul(Ctxt& acc, const Ctxt& ctxt, long dim, long idx) const;

};

//====================================

// Class used to multiply an encrypted row vector by a full block linear transformation.
class BlockMatMulFullExec : public MatMulExecBase {
public:

  const EncryptedArray& ea;
  bool minimal;
  std::vector<long> dims;
  std::vector<BlockMatMul1DExec> transforms;

  // The constructor encodes all the constants for a given
  // matrix in zzX format.
  // The mat argument defines the entries of the matrix.
  // Use the upgrade method (below) to convert to DoubleCRT format.
  // If the minimal flag is set to true, a strategy that relies
  // on a minimal number of key switching matrices will be used;
  // this is intended for use in conjunction with the 
  // addMinimal{1D,Frb}Matrices routines decalred in FHE.h.
  // If the minimal flag is false, it is best to use the
  // addSome{1D,Frb}Matrices routines declared in FHE.h.
  explicit
  BlockMatMulFullExec(const BlockMatMulFull& mat, bool minimal=false);

  // Replaces an encryption of row vector v by encryption of v*mat
  void mul(Ctxt& ctxt) const override;

  // Upgrades encoded constants from zzX to DoubleCRT.
  void upgrade() override { 
    for (auto& t: transforms) t.upgrade();
  }

  const EncryptedArray& getEA() const override { return ea; }

  // This really should be private.
  long rec_mul(Ctxt& acc, const Ctxt& ctxt, long dim, long idx) const;

};

//===================================

// ctxt = \sum_{i=0}^{d-1} \sigma^i(ctxt),
//   where d = order of p mod m, and \sigma is the Frobenius map

void traceMap(Ctxt& ctxt);

//====================================

// These routines apply linear transformation to plaintext arrays.
// Mainly for testing purposes.
void mul(NewPlaintextArray& pa, const MatMul1D& mat);
void mul(NewPlaintextArray& pa, const BlockMatMul1D& mat);
void mul(NewPlaintextArray& pa, const MatMulFull& mat);
void mul(NewPlaintextArray& pa, const BlockMatMulFull& mat);


// These are used mainly for performance evaluation.

extern int fhe_test_force_bsgs;
// Controls whether or not we use BSGS multiplication.
// 1 to force on, -1 to force off, 0 for default behaviour.


extern int fhe_test_force_hoist;
// Controls whether ot not we use hoisting.
// -1 to force off, 0 for default behaviour.

#endif
