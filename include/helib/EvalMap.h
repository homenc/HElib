/* Copyright (C) 2012-2020 IBM Corp.
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
#ifndef HELIB_EVALMAP_H
#define HELIB_EVALMAP_H
/** @file EvalMap.h
 *  @brief Implementing the recryption linear transformations
 */

#include <helib/EncryptedArray.h>
#include <helib/matmul.h>

namespace helib {

//! @class EvalMap
//! @brief Class that provides the functionality for the
//! linear transforms used in boostrapping.
//! The constructor is invoked with three arguments:
//!   - an EncryptedArray object ea
//!   - an integer vector mvec
//!   - a boolean flag invert
//! The mvec vector specifies the factorization of m
//! to use in the "powerful basis" decomposition.
//!
//! If the invert flag is false, the forward transformation is
//! used.  This transformation views the slots as being packed
//! with powerful-basis coefficients and performs a multi-point
//! polynomial evaluation.  This is the second transformation
//! used in bootstrapping.
//!
//! If invert flag is true, the inverse transformation is used.
//! In addition, the current implementation folds into the
//! inverse transformation a transformation that moves
//! the coefficients in each slot into a normal-basis representation,
//! which helps with the unpacking procedure.
//!
//! The constructor precomputes certain values, and the linear
//! transformation itself is effected using the apply method.
//!
//! Note that the factorization in mvec must correspond to the
//! generators used in PAlgebra.  The best way to ensure this is
//! to directly use the output of the program in  params.cpp: that
//! program computes values for mvec (to be used here), and gens
//! and ords (to be used in initialization of the Context).

class EvalMap
{
private:
  const EncryptedArray& ea;
  bool invert;   // apply transformation in inverse order?
  long nfactors; // how many factors of m
  std::unique_ptr<BlockMatMul1DExec> mat1;        // one block matrix
  NTL::Vec<std::unique_ptr<MatMul1DExec>> matvec; // regular matrices

public:
  EvalMap(const EncryptedArray& _ea,
          bool minimal,
          const NTL::Vec<long>& mvec,
          bool _invert,
          bool build_cache,
          bool normal_basis = true);

  // the normal_basis parameter indicates that we want the
  // normal basis transformation when invert == true.
  // On by default, off for testing

  void upgrade();
  void apply(Ctxt& ctxt) const;
};

//! @class ThinEvalMap
//! @brief Class that provides the functionality for the
//! linear transforms used in "thin" boostrapping,
//! where slots are assumed to contain constants.
//! The interface is exactly the same as for EvalMap,
//! except that the constructor does not have a normal_basis
//! parameter.

class ThinEvalMap
{
private:
  const EncryptedArray& ea;
  bool invert;   // apply transformation in inverse order?
  long nfactors; // how many factors of m
  NTL::Vec<std::unique_ptr<MatMulExecBase>> matvec; // regular matrices

public:
  ThinEvalMap(const EncryptedArray& _ea,
              bool minimal,
              const NTL::Vec<long>& mvec,
              bool _invert,
              bool build_cache);

  void upgrade();
  void apply(Ctxt& ctxt) const;
};

} // namespace helib

#endif // ifndef HELIB_EVALMAP_H
