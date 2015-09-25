#ifndef _EVAL_MAP_H_
#define _EVAL_MAP_H_
/* Copyright (C) 2012-2014 IBM Corp.
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
/** @file EvalMap.h
 *  @brief Implementing the recryption linear transformations
 */

#include "EncryptedArray.h"
#include "matrix.h"

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
//! The constructor precomputes certain values, but the linear
//! transformation itself is effected using the apply method.
//!
//! Note that the factorization in mvec must correspond to the
//! generators used in PAlgebra.  The best way to ensure this is
//! to used directly the output of the params program, which
//! will supply values for mvec (to be used here), and gens and ords
//! (to be used in initialize the FHEcontext).

class EvalMap {
private:
  const EncryptedArray& ea;
  bool invert;   // apply transformation in inverser order?
  long nfactors; // how many factors of m

#ifndef EVALMAP_CACHED    // no caching
  shared_ptr<PlaintextBlockMatrixBaseInterface>   mat1;   // one block matrix
  Vec< shared_ptr<PlaintextMatrixBaseInterface> > matvec; // regular matrices
#else
#if (EVALMAP_CACHED==0) // ZZX caching
  CachedPtxtBlockMatrix mat1;
  Vec<CachedPtxtMatrix> matvec;
#else               // DoubleCRT cashing
  CachedDCRTPtxtBlockMatrix mat1;
  Vec<CachedDCRTPtxtMatrix> matvec;
#endif
#endif
public:

  EvalMap(const EncryptedArray& _ea, const Vec<long>& mvec, bool _invert,
          bool normal_basis = true);

  // the normal_basis parameter indicates that we want the
  // normal basis transformation when invert == true.
  // On by default, off for testing

  void apply(Ctxt& ctxt) const;
};

#endif
