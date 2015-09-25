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
/** @file OldEvalMap.h
 *  @brief Old implementation of the recryption linear transformations
 */
#ifndef _OLD_EVAL_MAP_H_
#define _OLD_EVAL_MAP_H_

#include "FHE.h"
#include "Ctxt.h"
#include "EncryptedArray.h"
#include "matrix.h"
#include "permutations.h"

//! \cond FALSE (make doxygen ignore these classes)
class Step2aShuffleBase;
class TowerBase;


// This is an older version of EvalMap.
// It is much less depth efficient than the new one.

class OldEvalMap {
private:
  const EncryptedArray& ea;
  bool invert; 

  bool easy;  // easy => d1 == d, 
              // !ease => d1 != d (but we d1 * d2 == d)

  long nfactors;

  shared_ptr<PlaintextBlockMatrixBaseInterface> mat1;
    // use for both easy and !easy

  shared_ptr<Step2aShuffleBase> shuffle;
  shared_ptr<TowerBase> tower;
    // use only in the !easy case

  Vec< shared_ptr<PlaintextMatrixBaseInterface> > matvec;

  shared_ptr<PermNetwork> net;
  Vec<long> slot_rotate;
    // used for the initial/final inter- and intra-slot rotations

  
public:

  // mvec: the factorization of m
  // width: a bound on the width used for the permutation network
  // invert: false => forward eval [ coeffs packed in slots to coeffs ]
  //         true  => reverse eval [ coeffs to coeffs packed in slots ]

  // so for bootstrapping, we want the reverse eval, followed by digit
  // extraction, followed by forward eval

  OldEvalMap(const EncryptedArray& _ea, const Vec<long>& mvec, long width, 
          bool _invert);

  void apply(Ctxt& ctxt) const;
};

//! \endcond
#endif
