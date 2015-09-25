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
#ifndef _replicate_H_
#define _replicate_H_
/**
 * @file replicate.h
 * @brief Procedures for replicating a ciphertext slot across a full ciphertext
 *
 * This module implements a recursive, O(1)-amortized algorithm for
 * replications. On an input ciphertext that encrypts (x_1, ..., x_n), we
 * generate the n encrypted vectors (x_1, ..., x_1), ..., (x_n, ..., x_n),
 * in that order.
 *
 * To process the output vectors, a "call back" mechanism is used (so that we
 * do not need to generate them all, and instead can return them one by one).
 * For this purpose, the caller should pass a pointer to a class derived from
 * the purely abstract class ReplicateHandler.
 *
 * The replication procedures are meant to be used for linear algebra operation
 * where a matrix-vector multiplication can be implemented for example by
 * replicating each entry of the vector as a stand-alone ciphertext, then use
 * the SIMD operations on these ciphertexts.
 **/

#include "FHE.h"
#include "EncryptedArray.h"


// set to true to see some more info... 
NTL_THREAD_LOCAL
extern bool replicateVerboseFlag;


//! @brief The value in slot #pos is replicated in all other slots.
//! On an n-slot ciphertext, this algorithm performs O(log n) 1D rotations.  
void replicate(const EncryptedArray& ea, Ctxt& ctx, long pos);

//! @brief A lower-level routine. Same as replicate, but assumes
//! all slots are zero except slot #pos.
void replicate0(const EncryptedArray& ea, Ctxt& ctxt, long pos);

//! A virtual class to handle call-backs to get the output of replicate
class ReplicateHandler {
public:
  virtual void handle(const Ctxt& ctxt) = 0;
  virtual ~ReplicateHandler() {};
};

/**
 * replicateAll uses a hybrid strategy, combining the O(log n) strategy of the
 * replicate method, with an O(1) strategy, which is faster but introduces more
 * noise. This tradeoff is controlled by the parameter recBound:
 *
 * \li recBound < 0: recursion to depth |recBound| (faster, noisier)
 * \li recBound ==0: no recursion (slower, less noise)
 * \li recBound > 0: the recursion depth is chosen heuristically,
 *     but is capped at recBound
 * 
 * The default value for recBound is 64, this ensures that the choice is
 * based only on the heuristic, which will introduce noise corresponding to
 * O(log log n) levels of recursion, but still gives an algorithm that
 * theoretically runs in time O(n).
 **/
void replicateAll(const EncryptedArray& ea, const Ctxt& ctxt, 
                         ReplicateHandler *handler, long recBound = 64);


//! This function is obsolete, and is kept for historical purposes only. It
//! was a first attempt at implementing the O(1)-amortized algorithm, but is
//! less efficient than the function above.
void replicateAllOrig(const EncryptedArray& ea, const Ctxt& ctxt,
                      ReplicateHandler *handler);



void replicate(const EncryptedArray& ea, NewPlaintextArray& pa, long i);

#endif
