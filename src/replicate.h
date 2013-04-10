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

#include "FHE.h"
#include "EncryptedArray.h"


// set to true to see some more info... 
extern bool replicateVerboseFlag;


// The value in slot #pos is replicated in all
// other slots.  If there are n slots, this algorithm performs
// O(log n) 1D rotations.  

void replicate(const EncryptedArray& ea, Ctxt& ctx, long pos);



// A lower-level routine. Same as above, but assumes
// all slots are zero except slot #pos.

void replicate0(const EncryptedArray& ea, Ctxt& ctxt, long pos);



// The following code implements a recursive, O(1)-amortized
// algorithm for replications.  

// If the given ciphertext
// encrypts (x_1, ..., x_n), the algorithm will generate
// the n encrypted vectors (x_1, ..., x_1), ..., (x_n, ..., x_n),
// in that order.

// To process each such vector, a "call back" mechanism is used.
// For this purpose, the caller should pass a pointer to
// a class derived from the purely abstract class ReplicateHandler.


class ReplicateHandler {
public:
  virtual void handle(const Ctxt& ctxt) = 0;
  virtual ~ReplicateHandler() {};
};



// The algorithm, in fact, uses a hybrid strategy, combining
// the O(log n) strategy of the replicate method, with
// the O(1) strategy, which is faster, but introduces more noise. 

// The parameter recBound controls this tradeoff:

// recBound < 0 => recursion to depth |recBound| (faster, noisier)
// recBound == 0 => no recursion (slower, less noise)
// otherwise, a recursion depth is chosen heuristically,
//   but is capped at recBound

// The default value of 64 ensures that the the choice
// is based only on the heuristic, which will introduce
// noise corresponding to O(log log n) levels of recursion,
// but still gives an algorithm that theoretically runs in
// time O(n).



void replicateAll(const EncryptedArray& ea, const Ctxt& ctxt, 
                         ReplicateHandler *handler, long recBound = 64);


// This is really only here for historical purposes.  It was a 
// first attempt at implementing the O(1)-amortized algorithm,
// but is less efficient than the function above.

void replicateAllOrig(const EncryptedArray& ea, const Ctxt& ctxt,
                      ReplicateHandler *handler);


#endif
