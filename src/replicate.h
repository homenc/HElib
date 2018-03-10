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

class RepAux; // forward declerations
class RepAuxDim;

//! @brief The value in slot #pos is replicated in all other slots.
//! On an n-slot ciphertext, this algorithm performs O(log n) 1D rotations.  
void replicate(const EncryptedArray& ea, Ctxt& ctx, long pos);

//! @brief A lower-level routine. Same as replicate, but assumes
//! all slots are zero except slot #pos.
void replicate0(const EncryptedArray& ea, Ctxt& ctxt, long pos);

//! A virtual class to handle call-backs to get the output of replicate.
class ReplicateHandler {
public:
  virtual void handle(const Ctxt& ctxt) {}
  virtual ~ReplicateHandler() {}

  // The earlyStop call can be used to quit the replication mid-way, leaving
  // a ciphertext with (e.g.) two different entries, each replicated n/2 times
  virtual bool earlyStop(long d, long k, long prodDim) { return false; }
};
// Applications will derive from this class a handler that actually
// does something with the replicated cipehrtexts. But it can be used
// by itself as a do-nothing replicator for debugging, or to calculate
// the required automorphisms (see automorphVals in numbTh.h)


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
		  ReplicateHandler *handler, long recBound = 64,
		  RepAuxDim* repAuxPtr=NULL);

//! return the result as a vector of ciphertexts, mostly useful for
//! debugging purposes (for real parameters would take a lot of memory)
void replicateAll(std::vector<Ctxt>& v, const EncryptedArray& ea,
       	          const Ctxt& ctxt, long recBound = 64,
		  RepAuxDim* repAuxPtr=NULL);

//! This function is obsolete, and is kept for historical purposes only. It
//! was a first attempt at implementing the O(1)-amortized algorithm, but is
//! less efficient than the function above.
void replicateAllOrig(const EncryptedArray& ea, const Ctxt& ctxt,
                      ReplicateHandler *handler, RepAux* repAuxPtr=NULL);

void replicate(const EncryptedArray& ea, NewPlaintextArray& pa, long i);

// Structures to keep tables of maskign constants that are used in
// replication. A calling application can either supply this structure
// itself (if it is going to replicate the same tables in multiple
// replicate operations), or is cal let the repliction code use its
// own tables (in which case they are destroyed at the end of the
// replication process)

//! @cond FALSE (make doxygen ignore this class)
class RepAux { // one table for the whole thing
private:
  vector< copied_ptr<DoubleCRT> > _tab;

public:
  copied_ptr<DoubleCRT>& tab(long i) {
    if (i >= lsize(_tab)) _tab.resize(i+1);
    return _tab[i];
  }
};

class RepAuxDim { // two tables per dimension
private:
  vector< vector< copied_ptr<DoubleCRT> > > _tab, _tab1;

public:
  copied_ptr<DoubleCRT>& tab(long d, long i) {
    if (d >= lsize(_tab)) _tab.resize(d+1);
    if (i >= lsize(_tab[d])) _tab[d].resize(i+1);
    return _tab[d][i];
  }

  copied_ptr<DoubleCRT>& tab1(long d, long i) {
    if (d >= lsize(_tab1)) _tab1.resize(d+1);
    if (i >= lsize(_tab1[d])) _tab1[d].resize(i+1);
    return _tab1[d][i];
  }
};
//! @endcond

#endif
