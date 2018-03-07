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

#include "replicate.h"
#include "timing.h"
#include "cloned_ptr.h"




NTL_THREAD_LOCAL 
bool replicateVerboseFlag = false;

// The value in slot #pos is replicated in all
// other slots.  If there are n slots, this algorithm performs
// O(log n) 1D rotations.  

void replicate(const EncryptedArray& ea, Ctxt& ctxt, long pos)
{
  long nSlots = ea.size();
  assert(pos >= 0 && pos < nSlots); 

  ZZX mask;
  ea.encodeUnitSelector(mask, pos);
  ctxt.multByConstant(mask);
  replicate0(ea, ctxt, pos);
}

// Assumes all slots are zero except slot #pos,
// which is duplicated in all other slots

void replicate0(const EncryptedArray& ea, Ctxt& ctxt, long pos)
{
  long dim = ea.dimension();

  for (long d = 0; d < dim; d++) {
    if (!ea.nativeDimension(d)) {
      long shamt = -ea.coordinate(d, pos);
      ea.rotate1D(ctxt, d, shamt, true); // "don't care"
    }

    Ctxt ctxt_orig = ctxt; 

    long sz = ea.sizeOfDimension(d);
    long k = NumBits(sz);
    long e = 1;

    // now process bits k-2 down to 0
    for (long j = k-2; j >= 0; j--) {
      // e -> 2*e
      Ctxt tmp = ctxt;
      ea.rotate1D(tmp, d, e, true); // "don't care"
      ctxt += tmp;
      e = 2*e;
      
      long b = bit(sz, j); // bit j of sz
      // e -> e+b
      if (b) {
        ea.rotate1D(ctxt, d, 1, true); // "don't care"
        ctxt += ctxt_orig;
        e++;
      }
    }
  }
}


// The following code implements a recursive, O(1)-amortized
// algorithm for replications


// returns greatest integer k such that 2^k <= n
static
long GreatestPowerOfTwo(long n)
{
  assert(n >0);

  long k;

  k = 0;
  while ((1L << k) <= n/2) k++;

  return k;
}

// selects range of slots [lo..hi)
static
void SelectRange(const EncryptedArray& ea, ZZX& mask, long lo, long hi)
{
  long nSlots = ea.size();

  assert(lo >= 0 && lo <= hi && hi <= nSlots);

  vector<long> maskArray;
  maskArray.resize(nSlots);
  for (long i = 0; i < nSlots; i++) maskArray[i] = 0;
  for (long i = lo; i < hi; i++) maskArray[i] = 1;
  
  ea.encode(mask, maskArray);
}

// selects range of slots [lo..hi)
static
void SelectRange(const EncryptedArray& ea, Ctxt& ctxt, long lo, long hi)
{
  ZZX mask;
  SelectRange(ea, mask, lo, hi);
  ctxt.multByConstant(mask);
}

// recursiveReplicate:
//   n = GreatestPowerOfTwo(ea.size())
//   0 <= k <= n: size of current interval
//   0 <= pos < ea.size(): relative position of first vector
//   0 <= limit < ea.size(): max # of positions to process

static
void recursiveReplicate(const EncryptedArray& ea, const Ctxt& ctxt, 
                        long n, long k, long pos, long limit,  
                        RepAux& repAux,
                        ReplicateHandler *handler)
{
  if (pos >= limit) return;

  if (replicateVerboseFlag) {
    // DEBUG code
    cerr << "check: " << k; CheckCtxt(ctxt, "");
  }

  long nSlots = ea.size();

  if (k == 0) {

    if ( (1L << n) >= nSlots) {
      handler->handle(ctxt);
      return;
    }

    // need to replicate to fill positions [ (1L << n) .. nSlots )
    if (repAux.tab(0).null()) {
      // need to generate mask
      ZZX mask;
      SelectRange(ea, mask, 0, nSlots - (1L << n));
      repAux.tab(0).set_ptr(new DoubleCRT(mask, ea.getContext()));
    }


    Ctxt ctxt_tmp = ctxt;
    ctxt_tmp.multByConstant(*repAux.tab(0));

    ea.rotate(ctxt_tmp, 1L << n);
    ctxt_tmp += ctxt;
    handler->handle(ctxt_tmp);
    return;
  }


  k--;

  Ctxt ctxt_masked = ctxt;

  { // artificial scope to miminize storage in
    // the recursion


    { // another artificial scope

      // mask should be at index k+1

      if (repAux.tab(k+1).null()) {
        // need to generate mask

        vector< long > maskArray;
        maskArray.resize(nSlots);
        for (long i = 0; i < (1L << n); i++)
          maskArray[i] = 1- bit(i, k); // the reverse of bit k of i
        for (long i = (1L << n); i < nSlots; i++)
          maskArray[i] = 0;

        ZZX mask;
        ea.encode(mask, maskArray);
        repAux.tab(k+1).set_ptr(new DoubleCRT(mask, ea.getContext()));
      }

      ctxt_masked.multByConstant(*repAux.tab(k+1));
    }

    Ctxt ctxt_left = ctxt_masked;
    ea.rotate(ctxt_left, 1L << k);
    ctxt_left += ctxt_masked;

    recursiveReplicate(ea, ctxt_left, n, k, pos, limit, repAux, handler);
  
  }
 
  pos += (1L << k);
  if (pos >= limit)
    return;

  Ctxt ctxt_right = ctxt;
  ctxt_right -= ctxt_masked; 
  ctxt_masked = ctxt_right; // reuse ctxt_masked as a temp
  ea.rotate(ctxt_masked, -(1L << k));
  ctxt_right += ctxt_masked;

  recursiveReplicate(ea, ctxt_right, n, k, pos, limit, repAux, handler);
}

void replicateAllOrig(const EncryptedArray& ea, const Ctxt& ctxt,
                      ReplicateHandler *handler, RepAux* repAuxPtr)
{
  long nSlots = ea.size();
  long n = GreatestPowerOfTwo(nSlots); // 2^n <= nSlots

  Ctxt ctxt1 = ctxt;

  if ((1L << n) < nSlots)
    SelectRange(ea, ctxt1, 0, 1L << n);

  RepAux repAux;
  if (repAuxPtr==NULL) repAuxPtr = &repAux;

  recursiveReplicate(ea, ctxt1, n, n, 0, 1L << n, 
                     *repAuxPtr, handler);

  if ((1L << n) < nSlots) {
    ctxt1 = ctxt;
    SelectRange(ea, ctxt1, 1L << n, nSlots);
    ea.rotate(ctxt1, -(1L << n));
    recursiveReplicate(ea, ctxt1, n, n, 1L << n, nSlots, *repAuxPtr, handler);
  }
    
}

// The following code is based on the same logic as the
// recursive, O(1)-amortized algorithm, but works one
// dimension at a time, which allows us to use "native"
// rotations

// selects range of slots [lo..hi) in dimension d
static
void SelectRangeDim(const EncryptedArray& ea, ZZX& mask, long lo, long hi,
                    long d)
{
  long nSlots = ea.size();

  assert(d >= 0 && d < ea.dimension());
  assert(lo >= 0 && lo <= hi && hi <= ea.sizeOfDimension(d));

  vector<long> maskArray;
  maskArray.resize(nSlots);
  for (long i = 0; i < nSlots; i++) {
    long c = ea.coordinate(d, i);
    if (c >= lo && c < hi) 
      maskArray[i] = 1;
    else
      maskArray[i] = 0;
  }
  
  ea.encode(mask, maskArray);
}

// selects range of slots [lo..hi)
static
void SelectRangeDim(const EncryptedArray& ea, Ctxt& ctxt, long lo, long hi,
                    long d)
{
  ZZX mask;
  SelectRangeDim(ea, mask, lo, hi, d);
  ctxt.multByConstant(mask);
}

// replicateOneBlock: assumes that all slots are zero, except for one
// "block" whose coordinates in dimension d lie in the interval
//            [ pos*blockSize .. pos*(blockSize+1) -1 ]
// This block is then replicated throught the range 
//            [ 0.. floor(dSize/blockSize)*blockSize -1 ]

static void replicateOneBlock(const EncryptedArray& ea, Ctxt& ctxt, 
			      long pos, long blockSize, long d)
{
  long dSize = ea.sizeOfDimension(d);

  // Move this block to position 0. We can skip this step in
  // "good dimensions" of size divisible by the block size
  if (pos != 0 &&  (!ea.nativeDimension(d) || dSize % blockSize != 0) ) {
    ea.rotate1D(ctxt, d, -pos*blockSize, true);
  }

  long sz = dSize/blockSize; // how many blocks fit in this dimension

  if (sz == 1) return; // nothing to do, only one block in this dimension

  // do the actual replication using "shift and add"

  long k = NumBits(sz);
  long e = 1;
  Ctxt ctxt_orig = ctxt;

  // now process bits k-2 down to 0
  for (long j = k-2; j >= 0; j--) {
    // e -> 2*e
    Ctxt tmp = ctxt;
    ea.rotate1D(tmp, d, e*blockSize, /*don't-care-flag=*/true);
    ctxt += tmp;
    e = 2*e;
    
    long b = bit(sz, j); // bit j of sz
    // e -> e+b
    if (b) {
      ea.rotate1D(ctxt, d, 1*blockSize, /*don't-care-flag=*/true);
      ctxt += ctxt_orig;
      e++;
    }
  }
}

// forward declaration...mutual recursion
static
void replicateAllNextDim(const EncryptedArray& ea, const Ctxt& ctxt,
                         long d, long dimProd, long recBound,
                         RepAuxDim& repAux, ReplicateHandler *handler);



// recursiveReplicateDim:
//   d = dimension
//   ea.sizeOfDimension(d)/2 <= extent <= ea.sizeOfDimension(d),
//     only positions [0..extent) are non-zero
//   1 <= 2^k <= extent: size of current interval
//   0 <= pos < ea.sizeOfDimension(d): relative position of first vector
//   0 <= limit < ea.sizeOfDimension(): max # of positions to process
//   dimProd: product of dimensions 0..d
//   recBound: recursion bound (controls noise) 
//
// SHAI: limit and extent are always the same, it seems
static
void recursiveReplicateDim(const EncryptedArray& ea, const Ctxt& ctxt, 
                           long d, long extent, long k, long pos, long limit,  
                           long dimProd, long recBound,
                           RepAuxDim& repAux,
                           ReplicateHandler *handler)
{
  if (pos >= limit) return;

  if (replicateVerboseFlag) { // DEBUG code
    cerr << "check: " << k; CheckCtxt(ctxt, "");
  }
  
  long dSize = ea.sizeOfDimension(d);
  long nSlots = ea.size();

  if (k == 0) { // last level in this dimension: blocks of size 2^k=1

    if ( extent >= dSize) { // nothing to do in this dimension
      replicateAllNextDim(ea, ctxt, d+1, dimProd, recBound, repAux, handler);
      return;
    } // SHAI: Will we ever have extent > dSize??

    // need to replicate to fill positions [ (1L << n) .. dSize-1 ]

    if (repAux.tab(d,0).null()) { // generate mask if not there already
      ZZX mask;
      SelectRangeDim(ea, mask, 0, dSize - extent, d);
      repAux.tab(d, 0).set_ptr(new DoubleCRT(mask, ea.getContext()));
    }

    Ctxt ctxt_tmp = ctxt;
    ctxt_tmp.multByConstant(*repAux.tab(d, 0));

    ea.rotate1D(ctxt_tmp, d, extent, /*don't-care-flag=*/true);
    ctxt_tmp += ctxt;
    replicateAllNextDim(ea, ctxt_tmp, d+1, dimProd, recBound, repAux, handler);
    return;
  }

  // If we need to stop early, call the handler
  if (handler->earlyStop(d, k, dimProd)) {
    handler->handle(ctxt);
    return;
  }

  k--;
  Ctxt ctxt_masked = ctxt;

  {   // artificial scope to miminize storage in the recursion
    { // another artificial scope (SHAI: this seems redundant)

      // generate mask at index k+1, if not there yet

      if (repAux.tab(d, k+1).null()) { // need to generate
        vector< long > maskArray(nSlots,0);
        for (long i = 0; i < nSlots; i++) {
          long c = ea.coordinate(d, i);
          if (c < extent && bit(c, k) == 0)
            maskArray[i] = 1;
        }
	// store this mask in the repAux table
        ZZX mask;
        ea.encode(mask, maskArray);
        repAux.tab(d, k+1).set_ptr(new DoubleCRT(mask, ea.getContext()));
      }

      // Apply mask to zero out slots in ctxt
      ctxt_masked.multByConstant(*repAux.tab(d, k+1));
    }

    Ctxt ctxt_left = ctxt_masked;
    ea.rotate1D(ctxt_left, d, 1L << k, /*don't-care-flag=*/true);
    ctxt_left += ctxt_masked;

    recursiveReplicateDim(ea, ctxt_left, d, extent, k, pos, limit, 
                          dimProd, recBound, repAux, handler);
  }

  pos += (1L << k);
  if (pos >= limit)
    return;

  Ctxt ctxt_right = ctxt;
  ctxt_right -= ctxt_masked; 
  ctxt_masked = ctxt_right; // reuse ctxt_masked as a temp
  ea.rotate1D(ctxt_masked, d, -(1L << k), /*don't-care-flag=*/true);
  ctxt_right += ctxt_masked;

  recursiveReplicateDim(ea, ctxt_right, d, extent, k, pos, limit, 
                        dimProd, recBound, repAux, handler);
}

void replicateAllNextDim(const EncryptedArray& ea, const Ctxt& ctxt,
                         long d, long dimProd, long recBound,
                         RepAuxDim& repAux, ReplicateHandler *handler)

{
  assert(d >= 0);

  // If already fully replicated (or we need to stop early), call the handler
  if (d >= ea.dimension() || handler->earlyStop(d,/*k=*/-1,dimProd)) {
    handler->handle(ctxt);
    return;
  }
  
  long dSize = ea.sizeOfDimension(d);
  dimProd *= dSize; // product of all dimensions including this one

  long n = GreatestPowerOfTwo(dSize); // 2^n <= dSize
  long k = n;

  // We replicate 2^k-size blocks along this dimension, then call the
  // recursive procedure to handle the smaller subblocks. Consider for
  // example a 2D 5x2 cube, so the original slots are
  //
  //    ( s0 s2 s4 s6 s8 )
  //    ( s1 s3 s5 s7 s9 )
  //
  // Say that we start with k=2 in the 1st dimension (of size 5), we
  // will prepare floor(5/2)=2 ciphertexts as follows:
  //
  //    ( s0 s2 s0 s2 0 )   ( s4 s6 s4 s6 0 )
  //    ( s1 s3 s1 s3 0 )   ( s5 s7 s5 s7 0 )
  //
  // The call to recursiveReplicateDim (still with k=2) will first copy
  // s0/s1 and s4/s5 to the zero column at the end, then make a recursive
  // call with k=1 that will complete the replication along the current
  // dimension, resulting in the 4 ciphertexts
  // 
  //  (s0 s0 s0 s0 s0) (s2 s2 s2 s2 s2) (s4 s4 s4 s4 s4) (s6 s6 s6 s6 s6)
  //  (s1 s1 s1 s1 s1) (s3 s3 s3 s3 s3) (s5 s5 s5 s5 s5) (s7 s7 s7 s7 s7)
  //
  // Then a recursive call for the next dimension will complete the
  // replication of these entries, and a final step will deal with the
  // "leftover" positions s8 s9
  


  // The logic below cut the recursion depth by starting from smaller
  // blocks (by default size approx n rather than 2^n).
  // The inital block size is controlled by the recBound parameter:
  //   + recBound>0: blocks of size min(~n, 2^recBound). this ensures
  //     recursion depth <= recBound, and typically much smaller (~log n)
  //   + recBound=0: blocks of size 1 (no recursion)
  //   + recBound<0: blocks of size 2^n (full recursion)

  if (recBound >= 0) { // use heuristic recursion bound
    k = 0;
    if (dSize > 2 && dimProd*NumBits(dSize) > ea.size() / 8) {
      k = NumBits(NumBits(dSize))-1;
      if (k > n) k = n;
      if (k > recBound) k = recBound;
    }
  }
  else { // SHAI: I don't understand this else case
    k = -recBound;
    if (k > n) k = n;
  }

  long blockSize = 1L << k;        // blocks of size 2^k
  long numBlocks = dSize/blockSize;
  long extent = numBlocks * blockSize;

  // extent is an integral multiple of the block size, the recursive
  // call replicates only these slots, and then we have a separate
  // call for the leftover slots.

  Ctxt ctxt1 = ctxt;

  if (extent < dSize) { // select only the slots 0..extent-1 in this dimension
    if (repAux.tab1(d, 0).null()) { // generate mask if not already there
      ZZX mask;
      SelectRangeDim(ea, mask, 0, extent, d);
      repAux.tab1(d, 0).set_ptr(new DoubleCRT(mask, ea.getContext()));
      // store mask in 2nd table (tab1)
    }
    ctxt1.multByConstant(*repAux.tab1(d, 0)); // mult by mask to zero out slots
  }

  if (numBlocks == 1) { // just one block, call the recursive replication
    recursiveReplicateDim(ea, ctxt1, d, extent, k, 0, extent, 
                          dimProd, recBound, repAux, handler);
  }
  else { // replicate the slots in each block separately
    for (long pos = 0; pos < numBlocks; pos++) {
      Ctxt ctxt2 = ctxt1;
      // zero-out all the slots outside the current block
      SelectRangeDim(ea, ctxt2, pos*blockSize, (pos+1)*blockSize, d);

      // replicate the current block across this dimenssion using a
      // simple shift-and-add procedure.
      replicateOneBlock(ea, ctxt2, pos, blockSize, d);

      // now call the recursive replication to do the rest of the work
      recursiveReplicateDim(ea, ctxt2, d, extent, k, 0, extent, 
                            dimProd, recBound, repAux, handler);
    }
  }

  // If dSize is not an integral number of blocks, then we still need
  // to deal with the leftover slots.
  if (extent < dSize) {
    // zero-out the slots from before, leaving only the leftover slots
    ctxt1 = ctxt;
    if (repAux.tab1(d, 1).null()) { // generate mask if not already there
      ZZX mask;
      SelectRangeDim(ea, mask, extent, dSize, d);
      repAux.tab1(d, 1).set_ptr(new DoubleCRT(mask, ea.getContext()));
    }
    ctxt1.multByConstant(*repAux.tab1(d,1)); // mult by mask to zero out slots

    // move relevant slots to the beginning
    ea.rotate1D(ctxt1, d, -extent, /*don't-care-flag=*/true);

    // replicate the leftover block across this dimenssion using a
    // simple shift-and-add procedure.
    replicateOneBlock(ea, ctxt1, 0, blockSize, d);

    // now call the recursive replication to do the rest of the work
    recursiveReplicateDim(ea, ctxt1, d, extent, k, extent, dSize, 
                          dimProd, recBound, repAux, handler);
  }
}

// recBound < 0 => pure recursion
// recBound == 0 => no recursion
// otherwise, a recursion depth is chosen heuristically,
//   but is capped at recBound
void
replicateAll(const EncryptedArray& ea, const Ctxt& ctxt, 
	     ReplicateHandler *handler, long recBound, RepAuxDim* repAuxPtr)
{
  FHE_TIMER_START;
  RepAuxDim repAux;
  if (repAuxPtr==NULL) repAuxPtr = &repAux;
  replicateAllNextDim(ea, ctxt, 0, 1, recBound, *repAuxPtr, handler);
}


//! @brief An implementation of ReplicateHandler that explicitly returns
//!   all the replicated cipehrtexts in one big vector.
//!
//! This is useful mostly for debugging purposes, for real parameters
//! it would take a lot of memory.
class ExplicitReplicator : public ReplicateHandler {
  std::vector<Ctxt>& v; // space to store all cipehrtexts
  long slot;
public:
  // _v must already be of the right size (=number-of-slots)
  ExplicitReplicator(std::vector<Ctxt>& _v): v(_v),slot(0) {}
  virtual void handle(const Ctxt& ctxt) { v[slot++] = ctxt; }
};

// Returns the result as a vector of ciphertexts
void replicateAll(std::vector<Ctxt>& v, const EncryptedArray& ea,
       	          const Ctxt& ctxt, long recBound, RepAuxDim* repAuxPtr)
{
  v.resize(ea.size(), ctxt);
  ExplicitReplicator handler(v);
  replicateAll(ea, ctxt, &handler, recBound, repAuxPtr);
}


//=======================================================================================

// procedures for replicating plaintext-arrays, useful for debugging

template<class type>
class replicate_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, long i)
  {
    PA_BOILER

    assert(i >= 0 && i < n);
    for (long j = 0; j < n; j++) {
      if (j != i) data[j] = data[i];
    }
  }
};


void replicate(const EncryptedArray& ea, NewPlaintextArray& pa, long i)
{
  ea.dispatch<replicate_pa_impl>(pa, i); 
}




