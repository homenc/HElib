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

//! @cond FALSE (make doxygen ignore this class)
class RepAux {
private:
  vector< copied_ptr<DoubleCRT> > _tab;

public:
  copied_ptr<DoubleCRT>& tab(long i) {
    if (i >= lsize(_tab)) _tab.resize(i+1);
    return _tab[i];
  }
};
//! @endcond

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
                      ReplicateHandler *handler)

{
  long nSlots = ea.size();
  long n = GreatestPowerOfTwo(nSlots);

  Ctxt ctxt1 = ctxt;

  if ((1L << n) < nSlots)
    SelectRange(ea, ctxt1, 0, 1L << n);

  RepAux repAux;

  recursiveReplicate(ea, ctxt1, n, n, 0, 1L << n, 
                     repAux, handler);

  if ((1L << n) < nSlots) {
    ctxt1 = ctxt;
    SelectRange(ea, ctxt1, 1L << n, nSlots);
    ea.rotate(ctxt1, -(1L << n));
    recursiveReplicate(ea, ctxt1, n, n, 1L << n, nSlots, repAux, handler);
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

//! @cond FALSE (make doxygen ignore this class)
class RepAuxDim {
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

// replicateOneBlock: assumes that all slots are zero,
// except for those in the "block" whose coordinate
// in dimension d lies in the interval [ pos*blockSize .. pos*(blockSize+1) )
// This block is then replicated throught the range 
// [ 0.. floor(dSize/blockSize)*blockSize )

static
void replicateOneBlock(const EncryptedArray& ea, Ctxt& ctxt, 
                       long pos, long blockSize, long d)
{
  long dSize = ea.sizeOfDimension(d);

  if (pos != 0 &&  (!ea.nativeDimension(d) || dSize % blockSize != 0) ) {
    ea.rotate1D(ctxt, d, -pos*blockSize, true);
  }

  long sz = dSize/blockSize;

  if (sz == 1) return;

  long k = NumBits(sz);
  long e = 1;

  Ctxt ctxt_orig = ctxt;

  // now process bits k-2 down to 0
  for (long j = k-2; j >= 0; j--) {
    // e -> 2*e
    Ctxt tmp = ctxt;
    ea.rotate1D(tmp, d, e*blockSize, true); // "don't care"
    ctxt += tmp;
    e = 2*e;
    
    long b = bit(sz, j); // bit j of sz
    // e -> e+b
    if (b) {
      ea.rotate1D(ctxt, d, 1*blockSize, true); // "don't care"
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

static
void recursiveReplicateDim(const EncryptedArray& ea, const Ctxt& ctxt, 
                           long d, long extent, long k, long pos, long limit,  
                           long dimProd, long recBound,
                           RepAuxDim& repAux,
                           ReplicateHandler *handler)
{
  if (pos >= limit) return;

  if (replicateVerboseFlag) {
    // DEBUG code
    cerr << "check: " << k; CheckCtxt(ctxt, "");
  }

  long dSize = ea.sizeOfDimension(d);
  long nSlots = ea.size();

  if (k == 0) {

    if ( extent >= dSize) {
      replicateAllNextDim(ea, ctxt, d+1, dimProd, recBound, repAux, handler);
      return;
    }

    // need to replicate to fill positions [ (1L << n) .. dSize )
    if (repAux.tab(d, 0).null()) {
      // need to generate mask
      ZZX mask;
      SelectRangeDim(ea, mask, 0, dSize - extent, d);
      repAux.tab(d, 0).set_ptr(new DoubleCRT(mask, ea.getContext()));
    }

    Ctxt ctxt_tmp = ctxt;
    ctxt_tmp.multByConstant(*repAux.tab(d, 0));

    ea.rotate1D(ctxt_tmp, d, extent, true);
    ctxt_tmp += ctxt;
    replicateAllNextDim(ea, ctxt_tmp, d+1, dimProd, recBound, repAux, handler);
    return;
  }

  k--;

  Ctxt ctxt_masked = ctxt;

  { // artificial scope to miminize storage in
    // the recursion


    { // another artificial scope

      // mask should be at index k+1

      if (repAux.tab(d, k+1).null()) {
        // need to generate mask

        vector< long > maskArray;
        maskArray.resize(nSlots);

        for (long i = 0; i < nSlots; i++) {
          long c = ea.coordinate(d, i);
          if (c < extent && bit(c, k) == 0)
            maskArray[i] = 1;
          else
            maskArray[i] = 0;
        }

        ZZX mask;
        ea.encode(mask, maskArray);
        repAux.tab(d, k+1).set_ptr(new DoubleCRT(mask, ea.getContext()));
      }

      ctxt_masked.multByConstant(*repAux.tab(d, k+1));
    }

    Ctxt ctxt_left = ctxt_masked;
    ea.rotate1D(ctxt_left, d, 1L << k, true);
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
  ea.rotate1D(ctxt_masked, d, -(1L << k), true);
  ctxt_right += ctxt_masked;

  recursiveReplicateDim(ea, ctxt_right, d, extent, k, pos, limit, 
                        dimProd, recBound, repAux, handler);
}

void replicateAllNextDim(const EncryptedArray& ea, const Ctxt& ctxt,
                         long d, long dimProd, long recBound,
                         RepAuxDim& repAux, ReplicateHandler *handler)

{
  long dim = ea.dimension();

  assert(d >= 0);

  if (d >= dim) {
    handler->handle(ctxt);
    return;
  }
  
  long dSize = ea.sizeOfDimension(d);
  long n = GreatestPowerOfTwo(dSize);

  dimProd = dimProd * dSize;

  long k = n;

  // Note: recBound < 0 => |recBound| levels of recursion recursion,
  // recBound == 0 => no recursion

  if (recBound >= 0) {
    // heuristic recursion bound
    k = 0;
    if (dSize > 2 && dimProd*NumBits(dSize) > ea.size() / 8) {
      k = NumBits(NumBits(dSize))-1;
      if (k > n) k = n;
      if (k > recBound) k = recBound;
    }
  }
  else {
    k = -recBound;
    if (k > n) k = n;
  }

  long blockSize = 1L << k;
  long numBlocks = dSize/blockSize;
  long extent = numBlocks * blockSize;

  Ctxt ctxt1 = ctxt;

  if (extent < dSize) {
    if (repAux.tab1(d, 0).null()) {
      ZZX mask;
      SelectRangeDim(ea, mask, 0, extent, d);
      repAux.tab1(d, 0).set_ptr(new DoubleCRT(mask, ea.getContext()));
    }
    ctxt1.multByConstant(*repAux.tab1(d, 0));
  }

  if (numBlocks == 1) {
    recursiveReplicateDim(ea, ctxt1, d, extent, k, 0, extent, 
                          dimProd, recBound, repAux, handler);
  }
  else {
    for (long pos = 0; pos < numBlocks; pos++) {
      Ctxt ctxt2 = ctxt1;
      SelectRangeDim(ea, ctxt2, pos*blockSize, (pos+1)*blockSize, d);
      replicateOneBlock(ea, ctxt2, pos, blockSize, d);
      recursiveReplicateDim(ea, ctxt2, d, extent, k, 0, extent, 
                            dimProd, recBound, repAux, handler);
    }
  }

  if (extent < dSize) {
    ctxt1 = ctxt;

    if (repAux.tab1(d, 1).null()) {
      ZZX mask;
      SelectRangeDim(ea, mask, extent, dSize, d);
      repAux.tab1(d, 1).set_ptr(new DoubleCRT(mask, ea.getContext()));
    }
    ctxt1.multByConstant(*repAux.tab1(d, 1));

    ea.rotate1D(ctxt1, d, -extent, true);

    replicateOneBlock(ea, ctxt1, 0, blockSize, d);
    recursiveReplicateDim(ea, ctxt1, d, extent, k, extent, dSize, 
                            dimProd, recBound, repAux, handler);
  }
    
}

// recBound < 0 => pure recursion
// recBound == 0 => no recursion
// otherwise, a recursion depth is chosen heuristically,
//   but is capped at recBound

void replicateAll(const EncryptedArray& ea, const Ctxt& ctxt, 
                         ReplicateHandler *handler, long recBound)
{
  RepAuxDim repAux;
  replicateAllNextDim(ea, ctxt, 0, 1, recBound, repAux, handler);
}






//=======================================================================================



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
  ea.dispatch<replicate_pa_impl>(Fwd(pa), i); 
}




