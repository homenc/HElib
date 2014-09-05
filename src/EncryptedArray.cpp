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
/* EncryptedArray.cpp - Data-movement operations on arrays of slots
 */
#include <algorithm>
#include <NTL/ZZ.h>
NTL_CLIENT
#include "EncryptedArray.h"
#include "timing.h"
#include "cloned_ptr.h"
#include "NumbTh.h"


EncryptedArrayBase* buildEncryptedArray(const FHEcontext& context, const ZZX& G,
					const PAlgebraMod& alMod)
{
  switch (alMod.getTag()) {
    case PA_GF2_tag: {
      return new EncryptedArrayDerived<PA_GF2>(context, conv<GF2X>(G), alMod);
    }

    case PA_zz_p_tag: {
      zz_pBak bak; bak.save(); alMod.restoreContext();
      return new EncryptedArrayDerived<PA_zz_p>(context, conv<zz_pX>(G), alMod);
    }

    default: return NULL;
  }
}


template<class type>
EncryptedArrayDerived<type>::EncryptedArrayDerived(
   const FHEcontext& _context, const RX& _G, const PAlgebraMod& alMod)
  : context(_context), tab(alMod.getDerived(type()))
{
  tab.genMaskTable();
  tab.mapToSlots(mappingData, _G); // Compute the base-G representation maps
}

// rotate ciphertext in dimension i by amt
template<class type>
void EncryptedArrayDerived<type>::rotate1D(Ctxt& ctxt, long i, long amt, bool dc) const
{
  FHE_TIMER_START;
  const PAlgebra& al = context.zMStar;

  const vector< vector< RX > >& maskTable = tab.getMaskTable();

  RBak bak; bak.save(); tab.restoreContext();

  assert(&context == &ctxt.getContext());
  assert(i >= 0 && i < (long)al.numOfGens());

  // Make sure amt is in the range [1,ord-1]
  long ord = al.OrderOf(i);
  amt %= ord;
  if (amt == 0) return;
  long signed_amt = amt;
  if (amt < 0) amt += ord;

  // DIRT: the above assumes division with remainder
  // follows C++11 and C99 rules

  if (al.SameOrd(i)) { // a "native" rotation
    long val = PowerMod(al.ZmStarGen(i), amt, al.getM());
    ctxt.smartAutomorph(val);
  }
  else if (dc) { 
    // the "don't care" case...it is presumed that any shifts
    // "off the end" are zero.  For this, we have to use 
    // the "signed" version of amt.
    long val = PowerMod(al.ZmStarGen(i), signed_amt, al.getM());
    ctxt.smartAutomorph(val);
  }
  else {
    // more expensive "non-native" rotation
    assert(maskTable[i].size() > 0);
    long val = PowerMod(al.ZmStarGen(i), amt, al.getM());
    long ival = PowerMod(al.ZmStarGen(i), amt-ord, al.getM());

    const RX& mask = maskTable[i][ord-amt];
    DoubleCRT m1(conv<ZZX>(mask), context, ctxt.getPrimeSet());
    Ctxt tmp(ctxt); // a copy of the ciphertext

    tmp.multByConstant(m1);    // only the slots in which m1=1
    ctxt -= tmp;               // only the slots in which m1=0
    ctxt.smartAutomorph(val);  // shift left by val
    tmp.smartAutomorph(ival);  // shift right by ord-val
    ctxt += tmp;               // combine the two parts
  }
  FHE_TIMER_STOP;
}

// Shift k positions along the i'th dimension with zero fill.
// Negative shift amount denotes shift in the opposite direction.
template<class type>
void EncryptedArrayDerived<type>::shift1D(Ctxt& ctxt, long i, long k) const
{
  FHE_TIMER_START;
  const PAlgebra& al = context.zMStar;

  const vector< vector< RX > >& maskTable = tab.getMaskTable();

  RBak bak; bak.save(); tab.restoreContext();

  assert(&context == &ctxt.getContext());
  assert(i >= 0 && i < (long)al.numOfGens());

  long ord = al.OrderOf(i);

  if (k <= -ord || k >= ord) {
    ctxt.multByConstant(to_ZZX(0));
    return;
  }

  // Make sure amt is in the range [1,ord-1]
  long amt = k % ord;
  if (amt == 0) return;
  if (amt < 0) amt += ord;

  RX mask = maskTable[i][ord-amt];

  long val;
  if (k < 0)
    val = PowerMod(al.ZmStarGen(i), amt-ord, al.getM());
  else {
    mask = 1 - mask;
    val = PowerMod(al.ZmStarGen(i), amt, al.getM());
  }
  DoubleCRT m1(conv<ZZX>(mask), context, ctxt.getPrimeSet());
  ctxt.multByConstant(m1);   // zero out slots where mask=0
  ctxt.smartAutomorph(val);  // shift left by val
  FHE_TIMER_STOP;
}


// NOTE: masking depth: if there are N dimensions, and if for i = 1..N
// we define c_i = 1 if dimension i is bad and 0 o/w, then the masking
// depth is N - 1 + \sum_{i=1} c_i.  

template<class type>
void EncryptedArrayDerived<type>::rotate(Ctxt& ctxt, long amt) const
{
  FHE_TIMER_START;

  const PAlgebra& al = context.zMStar;

  const vector< vector< RX > >& maskTable = tab.getMaskTable();

  RBak bak; bak.save(); tab.restoreContext();

  assert(&context == &ctxt.getContext());

  // Simple case: just one generator
  if (al.numOfGens()==1) { // VJS: bug fix: <= must be ==
    rotate1D(ctxt, 0, amt);
    return;
  }

  // Make sure that amt is in [1,nslots-1]
  amt %= (long) al.getNSlots();
  if (amt == 0) { return; }
  if (amt < 0) amt += al.getNSlots();

  // rotate the ciphertext, one dimension at a time
  long i = al.numOfGens()-1;
  long v = al.coordinate(i, amt);
  RX mask = maskTable[i][v];
  Ctxt tmp(ctxt.getPubKey());
  const RXModulus& PhimXmod = tab.getPhimXMod();

  // optimize for the common case where the last generator has order in
  // Zm*/(p) different than its order in Zm*. In this case we can combine
  // the rotate1D relative to this generator with the masking after the
  // rotation. This saves one mult-by-constant, since we use the same mask
  // inside rotate1D as in the loop below.

  if (al.SameOrd(i) || v==0) rotate1D(ctxt, i, v); // no need to optimize
  else {
    long ord = al.OrderOf(i);
    long val = PowerMod(al.ZmStarGen(i), v, al.getM());
    long ival = PowerMod(al.ZmStarGen(i), v-ord, al.getM());

    DoubleCRT m1(conv<ZZX>(maskTable[i][ord-v]), context, ctxt.getPrimeSet());
    tmp = ctxt;  // a copy of the ciphertext

    tmp.multByConstant(m1);    // only the slots in which m1=1
    ctxt -= tmp;               // only the slots in which m1=0
    ctxt.smartAutomorph(val);  // shift left by val
    tmp.smartAutomorph(ival);  // shift right by ord-val

    // apply rotation relative to next generator before combining the parts
    --i;
    v = al.coordinate(i, amt);
    rotate1D(ctxt, i, v); 
    rotate1D(tmp, i, v+1);
    ctxt += tmp;         // combine the two parts

    if (i <= 0) { return; }  // no more generators

    mask = ((mask * (maskTable[i][v] - maskTable[i][v+1])) % PhimXmod)
             + maskTable[i][v+1];  // update the mask for next iteration
  }

  // Handle rotation relative to all the other generators (if any)
  for (i--; i >= 0; i--) {
    v = al.coordinate(i, amt);

    DoubleCRT m1(conv<ZZX>(mask), context, ctxt.getPrimeSet());
    tmp = ctxt;
    tmp.multByConstant(m1); // only the slots in which mask=1
    ctxt -= tmp;            // only the slots in which mask=0

    rotate1D(tmp, i, v); 
    rotate1D(ctxt, i, v+1);
    ctxt += tmp;
    if (i>0) {
      mask = ((mask * (maskTable[i][v] - maskTable[i][v+1])) % PhimXmod)
             + maskTable[i][v+1];  // update the mask for next iteration
    }
  }
  FHE_TIMER_STOP;
}

template<class type>
void EncryptedArrayDerived<type>::shift(Ctxt& ctxt, long k) const
{
  FHE_TIMER_START;


  const PAlgebra& al = context.zMStar;

  const vector< vector< RX > >& maskTable = tab.getMaskTable();

  RBak bak; bak.save(); tab.restoreContext();

  assert(&context == &ctxt.getContext());

  // Simple case: just one generator
  if (al.numOfGens()==1) {
    shift1D(ctxt, 0, k);
    return;
  }

  long nSlots = al.getNSlots();

  // Shifting by more than the number of slots gives an all-zero cipehrtext
  if (k <= -nSlots || k >= nSlots) {
    ctxt.multByConstant(to_ZZX(0));
    return;
  }

  // Make sure that amt is in [1,nslots-1]
  long amt = k % nSlots;
  if (amt == 0) return;
  if (amt < 0) amt += nSlots;

  // rotate the ciphertext, one dimension at a time
  long i = al.numOfGens()-1;
  long v = al.coordinate(i, amt);
  RX mask = maskTable[i][v];
  Ctxt tmp(ctxt.getPubKey());
  const RXModulus& PhimXmod = tab.getPhimXMod();

  rotate1D(ctxt, i, v);
  for (i--; i >= 0; i--) {
    v = al.coordinate(i, amt);

    DoubleCRT m1(conv<ZZX>(mask), context, ctxt.getPrimeSet());
    tmp = ctxt;
    tmp.multByConstant(m1); // only the slots in which mask=1
    ctxt -= tmp;            // only the slots in which mask=0
    if (i>0) {
      rotate1D(ctxt, i, v+1);
      rotate1D(tmp, i, v); 
      ctxt += tmp;                    // combine the two parts

      mask = ((mask * (maskTable[i][v] - maskTable[i][v+1])) % PhimXmod)
             + maskTable[i][v+1];  // update the mask before next iteration
    }
    else { // i == 0
      if (k < 0) v -= al.OrderOf(0);
      shift1D(tmp, 0, v);
      shift1D(ctxt, 0, v+1);
      ctxt += tmp;
    } 
  }
  FHE_TIMER_STOP;
}


template<class type>
void EncryptedArrayDerived<type>::encodeUnitSelector(ZZX& ptxt, long i) const
{
  assert(i >= 0 && i < (long)context.zMStar.getNSlots());
  RBak bak; bak.save(); tab.restoreContext();
  RX res;
  div(res, tab.getPhimXMod(), tab.getFactors()[i]); 
  mul(res, res, tab.getCrtCoeffs()[i]);
  conv(ptxt, res);
}

template<class type>
void EncryptedArrayDerived<type>::encode(ZZX& ptxt, const vector< RX >& array) const
{
  RX pp;
  tab.embedInSlots(pp, array, mappingData); 
  ptxt = conv<ZZX>(pp); 
}

template<class type>
void EncryptedArrayDerived<type>::decode(vector< RX >& array, const ZZX& ptxt) const
{
  FHE_TIMER_START;
  RX pp;
  conv(pp, ptxt);
  tab.decodePlaintext(array, pp, mappingData); 
  FHE_TIMER_STOP;
}

template<class type>
void EncryptedArrayDerived<type>::encode(ZZX& ptxt, const PlaintextArray& array) const
{
  assert(this == &(array.getEA().getDerived(type())));

  const PlaintextArrayDerived<type>& arr = array.getDerived(type());

  RBak bak; bak.save(); tab.restoreContext();
  encode(ptxt, arr.getData());
}

template<class type>
void EncryptedArrayDerived<type>::decode(PlaintextArray& array, const ZZX& ptxt) const
{
  assert(this == &(array.getEA().getDerived(type())));

  PlaintextArrayDerived<type>& arr = array.getDerived(type());

  RBak bak; bak.save(); tab.restoreContext();
  vector< RX > array1;
  decode(array1, ptxt);
  arr.setData(array1);
}

PlaintextArrayBase* buildPlaintextArray(const EncryptedArray& ea)
{
  switch (ea.getAlMod().getTag()) {
    case PA_GF2_tag: 
      return new PlaintextArrayDerived<PA_GF2>(ea);

    case PA_zz_p_tag: 
      return new PlaintextArrayDerived<PA_zz_p>(ea);

    default: return 0;
  }
}

// Other functions...

void runningSums(const EncryptedArray& ea, Ctxt& ctxt)
{
  long n = ea.size();

  long shamt = 1;
  while (shamt < n) {
    Ctxt tmp = ctxt;
    ea.shift(tmp, shamt);
    ctxt += tmp; // ctxt = ctxt + (ctxt >> shamt)
    shamt = 2*shamt;
  }
}

void totalSums(const EncryptedArray& ea, Ctxt& ctxt)
{
  long n = ea.size();

  if (n == 1) return;

  Ctxt orig = ctxt;

  long k = NumBits(n);
  long e = 1;

  for (long i = k-2; i >= 0; i--) {
    Ctxt tmp1 = ctxt;
    ea.rotate(tmp1, e);
    ctxt += tmp1; // ctxt = ctxt + (ctxt >>> e)
    e = 2*e;

    if (bit(n, i)) {
      Ctxt tmp2 = orig;
      ea.rotate(tmp2, e);
      ctxt += tmp2; // ctxt = ctxt + (orig >>> e)
                    // NOTE: we could have also computed
                    // ctxt =  (ctxt >>> e) + orig, however,
                    // this would give us greater depth/noise
      e += 1;
    }
  }
}

// Explicit instantiation

template class EncryptedArrayDerived<PA_GF2>;
template class EncryptedArrayDerived<PA_zz_p>;

template class PlaintextArrayDerived<PA_GF2>;
template class PlaintextArrayDerived<PA_zz_p>;
