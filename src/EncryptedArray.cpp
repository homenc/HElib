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

#include "EncryptedArray.h"

#include <algorithm>
#include "timing.h"
#include "cloned_ptr.h"


#ifdef FHE_BOOT_THREADS
NTL_THREAD_LOCAL MultiTask *bootTask = 0;
#endif


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

// this routine generates a random normal element
// and initializes a matrix mapping from polynomial to 
// normal basis and its inverse.  

template<class type>
void EncryptedArrayDerived<type>::initNormalBasisMatrix() const
{
  do {
    typename Lazy< Pair< Mat<R>, Mat<R> > >::Builder 
      builder(normalBasisMatrices); 

    if (!builder()) break;

    RBak bak; bak.save(); restoreContext();
    REBak ebak; ebak.save(); restoreContextForG();

    long d = RE::degree();
    long p = tab.getZMStar().getP();
    long r = tab.getR();

    // compute change of basis matrix CB
    mat_R CB;
    CB.SetDims(d, d);
    RE normal_element;
    RE H;
    bool got_it = false;

    H = power(conv<RE>(RX(1, 1)), p);
    

    do {
      NTL::random(normal_element);
   
      RE pow;
      pow = normal_element; 
      VectorCopy(CB[0], rep(pow), d);
      for (long i = 1; i < d; i++) {
        pow = eval(rep(pow), H);
        VectorCopy(CB[i], rep(pow), d);
      }

      Mat<ZZ> CB1;
      conv(CB1, CB);

      {
         zz_pBak bak1; bak1.save(); zz_p::init(p);
         Mat<zz_p> CB2;
         conv(CB2, CB1);
         got_it = determinant(CB2) != 0;
      }
    } while (!got_it);

    Mat<R> CBi;
    ppInvert(CBi, CB, p, r);

    UniquePtr< Pair< Mat<R>, Mat<R> > > ptr;
    ptr.make(CB, CBi);
    builder.move(ptr);
  } while(0);
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



/*****************************************************************/
/****************** Linear transformation code *******************/

// plaintextAutomorph: an auxilliary routine...maybe palce in NumbTh?
// Compute b(X) = a(X^k) mod Phi_m(X). Result is calclated in the output b
// "in place", so a should not alias b.

template <class RX, class RXModulus>
void plaintextAutomorph(RX& b, const RX& a, long k, const PAlgebra& zMStar, 
                       const RXModulus& PhimX)
{
  long m  = zMStar.getM();

  assert(zMStar.inZmStar(k));

  b.SetLength(m);
  for (long j = 0; j < m; j++) b[j] = 0;

  long d = deg(a);

  // compute b(X) = a(X^k) mod (X^m-1)
  mulmod_precon_t precon = PrepMulModPrecon(k, m);
  for (long j = 0; j <= d; j++) 
    b[MulModPrecon(j, k, m, precon)] = a[j]; // b[j*k mod m] = a[j]
  b.normalize();

  rem(b, b, PhimX); // reduce modulo the m'th cyclotomic
}

// A recursive matrix-by-vector multiply, used by the dense matrix code.
// This routine is optimized to use only the rotate1D routine rather
// than the more expensive linear-array rotations.
template<class type>
void EncryptedArrayDerived<type>::rec_mul(long dim, 
          Ctxt& res, const Ctxt& pdata, const vector<long>& idx,
          const PlaintextMatrixInterface<type>& mat,
          const vector<long>& dimx) const
{
  long ndims = dimension();
  long nslots = size();

  if (dim >= ndims) { // Last dimension (recursion edge condition)
    vector<RX> pmat;  // the plaintext diagonal
    pmat.resize(nslots);
    bool zDiag = true; // is this a zero diagonal
    for (long j = 0; j < nslots; j++) {
      long i = idx[j];
      RX val;
      if (mat.get(val, i, j)) // returns true if the entry is zero
        clear(pmat[j]);
      else {           // not a zero entry
        pmat[j] = val;
	zDiag = false; // not a zero diagonal
      }
    }
    if (zDiag) return; // zero diagonal, nothing to do

    // Now we have the constants for all the diagonal entries, encode the
    // diagonal as a single polynomial with these constants in the slots
    ZZX epmat;
    encode(epmat, pmat);

    // multiply by the polynomial, then add to the result
    Ctxt tmp = pdata;
    tmp.multByConstant(epmat);
    res += tmp;
  }
  else { // not the last dimension, make a recursive call
    long sdim = sizeOfDimension(dimx[dim]);

    // compute "in spirit" sum_i (pdata >> i) * i'th-diagonal, but
    // adjust the indexes so that we only need to rotate the cipehrtext
    // along the different dimensions separately
    for (long offset = 0; offset < sdim; offset++) {
      Ctxt pdata1 = pdata;
      vector<long> idx1;
      rotate1D(pdata1, dimx[dim], offset);
      this->EncryptedArrayBase::rotate1D(idx1, idx, dimx[dim], offset);
      // indexes adjusted, make the recursive call
      rec_mul(dim+1, res, pdata1, idx1, mat, dimx);
    }
  }
}


// helper class to sort dimensions, so that
//    - bad dimensions come before good dimensions (primary sort key)
//    - small dimensions come before large dimesnions (secondary sort key)
// this is a good order to process the dimensions in the recursive mat_mul_dense
// routine: it ensures that the work done at the leaves of the recursion is
// minimized, and that the work done at the non-leaves is dominated by the
// work done at the leaves.

//! \cond FALSE (make doxygen ignore these classes)
template<class type>
struct MatMulDimComp {
  const EncryptedArrayDerived<type> *ea;
  MatMulDimComp(const EncryptedArrayDerived<type> *_ea) : ea(_ea) {}

  bool operator()(long i, long j) { 
    return (!ea->nativeDimension(i) && ea->nativeDimension(j)) ||
           (  (ea->nativeDimension(i) == ea->nativeDimension(j)) &&
              (ea->sizeOfDimension(i) < ea->sizeOfDimension(j))  );
  }
};
//! \endcond

// Multiply a ciphertext vector by a plaintext dense matrix
template<class type>
void EncryptedArrayDerived<type>::mat_mul_dense(Ctxt& ctxt, const PlaintextMatrixBaseInterface& mat) const
{
  FHE_TIMER_START;
  assert(this == &mat.getEA().getDerived(type()));
  assert(&context == &ctxt.getContext());

  RBak bak; bak.save(); tab.restoreContext();

  // Get the derived type
  const PlaintextMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

  ctxt.cleanUp(); // not sure, but this may be a good idea

  Ctxt res(ctxt.getPubKey(), ctxt.getPtxtSpace()); // fresh encryption of zero
  vector<long> idx;
  idx.resize(size());
  for (long i = 0; i < size(); i++)
     idx[i] = i;

  vector<long> dimx;
  dimx.resize(dimension());
  for (long i = 0; i < dimension(); i++)
    dimx[i] = i;

  sort(dimx.begin(), dimx.end(), MatMulDimComp<type>(this));
  // sort the dimenesions so that bad ones come before good,
  // and then small ones come before large

  // call the recursive procedure to do the actual work
  rec_mul(0, res, ctxt, idx, mat1, dimx);

  ctxt = res; // copy the result back to ctxt
}

template<class type>
void EncryptedArrayDerived<type>::compMat_dense(CachedPtxtMatrix& cmat,
                               const PlaintextMatrixBaseInterface& mat) const
{
  FHE_TIMER_START;
  NTL::Error("cached compMat_dense not implemented yet");
}

template<class type>
void EncryptedArrayDerived<type>::compMat_dense(CachedDCRTPtxtMatrix& cmat,
                                const PlaintextMatrixBaseInterface& mat) const
{
  FHE_TIMER_START;
  CachedPtxtMatrix zzxMat;
  compMat_dense(zzxMat, mat);
  long n = zzxMat.length();
  cmat.SetLength(n);
  for (long i=0; i<n; i++) if (zzxMat[i])
    cmat[i] = DCRTptr(new DoubleCRT(*zzxMat[i], context));  
    // DoubleCRT defined relative to all primes, even the "special" ones
}

template<class CachedMatrix>
static void mat_mul_dense_tmpl(Ctxt& ctxt, const CachedMatrix& cmat,
			       const EncryptedArray& ea)
{
  FHE_TIMER_START;
  NTL::Error("cached mat_mul_dense not implemented yet");
}
void mat_mul_dense(Ctxt& ctxt, const CachedPtxtMatrix& cmat,
		   const EncryptedArray& ea)
{
  mat_mul_dense_tmpl(ctxt, cmat, ea);
}
void mat_mul_dense(Ctxt& ctxt, const CachedDCRTPtxtMatrix& cmat,
		   const EncryptedArray& ea)
{
  mat_mul_dense_tmpl(ctxt, cmat, ea);
}


// this mat_mul is optimized for diagonally sparse matrices

template<class type>
void EncryptedArrayDerived<type>::mat_mul(Ctxt& ctxt, const PlaintextMatrixBaseInterface& mat) const
{
  FHE_TIMER_START;
  assert(this == &mat.getEA().getDerived(type()));
  assert(&context == &ctxt.getContext());

  RBak bak; bak.save(); tab.restoreContext();

  // Get the derived type
  const PlaintextMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

  ctxt.cleanUp(); // not sure, but this may be a good idea

  Ctxt res(ctxt.getPubKey(), ctxt.getPtxtSpace()); // fresh encryption of zero

  long nslots = size();
  long d = getDegree();

  RX entry;
  vector<RX> diag;
  diag.resize(nslots);

  // Process the diagonals one at a time
  for (long i = 0; i < nslots; i++) {  // process diagonal i
    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry on this diagonal

    // Compute constants for each entry on this diagonal
    for (long j = 0; j < nslots; j++) { // process entry j
      bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j); // callback
      assert(zEntry || deg(entry) < d);

      if (!zEntry && IsZero(entry)) zEntry = true; // check for zero

      if (!zEntry) { // non-zero diagonal entry

        zDiag = false; // diagonal is non-zero

        // clear entries between last nonzero entry and this one
        for (long jj = nzLast+1; jj < j; jj++) clear(diag[jj]);
        nzLast = j;

        diag[j] = entry;
      }
    }
    
    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries
    for (long jj = nzLast+1; jj < nslots; jj++) clear(diag[jj]);

    // Now we have the constants for all the diagonal entries, encode the
    // diagonal as a single polynomial with these constants in the slots
    ZZX cpoly;
    encode(cpoly, diag);

    // rotate by i, multiply by the polynomial, then add to the result
    Ctxt shCtxt = ctxt;
    rotate(shCtxt, i); // rotate by i
    shCtxt.multByConstant(cpoly);
    res += shCtxt;
  }
  ctxt = res;
}

template<class type>
void EncryptedArrayDerived<type>::compMat(CachedPtxtMatrix& cmat,
                          const PlaintextMatrixBaseInterface& mat) const
{
  FHE_TIMER_START;
  assert(this == &mat.getEA().getDerived(type()));

  RBak bak; bak.save(); tab.restoreContext();

  // Get the derived type
  const PlaintextMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

  long nslots = size();
  long d = getDegree();
  RX entry;
  vector<RX> diag;
  diag.resize(nslots);
  cmat.SetLength(nslots);

  // Process the diagonals one at a time
  for (long i = 0; i < nslots; i++) {  // process diagonal i
    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry on this diagonal

    // Compute constants for each entry on this diagonal
    for (long j = 0; j < nslots; j++) { // process entry j
      bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j); // callback
      assert(zEntry || deg(entry) < d);

      if (!zEntry && IsZero(entry)) zEntry = true; // check for zero

      if (!zEntry) { // non-zero diagonal entry

        zDiag = false; // diagonal is non-zero

        // clear entries between last nonzero entry and this one
        for (long jj = nzLast+1; jj < j; jj++) clear(diag[jj]);
        nzLast = j;

        diag[j] = entry;
      }
    }
    
    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries
    for (long jj = nzLast+1; jj < nslots; jj++) clear(diag[jj]);

    // Now we have the constants for all the diagonal entries, encode the
    // diagonal as a single polynomial with these constants in the slots
    ZZX cpoly;
    encode(cpoly, diag);
    cmat[i] = ZZXptr(new ZZX(cpoly));
  }
}

template<class type>
void EncryptedArrayDerived<type>::compMat(CachedDCRTPtxtMatrix& cmat,
                          const PlaintextMatrixBaseInterface& mat) const
{
  FHE_TIMER_START;
  CachedPtxtMatrix zzxMat;
  compMat(zzxMat, mat);
  long n = zzxMat.length();
  cmat.SetLength(n);
  for (long i=0; i<n; i++) if (zzxMat[i])
    cmat[i] = DCRTptr(new DoubleCRT(*zzxMat[i], context));
    // DoubleCRT defined relative to all primes, even the "special" ones
}

template<class CachedMatrix>
void mat_mul_tmpl(Ctxt& ctxt, const CachedMatrix& cmat,
		  const EncryptedArray& ea)
{
  FHE_TIMER_START;
  ctxt.cleanUp(); // not sure, but this may be a good idea
  Ctxt res(ctxt.getPubKey(), ctxt.getPtxtSpace()); // fresh encryption of zero

  // Process the diagonals one at a time
  long nslots = ea.size();
  for (long i = 0; i < nslots; i++) {  // process diagonal i
    if (!cmat[i]) continue; // a zero diagonal

    // rotate by i, multiply, and add to the result
    Ctxt shCtxt = ctxt;
    ea.rotate(shCtxt, i); // rotate by i
    shCtxt.multByConstant(*cmat[i]);
    res += shCtxt;
  }
  ctxt = res;
}
void mat_mul(Ctxt& ctxt, const CachedPtxtMatrix& cmat,
	     const EncryptedArray& ea)
{
  mat_mul_tmpl(ctxt, cmat, ea);
}
void mat_mul(Ctxt& ctxt, const CachedDCRTPtxtMatrix& cmat,
	     const EncryptedArray& ea)
{
  mat_mul_tmpl(ctxt, cmat, ea);
}



template<class type, class RX>
static bool processDiagonal(vector<RX>& diag, long dim, long i,
       vector<RX>& tmpDiag, const PlaintextMatrixInterface<type>& mat,
       const PAlgebra& zMStar, long d, bool special)
{
  long D = tmpDiag.size();
  long nslots = diag.size();
  bool zDiag = true; // is this a zero diagonal?
  long nzLast = -1;  // index of last non-zero entry
  RX entry;

  // Process the entries in this diagonal one at a time
  for (long j = 0; j < D; j++) { // process entry j
    bool zEntry = mat.get(entry, mcMod(j-i, D), j); // entry [i,j-i mod D]
    assert(zEntry || deg(entry) < d);
    // get(...) returns true if the entry is empty, false otherwise

    if (!zEntry && IsZero(entry)) zEntry = true; // zero is an empty entry too

    if (!zEntry) {   // not a zero entry
      zDiag = false; // mark diagonal as non-empty

      // clear entries between last nonzero entry and this one
      for (long jj = nzLast+1; jj < j; jj++) clear(tmpDiag[jj]);
      nzLast = j;

      tmpDiag[j] = entry;
    }
  }    
  if (zDiag) return true; // zero diagonal, nothing to do

  // clear trailing zero entries
  for (long jj = nzLast+1; jj < D; jj++) clear(tmpDiag[jj]);

  if (special) diag.assign(nslots, tmpDiag[0]); // order-1 dimension
  else for (long j = 0; j < nslots; j++)
    diag[j] = tmpDiag[ zMStar.coordinate(dim,j) ];
    // rearrange the indexes based on the current dimension

  return false; // a nonzero diagonal
}


// Multiply ctx by plaintext matrix over the extension field/ring.
// Ctxt is treated as a row matrix v, and replaced by en encryption of
// v * mat' where mat' is the block-diagonal matrix defined by mat in
// dimension dim. Here, mat should represent a D x D matrix, where D is
// the order of generator dim.
// We also allow dim to be one greater than the number of generators in
// zMStar, as if there were an implicit generator of order 1, this is
// convenient in some applications.
template<class type>
void EncryptedArrayDerived<type>::mat_mul1D(Ctxt& ctxt,
     const PlaintextMatrixBaseInterface& mat, long dim) const
{
  FHE_TIMER_START;
  const PAlgebra& zMStar = context.zMStar;

  assert(this == &mat.getEA().getDerived(type()));
  assert(&context == &ctxt.getContext());
  assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator

  RBak bak; bak.save(); tab.restoreContext(); // backup the NTL modulus

  // Get the derived type
  const PlaintextMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

  ctxt.cleanUp();  // not sure, but this may be a good idea
  Ctxt res(ZeroCtxtLike, ctxt); // fresh encryption of zero

  long nslots = size();
  long d = getDegree();
  vector<RX> diag, diag1;
  diag.resize(D);
  diag1.resize(nslots);
  ZZX cpoly;

  // Process the diagonals one at a time
  for (long i = 0; i < D; i++) { // process diagonal i
    if (processDiagonal(diag1, dim, i, diag, mat1, zMStar, d, special))
      continue; // zero diagonal

    // encode as a polynomial, then multiply and add
    encode(cpoly, diag1);
    Ctxt shCtxt = ctxt;
    if (i != 0) rotate1D(shCtxt, dim, i);   
    shCtxt.multByConstant(cpoly);
    res += shCtxt;
  }

  ctxt = res;
}

template<class type>
void EncryptedArrayDerived<type>::compMat1D(CachedPtxtMatrix& cmat,
                 const PlaintextMatrixBaseInterface& mat, long dim) const
{
  FHE_TIMER_START;
  const PAlgebra& zMStar = context.zMStar;

  assert(this == &mat.getEA().getDerived(type()));
  assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator

  RBak bak; bak.save(); tab.restoreContext(); // backup the NTL modulus

  // Get the derived type
  const PlaintextMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextMatrixInterface<type>& >( mat );

  long nslots = size();
  long d = getDegree();
  vector<RX> diag, diag1;
  diag.resize(D);
  diag1.resize(nslots);
  cmat.SetLength(D);
  ZZX cpoly;

  // Process the diagonals one at a time
  for (long i = 0; i < D; i++) { // process diagonal i
    if (processDiagonal(diag1, dim, i, diag, mat1, zMStar, d, special))
      continue; // zero diagonal

    // encode as a polynomial, then multiply and add
    encode(cpoly, diag1);
    cmat[i] = ZZXptr(new ZZX(cpoly));
  }
}

template<class type>
void EncryptedArrayDerived<type>::compMat1D(CachedDCRTPtxtMatrix& cmat,
                     const PlaintextMatrixBaseInterface& mat, long dim) const
{
  FHE_TIMER_START;
  CachedPtxtMatrix zzxMat;
  compMat1D(zzxMat, mat, dim);
  long n = zzxMat.length();
  cmat.SetLength(n);
  for (long i=0; i<n; i++) if (zzxMat[i])
    cmat[i] = DCRTptr(new DoubleCRT(*zzxMat[i], context));
    // DoubleCRT defined relative to all primes, even the "special" ones
}

#ifdef FHE_BOOT_THREADS

template<class Matrix, class EA>
static void mat_mul1D_tmpl(Ctxt& ctxt, const Matrix& cmat, long dim,
			   const EA& ea)
{
  FHE_TIMER_START;
  const FHEcontext& context = ctxt.getContext();
  const PAlgebra& zMStar = context.zMStar;
  assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator

  ctxt.cleanUp();  // not sure, but this may be a good idea
  Ctxt res(ZeroCtxtLike, ctxt); // fresh encryption of zero

  Vec< shared_ptr<Ctxt> > tvec;
  tvec.SetLength(D);
  for (long i = 0; i < D; i++)
    tvec[i] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

  bootTask->exec1(D,
    [&](long first, long last) {
      for (long i = first; i < last; i++) { // process diagonal i
        if (!cmat[i]) continue;      // zero diagonal
    
        // rotate and multiply
        (*tvec[i]) = ctxt;
        if (i != 0) ea.rotate1D(*tvec[i], dim, i);   
        tvec[i]->multByConstant(*cmat[i]);
      }
    }
  );

  for (long i = 0; i < D; i++)
    res += *tvec[i];

  ctxt = res;
}

#else

template<class Matrix, class EA>
static void mat_mul1D_tmpl(Ctxt& ctxt, const Matrix& cmat, long dim,
			   const EA& ea)
{
  FHE_TIMER_START;
  const FHEcontext& context = ctxt.getContext();
  const PAlgebra& zMStar = context.zMStar;
  assert(dim >= 0 && dim <= LONG(zMStar.numOfGens()));

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator

  ctxt.cleanUp();  // not sure, but this may be a good idea
  Ctxt res(ZeroCtxtLike, ctxt); // fresh encryption of zero

  // Process the diagonals one at a time
  for (long i = 0; i < D; i++) { // process diagonal i
    if (!cmat[i]) continue;      // zero diagonal

    // rotate, multiply and add
    Ctxt shCtxt = ctxt;
    if (i != 0) ea.rotate1D(shCtxt, dim, i);   
    shCtxt.multByConstant(*cmat[i]);
    res += shCtxt;
  }
  ctxt = res;
}


#endif



void mat_mul1D(Ctxt& ctxt, const CachedPtxtMatrix& cmat, long dim,
	       const EncryptedArray& ea)
{ mat_mul1D_tmpl(ctxt, cmat, dim, ea); }
void mat_mul1D(Ctxt& ctxt, const CachedPtxtMatrix& cmat, long dim,
	       const EncryptedArrayDerived<PA_GF2>& ea)
{ mat_mul1D_tmpl(ctxt, cmat, dim, ea); }
void mat_mul1D(Ctxt& ctxt, const CachedPtxtMatrix& cmat, long dim,
	       const EncryptedArrayDerived<PA_zz_p>& ea)
{ mat_mul1D_tmpl(ctxt, cmat, dim, ea); }

void mat_mul1D(Ctxt& ctxt, const CachedDCRTPtxtMatrix& cmat, long dim,
	       const EncryptedArray& ea)
{ mat_mul1D_tmpl(ctxt, cmat, dim, ea); }
void mat_mul1D(Ctxt& ctxt, const CachedDCRTPtxtMatrix& cmat, long dim,
	       const EncryptedArrayDerived<PA_GF2>& ea)
{ mat_mul1D_tmpl(ctxt, cmat, dim, ea); }
void mat_mul1D(Ctxt& ctxt, const CachedDCRTPtxtMatrix& cmat, long dim,
	       const EncryptedArrayDerived<PA_zz_p>& ea)
{ mat_mul1D_tmpl(ctxt, cmat, dim, ea); }


// This code has a complexity of N+d (instead of N*d) where N is the number of
// nonzero diagonal blocks. However, it requires space for d extra ciphertexts
template<class type>
void EncryptedArrayDerived<type>::mat_mul(Ctxt& ctxt, const PlaintextBlockMatrixBaseInterface& mat) const
{
  FHE_TIMER_START;
  assert(this == &mat.getEA().getDerived(type()));
  assert(&context == &ctxt.getContext());

  const PAlgebra& zMStar = context.zMStar;
  long p = zMStar.getP(); 
  long m = zMStar.getM();
  const RXModulus& F = tab.getPhimXMod();

  RBak bak; bak.save(); tab.restoreContext();

  const PlaintextBlockMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

  ctxt.cleanUp(); // not sure, but this may be a good idea

  long nslots = size();
  long d = getDegree();

  Vec< shared_ptr<Ctxt> > acc;
  acc.SetLength(d);
  for (long k = 0; k < d; k++)
    acc[k] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

  mat_R entry;
  entry.SetDims(d, d);

  vector<RX> entry1;
  entry1.resize(d);
  
  vector< vector<RX> > diag;
  diag.resize(nslots);
  for (long j = 0; j < nslots; j++) diag[j].resize(d);

  for (long i = 0; i < nslots; i++) { // process diagonal i
    bool zDiag = true;
    long nzLast = -1;

    for (long j = 0; j < nslots; j++) {
      bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j);
      assert(zEntry || (entry.NumRows() == d && entry.NumCols() == d));
      // get(...) returns true if the entry is empty, false otherwise

      if (!zEntry && IsZero(entry)) zEntry=true; // zero is an empty entry too

      if (!zEntry) {    // non-empty entry
        zDiag = false;  // mark diagonal as non-empty

        // clear entries between last nonzero entry and this one
        for (long jj = nzLast+1; jj < j; jj++) {
          for (long k = 0; k < d; k++)
            clear(diag[jj][k]);
        }
        nzLast = j;

        // recode entry as a vector of polynomials
        for (long k = 0; k < d; k++) conv(entry1[k], entry[k]);

        // compute the lin poly coeffs
        buildLinPolyCoeffs(diag[j], entry1);
      }
    }
    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries    
    for (long jj = nzLast+1; jj < nslots; jj++) {
      for (long k = 0; k < d; k++)
        clear(diag[jj][k]);
    }

    // now diag[j] contains the lin poly coeffs

    Ctxt shCtxt = ctxt;
    rotate(shCtxt, i); 
    shCtxt.cleanUp();

    RX cpoly1, cpoly2;
    ZZX cpoly;

    // apply the linearlized polynomial
    for (long k = 0; k < d; k++) {

      // compute the constant
      bool zConst = true;
      vector<RX> cvec;
      cvec.resize(nslots);
      for (long j = 0; j < nslots; j++) {
        cvec[j] = diag[j][k];
        if (!IsZero(cvec[j])) zConst = false;
      }

      if (zConst) continue;

      encode(cpoly, cvec);
      conv(cpoly1, cpoly);

      // apply inverse automorphism to constant
      plaintextAutomorph(cpoly2, cpoly1, PowerMod(p, mcMod(-k, d), m), zMStar, F);
      conv(cpoly, cpoly2);
      Ctxt shCtxt1 = shCtxt;
      shCtxt1.multByConstant(cpoly);
      *acc[k] += shCtxt1;
    }
  }

  Ctxt res(ZeroCtxtLike, ctxt);

  for (long k = 0; k < d; k++) {
    acc[k]->frobeniusAutomorph(k);
    res += *acc[k];
  }

  ctxt = res;
}

template<class type>
void EncryptedArrayDerived<type>::compMat(CachedPtxtBlockMatrix& cmat,
			  const PlaintextBlockMatrixBaseInterface& mat) const
{
  FHE_TIMER_START;
  assert(this == &mat.getEA().getDerived(type()));
  const PAlgebra& zMStar = context.zMStar;
  long p = zMStar.getP(); 
  long m = zMStar.getM();
  const RXModulus& F = tab.getPhimXMod();

  RBak bak; bak.save(); tab.restoreContext();

  // Get the derived type
  const PlaintextBlockMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

  long nslots = size();
  long d = getDegree();

  mat_R entry;
  entry.SetDims(d, d);
  vector<RX> entry1;
  entry1.resize(d);
  
  vector< vector<RX> > diag;
  diag.resize(nslots);
  for (long j = 0; j < nslots; j++) diag[j].resize(d);
  cmat.SetDims(nslots, d);

  for (long i = 0; i < nslots; i++) { // process diagonal i
    bool zDiag = true;
    long nzLast = -1;

    for (long j = 0; j < nslots; j++) {
      bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j);
      assert(zEntry || (entry.NumRows() == d && entry.NumCols() == d));
      // get(...) returns true if the entry is empty, false otherwise

      if (!zEntry && IsZero(entry)) zEntry=true; // zero is an empty entry too

      if (!zEntry) {    // non-empty entry
        zDiag = false;  // mark diagonal as non-empty

        // clear entries between last nonzero entry and this one
        for (long jj = nzLast+1; jj < j; jj++) {
          for (long k = 0; k < d; k++)
            clear(diag[jj][k]);
        }
        nzLast = j;

        // recode entry as a vector of polynomials
        for (long k = 0; k < d; k++) conv(entry1[k], entry[k]);

        // compute the lin poly coeffs
        buildLinPolyCoeffs(diag[j], entry1);
      }
    }
    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries    
    for (long jj = nzLast+1; jj < nslots; jj++) {
      for (long k = 0; k < d; k++)
        clear(diag[jj][k]);
    }
    // now diag[j] contains the lin poly coeffs

    RX cpoly1, cpoly2;
    ZZX cpoly;

    // apply the linearlized polynomial
    for (long k = 0; k < d; k++) {

      // compute the constant
      bool zConst = true;
      vector<RX> cvec;
      cvec.resize(nslots);
      for (long j = 0; j < nslots; j++) {
        cvec[j] = diag[j][k];
        if (!IsZero(cvec[j])) zConst = false;
      }
      if (zConst) continue;

      encode(cpoly, cvec);
      conv(cpoly1, cpoly);

      // apply inverse automorphism to constant
      plaintextAutomorph(cpoly2, cpoly1, PowerMod(p, mcMod(-k, d), m), zMStar, F);
      conv(cpoly, cpoly2);
      cmat[i][k] = ZZXptr(new ZZX(cpoly));
    }
  }
}

template<class type>
void EncryptedArrayDerived<type>::compMat(CachedDCRTPtxtBlockMatrix& cmat,
                          const PlaintextBlockMatrixBaseInterface& mat) const
{
  FHE_TIMER_START;
  CachedPtxtBlockMatrix zzxMat;
  compMat(zzxMat, mat);
  long m = zzxMat.NumRows();
  long n = zzxMat.NumCols();
  cmat.SetDims(m,n);
  for (long i=0; i<m; i++) for (long j=0; j<n; j++) if (zzxMat[i][j])
    cmat[i][j] = DCRTptr(new DoubleCRT(*zzxMat[i][j], context, context.ctxtPrimes));
    // DoubleCRT defined relative only to the ciphertxt primes, not the "special" ones

}

template<class CachedMatrix>
void blockMat_mul_tmpl(Ctxt& ctxt, const CachedMatrix& cmat,
		       const EncryptedArray& ea)
{
  FHE_TIMER_START;
  ctxt.cleanUp(); // not sure, but this may be a good idea

  long nslots = ea.size();
  long d = ea.getDegree();

  Vec< shared_ptr<Ctxt> > acc;
  acc.SetLength(d);
  for (long k = 0; k < d; k++)
    acc[k] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

  for (long i = 0; i < nslots; i++) { // process diagonal i
    // apply the linearlized polynomial
    Ctxt shCtxt = ctxt;
    ea.rotate(shCtxt, i); 
    shCtxt.cleanUp();

    for (long k = 0; k < d; k++) {
      if (!cmat[i][k]) continue; // a zero constant
      Ctxt shCtxt1 = shCtxt;
      shCtxt1.multByConstant(*cmat[i][k]);
      *acc[k] += shCtxt1;
    }
  }

  Ctxt res(ZeroCtxtLike, ctxt);
  for (long k = 0; k < d; k++) {
    acc[k]->frobeniusAutomorph(k);
    res += *acc[k];
  }
  ctxt = res;
}
void mat_mul(Ctxt& ctxt, const CachedPtxtBlockMatrix& cmat,
	     const EncryptedArray& ea)
{
  blockMat_mul_tmpl(ctxt, cmat, ea);
}
void mat_mul(Ctxt& ctxt, const CachedDCRTPtxtBlockMatrix& cmat,
	     const EncryptedArray& ea)
{
  blockMat_mul_tmpl(ctxt, cmat, ea);
}

// Multiply ctx by plaintext matrix over the base field/ring.
// Ctxt is treated as a row matrix v, and replaced by en encryption of
// v * mat' where mat' is the block-diagonal matrix defined by mat in
// dimension dim. Here, mat should represent a D x D matrix, where D is
// the order of generator dim.
// We also allow dim to be one greater than the number of generators in
// zMStar, as if there were an implicit generator of order 1, this is
// convenient in some applications.
template<class type> void EncryptedArrayDerived<type>::mat_mul1D(Ctxt& ctxt,
               const PlaintextBlockMatrixBaseInterface& mat, long dim) const
{
  FHE_TIMER_START;
  const PAlgebra& zMStar = context.zMStar;

  assert(this == &mat.getEA().getDerived(type()));
  assert(&context == &ctxt.getContext());
  assert(dim >= 0 && dim <=  LONG(zMStar.numOfGens()));

  long p = zMStar.getP(); 
  long m = zMStar.getM();
  const RXModulus& F = tab.getPhimXMod();

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator
  bool bad = !special && !zMStar.SameOrd(dim);

  RBak bak; bak.save(); tab.restoreContext(); // backup the NTL modulus

  // Get the derived type
  const PlaintextBlockMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

  ctxt.cleanUp(); // not sure, but this may be a good idea

  long nslots = size();
  long d = getDegree();

  Vec< shared_ptr<Ctxt> > acc;
  acc.SetLength(d);
  for (long k = 0; k < d; k++)
    acc[k] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

  mat_R entry;
  entry.SetDims(d, d);

  vector<RX> entry1;
  entry1.resize(d);
  
  vector< vector<RX> > diag;
  diag.resize(D);
  for (long j = 0; j < D; j++) diag[j].resize(d);

  // Process the diagonals one at a time
  for (long i = 0; i < D; i++) { // process diagonal i
    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry

    // Process the entries in this diagonal one at a time
    for (long j = 0; j < D; j++) { // process entry j
      bool zEntry = mat1.get(entry, mcMod(j-i, D), j); // entry [i,j-i mod D]
      assert(zEntry || (entry.NumRows() == d && entry.NumCols() == d));
        // get(...) returns true if the entry is empty, false otherwise

      if (!zEntry && IsZero(entry)) zEntry=true; // zero is an empty entry too

      if (!zEntry) {    // non-empty entry
        zDiag = false;  // mark diagonal as non-empty

        // clear entries between last nonzero entry and this one
        for (long jj = nzLast+1; jj < j; jj++) {
          for (long k = 0; k < d; k++)
            clear(diag[jj][k]);
        }
        nzLast = j;

        // recode entry as a vector of polynomials
        for (long k = 0; k < d; k++) conv(entry1[k], entry[k]);

        // compute the lin poly coeffs
        buildLinPolyCoeffs(diag[j], entry1);
      }
    }
    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries    
    for (long jj = nzLast+1; jj < D; jj++) {
      for (long k = 0; k < d; k++)
        clear(diag[jj][k]);
    }

    // now diag[j] contains the lin poly coeffs

    vector<Ctxt> shCtxt;
    vector<RX> shMask;
    if (i == 0) {
      shCtxt.resize(1, ctxt);
      shMask.resize(1, conv<RX>(1));
    }
    else if (!bad) {
      shCtxt.resize(1, ctxt);
      shMask.resize(1, conv<RX>(1));
      rotate1D(shCtxt[0], dim, i);
      shCtxt[0].cleanUp();
    }
    else {
      // we fold the masking constants into the linearized polynomial
      // constants to save a level. We lift some code out of rotate1D
      // to do this.

      shCtxt.resize(2, ctxt);
      shMask.resize(2);

      long val = PowerMod(zMStar.ZmStarGen(dim), i, m);
      long ival = PowerMod(zMStar.ZmStarGen(dim), i-D, m);
      const RX& mask = tab.getMaskTable()[dim][D-i];

      shCtxt[0].smartAutomorph(val); 
      shCtxt[0].cleanUp();

      shCtxt[1].smartAutomorph(ival); 
      shCtxt[1].cleanUp();

      plaintextAutomorph(shMask[0], 1 - mask, val, zMStar, F);
      plaintextAutomorph(shMask[1], mask, ival, zMStar, F);
    }
    
    RX cpoly1, cpoly2, cpoly3;
    ZZX cpoly;

    // apply the linearlized polynomial
    for (long k = 0; k < d; k++) {

      // compute the constant
      bool zConst = true;
      vector<RX> cvec;
      cvec.resize(nslots);
      for (long j = 0; j < nslots; j++) {
        cvec[j] = diag[ special ? 0 : zMStar.coordinate(dim, j) ][k];
        if (!IsZero(cvec[j])) zConst = false;
      }

      if (zConst) continue;

      encode(cpoly, cvec);
      conv(cpoly1, cpoly);

      // apply inverse automorphism to constant
      plaintextAutomorph(cpoly2,cpoly1, PowerMod(p, mcMod(-k,d), m), zMStar, F);

      for (long j = 0; j < LONG(shCtxt.size()); j++) {
        MulMod(cpoly3, cpoly2, shMask[j], F);
        conv(cpoly, cpoly3);
        Ctxt shCtxt1 = shCtxt[j];;
        shCtxt1.multByConstant(cpoly);
        *acc[k] += shCtxt1;
      }
    }
  }

  Ctxt res(ZeroCtxtLike, ctxt);

  for (long k = 0; k < d; k++) {
    acc[k]->frobeniusAutomorph(k);
    res += *acc[k];
  }

  ctxt = res;
}


template<class type>
void EncryptedArrayDerived<type>::compMat1D(CachedPtxtBlockMatrix& cmat,
                 const PlaintextBlockMatrixBaseInterface& mat, long dim) const
{
  FHE_TIMER_START;
  const PAlgebra& zMStar = context.zMStar;

  assert(this == &mat.getEA().getDerived(type()));
  assert(dim >= 0 && dim <=  LONG(zMStar.numOfGens()));

  long p = zMStar.getP(); 
  long m = zMStar.getM();
  const RXModulus& F = tab.getPhimXMod();

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator
  bool bad = !special && !zMStar.SameOrd(dim);

  long nslots = size();
  long d = getDegree();

  RBak bak; bak.save(); tab.restoreContext(); // backup the NTL modulus

  // Get the derived type
  const PlaintextBlockMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

  mat_R entry;
  entry.SetDims(d, d);

  vector<RX> entry1;
  entry1.resize(d);
  
  if (bad)
    cmat.SetDims(D, 2*d);
  else
    cmat.SetDims(D, d);

  vector< vector<RX> > diag;
  diag.resize(D);
  for (long j = 0; j < D; j++) diag[j].resize(d);

  // Process the diagonals one at a time
  for (long i = 0; i < D; i++) { // process diagonal i
    bool zDiag = true; // is this a zero diagonal?
    long nzLast = -1;  // index of last non-zero entry

    // Process the entries in this diagonal one at a time
    for (long j = 0; j < D; j++) { // process entry j
      bool zEntry = mat1.get(entry, mcMod(j-i, D), j); // entry [i,j-i mod D]
      assert(zEntry || (entry.NumRows() == d && entry.NumCols() == d));
      // get(...) returns true if the entry is empty, false otherwise

      if (!zEntry && IsZero(entry)) zEntry=true; // zero is an empty entry too

      if (!zEntry) {    // non-empty entry
        zDiag = false;  // mark diagonal as non-empty

        // clear entries between last nonzero entry and this one
        for (long jj = nzLast+1; jj < j; jj++) {
          for (long k = 0; k < d; k++)
            clear(diag[jj][k]);
        }
        nzLast = j;

        // recode entry as a vector of polynomials
        for (long k = 0; k < d; k++) conv(entry1[k], entry[k]);

        // compute the lin poly coeffs
        buildLinPolyCoeffs(diag[j], entry1);
      }
    }
    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries    
    for (long jj = nzLast+1; jj < D; jj++) {
      for (long k = 0; k < d; k++)
        clear(diag[jj][k]);
    }

    // now diag[j] contains the lin poly coeffs

    vector<RX> shMask;
    if (i == 0) {
      shMask.resize(1, conv<RX>(1));
    }
    else if (!bad) {
      shMask.resize(1, conv<RX>(1));
    }
    else {
      // we fold the masking constants into the linearized polynomial
      // constants to save a level. We lift some code out of rotate1D
      // to do this.

      shMask.resize(2);

      long val = PowerMod(zMStar.ZmStarGen(dim), i, m);
      long ival = PowerMod(zMStar.ZmStarGen(dim), i-D, m);
      const RX& mask = tab.getMaskTable()[dim][D-i];

      plaintextAutomorph(shMask[0], 1 - mask, val, zMStar, F);
      plaintextAutomorph(shMask[1], mask, ival, zMStar, F);
    }
    
    RX cpoly1, cpoly2, cpoly3;
    ZZX cpoly;

    // apply the linearlized polynomial
    for (long k = 0; k < d; k++) {

      // compute the constant
      bool zConst = true;
      vector<RX> cvec;
      cvec.resize(nslots);
      for (long j = 0; j < nslots; j++) {
        cvec[j] = diag[ special ? 0 : zMStar.coordinate(dim, j) ][k];
        if (!IsZero(cvec[j])) zConst = false;
      }

      if (zConst) continue;

      encode(cpoly, cvec);
      conv(cpoly1, cpoly);

      // apply inverse automorphism to constant
      plaintextAutomorph(cpoly2,cpoly1, PowerMod(p, mcMod(-k,d), m), zMStar, F);

      for (long j = 0; j < LONG(shMask.size()); j++) {
        MulMod(cpoly3, cpoly2, shMask[j], F);
        conv(cpoly, cpoly3);
        cmat[i][k + d*j] = ZZXptr(new ZZX(cpoly));
      }
    }
  }
}

template<class type>
void EncryptedArrayDerived<type>::compMat1D(CachedDCRTPtxtBlockMatrix& cmat,
              const PlaintextBlockMatrixBaseInterface& mat, long dim) const
{
  FHE_TIMER_START;
  CachedPtxtBlockMatrix zzxMat;
  compMat1D(zzxMat, mat, dim);
  long m = zzxMat.NumRows();
  long n = zzxMat.NumCols();
  cmat.SetDims(m,n);
  for (long i=0; i<m; i++) for (long j=0; j<n; j++) if (zzxMat[i][j])
    cmat[i][j] = DCRTptr(new DoubleCRT(*zzxMat[i][j], context, context.ctxtPrimes));
    // DoubleCRT defined relative only to the ciphertxt primes, not the "special" ones
}

#ifdef FHE_BOOT_THREADS

template<class Matrix, class EA>
static void blockMat_mul1D_tmpl(Ctxt& ctxt, const Matrix& cmat, long dim,
				const EA& ea)
{
  FHE_TIMER_START;
  const FHEcontext& context = ctxt.getContext();
  const PAlgebra& zMStar = context.zMStar;
  assert(dim >= 0 && dim <=  LONG(zMStar.numOfGens()));

  long m = zMStar.getM();

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator
  bool bad = !special && !zMStar.SameOrd(dim);
  long d = ea.getDegree();

  Vec< shared_ptr<Ctxt> > acc;
  acc.SetLength(d);
  ctxt.cleanUp(); // not sure, but this may be a good idea
  for (long k = 0; k < d; k++)
    acc[k] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

  FHE_NTIMER_START(blockMat1);
  vector< vector<Ctxt> > shCtxt;

  shCtxt.resize(D);

  // Process the diagonals one at a time
  bootTask->exec1(D,
    [&](long first, long last) {
      for (long i = first; i < last; i++) { // process diagonal i
        if (i == 0) {
          shCtxt[i].resize(1, ctxt);
        }
        else if (!bad) {
          shCtxt[i].resize(1, ctxt);
          ea.rotate1D(shCtxt[i][0], dim, i);
          shCtxt[i][0].cleanUp();
        }
        else {
          // we fold the masking constants into the linearized polynomial
          // constants to save a level. We lift some code out of rotate1D
          // to do this.
    
          shCtxt[i].resize(2, ctxt);
    
          long val = PowerMod(zMStar.ZmStarGen(dim), i, m);
          long ival = PowerMod(zMStar.ZmStarGen(dim), i-D, m);
    
          shCtxt[i][0].smartAutomorph(val); 
          shCtxt[i][0].cleanUp();
    
          shCtxt[i][1].smartAutomorph(ival); 
          shCtxt[i][1].cleanUp();
        }
      }
    }
  );
  FHE_NTIMER_STOP(blockMat1);


  FHE_NTIMER_START(blockMat3);
  bootTask->exec1(d,
    [&](long first, long last) {
      for (long k = first; k < last; k++) {
        for (long i = 0; i < D; i++) {
          for (long j = 0; j < LONG(shCtxt[i].size()); j++) {
            if (!cmat[i][k + d*j]) continue; // zero constant
    
            Ctxt shCtxt1 = shCtxt[i][j];
            shCtxt1.multByConstant(*cmat[i][k + d*j]);
            *acc[k] += shCtxt1;
          }
        }
        acc[k]->frobeniusAutomorph(k);
      }
    }
  );
  FHE_NTIMER_STOP(blockMat3);

  FHE_NTIMER_START(blockMat4);
  Ctxt res(ZeroCtxtLike, ctxt);
  for (long k = 0; k < d; k++) {
    res += *acc[k];
  }
  FHE_NTIMER_STOP(blockMat4);

  ctxt = res;
}


#else


template<class Matrix, class EA>
static void blockMat_mul1D_tmpl(Ctxt& ctxt, const Matrix& cmat, long dim,
				const EA& ea)
{
  FHE_TIMER_START;
  const FHEcontext& context = ctxt.getContext();
  const PAlgebra& zMStar = context.zMStar;
  assert(dim >= 0 && dim <=  LONG(zMStar.numOfGens()));

  long m = zMStar.getM();

  // special case fo the extra dimension
  bool special = (dim == LONG(zMStar.numOfGens()));
  long D = special ? 1 : zMStar.OrderOf(dim); // order of current generator
  bool bad = !special && !zMStar.SameOrd(dim);
  long d = ea.getDegree();

  Vec< shared_ptr<Ctxt> > acc;
  acc.SetLength(d);
  ctxt.cleanUp(); // not sure, but this may be a good idea
  for (long k = 0; k < d; k++)
    acc[k] = shared_ptr<Ctxt>(new Ctxt(ZeroCtxtLike, ctxt));

  FHE_NTIMER_START(blockMat1);
  // Process the diagonals one at a time
  for (long i = 0; i < D; i++) { // process diagonal i
    vector<Ctxt> shCtxt;
    if (i == 0) {
      shCtxt.resize(1, ctxt);
    }
    else if (!bad) {
      shCtxt.resize(1, ctxt);
      ea.rotate1D(shCtxt[0], dim, i);
      shCtxt[0].cleanUp();
    }
    else {
      // we fold the masking constants into the linearized polynomial
      // constants to save a level. We lift some code out of rotate1D
      // to do this.

      shCtxt.resize(2, ctxt);

      long val = PowerMod(zMStar.ZmStarGen(dim), i, m);
      long ival = PowerMod(zMStar.ZmStarGen(dim), i-D, m);

      shCtxt[0].smartAutomorph(val); 
      shCtxt[0].cleanUp();

      shCtxt[1].smartAutomorph(ival); 
      shCtxt[1].cleanUp();
    }

    // apply the linearlized polynomial
    for (long k = 0; k < d; k++) {
      for (long j = 0; j < LONG(shCtxt.size()); j++) {
        if (!cmat[i][k + d*j]) continue; // zero constant

        Ctxt shCtxt1 = shCtxt[j];
        shCtxt1.multByConstant(*cmat[i][k + d*j]);
        *acc[k] += shCtxt1;
      }
    }
  }
  FHE_NTIMER_STOP(blockMat1);

  FHE_NTIMER_START(blockMat2);
  Ctxt res(ZeroCtxtLike, ctxt);
  for (long k = 0; k < d; k++) {
    acc[k]->frobeniusAutomorph(k);
    res += *acc[k];
  }
  FHE_NTIMER_STOP(blockMat2);

  ctxt = res;
}


#endif



void mat_mul1D(Ctxt& ctxt, const CachedPtxtBlockMatrix& cmat, long dim,
	       const EncryptedArray& ea)
{ blockMat_mul1D_tmpl(ctxt, cmat, dim, ea); }
void mat_mul1D(Ctxt& ctxt, const CachedPtxtBlockMatrix& cmat, long dim,
	       const EncryptedArrayDerived<PA_GF2>& ea)
{ blockMat_mul1D_tmpl(ctxt, cmat, dim, ea); }
void mat_mul1D(Ctxt& ctxt, const CachedPtxtBlockMatrix& cmat, long dim,
	       const EncryptedArrayDerived<PA_zz_p>& ea)
{ blockMat_mul1D_tmpl(ctxt, cmat, dim, ea); }

void mat_mul1D(Ctxt& ctxt, const CachedDCRTPtxtBlockMatrix& cmat, long dim,
	       const EncryptedArray& ea)
{ blockMat_mul1D_tmpl(ctxt, cmat, dim, ea); }
void mat_mul1D(Ctxt& ctxt, const CachedDCRTPtxtBlockMatrix& cmat, long dim,
	       const EncryptedArrayDerived<PA_GF2>& ea)
{ blockMat_mul1D_tmpl(ctxt, cmat, dim, ea); }
void mat_mul1D(Ctxt& ctxt, const CachedDCRTPtxtBlockMatrix& cmat, long dim,
	       const EncryptedArrayDerived<PA_zz_p>& ea)
{ blockMat_mul1D_tmpl(ctxt, cmat, dim, ea); }


// Linearized polynomials.
// L describes a linear map M by describing its action on the standard
// power basis: M(x^j mod G) = (L[j] mod G), for j = 0..d-1.  
// The result is a coefficient vector C for the linearized polynomial
// representing M: a polynoamial h in Z/(p^r)[X] of degree < d is sent to
//
//    M(h(X) \bmod G)= \sum_{i=0}^{d-1}(C[j] \cdot h(X^{p^j}))\bmod G).
template<class type> void
EncryptedArrayDerived<type>::buildLinPolyCoeffs(vector<ZZX>& C, 
						const vector<ZZX>& L) const
{
  RBak bak; bak.save(); restoreContext();
  vector<RX> CC, LL;
  convert(LL, L);
  buildLinPolyCoeffs(CC, LL);
  convert(C, CC);
}

template<class type> void
EncryptedArrayDerived<type>::buildLinPolyCoeffs(vector<RX>& C, 
						const vector<RX>& L) const
{
  FHE_TIMER_START;

  RBak bak; bak.save(); restoreContext();  // the NTL context for mod p^r
  REBak ebak; ebak.save(); restoreContextForG(); // The NTL context for mod G

  do {
    typename Lazy< Mat<RE> >::Builder builder(linPolyMatrix);
    if (!builder()) break;

   
    long p = tab.getZMStar().getP();
    long r = tab.getR();

    Mat<RE> M1;
    // build d x d matrix, d is taken from the surrent NTL context for G
    buildLinPolyMatrix(M1, p);
    Mat<RE> M2;
    ppInvert(M2, M1, p, r); // invert modulo prime-power p^r

    UniquePtr< Mat<RE> > ptr;
    ptr.make(M2);
    builder.move(ptr);
  } while (0);

  Vec<RE> CC, LL;
  convert(LL, L);
  mul(CC, LL, *linPolyMatrix);
  convert(C, CC);
}


// Apply the same linear transformation to all the slots.
// C is the output of ea.buildLinPolyCoeffs
void applyLinPoly1(const EncryptedArray& ea, Ctxt& ctxt, const vector<ZZX>& C)
{
  assert(&ea.getContext() == &ctxt.getContext());
  long d = ea.getDegree();
  assert(d == lsize(C));

  long nslots = ea.size();

  vector<ZZX> encodedC(d);
  for (long j = 0; j < d; j++) {
    vector<ZZX> v(nslots);
    for (long i = 0; i < nslots; i++) v[i] = C[j];
    ea.encode(encodedC[j], v);
  }

  applyLinPolyLL(ctxt, encodedC, ea.getDegree());
}


// Apply different transformations to different slots. Cvec is a vector of
// length ea.size(), with each entry the output of ea.buildLinPolyCoeffs
void applyLinPolyMany(const EncryptedArray& ea, Ctxt& ctxt, 
                      const vector< vector<ZZX> >& Cvec)
{
  assert(&ea.getContext() == &ctxt.getContext());
  long d = ea.getDegree();
  long nslots = ea.size();

  assert(nslots == lsize(Cvec));
  for (long i = 0; i < nslots; i++)
    assert(d == lsize(Cvec[i]));

  vector<ZZX> encodedC(d);
  for (long j = 0; j < d; j++) {
    vector<ZZX> v(nslots);
    for (long i = 0; i < nslots; i++) v[i] = Cvec[i][j];
    ea.encode(encodedC[j], v);
  }

  applyLinPolyLL(ctxt, encodedC, ea.getDegree());
}

// A low-level variant: encodedCoeffs has all the linPoly coeffs encoded
// in slots; different transformations can be encoded in different slots
template<class P>
void applyLinPolyLL(Ctxt& ctxt, const vector<P>& encodedC, long d)
{
  assert(d == lsize(encodedC));

  ctxt.cleanUp();  // not sure, but this may be a good idea

  Ctxt tmp(ctxt);

  ctxt.multByConstant(encodedC[0]);
  for (long j = 1; j < d; j++) {
    Ctxt tmp1(tmp);
    tmp1.frobeniusAutomorph(j);
    tmp1.multByConstant(encodedC[j]);
    ctxt += tmp1;
  }
}
template void applyLinPolyLL(Ctxt& ctxt, const vector<ZZX>& encodedC, long d);
template void applyLinPolyLL(Ctxt& ctxt, const vector<DoubleCRT>& encodedC, long d);

#if 0
/************************* OLD UNUSED CODE *************************/
template<class type>
void EncryptedArrayDerived<type>::mat_mul(Ctxt& ctxt, const PlaintextBlockMatrixBaseInterface& mat) const
{
  FHE_TIMER_START;
  assert(this == &mat.getEA().getDerived(type()));
  assert(&context == &ctxt.getContext());

  RBak bak; bak.save(); tab.restoreContext();

  const PlaintextBlockMatrixInterface<type>& mat1 = 
    dynamic_cast< const PlaintextBlockMatrixInterface<type>& >( mat );

  ctxt.cleanUp(); // not sure, but this may be a good idea

  Ctxt res(ctxt.getPubKey(), ctxt.getPtxtSpace());
  // a new ciphertext, encrypting zero
  

  long nslots = size();
  long d = getDegree();

  mat_R entry;
  entry.SetDims(d, d);

  vector<RX> entry1;
  entry1.resize(d);
  
  vector< vector<RX> > diag;
  diag.resize(nslots);
  for (long j = 0; j < nslots; j++) diag[j].resize(d);

  for (long i = 0; i < nslots; i++) {
    // process diagonal i


    bool zDiag = true;
    long nzLast = -1;

    for (long j = 0; j < nslots; j++) {
      bool zEntry = mat1.get(entry, mcMod(j-i, nslots), j);
      assert(zEntry || (entry.NumRows() == d && entry.NumCols() == d));
        // get(...) returns true if the entry is empty, false otherwise

      if (!zEntry && IsZero(entry)) zEntry=true; // zero is an empty entry too

      if (!zEntry) {    // non-empty entry

        zDiag = false;  // mark diagonal as non-empty

        // clear entries between last nonzero entry and this one

        for (long jj = nzLast+1; jj < j; jj++) {
          for (long k = 0; k < d; k++)
            clear(diag[jj][k]);
        }

        nzLast = j;

        // recode entry as a vector of polynomials
        for (long k = 0; k < d; k++) conv(entry1[k], entry[k]);

        // compute the lin poly coeffs
        buildLinPolyCoeffs(diag[j], entry1);
      }
    }

    if (zDiag) continue; // zero diagonal, continue

    // clear trailing zero entries    
    for (long jj = nzLast+1; jj < nslots; jj++) {
      for (long k = 0; k < d; k++)
        clear(diag[jj][k]);
    }

    // now diag[j] contains the lin poly coeffs

    Ctxt shCtxt = ctxt;
    rotate(shCtxt, i); 

    // apply the linearlized polynomial
    for (long k = 0; k < d; k++) {

      // compute the constant
      bool zConst = true;
      vector<RX> cvec;
      cvec.resize(nslots);
      for (long j = 0; j < nslots; j++) {
        cvec[j] = diag[j][k];
        if (!IsZero(cvec[j])) zConst = false;
      }

      if (zConst) continue;

      ZZX cpoly;
      encode(cpoly, cvec);
      // FIXME: record the encoded polynomial for future use

      Ctxt shCtxt1 = shCtxt;
      shCtxt1.frobeniusAutomorph(k);
      shCtxt1.multByConstant(cpoly);
      res += shCtxt1;
    }
  }
  ctxt = res;
}
#endif
/****************** End linear transformation code ******************/
/********************************************************************/


// Explicit instantiation

template class EncryptedArrayDerived<PA_GF2>;
template class EncryptedArrayDerived<PA_zz_p>;

template class PlaintextArrayDerived<PA_GF2>;
template class PlaintextArrayDerived<PA_zz_p>;
