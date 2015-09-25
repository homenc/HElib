
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
void EncryptedArrayDerived<type>::encode(ZZX& ptxt, const NewPlaintextArray& array) const
{
  RBak bak; bak.save(); tab.restoreContext();
  encode(ptxt, array.getData<type>());
}




template<class type>
void EncryptedArrayDerived<type>::decode(NewPlaintextArray& array, const ZZX& ptxt) const
{
  RBak bak; bak.save(); tab.restoreContext();
  decode(array.getData<type>(), ptxt);
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


// NewPlaintextArray


template<class type>
class rotate_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, long k)
  {
    PA_BOILER

    vector<RX> tmp(n); 

    for (long i = 0; i < n; i++)
      tmp[((i+k)%n + n)%n] = data[i];

    data = tmp;
  }
};

void rotate(const EncryptedArray& ea, NewPlaintextArray& pa, long k)
{
  ea.dispatch<rotate_pa_impl>(Fwd(pa), k); 
}


//=======================================================================================


template<class type>
class shift_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, long k)
  {
    PA_BOILER

    for (long i = 0; i < n; i++)
      if (i + k >= n || i + k < 0)
        clear(data[i]);

    rotate_pa_impl<type>::apply(ea, pa, k); 
  }
};

void shift(const EncryptedArray& ea, NewPlaintextArray& pa, long k)
{
  ea.dispatch<shift_pa_impl>(Fwd(pa), k); 
}




//=======================================================================================



template<class type>
class encode_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, 
    const vector<long>& array)
  {
    PA_BOILER

    assert(lsize(array) == n);
    convert(data, array);
  }

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, 
    const vector<ZZX>& array)
  {
    PA_BOILER

    assert(lsize(array) == n);
    convert(data, array);
    for (long i = 0; i < n; i++) assert(deg(data[i]) < d);
  }

};



void encode(const EncryptedArray& ea, NewPlaintextArray& pa, const vector<long>& array)
{
  ea.dispatch<encode_pa_impl>(Fwd(pa), array); 
}

void encode(const EncryptedArray& ea, NewPlaintextArray& pa, const vector<ZZX>& array)
{
  ea.dispatch<encode_pa_impl>(Fwd(pa), array); 
}

void encode(const EncryptedArray& ea, NewPlaintextArray& pa, long val)
{
   long n = ea.size();
   vector<long> array;
   array.resize(n);
   for (long i = 0; i < n; i++) array[i] = val;
   encode(ea, pa, array);
}

void encode(const EncryptedArray& ea, NewPlaintextArray& pa, const ZZX& val)
{
   long n = ea.size();
   vector<ZZX> array;
   array.resize(n);
   for (long i = 0; i < n; i++) array[i] = val;
   encode(ea, pa, array);
}



//=======================================================================================




template<class type>
class random_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa)
  {
    PA_BOILER

    for (long i = 0; i < n; i++)
      random(data[i], d);
  }
}; 


void random(const EncryptedArray& ea, NewPlaintextArray& pa)
{
  ea.dispatch<random_pa_impl>(Fwd(pa)); 
}



//=======================================================================================




template<class type>
class decode_pa_impl {
public:
  PA_INJECT(type)

  template<class T>
  static void apply(const EncryptedArrayDerived<type>& ea, 
    vector<T>& array, const NewPlaintextArray& pa)
  {
    CPA_BOILER

    convert(array, data);
  }

}; 


void decode(const EncryptedArray& ea, vector<long>& array, const NewPlaintextArray& pa)
{
  ea.dispatch<decode_pa_impl>(Fwd(array), pa); 
}


void decode(const EncryptedArray& ea, vector<ZZX>& array, const NewPlaintextArray& pa)
{
  ea.dispatch<decode_pa_impl>(Fwd(array), pa); 
}




//=======================================================================================



template<class type>
class equals_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, bool& res, 
    const NewPlaintextArray& pa, const  NewPlaintextArray& other)
  {
    CPA_BOILER

    const vector<RX>& odata = other.getData<type>(); 
    res = (data == odata);
  }


  static void apply(const EncryptedArrayDerived<type>& ea, bool& res, 
    const NewPlaintextArray& pa, const vector<long>& other)
  {
    CPA_BOILER

    vector<RX> odata;
    convert(odata, other);
    res = (data == odata);
  }


  static void apply(const EncryptedArrayDerived<type>& ea, bool& res,
    const NewPlaintextArray& pa, const vector<ZZX>& other)
  {
    CPA_BOILER

    vector<RX> odata;
    convert(odata, other);
    res = (data == odata);
  }

};



bool equals(const EncryptedArray& ea, const NewPlaintextArray& pa, const NewPlaintextArray& other)
{
  bool res;
  ea.dispatch<equals_pa_impl>(Fwd(res), pa, other); 
  return res;
}

bool equals(const EncryptedArray& ea, const NewPlaintextArray& pa, const vector<long>& other)
{
  bool res;
  ea.dispatch<equals_pa_impl>(Fwd(res), pa, other); 
  return res;
}


bool equals(const EncryptedArray& ea, const NewPlaintextArray& pa, const vector<ZZX>& other)
{
  bool res;
  ea.dispatch<equals_pa_impl>(Fwd(res), pa, other); 
  return res;
}





//=======================================================================================




template<class type>
class add_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, 
    const NewPlaintextArray& other)
  {
    PA_BOILER

    const vector<RX>& odata = other.getData<type>(); 

    for (long i = 0; i < n; i++)
      data[i] += odata[i];
  }
}; 


void add(const EncryptedArray& ea, NewPlaintextArray& pa, const NewPlaintextArray& other)
{
  ea.dispatch<add_pa_impl>(Fwd(pa), other); 
}





//=======================================================================================




template<class type>
class sub_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, 
    const NewPlaintextArray& other)
  {
    PA_BOILER

    const vector<RX>& odata = other.getData<type>(); 

    for (long i = 0; i < n; i++)
      data[i] -= odata[i];
  }
}; 


void sub(const EncryptedArray& ea, NewPlaintextArray& pa, const NewPlaintextArray& other)
{
  ea.dispatch<sub_pa_impl>(Fwd(pa), other); 
}



//=======================================================================================




template<class type>
class mul_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, 
    const NewPlaintextArray& other)
  {
    PA_BOILER

    const vector<RX>& odata = other.getData<type>(); 

    for (long i = 0; i < n; i++)
      data[i] = (data[i] * odata[i]) % G;
  }
}; 


void mul(const EncryptedArray& ea, NewPlaintextArray& pa, const NewPlaintextArray& other)
{
  ea.dispatch<mul_pa_impl>(Fwd(pa), other); 
}





//=======================================================================================




template<class type>
class negate_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa) 
  {
    PA_BOILER

    for (long i = 0; i < n; i++)
      NTL::negate(data[i], data[i]);
  }
}; 


void negate(const EncryptedArray& ea, NewPlaintextArray& pa)
{
  ea.dispatch<negate_pa_impl>(Fwd(pa)); 
}






//=======================================================================================



template<class type>
class frobeniusAutomorph_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa, long j)
  {
    PA_BOILER

    long p = ea.getTab().getZMStar().getP();

    j = mcMod(j, d);
    RX H = PowerMod(RX(1, 1), power_ZZ(p, j), G);

    for (long i = 0; i < n; i++)
      data[i] = CompMod(data[i], H, G);
  }

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa,
    const Vec<long>& vec) 
  {
    PA_BOILER

    assert(vec.length() == n);

    long p = ea.getTab().getZMStar().getP();

    for (long i = 0; i < n; i++) {
      long j = mcMod(vec[i], d);
      RX H = PowerMod(RX(1, 1), power_ZZ(p, j), G);
      data[i] = CompMod(data[i], H, G);
    }
  }
};




void frobeniusAutomorph(const EncryptedArray& ea, NewPlaintextArray& pa, long j)
{
  ea.dispatch<frobeniusAutomorph_pa_impl>(Fwd(pa), j); 
}


void frobeniusAutomorph(const EncryptedArray& ea, NewPlaintextArray& pa, const Vec<long>& vec)
{
  ea.dispatch<frobeniusAutomorph_pa_impl>(Fwd(pa), vec); 
}




//=======================================================================================



template<class type>
class applyPerm_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa,
    const Vec<long>& pi) 
  {
    PA_BOILER

    assert(pi.length() == n);

    vector<RX> tmp;
    tmp.resize(n);
    for (long i = 0; i < n; i++)
      tmp[i] = data[pi[i]];

    data = tmp;
  }
};




void applyPerm(const EncryptedArray& ea, NewPlaintextArray& pa, const Vec<long>& pi)
{
  ea.dispatch<applyPerm_pa_impl>(Fwd(pa), pi); 
}



//=======================================================================================




template<class type>
class print_pa_impl {
public:
  PA_INJECT(type)

  static void apply(const EncryptedArrayDerived<type>& ea, 
    ostream& s, const NewPlaintextArray& pa)
  {
    CPA_BOILER


    if (n == 0) 
      s << "[]";
    else {
      if (IsZero(data[0])) s << "[[0]";
      else                 s << "[" << data[0];
      for (long i = 1; i < lsize(data); i++)
        if (IsZero(data[i])) s << " [0]";
	else                 s << " " << data[i];
      s << "]";
    }
  }

}; 


void print(const EncryptedArray& ea, ostream& s, const NewPlaintextArray& pa)
{
  ea.dispatch<print_pa_impl>(Fwd(s), pa); 
}



















// Explicit instantiation

template class EncryptedArrayDerived<PA_GF2>;
template class EncryptedArrayDerived<PA_zz_p>;


template class NewPlaintextArrayDerived<PA_GF2>;
template class NewPlaintextArrayDerived<PA_zz_p>;

