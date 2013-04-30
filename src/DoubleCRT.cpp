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
/* DoubleCRT.cpp - This class holds an integer polynomial in double-CRT form
 *
 * Double-CRT form is a matrix of L rows and phi(m) columns. The i'th row
 * contains the FFT of the element wrt the ith prime, i.e. the evaluations of
 * the polynomial at the primitive mth roots of unity mod the ith prime. The
 * polynomial thus represented is defined modulo the product of all the primes
 * in use. The list of primes is defined by the data member modChain, which is
 * a vector of Cmodulus objects. 
 */


#if 0

#define DoubleCRT AltCRT
#include "AltCRT.cpp"

#else


#include <NTL/ZZX.h>
//NTL_CLIENT
#include "NumbTh.h"
#include "PAlgebra.h"
#include "CModulus.h"
#include "DoubleCRT.h"
#include "SingleCRT.h"
#include "timing.h"

// NTL implementation of mat_long

//NTL_matrix_impl(long,vec_long,vec_vec_long,mat_long)
//NTL_io_matrix_impl(long,vec_long,vec_vec_long,mat_long)
//NTL_eq_matrix_impl(long,vec_long,vec_vec_long,mat_long)


bool DoubleCRT::dryRun = false;

// representing an integer polynomial as DoubleCRT. If the number of moduli
// to use is not specified, the resulting object uses all the moduli in
// the context. If the coefficients of poly are larger than the product of
// the used moduli, they are effectively reduced modulo that product


// a "sanity check" function, verifies consistency of matrix with current
// moduli chain an error is raised if they are not consistent
void DoubleCRT::verify()
{
  assert(map.getIndexSet() <= (context.specialPrimes | context.ctxtPrimes));
  const IndexSet& s = map.getIndexSet();

  long phim = context.zMStar.getPhiM();

  // check that the content of i'th row is in [0,pi) for all i
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];

    if (row.length() != phim) 
      Error("DoubleCRT object has bad row length");

    long pi = context.ithPrime(i); // the i'th modulus
    for (long j=0; j<phim; j++)
      if (row[j]<0 || row[j]>= pi) 
	Error("DoubleCRT object has inconsistent data");
  }
}

// Arithmetic operations. Only the "destructive" versions are used,
// i.e., a += b is implemented but not a + b.

// Generic operation, Fnc is AddMod, SubMod, or MulMod (from NTL's ZZ module)
template<class Fun>
DoubleCRT& DoubleCRT::Op(const DoubleCRT &other, Fun fun,
			 bool matchIndexSets)
{
  if (dryRun) return *this;

  if (&context != &other.context)
    Error("DoubleCRT::Op: incompatible objects");

  // Match the index sets, if needed
  if (matchIndexSets && !(map.getIndexSet() >= other.map.getIndexSet()))
    addPrimes(other.map.getIndexSet() / map.getIndexSet()); // This is expensive

  // If you need to mod-up the other, do it on a temporary scratch copy
  DoubleCRT tmp(context, IndexSet()); 
  const IndexMap<vec_long>* other_map = &other.map;
  if (!(map.getIndexSet() <= other.map.getIndexSet())){ // Even more expensive
    tmp = other;
    tmp.addPrimes(map.getIndexSet() / other.map.getIndexSet());
    other_map = &tmp.map;
  }

  const IndexSet& s = map.getIndexSet();
  long phim = context.zMStar.getPhiM();

  // add/sub/mul the data, element by element, modulo the respective primes
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    vec_long& row = map[i];
    const vec_long& other_row = (*other_map)[i];
    
    for (long j = 0; j < phim; j++)
      row[j] = fun.apply(row[j], other_row[j], pi);
  }
  return *this;
}

template
DoubleCRT& DoubleCRT::Op<DoubleCRT::MulFun>(const DoubleCRT &other, MulFun fun,
			 bool matchIndexSets);

template
DoubleCRT& DoubleCRT::Op<DoubleCRT::AddFun>(const DoubleCRT &other, AddFun fun,
			 bool matchIndexSets);

template
DoubleCRT& DoubleCRT::Op<DoubleCRT::SubFun>(const DoubleCRT &other, SubFun fun,
			 bool matchIndexSets);

template<class Fun>
DoubleCRT& DoubleCRT::Op(const ZZ &num, Fun fun)
{
  if (dryRun) return *this;

  const IndexSet& s = map.getIndexSet();
  long phim = context.zMStar.getPhiM();
  
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    long n = rem(num, pi);  // n = num % pi
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++)
      row[j] = fun.apply(row[j], n, pi);
  }
  return *this;
}

template
DoubleCRT& DoubleCRT::Op<DoubleCRT::MulFun>(const ZZ &num, MulFun fun);

template
DoubleCRT& DoubleCRT::Op<DoubleCRT::AddFun>(const ZZ &num, AddFun fun);

template
DoubleCRT& DoubleCRT::Op<DoubleCRT::SubFun>(const ZZ &num, SubFun fun);

DoubleCRT& DoubleCRT::Negate(const DoubleCRT& other)
{
  if (dryRun) return *this;

  if (&context != &other.context) 
    Error("DoubleCRT Negate: incompatible contexts");

  if (map.getIndexSet() != other.map.getIndexSet()) {
    map = other.map; // copy the data
  }
  const IndexSet& s = map.getIndexSet();
  long phim = context.zMStar.getPhiM();
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    vec_long& row = map[i];
    const vec_long& other_row = other.map[i];
    for (long j = 0; j < phim; j++)
      row[j] = NegateMod(other_row[j], pi);
  }
  return *this;
}

template<class Fun>
DoubleCRT& DoubleCRT::Op(const ZZX &poly, Fun fun)
{
  if (dryRun) return *this;

  const IndexSet& s = map.getIndexSet();
  DoubleCRT other(poly, context, s); // other defined wrt same primes as *this

  return Op(other, fun);
}

template
DoubleCRT& DoubleCRT::Op<DoubleCRT::MulFun>(const ZZX &poly, MulFun fun);

template
DoubleCRT& DoubleCRT::Op<DoubleCRT::AddFun>(const ZZX &poly, AddFun fun);

template
DoubleCRT& DoubleCRT::Op<DoubleCRT::SubFun>(const ZZX &poly, SubFun fun);

// break *this into n digits,according to the primeSets in context.digits
void DoubleCRT::breakIntoDigits(vector<DoubleCRT>& digits, long n) const
{
  FHE_TIMER_START;
  IndexSet allPrimes = getIndexSet() | context.specialPrimes;
  assert(n <= (long)context.digits.size());

  digits.resize(n, DoubleCRT(context, IndexSet::emptySet()));
  if (dryRun) return;

  for (long i=0; i<(long)digits.size(); i++) {
    digits[i]=*this;
    IndexSet notInDigit = digits[i].getIndexSet()/context.digits[i];
    digits[i].removePrimes(notInDigit); // reduce modulo the digit primes
  }
  
  for (long i=0; i<(long)digits.size(); i++) {
    IndexSet notInDigit = allPrimes / digits[i].getIndexSet();
    digits[i].addPrimes(notInDigit); // add back all the primes

    // subtract this digits from all the others, then divide by pi
    ZZ pi = context.productOfPrimes(context.digits[i]);
    for (long j=i+1; j<(long)digits.size(); j++) {
      digits[j].Sub(digits[i], /*matchIndexSets=*/false);
      digits[j] /= pi;
    }
  }
#if 0
  dgts.resize(n, DoubleCRT(context, IndexSet::emptySet()));
  for (long i=0; i<n; i++) // copy only the primes for this digit
    dgts[i].partialCopy(*this, context.digits[i]);

  IndexSet allPrimes = getIndexSet() | context.specialPrimes;
  for (long i=0; i<n; i++) {
    IndexSet notInDigit = allPrimes / dgts[i].getIndexSet();
    dgts[i].addPrimes(notInDigit); // add back all the primes

    // subtract this digits from all the others, then divide by pi
    ZZ pi = context.productOfPrimes(context.digits[i]);
    for (long j=i+1; j<n; j++) {
      dgts[j].Sub(dgts[i], /*matchIndexSets=*/false);
      dgts[j] /= pi;
    }
  }
#endif
  FHE_TIMER_STOP;
}

// expand index set by s1.
// it is assumed that s1 is disjoint from the current index set.
void DoubleCRT::addPrimes(const IndexSet& s1)
{
  if (empty(s1)) return; // nothing to do
  assert( disjoint(s1,map.getIndexSet()) ); // s1 is disjoint from *this

  ZZX poly;
  toPoly(poly); // recover in coefficient representation

  map.insert(s1);  // add new rows to the map
  if (dryRun) return;

  // fill in new rows
  for (long i = s1.first(); i <= s1.last(); i = s1.next(i)) {
    context.ithModulus(i).FFT(map[i], poly); // reduce mod p_i and store FFT image
  }
}

// Expand index set by s1, and multiply by \prod{q \in s1}. s1 is assumed to
// be disjoint from the current index set. Returns the logarithm of product.
double DoubleCRT::addPrimesAndScale(const IndexSet& s1)
{
  if (empty(s1)) return 0.0; // nothing to do
  assert(empty(s1 & map.getIndexSet())); // s1 is disjoint from *this

  // compute factor to scale existing rows
  ZZ factor = to_ZZ(1);
  double logFactor = 0.0;
  for (long i = s1.first(); i <= s1.last(); i = s1.next(i)) {
    long qi = context.ithPrime(i);
    factor *= qi;
    logFactor += log((double)qi);
  }

  // scale existing rows
  long phim = context.zMStar.getPhiM();
  const IndexSet& iSet = map.getIndexSet();
  for (long i = iSet.first(); i <= iSet.last(); i = iSet.next(i)) {
    long qi = context.ithPrime(i);
    long f = rem(factor, qi);     // f = factor % qi
    vec_long& row = map[i];
    // scale row by a factor of f modulo qi
    mulmod_precon_t bninv = PrepMulModPrecon(f, qi, 1.0/(double)qi);
    for (long j=0; j<phim; j++) 
      row[j] = MulModPrecon(row[j], f, qi, bninv);
  }

  // insert new rows and fill them with zeros
  map.insert(s1);  // add new rows to the map
  for (long i = s1.first(); i <= s1.last(); i = s1.next(i)) {
    vec_long& row = map[i];
    for (long j=0; j<phim; j++) row[j] = 0;
  }

  return logFactor;
}

DoubleCRT::DoubleCRT(const ZZX& poly, const FHEcontext &_context, const IndexSet& s)
: context(_context), map(new DoubleCRTHelper(_context))
{
  assert(s.last() < context.numPrimes());

  map.insert(s);
  if (dryRun) return;

  // convert the integer polynomial to FFT representation modulo the primes
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    const Cmodulus &pi = context.ithModulus(i);
    pi.FFT(map[i], poly); // reduce mod pi and store FFT image
  }
}

DoubleCRT::DoubleCRT(const ZZX& poly, const FHEcontext &_context)
: context(_context), map(new DoubleCRTHelper(_context))
{
  IndexSet s = IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);
  if (dryRun) return;

  // convert the integer polynomial to FFT representation modulo the primes
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    const Cmodulus &pi = context.ithModulus(i);
    pi.FFT(map[i], poly); // reduce mod pi and store FFT image
  }
}

DoubleCRT::DoubleCRT(const ZZX& poly)
: context(*activeContext), map(new DoubleCRTHelper(*activeContext))
{
  IndexSet s = IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);
  if (dryRun) return;

  // convert the integer polynomial to FFT representation modulo the primes
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    const Cmodulus &pi = context.ithModulus(i);
    pi.FFT(map[i], poly); // reduce mod pi and store FFT image
  }
}

DoubleCRT::DoubleCRT(const FHEcontext &_context, const IndexSet& s)
: context(_context), map(new DoubleCRTHelper(_context))
{
  assert(s.last() < context.numPrimes());

  map.insert(s);
  if (dryRun) return;

  long phim = context.zMStar.getPhiM();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++) row[j] = 0;
  }
}

DoubleCRT::DoubleCRT(const FHEcontext &_context)
: context(_context), map(new DoubleCRTHelper(_context))
{
  IndexSet s = IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);
  if (dryRun) return;

  long phim = context.zMStar.getPhiM();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++) row[j] = 0;
  }
}

DoubleCRT& DoubleCRT::operator=(const DoubleCRT& other)
// optimized for the case of matching index sets
{
   if (this == &other) return *this;

   if (&context != &other.context) 
      Error("DoubleCRT assignment: incompatible contexts");

   if (map.getIndexSet() != other.map.getIndexSet()) {
      map = other.map; // copy the data
   }
   else {
      const IndexSet& s = map.getIndexSet();
      long phim = context.zMStar.getPhiM();
      for (long i = s.first(); i <= s.last(); i = s.next(i)) {
         vec_long& row = map[i];
         const vec_long& other_row = other.map[i];
         for (long j = 0; j < phim; j++)
            row[j] = other_row[j];
      }
   }
   return *this;
}

#if 0
// Copy only the primes in s \intersect other.getIndexSet()
void DoubleCRT::partialCopy(const DoubleCRT& other, const IndexSet& _s)
{
   if (&context != &other.context) 
      Error("DoubleCRT::partialCopy: incompatible contexts");

   // set the primes of *this to s \intersect other.getIndexSet()
   IndexSet s = _s;
   s.retain(other.getIndexSet());
   map.remove(getIndexSet() / s);
   map.insert(s / getIndexSet());

   long phim = context.zMStar.getPhiM();
   for (long i = s.first(); i <= s.last(); i = s.next(i)) {
     vec_long& row = map[i];
     const vec_long& other_row = other.map[i];
     for (long j = 0; j < phim; j++)
       row[j] = other_row[j];
   }
}
#endif

DoubleCRT& DoubleCRT::operator=(const ZZX&poly)
{
  if (dryRun) return *this;

  const IndexSet& s = map.getIndexSet();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) 
    context.ithModulus(i).FFT(map[i],poly); // reduce mod pi and store FFT image

  return *this;
}

DoubleCRT& DoubleCRT::operator=(const ZZ& num)
{
  const IndexSet& s = map.getIndexSet();
  if (dryRun) return *this;

  long phim = context.zMStar.getPhiM();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];
    long pi = context.ithPrime(i);
    long n = rem(num, pi);

    for (long j = 0; j < phim; j++) row[j] = n;
  }

  return *this;
}

void DoubleCRT::toPoly(ZZX& poly, const IndexSet& s,
		       bool positive) const
{
FHE_TIMER_START
  if (dryRun) return;

  IndexSet s1 = map.getIndexSet() & s;

  if (empty(s1)) {
    clear(poly);
    return;
  }

  clear(poly);
  ZZ prod;
  prod = 1;

  zz_pBak bak; bak.save();

  for (long i = s1.first(); i <= s1.last(); i = s1.next(i)) {
    context.ithModulus(i).restoreModulus();
    zz_pX& tmp = context.ithModulus(i).getScratch();
    context.ithModulus(i).iFFT(tmp, map[i]); 
    CRT(poly, prod, tmp);  // NTL :-)
  }

  if (positive) {
    long d = deg(poly);
    for (long j = 0; j <= d; j++) 
      if (poly.rep[j] < 0)
        poly.rep[j] += prod;   

    // no need to normalize poly here
  }
FHE_TIMER_STOP
}

#if 0
{
FHE_TIMER_START
  if (dryRun) return;

  IndexSet s1 = map.getIndexSet() & s;

  if (empty(s1)) {
    clear(poly);
    return;
  }

  ZZ p = to_ZZ(context.ithPrime(s1.first()));  // the first modulus

  // Get poly modulo the first prime in coefficent form
  long i = s1.first();
  const Cmodulus& mod = context.ithModulus(i);
  mod.iFFT(poly, map[i]);

  vec_ZZ& vp = poly.rep;

  // ensure that vp is of size phi(m) with entries in [-p/2,p/2]
  long phim = context.zMStar.getPhiM();
  long vpLength = vp.length();
  if (vpLength < phim) { // just in case of leading zeros in poly
    vp.SetLength(phim);
    for (long j = vpLength; j < phim; j++) vp[j]=0;
  }
  ZZ p_over_2 = p/2;
  for (long j = 0; j < phim; j++) if (vp[j] > p_over_2) vp[j] -= p;

  // do incremental integer CRT for other levels

  ZZX current;
  for (i = s1.next(i); i <= s1.last(); i = s1.next(i)) {
    long q = context.ithPrime(i);          // the next modulus
    context.ithModulus(i).iFFT(current, map[i]); // Poly mod q in coeff form

    // CRT the coefficient vectors of poly and current
    intVecCRT(vp, p, current.rep, q);    // defined in the module NumbTh
    p *= q;     // update the modulus
  }

  // The above yeilds polynomial with coefficients in [-p/2,p/2]
  // If we need positive, just add p to all the negative coefficients
  if (positive) 
    for (long j=0; j<poly.rep.length(); j++) {
      if (poly.rep[j] < 0) poly.rep[j] += p;
    }

  poly.normalize(); // need to call this after we work on the coeffs
FHE_TIMER_STOP
}
#endif

void DoubleCRT::toPoly(ZZX& p, bool positive) const
{
  const IndexSet& s = map.getIndexSet();
  toPoly(p, s, positive);
}

// Division by constant
DoubleCRT& DoubleCRT::operator/=(const ZZ &num)
{
  if (dryRun) return *this;

  const IndexSet& s = map.getIndexSet();
  long phim = context.zMStar.getPhiM();
  
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    long n = InvMod(rem(num, pi),pi);  // n = num^{-1} mod pi
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++)
      row[j] = MulMod(row[j], n, pi);
  }
  return *this;
}

// Small-exponent polynomial exponentiation
void DoubleCRT::Exp(long e)
{
  if (dryRun) return;

  const IndexSet& s = map.getIndexSet();
  long phim = context.zMStar.getPhiM();
  
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++)
      row[j] = PowerMod(row[j], e, pi);
  }
}

// Apply the automorphism F(X) --> F(X^k)  (with gcd(k,m)=1)
void DoubleCRT::automorph(long k)
{
  if (dryRun) return;

  const PAlgebra& zMStar = context.zMStar;
  if (!zMStar.inZmStar(k))
    Error("DoubleCRT::automorph: k not in Zm*");

  long m = zMStar.getM();
  vector<long> tmp(m);  // temporary array of size m

  const IndexSet& s = map.getIndexSet();

  // go over the rows, permute them one at a time
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];

    for (long j=1; j<m; j++) { // 1st pass: copy to temporary array
      long idx = zMStar.indexInZmstar(j); // returns -1 if j \notin (Z/mZ)*
      if (idx>=0) tmp[j] = row[idx];
    }

    for (long j=1; j<m; j++) { // 2nd pass: copy back from temporary array
      long idx = zMStar.indexInZmstar(j); // returns -1 if j \notin (Z/mZ)*
      if (idx>=0) row[idx] = tmp[MulMod(j,k,m)];
                                           // new[j] = old[j*k mod m]
    }    
  }
}

// fills each row i with random integers mod pi
void DoubleCRT::randomize(const ZZ* seed) 
{
  if (dryRun) return;

  if (seed != NULL) SetSeed(*seed);

  const IndexSet& s = map.getIndexSet();
  long phim = context.zMStar.getPhiM();
  
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++)
      row[j] = RandomBnd(pi);   // RandomBnd is defined in NTL's module ZZ
  }
}

DoubleCRT& DoubleCRT::operator=(const SingleCRT& scrt)
{
  if (&context != &scrt.getContext())
    Error("DoubleCRT=SingleCRT -- incompatible contexts");

  map.clear();  // empty the map
  const IndexSet& s = scrt.getMap().getIndexSet();
  map.insert(s);
  
  if (dryRun) return *this;

  for (long i = s.first(); i <= s.last(); i = s.next(i)) 
    context.ithModulus(i).FFT(map[i],scrt.getMap()[i]); // compute FFT image
  return *this;
}

void DoubleCRT::toSingleCRT(SingleCRT& scrt, const IndexSet& s) const 
{
  if (&context != &scrt.getContext())
    Error("DoubleCRT::toSingleCRT -- incompatible contexts");

  IndexSet s1 = s & map.getIndexSet();
  scrt.map.clear();
  scrt.map.insert(s1);

  for (long i = s1.first(); i <= s1.last(); i = s1.next(i)) {
    context.ithModulus(i).restoreModulus();
    zz_pX tmp;
    context.ithModulus(i).iFFT(tmp, map[i]); // inverse FFT
    conv(scrt.map[i], tmp);
  }
}

void DoubleCRT::toSingleCRT(SingleCRT& scrt) const 
{
  const IndexSet& s = map.getIndexSet();
  toSingleCRT(scrt, s);
}


void DoubleCRT::scaleDownToSet(const IndexSet& s, long ptxtSpace)
{
  assert(ptxtSpace >= 2);
  assert(!(getIndexSet()<=s)); // cannot mod-down to the empty set

  IndexSet diff = getIndexSet() / s;
  if (empty(diff)) return;     // nothing to do

  if (dryRun) {
    removePrimes(diff);// remove the primes from consideration
    return;
  }

  ZZX delta;
  ZZ diffProd = context.productOfPrimes(diff); // mod-down by this factor
  toPoly(delta, diff); // convert to coeff-representation modulo diffProd

  long delta_len = delta.rep.length();
  if (ptxtSpace == 2) { // simpler handling for plaintext space mod 2
    for (long i = 0; i < delta_len; i++) {
      if (IsOdd(delta.rep[i])) { // add or subtract diffProd to make it even
	if (sign(delta.rep[i]) < 0) delta.rep[i] += diffProd;
	else                        delta.rep[i] -= diffProd;
      }
    }
  }
  // The general case of plaintext space modulo some p > 2, we need to
  // subtract from each coefficient delta[i] the ineteger
  //               diffProd * (delta[i] * diffProd^{-1} mod ptxtSpace).
  // This does not change delta modulo diffProd, but makes it divisible
  // by ptxtSpace.
  else {
    long p_over_2 = ptxtSpace/2;
    long prodInv = InvMod(rem(diffProd,ptxtSpace), ptxtSpace);
    for (long i = 0; i < delta_len; i++) {
      long delta_i_modP = rem(delta.rep[i],ptxtSpace);
      if (delta_i_modP != 0) { // if not already 0 mod ptxtSpace
	delta_i_modP = MulMod(delta_i_modP, prodInv, ptxtSpace);
	if (delta_i_modP > p_over_2) delta_i_modP -= ptxtSpace;
	delta.rep[i] -= diffProd * delta_i_modP;
      }
    }
  }
  delta.normalize(); // need to normalize after working directly on then coeffs

  removePrimes(diff);// remove the primes from consideration
  *this -= delta;    // convert delta to DoubleCRT, then subtract
  *this /= diffProd; // *this is divisible by diffProd, so this operation actually scales it down
}

ostream& operator<< (ostream &str, const DoubleCRT &d)
{
  const IndexSet& set = d.map.getIndexSet();

  // check that the content of i'th row is in [0,pi) for all i
  str << "[" << set << endl;
  for (long i = set.first(); i <= set.last(); i = set.next(i))
    str << " " << d.map[i] << "\n";
  str << "]";
  return str;
}

istream& operator>> (istream &str, DoubleCRT &d)
{
  //  cerr << "DoubleCRT[";
  // Advance str beyond first '[' 
  seekPastChar(str, '[');  // this function is defined in NumbTh.cpp

  IndexSet set;
  const FHEcontext& context = d.context;
  long phim = context.zMStar.getPhiM();

  str >> set; // read in the indexSet
  assert(set <= (context.specialPrimes | context.ctxtPrimes));
  d.map.clear();
  d.map.insert(set); // fix the index set for the data

  for (long i = set.first(); i <= set.last(); i = set.next(i)) {
    str >> d.map[i]; // read the actual data

    // verify that the data is valid
    assert (d.map[i].length() == phim);
    for (long j=0; j<phim; j++)
      assert(d.map[i][j]>=0 && d.map[i][j]<context.ithPrime(i));
  }

  // Advance str beyond closing ']'
  seekPastChar(str, ']');
  //  cerr << "]";
  return str;
}
#endif

