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

#include "AltCRT.h"
#include "DoubleCRT.h"
#include "timing.h"


/*
 * Current incarnation of lazy strategy: 
 *  - results are generally left reduced mod X^m-1
 *  - right before two AltCRT objects are multiplied,
 *    they are reduced mod Phi_m(X), unless the "lazy"
 *    flag in the context is set, in which case, we are
 *    *really* lazy, and don't reduce the inputs to
 *    multiplication
 *  - objects are also reduced in toPoly
 *  - the reduce() method can be explicitly called
 *    to force a reduction mod Phi_m(X), and the reduce()
 *    method in Ctxt calls this method.
 *    Right now, this is invoked both before and after 
 *    we relinearize a ciphertext
 */

void AltCRT::verify()
{
  assert(map.getIndexSet() <= (context.specialPrimes | context.ctxtPrimes));
}

// reduce x mod X^m-1
static void reduce1(zz_pX& x, long m) 
{
  long d = deg(x);
  if (d >= m) {
    long j = 0;
    for (long i = m; i <= d; i++) {
      x[j] += x[i];
      j++;
      if (j >= m) j = 0;
    }

    x.SetLength(m);
    x.normalize();
  }
}

// reduce x mod f, assuming already reduced mod X^m-1
static void reduce2(zz_pX& x, const zz_pXModulus1& f)
{
   assert(deg(x) < f.m);
   if (deg(x) < f.n) return;
   rem(x, x, f);
}

static void reduce12(zz_pX& x, const zz_pXModulus1& f)
{
  reduce1(x, f.m);
  reduce2(x, f);
}



void MulMod1(zz_pX& x, const zz_pX& a, const zz_pX& b, const zz_pXModulus1& f)
{
  assert(deg(a) < f.m && deg(b) < f.m);

  zz_pX t;
  mul(t, a, b);

  reduce1(t, f.m);

  x = t;
}




// Arithmetic operations. Only the "destructive" versions are used,
// i.e., a += b is implemented but not a + b.

template<class Fun>
AltCRT& AltCRT::Op(const AltCRT &other, Fun fun,
			 bool matchIndexSets)
{
  if (isDryRun()) return *this;

  if (&context != &other.context)
    Error("AltCRT::Op: incompatible objects");

  if (fun.reduce()) {
    this->reduce();
    other.reduce();
  }


  // Match the index sets, if needed
  if (matchIndexSets && !(map.getIndexSet() >= other.map.getIndexSet()))
    addPrimes(other.map.getIndexSet() / map.getIndexSet()); // This is expensive

  // If you need to mod-up the other, do it on a temporary scratch copy
  AltCRT tmp(context, IndexSet()); 
  const IndexMap<zz_pX>* other_map = &other.map;
  if (!(map.getIndexSet() <= other.map.getIndexSet())){ // Even more expensive
    tmp = other;
    tmp.addPrimes(map.getIndexSet() / other.map.getIndexSet());
    other_map = &tmp.map;
  }

  const IndexSet& s = map.getIndexSet();

  zz_pBak bak; bak.save();


  // add/sub/mul the data, element by element, modulo the respective primes
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    context.ithModulus(i).restoreModulus();
    const zz_pXModulus1& phimx = context.ithModulus(i).getPhimX();
    fun.apply(map[i], (*other_map)[i], phimx);
  }

  return *this;
}

template
AltCRT& AltCRT::Op<AltCRT::MulFun>(const AltCRT &other, MulFun fun,
			 bool matchIndexSets);

template
AltCRT& AltCRT::Op<AltCRT::AddFun>(const AltCRT &other, AddFun fun,
			 bool matchIndexSets);

template
AltCRT& AltCRT::Op<AltCRT::SubFun>(const AltCRT &other, SubFun fun,
			 bool matchIndexSets);





template<class Fun>
AltCRT& AltCRT::Op(const ZZ &num, Fun fun)
{
  if (isDryRun()) return *this;

  const IndexSet& s = map.getIndexSet();

  zz_pBak bak; bak.save();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    context.ithModulus(i).restoreModulus();
    zz_p n;
    conv(n, num);
    fun.apply(map[i], n);
  }

  return *this;
}

template
AltCRT& AltCRT::Op<AltCRT::MulFun>(const ZZ &num, MulFun fun);

template
AltCRT& AltCRT::Op<AltCRT::AddFun>(const ZZ &num, AddFun fun);

template
AltCRT& AltCRT::Op<AltCRT::SubFun>(const ZZ &num, SubFun fun);



AltCRT& AltCRT::Negate(const AltCRT& other)
{
  if (isDryRun()) return *this;

  if (&context != &other.context) 
    Error("AltCRT Negate: incompatible contexts");

  if (map.getIndexSet() != other.map.getIndexSet()) {
    map = other.map; // DIRT: copies zz_p objects out of context
  }
  const IndexSet& s = map.getIndexSet();

  zz_pBak bak; bak.save();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    context.ithModulus(i).restoreModulus();
    NTL::negate(map[i], other.map[i]);
  }

  return *this;
}



// The following is identical to definition in DoubleCRT

template<class Fun>
AltCRT& AltCRT::Op(const ZZX &poly, Fun fun)
{
  if (isDryRun()) return *this;

  const IndexSet& s = map.getIndexSet();
  AltCRT other(poly, context, s); // other defined wrt same primes as *this

  return Op(other, fun);
}


template
AltCRT& AltCRT::Op<AltCRT::MulFun>(const ZZX &poly, MulFun fun);

template
AltCRT& AltCRT::Op<AltCRT::AddFun>(const ZZX &poly, AddFun fun);

template
AltCRT& AltCRT::Op<AltCRT::SubFun>(const ZZX &poly, SubFun fun);



void AltCRT::reduce() const
// logically but not really const
{
  FHE_TIMER_START;

  if (context.lazy) return;

  const IndexSet& s = map.getIndexSet();

  zz_pBak bak; bak.save();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    context.ithModulus(i).restoreModulus();
    const zz_pXModulus1& phimx = context.ithModulus(i).getPhimX();
    reduce12(const_cast<zz_pX&>(map[i]), phimx);
  }
}



// The following is identical to definition in DoubleCRT
// ...not quite...we first reduce *this

// break *this into n digits,according to the primeSets in context.digits
void AltCRT::breakIntoDigits(vector<AltCRT>& digits, long n) const
{
  FHE_TIMER_START;
  IndexSet allPrimes = getIndexSet() | context.specialPrimes;
  assert(n <= (long)context.digits.size());

  digits.resize(n, AltCRT(context, IndexSet::emptySet()));
  if (isDryRun()) return;

  for (long i=0; i<(long)digits.size(); i++) {
    digits[i]=*this;;
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
  dgts.resize(n, AltCRT(context, IndexSet::emptySet()));
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
void AltCRT::addPrimes(const IndexSet& s1)
{
  if (empty(s1)) return; // nothing to do
  assert( disjoint(s1,map.getIndexSet()) ); // s1 is disjoint from *this

  ZZX poly;
  toPoly(poly); // recover in coefficient representation

  map.insert(s1);  // add new rows to the map
  if (isDryRun()) return;


  zz_pBak bak; bak.save();

  // fill in new rows
  for (long i = s1.first(); i <= s1.last(); i = s1.next(i)) {
    context.ithModulus(i).restoreModulus();
    conv(map[i], poly);
  }
}






// Expand index set by s1, and multiply by \prod{q \in s1}. s1 is assumed to
// be disjoint from the current index set. Returns the logarithm of product.
double AltCRT::addPrimesAndScale(const IndexSet& s1)
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


  zz_pBak bak; bak.save();

  // scale existing rows
  const IndexSet& iSet = map.getIndexSet();
  for (long i = iSet.first(); i <= iSet.last(); i = iSet.next(i)) {
    context.ithModulus(i).restoreModulus();
    zz_p f;
    conv(f, factor);
    map[i] *= f;
  }

  // insert new rows and fill them with zeros
  map.insert(s1);  // add new rows to the map
  for (long i = s1.first(); i <= s1.last(); i = s1.next(i)) {
    context.ithModulus(i).restoreModulus();
    clear(map[i]);
  }

  return logFactor;
}





AltCRT::AltCRT(const ZZX& poly, const FHEcontext &_context, const IndexSet& s)
: context(_context), map(new AltCRTHelper(_context))
{
  assert(s.last() < context.numPrimes());

  map.insert(s);
  if (isDryRun()) return;

  zz_pBak bak; bak.save();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    context.ithModulus(i).restoreModulus();
    conv(map[i], poly);
  }
}





AltCRT::AltCRT(const ZZX& poly, const FHEcontext &_context)
: context(_context), map(new AltCRTHelper(_context))
{
  IndexSet s = IndexSet(0, context.numPrimes()-1);

  map.insert(s);
  if (isDryRun()) return;


  zz_pBak bak; bak.save();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    context.ithModulus(i).restoreModulus();
    conv(map[i], poly);
  }
}





AltCRT::AltCRT(const ZZX& poly)
: context(*activeContext), map(new AltCRTHelper(*activeContext))
{
  IndexSet s = IndexSet(0, context.numPrimes()-1);

  map.insert(s);
  if (isDryRun()) return;

  zz_pBak bak; bak.save();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    context.ithModulus(i).restoreModulus();
    conv(map[i], poly);
  }
}





AltCRT::AltCRT(const FHEcontext &_context, const IndexSet& s)
: context(_context), map(new AltCRTHelper(_context))
{
  assert(s.last() < context.numPrimes());

  map.insert(s);
  if (isDryRun()) return;

  zz_pBak bak; bak.save();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    context.ithModulus(i).restoreModulus();
    clear(map[i]);
  }
}





AltCRT::AltCRT(const FHEcontext &_context)
: context(_context), map(new AltCRTHelper(_context))
{
  IndexSet s = IndexSet(0, context.numPrimes()-1);

  map.insert(s);
  if (isDryRun()) return;

  zz_pBak bak; bak.save();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    context.ithModulus(i).restoreModulus();
    clear(map[i]);
  }
}


AltCRT::AltCRT(const AltCRT& other) 
: context(other.context), map(new AltCRTHelper(other.context))

{
   map = other.map;  // DIRT: copies zz_p objects out of context
}


AltCRT& AltCRT::operator=(const AltCRT& other)
// optimized for the case of matching index sets
{
   if (this == &other) return *this;

   if (&context != &other.context) 
      Error("AltCRT assignment: incompatible contexts");


   if (map.getIndexSet() != other.map.getIndexSet()) {
      map = other.map; // DIRT: copies zz_p objects out of context
   }
   else {
      const IndexSet& s = map.getIndexSet();
      zz_pBak bak; bak.save(); 

      for (long i = s.first(); i <= s.last(); i = s.next(i)) {
         context.ithModulus(i).restoreModulus();
         map[i] = other.map[i];
      }
   }
   return *this;
}




AltCRT& AltCRT::operator=(const ZZX&poly)
{
  if (isDryRun()) return *this;

  const IndexSet& s = map.getIndexSet();
  zz_pBak bak; bak.save();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) { 
    context.ithModulus(i).restoreModulus();
    conv(map[i], poly);
  }

  return *this;
}





AltCRT& AltCRT::operator=(const ZZ& num)
{
  if (isDryRun()) return *this;

  const IndexSet& s = map.getIndexSet();
  zz_pBak bak; bak.save();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) { 
    context.ithModulus(i).restoreModulus();
    conv(map[i], num);
  }

  return *this;
}

// DIRT: this method affect the NTL zz_p::modulus
long AltCRT::getOneRow(zz_pX& row, long idx) const
{
  if (!map.getIndexSet().contains(idx)) // idx not in the primeset
    return 0;

  context.ithModulus(idx).restoreModulus(); // recover NTL modulus for prime
  long q = context.ithPrime(idx);

  // LAZY: first reduce map[idx] mod phimx
  const zz_pXModulus1& phimx = context.ithModulus(idx).getPhimX();
  reduce12(const_cast<zz_pX&>(map[idx]), phimx);
  row = map[idx];

  return q;
}

// Get the row corresponding to the i'th moduli, in Vec<long> format.
// For conveience, returns the modulus that was used for this row.
// If idx is not in the current primesSet then do nothing and return 0;
long AltCRT::getOneRow(Vec<long>& row, long idx, bool positive) const
{
  if (!map.getIndexSet().contains(idx)) // idx not in the primeset
    return 0;

  zz_pBak bak; bak.save();  // backup NTL's current modulus
  context.ithModulus(idx).restoreModulus(); // recover NTL modulus for prime
  long q = context.ithPrime(idx);

  // LAZY: first reduce map[idx] mod phimx
  const zz_pXModulus1& phimx = context.ithModulus(idx).getPhimX();
  reduce12(const_cast<zz_pX&>(map[idx]), phimx);

  // copy the row to Vec<long> format
  conv(row, map[idx].rep); // does this work??

  // By default, integers are in [0,q).
  // If we need the symmetric interval then make it so.
  if (positive) {
    long phim = context.zMStar.getPhiM();
    for (long j = 0; j < phim; j++) if (row[j] > q/2) row[j] -= q;
  }
  return q;
} // NTL's current modulus restored upon exit



// DIRT: I am not sure if this function behaves the same
// as in DoubleCRT if the prime 2 is allowed: the endpoints
// of the interval [-P/2,P/2] may be handled differently.
// But all primes should be odd (in fact, it is required
// that p = 1 (mod 2m)), so the point is academic (VJS)

void AltCRT::toPoly(ZZX& poly, const IndexSet& s,
		       bool positive) const
{
  if (isDryRun()) return;

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

    // LAZY: first reduce map[i] mod phimx
    const zz_pXModulus1& phimx = context.ithModulus(i).getPhimX();
    reduce12(const_cast<zz_pX&>(map[i]), phimx);
    
    CRT(poly, prod, map[i]);  // NTL :-)
  }

  if (positive) {
    long d = deg(poly);
    for (long j = 0; j <= d; j++) 
      if (poly.rep[j] < 0)
        poly.rep[j] += prod;   

    // no need to normalize poly here
  }
}





// The following is identical to definition in DoubleCRT

void AltCRT::toPoly(ZZX& p, bool positive) const
{
  const IndexSet& s = map.getIndexSet();
  toPoly(p, s, positive);
}




// Division by constant
AltCRT& AltCRT::operator/=(const ZZ &num)
{
  if (isDryRun()) return *this;

  const IndexSet& s = map.getIndexSet();
  zz_pBak bak; bak.save();
  
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    context.ithModulus(i).restoreModulus();
    zz_p n;
    conv(n, num);
    map[i] /= n;
  }
  return *this;
}




// Small-exponent polynomial exponentiation
void AltCRT::Exp(long e)
{
  if (isDryRun()) return;

  const IndexSet& s = map.getIndexSet();
  zz_pBak bak; bak.save();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    context.ithModulus(i).restoreModulus();
    const zz_pXModulus1& phimx = context.ithModulus(i).getPhimX();

    // LAZY: first reduce map[i] mod phimx
    reduce12(map[i], phimx);

    PowerMod(map[i], map[i], e, phimx.upcast());
  }
}




// Apply the automorphism F(X) --> F(X^k)  (with gcd(k,m)=1)
void AltCRT::automorph(long k)
{
  if (isDryRun()) return;

  const PAlgebra& zMStar = context.zMStar;
  if (!zMStar.inZmStar(k))
    Error("AltCRT::automorph: k not in Zm*");

  long m = zMStar.getM();

  const IndexSet& s = map.getIndexSet();
  zz_pBak bak; bak.save();

  // go over the rows, permute them one at a time

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    context.ithModulus(i).restoreModulus();
    zz_pX& tmp = Cmodulus::getScratch_zz_pX();
    zz_pX& row = map[i];
    long d = deg(row);

    tmp.rep.SetLength(m);
    for (long j = 0; j < m; j++) tmp.rep[j] = 0;

    mulmod_precon_t precon = PrepMulModPrecon(k, m);
    for (long j = 0; j <= d; j++) 
      tmp.rep[MulModPrecon(j, k, m, precon)] = row.rep[j];

    tmp.normalize();

    row = tmp;
  }
}








// FIXME: there is a potential incompatibilty here
// with DoubleCRT -- starting from the same seed,
// we will get different polynomials.  This may lead
// to trouble if we start mixing DoubleCRT's and AltCRT's,
// especially with respect to the way we currently do 
// key switching (VJS)

// fills each row i with random numbers
void AltCRT::randomize(const ZZ* seed) 
{
  if (isDryRun()) return;

  if (seed != NULL) SetSeed(*seed);

  const IndexSet& s = map.getIndexSet();
  long phim = context.zMStar.getPhiM();
  zz_pBak bak; bak.save();
  
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    context.ithModulus(i).restoreModulus();
    random(map[i], phim);
  }
}


// The following is identical to definition in DoubleCRT


void AltCRT::scaleDownToSet(const IndexSet& s, long ptxtSpace)
{
  assert(ptxtSpace >= 2);

  IndexSet diff = getIndexSet() / s;
  assert(diff!=s);          // cannot mod-down to the empty set
  if (empty(diff)) return;  // nothing to do

  if (isDryRun()) {
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
    mulmod_precon_t precon = PrepMulModPrecon(prodInv, ptxtSpace); // optimization
    for (long i = 0; i < delta_len; i++) {
      long delta_i_modP = rem(delta.rep[i],ptxtSpace);
      if (delta_i_modP != 0) { // if not already 0 mod ptxtSpace
	delta_i_modP = MulModPrecon(delta_i_modP, prodInv, ptxtSpace, precon);
	if (delta_i_modP > p_over_2) delta_i_modP -= ptxtSpace;
	delta.rep[i] -= diffProd * delta_i_modP;
      }
    }
  }
  delta.normalize(); // need to normalize after working directly on then coeffs

  removePrimes(diff);// remove the primes from consideration
  *this -= delta;    // convert delta to AltCRT, then subtract
  *this /= diffProd; // *this is divisible by diffProd, so this operation actually scales it down
}








ostream& operator<< (ostream &str, const AltCRT &d)
{
  // output the set of "active" primes
  const IndexSet& set = d.map.getIndexSet();
  str << "[" << set << endl;

  // For each "active" prime, output the zz_pX object
  zz_pBak bak; bak.save();
  for (long i = set.first(); i <= set.last(); i = set.next(i)) {
    d.context.ithModulus(i).restoreModulus();

    // LAZY: first reduce d.map[i] mod phimx
    const zz_pXModulus1& phimx = d.context.ithModulus(i).getPhimX();
    reduce12(const_cast<zz_pX&>(d.map[i]), phimx);

    str << " " << d.map[i] << "\n";
  }
  str << "]";
  return str;
}

istream& operator>> (istream &str, AltCRT &d)
{
  // Advance str beyond first '[' 
  seekPastChar(str, '[');  // this function is defined in NumbTh.cpp

  IndexSet set;
  const FHEcontext& context = d.context;

  str >> set; // read in the indexSet, describing the "active" primes
  assert(set <= (context.specialPrimes | context.ctxtPrimes));
  d.map.clear();
  d.map.insert(set); // fix the index set for the data

  zz_pBak bak; bak.save();
  for (long i = set.first(); i <= set.last(); i = set.next(i)) {
    d.context.ithModulus(i).restoreModulus();
    str >> d.map[i]; // read the actual data
  }
  // Advance str beyond closing ']'
  seekPastChar(str, ']');
  return str;
}



