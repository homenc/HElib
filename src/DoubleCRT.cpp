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
/* DoubleCRT.cpp - This class holds an integer polynomial in double-CRT form
 *
 * Double-CRT form is a matrix of L rows and phi(m) columns. The i'th row
 * contains the FFT of the element wrt the ith prime, i.e. the evaluations of
 * the polynomial at the primitive mth roots of unity mod the ith prime. The
 * polynomial thus represented is defined modulo the product of all the primes
 * in use. The list of primes is defined by the data member modChain, which is
 * a vector of Cmodulus objects. 
 */
#include <NTL/ZZVec.h>
#include <NTL/BasicThreadPool.h>

#include "DoubleCRT.h"
#include "timing.h"


// A threaded implementation of DoubleCRT operations

static
long MakeIndexVector(const IndexSet& s, Vec<long>& v)
{
  long sz = s.card();
  v.SetLength(sz);
  for (long i = s.first(), j = 0; i <= s.last(); i = s.next(i), j++)
    v[j] = i;
  return sz;
}

// representing an integer polynomial as DoubleCRT. If the number of moduli
// to use is not specified, the resulting object uses all the moduli in
// the context. If the coefficients of poly are larger than the product of
// the used moduli, they are effectively reduced modulo that product


void DoubleCRT::FFT(const ZZX& poly, const IndexSet& s)
{
  FHE_TIMER_START;

  if (empty(s)) return;

  static thread_local Vec<long> tls_ivec;
  Vec<long>& ivec = tls_ivec;

  long icard = MakeIndexVector(s, ivec);
  NTL_EXEC_RANGE(icard, first, last)
      for (long j = first; j < last; j++) {
        long i = ivec[j];
        context.ithModulus(i).FFT(map[i], poly); 
      }
  NTL_EXEC_RANGE_END
}

// FIXME: "code bloat": this just replicates the above with ZZX -> zzX
void DoubleCRT::FFT(const zzX& poly, const IndexSet& s)
{
  FHE_TIMER_START;

  if (empty(s)) return;

  static thread_local Vec<long> tls_ivec;
  Vec<long>& ivec = tls_ivec;

  long icard = MakeIndexVector(s, ivec);
  NTL_EXEC_RANGE(icard, first, last)
      for (long j = first; j < last; j++) {
        long i = ivec[j];
        context.ithModulus(i).FFT(map[i], poly); 
      }
  NTL_EXEC_RANGE_END
}


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
  if (isDryRun()) return *this;

  if (&context != &other.context)
    Error("DoubleCRT::Op: incompatible objects");

  // Match the index sets, if needed
  if (matchIndexSets && !(map.getIndexSet() >= other.map.getIndexSet())) {
    FHE_NTIMER_START(addPrimes_1); 
    addPrimes(other.map.getIndexSet() / map.getIndexSet()); // This is expensive
  }

  // If you need to mod-up the other, do it on a temporary scratch copy
  DoubleCRT tmp(context, IndexSet()); 
  const IndexMap<vec_long>* other_map = &other.map;
  if (!(map.getIndexSet() <= other.map.getIndexSet())){ // Even more expensive
    FHE_NTIMER_START(addPrimes_2); 
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



// Victor says: I added this routine so I could look
// examine its performance more carefully

DoubleCRT& DoubleCRT::do_mul(const DoubleCRT &other, 
			     bool matchIndexSets)
{
  FHE_TIMER_START;

  if (isDryRun()) return *this;

  if (&context != &other.context)
    Error("DoubleCRT::Op: incompatible objects");

  // Match the index sets, if needed
  if (matchIndexSets && !(map.getIndexSet() >= other.map.getIndexSet())) {
    FHE_NTIMER_START(addPrimes_3);
    addPrimes(other.map.getIndexSet() / map.getIndexSet()); // This is expensive
  }

  // If you need to mod-up the other, do it on a temporary scratch copy
  DoubleCRT tmp(context, IndexSet()); 
  const IndexMap<vec_long>* other_map = &other.map;
  if (!(map.getIndexSet() <= other.map.getIndexSet())){ // Even more expensive
    FHE_NTIMER_START(addPrimes_4);
    tmp = other;
    tmp.addPrimes(map.getIndexSet() / other.map.getIndexSet());
    other_map = &tmp.map;
  }

  const IndexSet& s = map.getIndexSet();
  long phim = context.zMStar.getPhiM();

  // add/sub/mul the data, element by element, modulo the respective primes
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    mulmod_t pi_inv = context.ithModulus(i).getQInv(); 
    vec_long& row = map[i];
    const vec_long& other_row = (*other_map)[i];


    for (long j = 0; j < phim; j++)
      row[j] = MulMod(row[j], other_row[j], pi, pi_inv);

  }
  return *this;
}

#if 0
template
DoubleCRT& DoubleCRT::Op<DoubleCRT::MulFun>(const DoubleCRT &other, MulFun fun,
			 bool matchIndexSets);
#endif

template
DoubleCRT& DoubleCRT::Op<DoubleCRT::AddFun>(const DoubleCRT &other, AddFun fun,
			 bool matchIndexSets);

template
DoubleCRT& DoubleCRT::Op<DoubleCRT::SubFun>(const DoubleCRT &other, SubFun fun,
			 bool matchIndexSets);

template<class Fun>
DoubleCRT& DoubleCRT::Op(const ZZ &num, Fun fun)
{
  if (isDryRun()) return *this;

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
  if (isDryRun()) return *this;

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
  if (isDryRun()) return *this;

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
  if (isDryRun()) return;

  for (long i=0; i<(long)digits.size(); i++) {
    digits[i]=*this;
    IndexSet notInDigit = digits[i].getIndexSet()/context.digits[i];
    digits[i].removePrimes(notInDigit); // reduce modulo the digit primes
  }
  
  for (long i=0; i<(long)digits.size(); i++) {
    FHE_NTIMER_START(addPrimes_5);
    IndexSet notInDigit = allPrimes / digits[i].getIndexSet();
    digits[i].addPrimes(notInDigit); // add back all the primes

    ZZ pi = context.productOfPrimes(context.digits[i]);
    for (long j=i+1; j<(long)digits.size(); j++) {
      digits[j].Sub(digits[i], /*matchIndexSets=*/false);
      digits[j] /= pi;
    }
  }
  FHE_TIMER_STOP;
}

// expand index set by s1.
// it is assumed that s1 is disjoint from the current index set.
void DoubleCRT::addPrimes(const IndexSet& s1)
{
  FHE_TIMER_START;

  if (empty(s1)) return; // nothing to do
  assert( disjoint(s1,map.getIndexSet()) ); // s1 is disjoint from *this

  ZZX poly;
  toPoly(poly); // recover in coefficient representation

  map.insert(s1);  // add new rows to the map
  if (isDryRun()) return;

  // fill in new rows
  FFT(poly, s1);
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
    mulmod_precon_t bninv = PrepMulModPrecon(f, qi);
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


// *****************************************************

DoubleCRT::DoubleCRT(const ZZX& poly, const FHEcontext &_context, const IndexSet& s)
: context(_context), map(new DoubleCRTHelper(_context))
{
  FHE_TIMER_START;
  assert(s.last() < context.numPrimes());

  map.insert(s);
  if (isDryRun()) return;

  // convert the integer polynomial to FFT representation modulo the primes
  FFT(poly, s);
}

DoubleCRT::DoubleCRT(const ZZX& poly, const FHEcontext &_context)
: context(_context), map(new DoubleCRTHelper(_context))
{
  FHE_TIMER_START;
  IndexSet s = IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);
  if (isDryRun()) return;

  // convert the integer polynomial to FFT representation modulo the primes
  FFT(poly, s);
}

DoubleCRT::DoubleCRT(const ZZX& poly)
: context(*activeContext), map(new DoubleCRTHelper(*activeContext))
{
  FHE_TIMER_START;
  IndexSet s = IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);
  if (isDryRun()) return;

  // convert the integer polynomial to FFT representation modulo the primes
  FFT(poly, s);
}


// *****************************************************
// FIXME: "code bloat": this just replicates the above with ZZX -> zzX

DoubleCRT::DoubleCRT(const zzX& poly, const FHEcontext &_context, const IndexSet& s)
: context(_context), map(new DoubleCRTHelper(_context))
{
  FHE_TIMER_START;
  assert(s.last() < context.numPrimes());

  map.insert(s);
  if (isDryRun()) return;

  // convert the integer polynomial to FFT representation modulo the primes
  FFT(poly, s);
}

DoubleCRT::DoubleCRT(const zzX& poly, const FHEcontext &_context)
: context(_context), map(new DoubleCRTHelper(_context))
{
  FHE_TIMER_START;
  IndexSet s = IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);
  if (isDryRun()) return;

  // convert the integer polynomial to FFT representation modulo the primes
  FFT(poly, s);
}

DoubleCRT::DoubleCRT(const zzX& poly)
: context(*activeContext), map(new DoubleCRTHelper(*activeContext))
{
  FHE_TIMER_START;
  IndexSet s = IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);
  if (isDryRun()) return;

  // convert the integer polynomial to FFT representation modulo the primes
  FFT(poly, s);
}

DoubleCRT::DoubleCRT(const FHEcontext &_context, const IndexSet& s)
: context(_context), map(new DoubleCRTHelper(_context))
{
  assert(s.last() < context.numPrimes());

  map.insert(s);
  if (isDryRun()) return;

  long phim = context.zMStar.getPhiM();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++) row[j] = 0;
  }
}

// *****************************************************

DoubleCRT::DoubleCRT(const FHEcontext &_context)
: context(_context), map(new DoubleCRTHelper(_context))
{
  IndexSet s = IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);
  if (isDryRun()) return;

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


DoubleCRT& DoubleCRT::operator=(const ZZX&poly)
{
  if (isDryRun()) return *this;

  const IndexSet& s = map.getIndexSet();

  FFT(poly, s);

  return *this;
}

DoubleCRT& DoubleCRT::operator=(const ZZ& num)
{
  const IndexSet& s = map.getIndexSet();
  if (isDryRun()) return *this;

  long phim = context.zMStar.getPhiM();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];
    long pi = context.ithPrime(i);
    long n = rem(num, pi);

    for (long j = 0; j < phim; j++) row[j] = n;
  }

  return *this;
}

// DIRT: this method affect the NTL zz_p::modulus
long DoubleCRT::getOneRow(zz_pX& row, long idx) const
{
  if (!map.getIndexSet().contains(idx)) // idx not in the primeset
    return 0;

  // convert from evaluation to standard coefficient representation
  context.ithModulus(idx).restoreModulus(); // recover NTL modulus for prime
  context.ithModulus(idx).iFFT(row, map[idx]);
  return context.ithPrime(idx);
}

// Get the row corresponding to the i'th moduli, in Vec<long> format.
// For conveience, returns the modulus that was used for this row.
// If idx is not in the current primesSet then do nothing and return 0;
long DoubleCRT::getOneRow(Vec<long>& row, long idx, bool positive) const
{
  zz_pBak bak; bak.save();   // backup NTL's current modulus

  zz_pX& tmp = Cmodulus::getScratch_zz_pX();
  long q = getOneRow(tmp,idx);
  if (q==0) return 0;        // no such index

  conv(row, tmp.rep);        // copy the row to Vec<long> format

  // By default, integers are in [0,q).
  // If we need the symmetric interval then make it so.
  if (!positive) {
    long phim = context.zMStar.getPhiM();
    for (long j = 0; j < phim; j++) if (row[j] > q/2) row[j] -= q;
  }
  return q;
}


void DoubleCRT::toPoly(ZZX& poly, const IndexSet& s,
		       bool positive) const
{
  FHE_TIMER_START;
  if (isDryRun()) return;

  IndexSet s1 = map.getIndexSet() & s;

  if (empty(s1)) {
    clear(poly);
    return;
  }


  static thread_local Vec<long> tls_ivec;
  static thread_local Vec<long> tls_pvec;
  static thread_local Vec< Vec<long> > tls_remtab;
  static thread_local Vec<zz_pX> tls_tmpvec;

  Vec<long>& ivec = tls_ivec;
  Vec<long>& pvec = tls_pvec;
  Vec< Vec<long> >& remtab = tls_remtab;
  Vec<zz_pX>& tmpvec = tls_tmpvec;

  long phim = context.zMStar.getPhiM();
  long icard = MakeIndexVector(s1, ivec);

  PartitionInfo pinfo(icard);
  long cnt = pinfo.NumIntervals();

  remtab.SetLength(phim);
  for (long h = 0; h < phim; h++) remtab[h].SetLength(icard);

  tmpvec.SetLength(cnt);
  for (long i = 0; i < cnt; i++) tmpvec[i].SetMaxLength(phim);

  { FHE_NTIMER_START(toPoly_FFT);
  
  NTL_EXEC_INDEX(cnt, index)
      long first, last;
      pinfo.interval(first, last, index);

      zz_pX& tmp = tmpvec[index];
  
      for (long j = first; j < last; j++) {
        long i = ivec[j];
        context.ithModulus(i).iFFT(tmp, map[i]); 
  
        long d = deg(tmp);
        for (long h = 0; h <= d; h++) remtab[h][j] = rep(tmp.rep[h]);
        for (long h = d+1; h < phim; h++) remtab[h][j] = 0;
      }
  NTL_EXEC_INDEX_END

  }

  {FHE_NTIMER_START(toPoly_CRT);

  PartitionInfo pinfo1(phim);
  long cnt1 = pinfo1.NumIntervals();

  static thread_local ZZ tls_prod;
  static thread_local ZZ tls_prod_half;
  static thread_local Vec<long> tls_qvec;
  static thread_local Vec<double> tls_qrecipvec;
  static thread_local Vec<long> tls_tvec;
  static thread_local Vec<mulmod_precon_t> tls_tqinvvec;

  static thread_local ZZVec tls_prod1vec;
  static thread_local ZZVec tls_resvec;

  ZZ& prod = tls_prod;
  ZZ& prod_half = tls_prod_half;
  ZZVec& prod1vec = tls_prod1vec;
  Vec<long>& qvec = tls_qvec;
  Vec<double>& qrecipvec = tls_qrecipvec;
  Vec<long>& tvec = tls_tvec;
  Vec<mulmod_precon_t>& tqinvvec = tls_tqinvvec;
  ZZVec& resvec = tls_resvec;

  qvec.SetLength(icard);
  qrecipvec.SetLength(icard);
  tvec.SetLength(icard);
  tqinvvec.SetLength(icard);

  prod = 1;
  for (long j = 0; j < icard; j++) {
    long i = ivec[j];
    long q = context.ithModulus(i).getQ();
    qvec[j] = q;
    qrecipvec[j] = 1/double(q);
    mul(prod, prod, q);
  }

  long sz = prod.size();

  if (prod1vec.length() != icard || prod1vec.BaseSize() != sz+1) {
    prod1vec.kill();
    prod1vec.SetSize(icard, sz+1);
  }

  for (long j = 0; j < icard; j++) {
    long q = qvec[j];
    div(prod1vec[j], prod, q);
    long t = rem(prod1vec[j], q);
    t = InvMod(t, q);
    tvec[j] = t;
    tqinvvec[j] = PrepMulModPrecon(t, q);
  }

  if (resvec.length() != phim || resvec.BaseSize() != sz+1) {
    resvec.kill();
    resvec.SetSize(phim, sz+1);
  }

  if (!positive) {
    // prod_half = (prod+1)/2
    add(prod_half, prod, 1);
    div(prod_half, prod_half, 2);
  }
  
  NTL_EXEC_INDEX(cnt1, index)
  NTL_IMPORT(icard)
      long first, last;
      pinfo1.interval(first, last, index);

      long *qvecp = qvec.elts();
      double *qrecipvecp = qrecipvec.elts();
      long *tvecp = tvec.elts();
      mulmod_precon_t *tqinvvecp = tqinvvec.elts();
      ZZ *prod1vecp = prod1vec.elts();

      ZZ tmp;
      tmp.SetSize(sz+4);
  
      for (long h = first; h < last; h++) {
        clear(tmp);
        double quotient = 0;
        long *remvec = remtab[h].elts();
        
        for (long j = 0; j < icard; j++) {
          long q = qvecp[j];
          long t = tvecp[j];
          mulmod_precon_t tqinv = tqinvvecp[j];
          long r = remvec[j];
          double qrecip = qrecipvecp[j];
          r = MulModPrecon(r, t, q, tqinv);
          MulAddTo(tmp, prod1vecp[j], r);
          quotient += r*qrecip;
        }
        
        MulSubFrom(tmp, prod, long(quotient));
        while (tmp < 0) add(tmp, tmp, prod);
        while (tmp >= prod) sub(tmp, tmp, prod);
        if (!positive && tmp >= prod_half) 
          tmp -= prod;
        resvec[h] = tmp;
      }
  NTL_EXEC_INDEX_END

  poly.SetLength(phim);
  for (long j = 0; j < phim; j++) poly[j] = resvec[j];
  poly.normalize();

  // NOTE: assigning to poly[j] within the parallel loop
  // leads to horrible performance, as there apparently is
  // a lot of contention within malloc.
  }
}




void DoubleCRT::toPoly(ZZX& p, bool positive) const
{
  const IndexSet& s = map.getIndexSet();
  toPoly(p, s, positive);
}

// Division by constant
DoubleCRT& DoubleCRT::operator/=(const ZZ &num)
{
  if (isDryRun()) return *this;

  const IndexSet& s = map.getIndexSet();
  long phim = context.zMStar.getPhiM();
  
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    long n = InvMod(rem(num, pi),pi);  // n = num^{-1} mod pi
    vec_long& row = map[i];
    mulmod_precon_t precon = PrepMulModPrecon(n, pi);
    for (long j = 0; j < phim; j++)
      row[j] = MulModPrecon(row[j], n, pi, precon);
  }
  return *this;
}

// Small-exponent polynomial exponentiation
void DoubleCRT::Exp(long e)
{
  if (isDryRun()) return;

  const IndexSet& s = map.getIndexSet();
  long phim = context.zMStar.getPhiM();
  
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    vec_long& row = map[i];
    for (long j = 0; j < phim; j++)
      row[j] = PowerMod(row[j], e, pi);
  }
}

#if 1

// Apply the automorphism F(X) --> F(X^k)  (with gcd(k,m)=1)
void DoubleCRT::automorph(long k)
{
  if (isDryRun()) return;

  const PAlgebra& zMStar = context.zMStar;
  if (!zMStar.inZmStar(k))
    Error("DoubleCRT::automorph: k not in Zm*");

  long m = zMStar.getM();
  long phim = zMStar.getPhiM();
  vector<long> tmp(m);  // temporary array of size m
  mulmod_precon_t precon = PrepMulModPrecon(k, m);

  const IndexSet& s = map.getIndexSet();

  // go over the rows, permute them one at a time
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];

#if 0
    for (long j=1; j<m; j++) { // 1st pass: copy to temporary array
      long idx = zMStar.indexInZmstar_unchecked(j); // returns -1 if j \notin (Z/mZ)*
      if (idx>=0) tmp[j] = row[idx];
    }

    for (long j=1; j<m; j++) { // 2nd pass: copy back from temporary array
      long idx = zMStar.indexInZmstar_unchecked(j); // returns -1 if j \notin (Z/mZ)*
      if (idx>=0) row[idx] = tmp[MulModPrecon(j,k,m,precon)];
                                           // new[j] = old[j*k mod m]
    }    
#else
    // slightly faster...

    for (long j = 0; j < phim; j++) {
       tmp[zMStar.repInZmstar_unchecked(j)] = row[j];
    }

    for (long j = 0; j < phim; j++) {
       row[j] = tmp[ MulModPrecon( zMStar.repInZmstar_unchecked(j), k, m, precon) ];
    }

#endif
  }
}

#else

// VJS: I tried this as an alternative...it is slower :-(

// Apply the automorphism F(X) --> F(X^k)  (with gcd(k,m)=1)
void DoubleCRT::automorph(long k)
{
  if (isDryRun()) return;

  const PAlgebra& zMStar = context.zMStar;
  if (!zMStar.inZmStar(k))
    Error("DoubleCRT::automorph: k not in Zm*");

  long m = zMStar.getM();
  long phim = zMStar.getPhiM();
  vector<long> tmp(phim);  // temporary array of size m

  k = InvMod(k, m);
  mulmod_precon_t precon = PrepMulModPrecon(k, m);

  const IndexSet& s = map.getIndexSet();

  // go over the rows, permute them one at a time
  // new[j*k mod m] = old[j]
  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    vec_long& row = map[i];

    for (long j = 0; j < phim; j++) tmp[j] = row[j];

    for (long j = 0; j < phim; j++) {
       long rep = zMStar.repInZmstar_unchecked(j);
       rep = MulModPrecon(rep, k, m, precon);
       long idx = zMStar.indexInZmstar_unchecked(rep);
       row[idx] = tmp[j];
    }
  }

}

#endif

// fills each row i with random integers mod pi
void DoubleCRT::randomize(const ZZ* seed) 
{
  FHE_TIMER_START;

  if (isDryRun()) return;

  if (seed != NULL) SetSeed(*seed);

  const IndexSet& s = map.getIndexSet();
  long phim = context.zMStar.getPhiM();

  RandomStream& stream = GetCurrentRandomStream();
  const long bufsz = 2048;

  Vec<unsigned char> buf_storage;
  buf_storage.SetLength(bufsz);

  unsigned char *buf = buf_storage.elts();

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i);
    long k = NumBits(pi-1);
    long nb = (k+7)/8;
    unsigned long mask = (1UL << k) - 1UL;

    vec_long& row = map[i];
    long j = 0;
    
    for (;;) {
      { FHE_NTIMER_START(randomize_stream);
      stream.get(buf, bufsz);
      }

      for (long pos = 0; pos <= bufsz-nb; pos += nb) {
#if 0
        unsigned long utmp = 0;
        for (long cnt = nb-1;  cnt >= 0; cnt--)
          utmp = (utmp << 8) | buf[pos+cnt]; 
#elif 0

        // "Duff's device" to avoid loops
        // It's a bit faster...but not much
       
        unsigned long utmp = buf[pos+nb-1];
        switch (nb) {
        case 8: utmp = (utmp << 8) | buf[pos+6];
        case 7: utmp = (utmp << 8) | buf[pos+5];
        case 6: utmp = (utmp << 8) | buf[pos+4];
        case 5: utmp = (utmp << 8) | buf[pos+3];
        case 4: utmp = (utmp << 8) | buf[pos+2];
        case 3: utmp = (utmp << 8) | buf[pos+1];
        case 2: utmp = (utmp << 8) | buf[pos+0];
        }

#else
        unsigned long utmp = buf[pos+nb-1];

        {

        // This is gcc non-standard. Works also on clang and icc.
        
        static void *dispatch_table[] =
           { &&L0, &&L1, &&L2, &&L3, &&L4, &&L5, &&L6, &&L7, &&L8 };

        goto *dispatch_table[nb];
   

        L8: utmp = (utmp << 8) | buf[pos+6];
        L7: utmp = (utmp << 8) | buf[pos+5];
        L6: utmp = (utmp << 8) | buf[pos+4];
        L5: utmp = (utmp << 8) | buf[pos+3];
        L4: utmp = (utmp << 8) | buf[pos+2];
        L3: utmp = (utmp << 8) | buf[pos+1];
        L2: utmp = (utmp << 8) | buf[pos+0];
        L1: ;
        L0: ;
        }

#endif

        utmp = (utmp & mask);
        
        long tmp = utmp;

        row[j] = tmp;
        j += (tmp < pi);
        if (j >= phim) break;
      }
      if (j >= phim) break;
    }
  }
}


void DoubleCRT::scaleDownToSet(const IndexSet& s, long ptxtSpace)
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
  *this -= delta;    // convert delta to DoubleCRT, then subtract
  *this /= diffProd; // *this is divisible by diffProd, so this operation actually scales it down
}
#if 0 // By riku.
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
std::ostream& operator<< (std::ostream &str, const DoubleCRT &d)
{
  const IndexSet& set = d.map.getIndexSet();
  str << '[';
  str << set;
  str << ']';
  NTL::ZZX poly;
  FHE_NTIMER_START(TO_POLY_OUTPUT);
  d.toPoly(poly, true);
  FHE_NTIMER_STOP(TO_POLY_OUTPUT);
  const FHEcontext &context = d.context;
  double bits = context.logOfProduct(set);
  bits /= std::log(2.);
  long bytes = long(std::ceil(bits) + 7) >> 3;
  long phim = context.zMStar.getPhiM();
  std::vector<uint8_t> buff(bytes);
  for (long i = 0; i < phim; i++) {
      const NTL::ZZ &e = NTL::coeff(poly, i);
      BytesFromZZ(buff.data(), e, bytes);
      str.write(reinterpret_cast<char *>(buff.data()), bytes);
  }
  return str;
}

std::istream& operator>> (std::istream &str, DoubleCRT &d)
{
  IndexSet set;
  const FHEcontext& context = d.context;
  long phim = context.zMStar.getPhiM();
  seekPastChar(str, '[');
  str >> set; // read in the indexSet
  seekPastChar(str, ']');
  assert(set <= (context.specialPrimes | context.ctxtPrimes));
  d.map.clear();
  d.map.insert(set); // fix the index set for the data

  NTL::ZZX poly;
  poly.SetLength(phim + 1);
  double bits = context.logOfProduct(set);
  bits /= std::log(2.0);
  long bytes = long(std::ceil(bits) + 7) >> 3;
  std::vector<uint8_t> buff(bytes);
  NTL::ZZ e;
  for (long i = 0; i < phim; i++) {
      str.read(reinterpret_cast<char *>(buff.data()), bytes);
      NTL::ZZFromBytes(e, buff.data(), bytes);
      NTL::SetCoeff(poly, i, e);
  }
  FHE_NTIMER_START(FROM_POLY_OUTPUT);
  d = poly;
  FHE_NTIMER_STOP(FROM_POLY_OUTPUT);
  return str;
}
