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
/* SingleCRT.cpp - This class hold integer polynomials modulo many small primes
 *
 * SingleCRT form is a vector of (ZZX) polynomials, where the i'th polynomial
 * contains the coefficients wrt the ith prime in the chain of moduli. The
 * polynomial thus represented is defined modulo the product of all the 
 * primes in use.
 *
 * This is mostly a helper class for DoubleCRT, and all the rules that apply
 * to DoubleCRT wrt modulus chain and levels apply also to SingleCRT objects.
 * Although SingleCRT and DoubleCRT objects can interact in principle,
 * translation back and forth are expensive since they involve FFT/iFFT.
 * Hence support for interaction between them is limited to explicit
 * conversions.
 */
#include <NTL/ZZ.h>

NTL_CLIENT

#include "NumbTh.h"
#include "SingleCRT.h"
#include "DoubleCRT.h"

// a "sanity check" function, verifies consistency of polys with current
// moduli chain an error is raised if they are not consistent
void SingleCRT::verify()
{
  const IndexSet& s = map.getIndexSet();
  if (s.last() >= context.numPrimes())
    Error("SingleCRT object has too many rows");

  // check that the content of i'th row is in [0,pi) for all i
  for (long i = s.first(); i < s.last(); i = s.next(i)) {
    long pi = context.ithPrime(i); // the i'th modulus
    vec_ZZ& vp = map[i].rep;
    for (long j=0; j<vp.length(); j++)
      if (vp[j]<0 || vp[j]>= pi) 
	Error("SingleCRT object has inconsistent data");
  }
}

// Generic operators, Fnc is either AddMod or SubMod
// This should have been a template, but gcc refuses to cooperate
SingleCRT& SingleCRT::Op(const SingleCRT &other,
			 void (*Fnc)(ZZ&, const ZZ&, const ZZ&, const ZZ&),
			 bool matchIndexSets)
{
 if (&context != &other.context)
    Error("SingleCRT::Op: incomopatible objects");

  // Match the index sets, if needed
  if (matchIndexSets && !(map.getIndexSet() >= other.map.getIndexSet()))
    addPrimes(other.map.getIndexSet() / map.getIndexSet()); // This is expensive

  // If you need to mod-up the other, do it on a temporary scratch copy
  SingleCRT tmp(context);
  const IndexMap<ZZX>* other_map = &other.map;
  if (!(map.getIndexSet() <= other.map.getIndexSet())){ // Even more expensive
    tmp = other;
    tmp.addPrimes(map.getIndexSet() / other.map.getIndexSet());
    other_map = &tmp.map;
  }

  const IndexSet& s = map.getIndexSet();

  // add/sub polynomial, modulo the respective primes

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    ZZ pi = to_ZZ(context.ithPrime(i));
    vec_ZZ& vp1 = map[i].rep;
    const vec_ZZ& vp2 = (*other_map)[i].rep;

    long len1 = vp1.length();
    long len2 = vp2.length();
    long maxlen = max(len1, len2);
    vp1.SetLength(maxlen);
    for (long j=len1; j < maxlen; j++) clear(vp1[j]);
    for (long j=0; j<len2; j++) 
      Fnc(vp1[j], vp1[j], vp2[j], pi);
    map[i].normalize();
  }
  return *this;
}

// Implementation of scrt += poly, scrt -= poly, or scrt *= poly. This
// implementation is safe for "in place" operation, e.g., s += s.map[i]
SingleCRT& SingleCRT::Op(const ZZX &poly,
			 void (*Fnc)(ZZ&, const ZZ&, const ZZ&, const ZZ&))
{
  const IndexSet& s = map.getIndexSet();

  ZZX poly1, poly2;

  poly1 = poly;

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    ZZ pi = to_ZZ(context.ithPrime(i));
    poly2 = poly1;
    PolyRed(poly2,pi,/*abs=*/true); // abs=true means reduce to [0,pi-1)

    vec_ZZ& vp1 = map[i].rep;
    vec_ZZ& vp2 = poly2.rep;

    long len1 = vp1.length();
    long len2 = vp2.length();
    long maxlen = max(len1, len2);
    vp1.SetLength(maxlen);
    for (long j=len1; j < maxlen; j++) clear(vp1[j]);
    for (long j=0; j<len2; j++) 
      Fnc(vp1[j], vp1[j], vp2[j], pi);

    map[i].normalize();
  }
  return *this;
}

// Here Fnc is either Add(ZZX,ZZX,ZZ), Sub(ZZX,ZZX,ZZ), or Mul(ZZX,ZZX,ZZ)
// FIXME: this is not alias friendly
SingleCRT& SingleCRT::Op(const ZZ &num, void (*Fnc)(ZZX&, const ZZX&, const ZZ&))
{
  const IndexSet& s = map.getIndexSet();
  ZZ pi;
  ZZ n;
  ZZX poly1;

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    conv(pi, context.ithPrime(i));
    rem(n, num, pi);  // n = num % pi
    poly1 = map[i];
    Fnc(poly1,poly1,n);
    PolyRed(poly1,pi,/*abs=*/true); // abs=true means reduce to [0,pi-1)
    map[i] = poly1;
  }
  return *this;
}

// Constructors

SingleCRT::SingleCRT(const ZZX&poly, const FHEcontext& _context, const IndexSet& s) : context(_context)
{
  assert(s.last() < context.numPrimes());

  map.insert(s);
  *this = poly;           // convert polynomial to singleCRT representation
}

SingleCRT::SingleCRT(const ZZX&poly, const FHEcontext& _context)
: context(_context)
{
  IndexSet s = IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);
  *this = poly;           // convert polynomial to singleCRT representation
}

// Uses the "active context", run-time error if it is NULL
SingleCRT::SingleCRT(const ZZX&poly)
: context(*activeContext)
{
  IndexSet s = IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);
  *this = poly;           // convert polynomial to singleCRT representation
}



// Without specifying a ZZX, we get the zero polynomial
SingleCRT::SingleCRT(const FHEcontext &_context, const IndexSet& s)
: context(_context)
{
  assert(s.last() < context.numPrimes());

  map.insert(s);
  // default constructor for ZZX creates the zero polynomial
}

SingleCRT::SingleCRT(const FHEcontext &_context)
: context(_context)
{
  IndexSet s = IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);
  // default constructor for ZZX creates the zero polynomial
}


// Uses the "active context", run-time error if it is NULL
SingleCRT::SingleCRT(): context(*activeContext)
{
  IndexSet s = IndexSet(0, context.numPrimes()-1);
  // FIXME: maybe the default index set should be determined by context?

  map.insert(s);
  // default constructor for ZZX creates the zero polynomial
}


  // Assignment operators

SingleCRT& SingleCRT::operator=(const SingleCRT& other) 
{
  if (&context != &other.context) 
    Error("SingleCRT assignment: context mismatch");

  map = other.map;
  return *this;
}


SingleCRT& SingleCRT::operator=(const DoubleCRT& other)
{
  other.toSingleCRT(*this);
  return *this;
}

SingleCRT& SingleCRT::operator=(const ZZX& poly)
{
  const IndexSet& s = map.getIndexSet();
  ZZX poly1;

  for (long i = s.first(); i <= s.last(); i = s.next(i)) { 
    ZZ pi = to_ZZ(context.ithPrime(i));
    poly1 = poly;
    PolyRed(poly1,pi,true); // the flag true means reduce to [0,pi-1)
    map[i] = poly1;
  }
  return *this;
}

// explicitly changing the number of levels
void SingleCRT::addPrimes(const IndexSet& s1)
{
  assert(empty(s1 & map.getIndexSet()));

  ZZX poly, poly1;
  toPoly(poly); // recover in coefficient representation

  map.insert(s1);  // add new rows to the map

  // fill in new rows
  for (long i = s1.first(); i <= s1.last(); i = s1.next(i)) {
    ZZ pi = to_ZZ(context.ithPrime(i));

    poly1 = poly;
    PolyRed(poly1,pi,true); // the flag true means reduce to [0,pi-1)
    map[i] = poly;
  }
}


// Division by constant
// FIXME: this is not alias friendly
SingleCRT& SingleCRT::operator/=(const ZZ &num)
{
  const IndexSet& s = map.getIndexSet();
  ZZ pi, n;

  for (long i = s.first(); i <= s.last(); i = s.next(i)) {
    pi = to_ZZ(context.ithPrime(i));
    rem(n,num,pi);
    InvMod(n,n,pi);   // n = num^{-1} mod pi

    vec_ZZ& vp = map[i].rep;
    for (long j=0; j<vp.length(); j++) MulMod(vp[j], vp[j], n, pi);
    map[i].normalize();
  }
  return *this;
}

// Recovering the polynomial in coefficient representation. This yields an
// integer polynomial with coefficients in [-P/2,P/2] (P is the product of
// all moduli used). The polynomial is reduced modulo the product of only
// the primes in the IndexSet parameter.
void SingleCRT::toPoly(ZZX& poly, const IndexSet& s) const
{
  IndexSet s1 = map.getIndexSet() & s;

  if (empty(s1)) {
    clear(poly);
    return;
  }

  ZZ p = to_ZZ(context.ithPrime(s1.first()));  // the first modulus

  poly = map[s1.first()];  // Get poly modulo the first prime
  vec_ZZ& vp = poly.rep;

  // ensure that coeficient vector is of size phi(m) with entries in [-p/2,p/2]
  long phim = context.zMStar.getPhiM();
  long vpLength = vp.length();
  if (vpLength<phim) { // just in case of leading zeros in poly
    vp.SetLength(phim);
    for (long j=vpLength; j<phim; j++) vp[j]=0;
  }
  ZZ p_over_2 = p/2;
  for (long j=0; j<phim; j++) if (vp[j] > p_over_2) vp[j] -= p;

  // do incremental integer CRT for other levels  
  for (long i = s1.next(s1.first()); i <= s1.last(); i = s1.next(i)) {
    long q = context.ithPrime(i);       // the next modulus

    // CRT the coefficient vectors of poly and current
    intVecCRT(vp, p, map[i].rep, q);    // defined in the module NumbTh
    p *= q;     // update the modulus
  }
  poly.normalize(); // need to call this after we work on the coeffs
}


void SingleCRT::toPoly(ZZX& p) const
{
  const IndexSet& s = map.getIndexSet();
  toPoly(p, s);
}
