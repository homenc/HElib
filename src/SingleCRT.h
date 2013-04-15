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
/**
 * @file SingleCRT.h
 * @brief Decleration for the helper SingleCRT class
 **/
#ifndef _SingleCRT_H_
#define _SingleCRT_H_
/**
 * @class SingleCRT
 * @brief This class hold integer polynomials modulo many small primes
 *
 * SingleCRT form is a map from an index set to (ZZX) polynomials, where the
 * i'th polynomial contains the coefficients wrt the ith prime in the chain of
 * moduli. The polynomial thus represented is defined modulo the product of
 * all the primes in use.
 *
 * This is mostly a helper class for DoubleCRT, and all the rules that apply
 * to DoubleCRT wrt moduli sets apply also to SingleCRT objects. Although
 * SingleCRT and DoubleCRT objects can interact in principle, translation back
 * and forth are expensive since they involve FFT/iFFT. Hence support for
 * interaction between them is limited to explicit conversions.
 */
#include <vector>
#include <iostream>
#include <NTL/ZZX.h>

#include "FHEContext.h"
#include "IndexMap.h"
#include "DoubleCRT.h"

class SingleCRT {
  const FHEcontext& context;
  IndexMap<ZZX> map;           // the SingleCRT data

  friend class DoubleCRT;

  // a "sanity check" function, verifies consistency of polys with current
  // moduli chain an error is raised if they are not consistent
  void verify();

  // Generic operators, Fnc is either AddMod, SubMod, or MulMod. (This should
  // have been a template, but gcc refuses to cooperate.) 

  // The behavior when the *this and other have different index sets depends
  // on the matchIndexSets flag. When it is set to true then the union of the
  // two index sets is used; if false, then the index set of *this is used.
  SingleCRT& Op(const SingleCRT &other, 
		void (*Fnc)(ZZ&, const ZZ&, const ZZ&, const ZZ&),
		bool matchIndexSets=true);

  SingleCRT& Op(const ZZX &poly,
		void (*Fnc)(ZZ&, const ZZ&, const ZZ&, const ZZ&));
  SingleCRT& Op(const ZZ &num, void (*Fnc)(ZZX&, const ZZX&, const ZZ&));


 public:


  // Constructors and assignment operators

  // representing an integer polynomial as SingleCRT. If the index set
  // is not specified, then all the primes from the context are used.
  // If the coefficients of poly are larger than the product of
  // the used moduli, they are effectively reduced modulo that product

  SingleCRT(const ZZX& poly, const FHEcontext& _context, const IndexSet& s);
  SingleCRT(const ZZX& poly, const FHEcontext& _context);
  SingleCRT(const ZZX& poly); // uses active context
  

 // Without specifying a ZZX, we get the zero polynomial

  SingleCRT(const FHEcontext& _context, const IndexSet& s);
  SingleCRT(const FHEcontext& _context);

  SingleCRT(); // uses the "active context", run-time error if it is NULL

  // Assignment operators

  SingleCRT& operator=(const SingleCRT& other);
  SingleCRT& operator=(const DoubleCRT& dcrt);
  SingleCRT& operator=(const ZZX& poly);
  SingleCRT& operator=(const ZZ& num)  { *this = to_ZZX(num); return *this; }
  SingleCRT& operator=(const long num) { *this = to_ZZX(num); return *this; }

  bool operator==(const SingleCRT& other) const {
    return &context == &other.context && map == other.map;
  }

  bool operator!=(const SingleCRT& other) const { 
    return !(*this==other); 
  }

  // Set to zero, one

  SingleCRT& setZero() { *this = ZZX::zero(); return *this; }
  SingleCRT& setOne()  { *this = 1; return *this; }


  // expand the index set by s1.
  // it is assumed that s1 is disjoint from the current index set.
  void addPrimes(const IndexSet& s1);


  // remove s1 from the index set
  void removePrimes(const IndexSet& s1) {
    map.remove(s1);
  }

  // Arithmetic operations. Only the "destructive" versions are used,
  // i.e., a += b is implemented but not a + b. 

  // Addition, negation, subtraction
  SingleCRT& operator+=(const SingleCRT &other){ return Op(other,NTL::AddMod);}
  SingleCRT& operator+=(const ZZX &poly)       { return Op(poly, NTL::AddMod);}
  SingleCRT& operator+=(const ZZ &num)         { return Op(num,  NTL::add); }
  SingleCRT& operator+=(long num)       { return Op(to_ZZ(num),  NTL::add); }

  SingleCRT& operator-=(const SingleCRT &other){ return Op(other,NTL::SubMod);}
  SingleCRT& operator-=(const ZZX &poly)       { return Op(poly, NTL::SubMod); }
  SingleCRT& operator-=(const ZZ &num)         { return Op(num,  NTL::sub); }
  SingleCRT& operator-=(long num)       { return Op(to_ZZ(num),  NTL::sub); }

  // Procedural equivalents, supporting also the matchIndexSets flag
  void Add(const SingleCRT &other, bool matchIndexSet=true) 
  { Op(other, NTL::AddMod, matchIndexSet); }
  void Sub(const SingleCRT &other, bool matchIndexSet=true) 
  { Op(other, NTL::SubMod, matchIndexSet); }

  // These are the prefix versions, ++dcrt and --dcrt. 
  SingleCRT& operator++() { return (*this += 1); };
  SingleCRT& operator--() { return (*this -= 1); };

  // Postfix version...no return value...just for style
  void operator++(int) { *this += 1; };
  void operator--(int) { *this -= 1; };

  // Multiplication by constant
  SingleCRT& operator*=(const ZZ &num)          { return Op(num,NTL::mul); }
  SingleCRT& operator*=(long num)        { return Op(to_ZZ(num),NTL::mul); }

  // Division by constant
  SingleCRT& operator/=(const ZZ &num);
  SingleCRT& operator/=(long num) { return (*this /= to_ZZ(num)); }

  // Recovering the polynomial in coefficient representation. This yields an
  // integer polynomial with coefficients in [-P/2,P/2] (P is the product of
  // all moduli used). The polynomial is reduced modulo the product of only
  // the primes in the IndexSet parameter.

  void toPoly(ZZX& p, const IndexSet& s) const;
  void toPoly(ZZX& p) const;

#if 0
  // I/O: ONLY the ZZX vector is outputted/recovered, not the moduli chain!! An
  // error is raised on input if this is not consistent with the current chain

  friend ostream& operator<<(ostream &s, const SingleCRT &scrt) 
  { s << scrt.polys; return s; }

  friend istream& operator>> (istream &s, SingleCRT &scrt)
  { s >> scrt.polys; scrt.verify(); return s; }
#endif

  // Access methods
  const IndexMap<ZZX>& getMap() const { return map; } 
  const FHEcontext& getContext() const { return context; }
};
inline void conv(SingleCRT &s, const ZZX &p) { s=p; }
// Cannot implement to_SingleCRT(p), since modChain is not defined

inline void conv(ZZX &p, const SingleCRT &s) { s.toPoly(p); }
inline ZZX to_ZZX(const SingleCRT &s)  { ZZX p; s.toPoly(p); return p; }

inline void conv(SingleCRT &s, const DoubleCRT &d) { s=d; }
#endif // #ifndef _SingleCRT_H_
