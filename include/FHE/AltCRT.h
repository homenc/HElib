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
#ifndef _AltCRT_H_
#define _AltCRT_H_
/**
 * @file AltCRT.h
 * @brief Alternative implementation of integer polynomials
 **/


#include "NumbTh.h"
#include "IndexMap.h"
#include "FHEContext.h"


/**
* @class AltCRTHelper
* @brief A helper class to enforce consistency within an AltCRTHelper object
*
* See Section 2.6.2 of the design document (IndexMap)
*/
class AltCRTHelper : public IndexMapInit<zz_pX> {
private: 
  long val;

public:

  AltCRTHelper(const FHEcontext& context) { 
    val = context.zMStar.getPhiM(); 
  }

  virtual void init(zz_pX& v) { 
    v.rep.SetMaxLength(val); 
  }

  virtual IndexMapInit<zz_pX> * clone() const { 
    return new AltCRTHelper(*this); 
  }

private:
  AltCRTHelper(); // disable default constructor
};


// helper routine: computes (a*b mod X^m-1) ...experimental "lazy"
// version.... 
void MulMod1(zz_pX& x, const zz_pX& a, const zz_pX& b, const zz_pXModulus1& f);

//! @class AltCRT
//! @brief Alternative implementation of integer polynomials
class AltCRT {
  const FHEcontext& context; // the context
  IndexMap<zz_pX> map; // the data itself: if the i'th prime is in use then
                          // map[i] is the vector of evaluations wrt this prime

  // a "sanity check" function, verifies consistency of the map with current
  // moduli chain, an error is raised if they are not consistent
  void verify();


  // Generic operators. 
  // The behavior when *this and other use different primes depends on the flag
  // matchIndexSets. When it is set to true then the effective modulus is
  // determined by the union of the two index sets; otherwise, the index set
  // of *this.

  class AddFun {
  public:
    void apply(zz_pX& a, const zz_pX& b, const zz_pXModulus1& f) 
    { return add(a, a, b); }

    void apply(zz_pX& a, zz_p b)
    { return add(a, a, b); }

    bool reduce() { return false; }
  };

  class SubFun {
  public:
    void apply(zz_pX& a, const zz_pX& b, const zz_pXModulus1& f) 
    { return sub(a, a, b); }

    void apply(zz_pX& a, zz_p b)
    { return sub(a, a, b); }

    bool reduce() { return false; }
  };

  class MulFun {
  public:
    void apply(zz_pX& a, const zz_pX& b, const zz_pXModulus1& f) 
    { return MulMod1(a, a, b, f); }

    void apply(zz_pX& a, zz_p b)
    { return mul(a, a, b); }

    bool reduce() { return true; }
  };


  template<class Fun>
  AltCRT& Op(const AltCRT &other, Fun fun,
		bool matchIndexSets=true);

  template<class Fun>
  AltCRT& Op(const ZZ &num, Fun fun);

  template<class Fun>
  AltCRT& Op(const ZZX &poly, Fun fun);

public:

  // Constructors and assignment operators

  // representing an integer polynomial as AltCRT. If the set of primes
  // to use is not specified, the resulting object uses all the primes in
  // the context. If the coefficients of poly are larger than the product of
  // the used primes, they are effectively reduced modulo that product

  // copy constructor: default but with reduction

  AltCRT(const AltCRT& other);

  AltCRT(const ZZX&poly, const FHEcontext& _context, const IndexSet& indexSet);
  AltCRT(const ZZX&poly, const FHEcontext& _context);

  explicit AltCRT(const ZZX&poly); 
  // uses the "active context", run-time error if it is NULL
  // declare "explicit" to avoid implicit type conversion

 // Without specifying a ZZX, we get the zero polynomial

  AltCRT(const FHEcontext &_context, const IndexSet& indexSet);

  explicit AltCRT(const FHEcontext &_context);
  // declare "explicit" to avoid implicit type conversion

  //  AltCRT(); 
  // uses the "active context", run-time error if it is NULL


  // Assignment operator, the following two lines are equivalent:
  //    AltCRT dCRT(poly, context, indexSet);
  // or
  //    AltCRT dCRT(context, indexSet); dCRT = poly;

  AltCRT& operator=(const AltCRT& other);

  // Copy only the primes in s \intersect other.getIndexSet()
  //  void partialCopy(const AltCRT& other, const IndexSet& s);

  AltCRT& operator=(const ZZX& poly);
  AltCRT& operator=(const ZZ& num);
  AltCRT& operator=(const long num) { *this = to_ZZ(num); return *this; }

  //! Get one row of a polynomial
  long getOneRow(Vec<long>& row, long idx, bool positive=false) const;
  long getOneRow(zz_pX& row, long idx) const; // This affects NTL's modulus

  // Recovering the polynomial in coefficient representation. This yields an
  // integer polynomial with coefficients in [-P/2,P/2], unless the positive
  // flag is set to true, in which case we get coefficients in [0,P-1] (P is
  // the product of all moduli used). Using the optional IndexSet param
  // we compute the polynomial reduced modulo the product of only the ptimes
  // in that set.

  void toPoly(ZZX& p, const IndexSet& s, bool positive=false) const;
  void toPoly(ZZX& p, bool positive=false) const;

  // The variant toPolyMod has another argument, which is a modulus Q, and it
  // computes toPoly() mod Q. This is offerred as a separate function in the
  // hope that one day we will figure out a more efficient method of computing
  // this. Right now it is not implemented
  // 
  // void toPolyMod(ZZX& p, const ZZ &Q, const IndexSet& s) const;


  bool operator==(const AltCRT& other) const {
    assert(&context == &other.context);
    return map == other.map;
  }

  bool operator!=(const AltCRT& other) const { 
    return !(*this==other);
  }

  // Set to zero, one
  AltCRT& SetZero() { 
    *this = ZZ::zero(); 
    return *this; 
  }

  AltCRT& SetOne()  { 
    *this = 1; 
    return *this; 
  }

  // break into n digits,according to the primeSets in context.digits
  void breakIntoDigits(vector<AltCRT>& dgts, long n) const;

  // expand the index set by s1.
  // it is assumed that s1 is disjoint from the current index set.
  void addPrimes(const IndexSet& s1);

  // Expand index set by s1, and multiply by \prod{q \in s1}. s1 is assumed to
  // be disjoint from the current index set. Returns the logarithm of product.
  double addPrimesAndScale(const IndexSet& s1);

  // remove s1 from the index set
  void removePrimes(const IndexSet& s1) {
    map.remove(s1);
  }


  // Arithmetic operations. Only the "destructive" versions are used,
  // i.e., a += b is implemented but not a + b.

  // Addition, negation, subtraction
  AltCRT& Negate(const AltCRT& other);
  AltCRT& Negate() { return Negate(*this); }

  AltCRT& operator+=(const AltCRT &other) {
    return Op(other, AddFun());
  }

  AltCRT& operator+=(const ZZX &poly) {
    return Op(poly, AddFun());
  }

  AltCRT& operator+=(const ZZ &num) { 
    return Op(num, AddFun());
  }

  AltCRT& operator+=(long num) { 
    return Op(to_ZZ(num), AddFun());
  }

  AltCRT& operator-=(const AltCRT &other) {
    return Op(other,SubFun());
  }

  AltCRT& operator-=(const ZZX &poly) {
    return Op(poly,SubFun());
  }
  
  AltCRT& operator-=(const ZZ &num) { 
    return Op(num, SubFun());
  }

  AltCRT& operator-=(long num) { 
    return Op(to_ZZ(num), SubFun());
  }

  // These are the prefix versions, ++dcrt and --dcrt. 
  AltCRT& operator++() { return (*this += 1); };
  AltCRT& operator--() { return (*this -= 1); };

  // These are the postfix versions -- return type is void,
  // so it is offered just for style...
  void operator++(int) { *this += 1; };
  void operator--(int) { *this -= 1; };


  // Multiplication
  AltCRT& operator*=(const AltCRT &other) {
    return Op(other,MulFun());
  }

  AltCRT& operator*=(const ZZX &poly) {
    return Op(poly,MulFun());
  }

  AltCRT& operator*=(const ZZ &num) { 
    return Op(num,MulFun());
  }

  AltCRT& operator*=(long num) { 
    return Op(to_ZZ(num),MulFun());
  }


  // Procedural equivalents, supporting also the matchIndexSets flag
  void Add(const AltCRT &other, bool matchIndexSets=true) {
    Op(other, AddFun(), matchIndexSets); 
  }

  void Sub(const AltCRT &other, bool matchIndexSets=true) {
    Op(other, SubFun(), matchIndexSets); 
  }

  void Mul(const AltCRT &other, bool matchIndexSets=true) {
    Op(other, MulFun(), matchIndexSets); 
  }

  // Division by constant
  AltCRT& operator/=(const ZZ &num);
  AltCRT& operator/=(long num) { return (*this /= to_ZZ(num)); }


  // Small-exponent polynomial exponentiation
  void Exp(long k);


  // Apply the automorphism F(X) --> F(X^k)  (with gcd(k,m)=1)
  void automorph(long k);
  AltCRT& operator>>=(long k) { automorph(k); return *this; }


  // Utilities

  const FHEcontext& getContext() const { return context; }
  const IndexMap<zz_pX>& getMap() const { return map; }
  const IndexSet& getIndexSet() const { return map.getIndexSet(); }

  // Choose random AltCRT's, either at random or with small/Gaussian
  // coefficients. 

  // fills each row i w/ random ints mod pi, uses NTL's PRG
  void randomize(const ZZ* seed=NULL);

  // Coefficients are -1/0/1, Prob[0]=1/2
  void sampleSmall() {
    ZZX poly; 
    ::sampleSmall(poly,context.zMStar.getPhiM()); // degree-(phi(m)-1) polynomial
    *this = poly; // convert to AltCRT
  }

  // Coefficients are -1/0/1 with pre-specified number of nonzeros
  void sampleHWt(long Hwt) {
    ZZX poly; 
    ::sampleHWt(poly,Hwt,context.zMStar.getPhiM());
    *this = poly; // convert to AltCRT
  }

  // Coefficients are Gaussians
  void sampleGaussian(double stdev=0.0) {
    if (stdev==0.0) stdev=to_double(context.stdev); 
    ZZX poly; 
    ::sampleGaussian(poly, context.zMStar.getPhiM(), stdev);
    *this = poly; // convert to AltCRT
  }


  // Coefficients are uniform in [-B..B] 
  void sampleUniform(const ZZ& B) {
    ZZX poly;
    ::sampleUniform(poly, B, context.zMStar.getPhiM());
    *this = poly;
  }



  // used to implement modulus switching
  void scaleDownToSet(const IndexSet& s, long ptxtSpace);

  void reduce() const; // reduce rows mod phimx...remove laziness

  // I/O: ONLY the matrix is outputted/recovered, not the moduli chain!! An
  // error is raised on input if this is not consistent with the current chain

  friend ostream& operator<< (ostream &s, const AltCRT &d);
  friend istream& operator>> (istream &s, AltCRT &d);
};




inline void conv(AltCRT &d, const ZZX &p) { d=p; }

inline AltCRT to_AltCRT(const ZZX& p) {
  return AltCRT(p);
}

inline void conv(ZZX &p, const AltCRT &d) { d.toPoly(p); }

inline ZZX to_ZZX(const AltCRT &d)  { ZZX p; d.toPoly(p); return p; }

typedef shared_ptr<AltCRT> DCRTptr;
typedef shared_ptr<ZZX> ZZXptr;

#endif
