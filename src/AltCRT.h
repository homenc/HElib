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
* @brief An alternative representation of ring elements
*
* The AltCRT module offers a drop-in replacement to DoubleCRT, it exposes
* the same interface but internally uses a single-CRT representation. That
* is, polynomials are stored in coefficient representation, modulo each of
* the small primes in our chain. Currently this class is used only for
* testing and debugging purposes.
**/
#include <vector>
#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include "NumbTh.h"
#include "IndexMap.h"
#include "FHEContext.h"
NTL_CLIENT
class SingleCRT;

/**
* @class AltCRTHelper
* @brief A helper class to enforce consistency within an AltCRT object
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

  /** @brief the init method ensures that all rows have the same size */
  virtual void init(zz_pX& v) { 
    v.rep.SetMaxLength(val); 
  }

  /** @brief clone allocates a new object and copies the content */
  virtual IndexMapInit<zz_pX> * clone() const { 
    return new AltCRTHelper(*this); 
  }
private:
  AltCRTHelper(); // disable default constructor
};


/**
* @class AltCRT
* @brief A single-CRT representation of a ring element
* 
* AltCRT offers the same interface as DoubleCRT, but with a different
* internal representation. That is, polynomials are stored in
* coefficient representation, modulo each of the small primes in our
* chain. Currently this class is used only for testing and debugging
* purposes.
*/
class AltCRT {
  const FHEcontext& context; //! the context of this ring element
  IndexMap<zz_pX> map; //! the data itself: if the i'th prime is in use then
                       //! map[i] is the vector of evaluations wrt this prime

  //! a "sanity check" function, verifies consistency of the map with
  //! current moduli chain, an error is raised if they are not consistent
  void verify();

  // Generic operators. 
  // The behavior when *this and other use different primes depends on the flag
  // matchIndexSets. When it is set to true then the effective modulus is
  // determined by the union of the two index sets; otherwise, the index set
  // of *this.

  class AddFun {
  public:
    void apply(zz_pX& a, const zz_pX& b, const zz_pXModulus& f) 
    { return add(a, a, b); }

    void apply(zz_pX& a, zz_p b)
    { return add(a, a, b); }
  };

  class SubFun {
  public:
    void apply(zz_pX& a, const zz_pX& b, const zz_pXModulus& f) 
    { return sub(a, a, b); }

    void apply(zz_pX& a, zz_p b)
    { return sub(a, a, b); }
  };

  class MulFun {
  public:
    void apply(zz_pX& a, const zz_pX& b, const zz_pXModulus& f) 
    { return MulMod(a, a, b, f); }

    void apply(zz_pX& a, zz_p b)
    { return mul(a, a, b); }
  };


  template<class Fun>
  AltCRT& Op(const AltCRT &other, Fun fun,
		bool matchIndexSets=true);

  template<class Fun>
  AltCRT& Op(const ZZ &num, Fun fun);

  template<class Fun>
  AltCRT& Op(const ZZX &poly, Fun fun);

  static bool dryRun; // do not actually perform any of the operations

public:

  // Constructors and assignment operators

  // representing an integer polynomial as AltCRT. If the set of primes
  // to use is not specified, the resulting object uses all the primes in
  // the context. If the coefficients of poly are larger than the product of
  // the used primes, they are effectively reduced modulo that product

  // copy constructor: default

  //! @brief Initializing AltCRT from a ZZX polynomial
  AltCRT(const ZZX&poly, const FHEcontext& _context);

  //! @param poly The ring element itself, zero if not specified
  //! @param _context The context for this AltCRT object, use "current active context" if not specified
  //! @param indexSet Which primes to use for this object, if not specified then use all of them
  AltCRT(const ZZX&poly, const FHEcontext& _context, const IndexSet& indexSet);

  //! @brief Context is not specified, use the "active context"
  //  (run-time error if active context is NULL)
  //  declared "explicit" to avoid implicit type conversion
  explicit AltCRT(const ZZX&poly);

  //! @brief Without specifying a ZZX, we get the zero polynomial
  //  declared "explicit" to avoid implicit type conversion
  explicit AltCRT(const FHEcontext &_context);

  //! @brief Also specify the IndexSet explicitly
  AltCRT(const FHEcontext &_context, const IndexSet& indexSet);


  // Assignment operator, the following two lines are equivalent:
  //    AltCRT dCRT(poly, context, indexSet);
  // or
  //    AltCRT dCRT(context, indexSet); dCRT = poly;

  AltCRT& operator=(const AltCRT& other);
  AltCRT& operator=(const SingleCRT& other);
  AltCRT& operator=(const ZZX& poly);
  AltCRT& operator=(const ZZ& num);
  AltCRT& operator=(const long num) { *this = to_ZZ(num); return *this; }

  //! @brief Recovering the polynomial in coefficient representation.
  void toPoly(ZZX& p, bool positive=false) const;

  //! @brief Recovering the polynomial in coefficient representation.
  //! This yields an integer polynomial with coefficients in [-P/2,P/2],
  //! unless the positive flag is set to true, in which case we get
  //! coefficients in [0,P-1] (P is the product of all moduli used).
  //! Using the optional IndexSet param we compute the polynomial
  //! reduced modulo the product of only the ptimes in that set.
  void toPoly(ZZX& p, const IndexSet& s, bool positive=false) const;

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

  // @brief Set to zero
  AltCRT& SetZero() { 
    *this = ZZ::zero(); 
    return *this; 
  }
  // @brief Set to one
  AltCRT& SetOne()  { 
    *this = 1; 
    return *this; 
  }

  //! @brief Break into n digits,according to the primeSets in context.digits.
  //! See Section 3.1.6 of the design document (re-linearization)
  void breakIntoDigits(vector<AltCRT>& dgts, long n) const;

  //! @brief Expand the index set by s1.
  //! It is assumed that s1 is disjoint from the current index set.
  void addPrimes(const IndexSet& s1);

  //! @brief Expand index set by s1, and multiply by Prod_{q in s1}.
  //! s1 is disjoint from the current index set, returns log(product).
  double addPrimesAndScale(const IndexSet& s1);

  //! @brief Remove s1 from the index set
  void removePrimes(const IndexSet& s1) {
    map.remove(s1);
  }

  /**
   * @name Arithmetic operation
   * @brief Only the "destructive" versions are used,
   * i.e., a += b is implemented but not a + b.
   **/
  ///@{

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


  //! @brief Small-exponent polynomial exponentiation
  void Exp(long k);


  // Apply the automorphism F(X) --> F(X^k)  (with gcd(k,m)=1)
  void automorph(long k);
  AltCRT& operator>>=(long k) { automorph(k); return *this; }
  ///@}

  // Utilities

  const FHEcontext& getContext() const { return context; }
  const IndexMap<zz_pX>& getMap() const { return map; }
  const IndexSet& getIndexSet() const { return map.getIndexSet(); }

  // Choose random AltCRT's, either at random or with small/Gaussian
  // coefficients. 

  //! @brief Fills each row i with random ints mod pi, uses NTL's PRG
  void randomize(const ZZ* seed=NULL);

  //! @brief Coefficients are -1/0/1, Prob[0]=1/2
  void sampleSmall() {
    ZZX poly; 
    ::sampleSmall(poly,context.zMStar.getPhiM()); // degree-(phi(m)-1) polynomial
    *this = poly; // convert to AltCRT
  }

  //! @brief Coefficients are -1/0/1 with pre-specified number of nonzeros
  void sampleHWt(long Hwt) {
    ZZX poly; 
    ::sampleHWt(poly,Hwt,context.zMStar.getPhiM());
    *this = poly; // convert to AltCRT
  }

  //! @brief Coefficients are Gaussians
  void sampleGaussian(double stdev=0.0) {
    if (stdev==0.0) stdev=to_double(context.stdev); 
    ZZX poly; 
    ::sampleGaussian(poly, context.zMStar.getPhiM(), stdev);
    *this = poly; // convert to AltCRT
  }


  //! @brief Makes a corresponding SingleCRT object.
  // Restricted to the given index set, if specified
  void toSingleCRT(SingleCRT& scrt, const IndexSet& s) const;
  void toSingleCRT(SingleCRT& scrt) const;

  // used to implement modulus switching
  void scaleDownToSet(const IndexSet& s, long ptxtSpace);

  // I/O: ONLY the matrix is outputted/recovered, not the moduli chain!! An
  // error is raised on input if this is not consistent with the current chain

  friend ostream& operator<< (ostream &s, const AltCRT &d);
  friend istream& operator>> (istream &s, AltCRT &d);

  //! @brief Used for testing/debugging
  //! The dry-run option disables most operations, to save time. This lets
  //! us quickly go over the evaluation of a circuit and estimate the
  //! resulting noise magnitude, without having to actually compute anything. 
  static bool setDryRun(bool toWhat=true) { dryRun=toWhat; return dryRun; }
};


inline void conv(AltCRT &d, const ZZX &p) { d=p; }

inline AltCRT to_AltCRT(const ZZX& p) {
  return AltCRT(p);
}

inline void conv(ZZX &p, const AltCRT &d) { d.toPoly(p); }

inline ZZX to_ZZX(const AltCRT &d)  { ZZX p; d.toPoly(p); return p; }

inline void conv(AltCRT& d, const SingleCRT& s) { d=s; }





#endif
