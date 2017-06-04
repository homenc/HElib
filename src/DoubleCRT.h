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
#ifndef _DoubleCRT_H_
#define _DoubleCRT_H_
/**
 * @file DoubleCRT.h
 * @brief Integer polynomials (elements in the ring R_Q) in double-CRT form
 **/

#include "NumbTh.h"
#include "IndexMap.h"
#include "FHEContext.h"
#include "timing.h"

/**
* @class DoubleCRTHelper
* @brief A helper class to enforce consistency within an DoubleCRTHelper object
*
* See Section 2.6.2 of the design document (IndexMap)
*/
class DoubleCRTHelper : public IndexMapInit<vec_long> {
private: 
  long val;

public:
  DoubleCRTHelper(const FHEcontext& context) { 
    val = context.zMStar.getPhiM(); 
  }

  /** @brief the init method ensures that all rows have the same size */
  virtual void init(vec_long& v) { 
    v.FixLength(val); 
  }

  /** @brief clone allocates a new object and copies the content */
  virtual IndexMapInit<vec_long> * clone() const { 
    return new DoubleCRTHelper(*this); 
  }
private:
  DoubleCRTHelper(); // disable default constructor
};


/**
 * @class DoubleCRT
 * @brief Implementatigs polynomials (elements in the ring R_Q) in double-CRT form
 *
 * Double-CRT form is a matrix of L rows and phi(m) columns. The i'th row
 * contains the FFT of the element wrt the ith prime, i.e. the evaluations
 * of the polynomial at the primitive mth roots of unity mod the ith prime.
 * The polynomial thus represented is defined modulo the product of all the
 * primes in use.
 *
 * The list of primes is defined by the data member indexMap.
 * indexMap.getIndexSet() defines the set of indices of primes
 * associated with this DoubleCRT object: they index the
 * primes stored in the associated FHEContext.
 *
 * Arithmetic operations are computed modulo the product of the primes in use
 * and also modulo Phi_m(X). Arithmetic operations can only be applied to
 * DoubleCRT objects relative to the same context, trying to add/multiply
 * objects that have different FHEContext objects will raise an error.
 **/
class DoubleCRT {
  const FHEcontext& context; // the context
  IndexMap<vec_long> map; // the data itself: if the i'th prime is in use then
                          // map[i] is the vector of evaluations wrt this prime

  //! a "sanity check" method, verifies consistency of the map with
  //! current moduli chain, an error is raised if they are not consistent
  void verify();

  // Generic operators. 
  // The behavior when *this and other use different primes depends on the flag
  // matchIndexSets. When it is set to true then the effective modulus is
  // determined by the union of the two index sets; otherwise, the index set
  // of *this.

  class AddFun {
  public:
    long apply(long a, long b, long n) { return AddMod(a, b, n); }
  };

  class SubFun {
  public:
    long apply(long a, long b, long n) { return SubMod(a, b, n); }
  };

  class MulFun {
  public:
    long apply(long a, long b, long n) { return MulMod(a, b, n); }
  };


  template<class Fun>
  DoubleCRT& Op(const DoubleCRT &other, Fun fun,
		bool matchIndexSets=true);

  DoubleCRT& do_mul(const DoubleCRT &other, 
		bool matchIndexSets=true);

  template<class Fun>
  DoubleCRT& Op(const ZZ &num, Fun fun);

  template<class Fun>
  DoubleCRT& Op(const ZZX &poly, Fun fun);

public:

  // Constructors and assignment operators

  // representing an integer polynomial as DoubleCRT. If the set of primes
  // to use is not specified, the resulting object uses all the primes in
  // the context. If the coefficients of poly are larger than the product of
  // the used primes, they are effectively reduced modulo that product

  // copy constructor: default

  //! @brief Initializing DoubleCRT from a ZZX polynomial
  //! @param poly The ring element itself, zero if not specified
  //! @param _context The context for this DoubleCRT object, use "current active context" if not specified
  //! @param indexSet Which primes to use for this object, if not specified then use all of them
  DoubleCRT(const ZZX&poly, const FHEcontext& _context, const IndexSet& indexSet);
  DoubleCRT(const ZZX&poly, const FHEcontext& _context);

  //! @brief Context is not specified, use the "active context"
  //  (run-time error if active context is NULL)
  //  declared "explicit" to avoid implicit type conversion
  explicit DoubleCRT(const ZZX&poly); 

  //! @brief Same as above, but with zzX's
  DoubleCRT(const zzX&poly, const FHEcontext& _context, const IndexSet& indexSet);
  DoubleCRT(const zzX&poly, const FHEcontext& _context);
  explicit DoubleCRT(const zzX&poly); 

 // Without specifying a ZZX, we get the zero polynomial
  explicit DoubleCRT(const FHEcontext &_context);
  // declare "explicit" to avoid implicit type conversion

  //! @brief Also specify the IndexSet explicitly
  DoubleCRT(const FHEcontext &_context, const IndexSet& indexSet);

  // Assignment operator, the following two lines are equivalent:
  //    DoubleCRT dCRT(poly, context, indexSet);
  // or
  //    DoubleCRT dCRT(context, indexSet); dCRT = poly;

  DoubleCRT& operator=(const DoubleCRT& other);

  // Copy only the primes in s \intersect other.getIndexSet()
  //  void partialCopy(const DoubleCRT& other, const IndexSet& s);

  DoubleCRT& operator=(const ZZX& poly);
  DoubleCRT& operator=(const ZZ& num);
  DoubleCRT& operator=(const long num) { *this = to_ZZ(num); return *this; }

  //! Get one row of a polynomial
  long getOneRow(Vec<long>& row, long idx, bool positive=false) const;
  long getOneRow(zz_pX& row, long idx) const; // This affects NTL's modulus

  //! @brief Recovering the polynomial in coefficient representation.
  //! This yields an integer polynomial with coefficients in [-P/2,P/2],
  //! unless the positive flag is set to true, in which case we get
  //! coefficients in [0,P-1] (P is the product of all moduli used).
  //! Using the optional IndexSet param we compute the polynomial
  //! reduced modulo the product of only the ptimes in that set.
  void toPoly(ZZX& p, const IndexSet& s, bool positive=false) const;
  void toPoly(ZZX& p, bool positive=false) const;

  // The variant toPolyMod has another argument, which is a modulus Q, and it
  // computes toPoly() mod Q. This is offerred as a separate function in the
  // hope that one day we will figure out a more efficient method of computing
  // this. Right now it is not implemented
  // 
  // void toPolyMod(ZZX& p, const ZZ &Q, const IndexSet& s) const;


  bool operator==(const DoubleCRT& other) const {
    assert(&context == &other.context);
    return map == other.map;
  }

  bool operator!=(const DoubleCRT& other) const { 
    return !(*this==other);
  }

  // @brief Set to zero
  DoubleCRT& SetZero() { 
    *this = ZZ::zero(); 
    return *this; 
  }

  // @brief Set to one
  DoubleCRT& SetOne()  { 
    *this = 1; 
    return *this; 
  }

  //! @brief Break into n digits,according to the primeSets in context.digits.
  //! See Section 3.1.6 of the design document (re-linearization)
  void breakIntoDigits(vector<DoubleCRT>& dgts, long n) const;

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
  DoubleCRT& Negate(const DoubleCRT& other);
  DoubleCRT& Negate() { return Negate(*this); }

  DoubleCRT& operator+=(const DoubleCRT &other) {
    return Op(other, AddFun());
  }

  DoubleCRT& operator+=(const ZZX &poly) {
    return Op(poly, AddFun());
  }

  DoubleCRT& operator+=(const ZZ &num) { 
    return Op(num, AddFun());
  }

  DoubleCRT& operator+=(long num) { 
    return Op(to_ZZ(num), AddFun());
  }

  DoubleCRT& operator-=(const DoubleCRT &other) {
    return Op(other,SubFun());
  }

  DoubleCRT& operator-=(const ZZX &poly) {
    return Op(poly,SubFun());
  }
  
  DoubleCRT& operator-=(const ZZ &num) { 
    return Op(num, SubFun());
  }

  DoubleCRT& operator-=(long num) { 
    return Op(to_ZZ(num), SubFun());
  }

  // These are the prefix versions, ++dcrt and --dcrt. 
  DoubleCRT& operator++() { return (*this += 1); };
  DoubleCRT& operator--() { return (*this -= 1); };

  // These are the postfix versions -- return type is void,
  // so it is offered just for style...
  void operator++(int) { *this += 1; };
  void operator--(int) { *this -= 1; };


  // Multiplication
  DoubleCRT& operator*=(const DoubleCRT &other) {
    //return Op(other,MulFun());
    return do_mul(other);
  }

  DoubleCRT& operator*=(const ZZX &poly) {
    return Op(poly,MulFun());
  }

  DoubleCRT& operator*=(const ZZ &num) { 
    return Op(num,MulFun());
  }

  DoubleCRT& operator*=(long num) { 
    return Op(to_ZZ(num),MulFun());
  }


  // Procedural equivalents, supporting also the matchIndexSets flag
  void Add(const DoubleCRT &other, bool matchIndexSets=true) {
    Op(other, AddFun(), matchIndexSets); 
  }

  void Sub(const DoubleCRT &other, bool matchIndexSets=true) {
    Op(other, SubFun(), matchIndexSets); 
  }

  void Mul(const DoubleCRT &other, bool matchIndexSets=true) {
    // Op(other, MulFun(), matchIndexSets); 
    do_mul(other, matchIndexSets); 
  }

  // Division by constant
  DoubleCRT& operator/=(const ZZ &num);
  DoubleCRT& operator/=(long num) { return (*this /= to_ZZ(num)); }

  //! @brief Small-exponent polynomial exponentiation
  void Exp(long k);

  // Apply the automorphism F(X) --> F(X^k)  (with gcd(k,m)=1)
  void automorph(long k);
  DoubleCRT& operator>>=(long k) { automorph(k); return *this; }
  ///@}

  // Utilities

  const FHEcontext& getContext() const { return context; }
  const IndexMap<vec_long>& getMap() const { return map; }
  const IndexSet& getIndexSet() const { return map.getIndexSet(); }

  // Choose random DoubleCRT's, either at random or with small/Gaussian
  // coefficients. 

  //! @brief Fills each row i with random ints mod pi, uses NTL's PRG
  void randomize(const ZZ* seed=NULL);

  //! @brief Coefficients are -1/0/1, Prob[0]=1/2
  void sampleSmall() {
    ZZX poly; 
    ::sampleSmall(poly,context.zMStar.getPhiM()); // degree-(phi(m)-1) polynomial
    *this = poly; // convert to DoubleCRT
  }

  //! @brief Coefficients are -1/0/1 with pre-specified number of nonzeros
  void sampleHWt(long Hwt) {
    ZZX poly; 
    ::sampleHWt(poly,Hwt,context.zMStar.getPhiM());
    *this = poly; // convert to DoubleCRT
  }

  //! @brief Coefficients are Gaussians
  void sampleGaussian(double stdev=0.0) {
    if (stdev==0.0) stdev=to_double(context.stdev); 
    ZZX poly; 
    ::sampleGaussian(poly, context.zMStar.getPhiM(), stdev);
    *this = poly; // convert to DoubleCRT
  }

  //! @brief Coefficients are uniform in [-B..B]
  void sampleUniform(const ZZ& B) {
    ZZX poly;
    ::sampleUniform(poly, B, context.zMStar.getPhiM());
    *this = poly;
  }


  // used to implement modulus switching
  void scaleDownToSet(const IndexSet& s, long ptxtSpace);


  void FFT(const ZZX& poly, const IndexSet& s);
  void FFT(const zzX& poly, const IndexSet& s);
  // for internal use


  void reduce() const {} // place-holder for consistenct with AltCRT


  // I/O: ONLY the matrix is outputted/recovered, not the moduli chain!! An
  // error is raised on input if this is not consistent with the current chain

  friend ostream& operator<< (ostream &s, const DoubleCRT &d);
  friend istream& operator>> (istream &s, DoubleCRT &d);
};




inline void conv(DoubleCRT &d, const ZZX &p) { d=p; }

inline DoubleCRT to_DoubleCRT(const ZZX& p) {
  return DoubleCRT(p);
}

inline void conv(ZZX &p, const DoubleCRT &d) { d.toPoly(p); }

inline ZZX to_ZZX(const DoubleCRT &d)  { ZZX p; d.toPoly(p); return p; }

typedef shared_ptr<DoubleCRT> DCRTptr;
typedef shared_ptr<ZZX> ZZXptr;

#endif // #ifndef _DoubleCRT_H_
