/* Copyright (C) 2012-2020 IBM Corp.
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
#ifndef HELIB_DOUBLECRT_H
#define HELIB_DOUBLECRT_H
/**
 * @file DoubleCRT.h
 * @brief Integer polynomials (elements in the ring R_Q) in double-CRT form
 **/
#include <helib/zzX.h>
#include <helib/NumbTh.h>
#include <helib/IndexMap.h>
#include <helib/timing.h>

namespace helib {

class Context;

/**
 * @class DoubleCRTHelper
 * @brief A helper class to enforce consistency within an DoubleCRTHelper object
 *
 * See Section 2.6.2 of the design document (IndexMap)
 */
class DoubleCRTHelper : public IndexMapInit<NTL::vec_long>
{
private:
  long val;

public:
  DoubleCRTHelper(const Context& context);

  /** @brief the init method ensures that all rows have the same size */
  virtual void init(NTL::vec_long& v) { v.FixLength(val); }

  /** @brief clone allocates a new object and copies the content */
  virtual IndexMapInit<NTL::vec_long>* clone() const
  {
    return new DoubleCRTHelper(*this);
  }

private:
  DoubleCRTHelper(); // disable default constructor
};

/**
 * @class DoubleCRT
 * @brief Implementing polynomials (elements in the ring R_Q) in double-CRT
 * form
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
 * primes stored in the associated Context.
 *
 * Arithmetic operations are computed modulo the product of the primes in use
 * and also modulo Phi_m(X). Arithmetic operations can only be applied to
 * DoubleCRT objects relative to the same context, trying to add/multiply
 * objects that have different Context objects will raise an error.
 **/
class DoubleCRT
{
  const Context& context; // the context

  // the data itself: if the i'th prime is in use then map[i] is the std::vector
  // of evaluations wrt this prime
  IndexMap<NTL::vec_long> map;

  //! a "sanity check" method, verifies consistency of the map with
  //! current moduli chain, an error is raised if they are not consistent
  void verify();

  // Generic operators.
  // The behavior when *this and other use different primes depends on the flag
  // matchIndexSets. When it is set to true then the effective modulus is
  // determined by the union of the two index sets; otherwise, the index set
  // of *this.

  class AddFun
  {
  public:
    long apply(long a, long b, long n) { return NTL::AddMod(a, b, n); }
  };

  class SubFun
  {
  public:
    long apply(long a, long b, long n) { return NTL::SubMod(a, b, n); }
  };

  class MulFun
  {
  public:
    long apply(long a, long b, long n) { return NTL::MulMod(a, b, n); }
  };

  template <typename Fun>
  DoubleCRT& Op(const DoubleCRT& other, Fun fun, bool matchIndexSets = true);

  DoubleCRT& do_mul(const DoubleCRT& other, bool matchIndexSets = true);

  template <typename Fun>
  DoubleCRT& Op(const NTL::ZZ& num, Fun fun);

  template <typename Fun>
  DoubleCRT& Op(const NTL::ZZX& poly, Fun fun);

public:
  // Constructors and assignment operators

  // representing an integer polynomial as DoubleCRT. If the set of primes
  // to use is not specified, the resulting object uses all the primes in
  // the context. If the coefficients of poly are larger than the product of
  // the used primes, they are effectively reduced modulo that product

  // Default copy-constructor:
  DoubleCRT(const DoubleCRT& other) = default;

  //! @brief Initializing DoubleCRT from a ZZX polynomial
  //! @param poly The ring element itself, zero if not specified
  //! @param _context The context for this DoubleCRT object, use "current active
  //! context" if not specified
  //! @param indexSet Which primes to use for this object, if not specified then
  //! use all of them
  DoubleCRT(const NTL::ZZX& poly,
            const Context& _context,
            const IndexSet& indexSet);

// FIXME-IndexSet
#if 0
  DoubleCRT(const NTL::ZZX&poly, const Context& _context);
#endif

// FIXME-IndexSet
#if 0
  /**
   * @brief Context is not specified, use the "active context"
   * @note (run-time error if active context is nullptr)
   */
  //  declared "explicit" to avoid implicit type conversion
  explicit DoubleCRT(const NTL::ZZX&poly);
#endif

  //! @brief Same as above, but with zzX's
  DoubleCRT(const zzX& poly, const Context& _context, const IndexSet& indexSet);

// FIXME-IndexSet
#if 0
  DoubleCRT(const zzX&poly, const Context& _context);
  explicit DoubleCRT(const zzX&poly);
#endif

// FIXME-IndexSet
#if 0
  // Without specifying a ZZX, we get the zero polynomial
  explicit DoubleCRT(const Context &_context);
  // declare "explicit" to avoid implicit type conversion
#endif

  //! @brief Also specify the IndexSet explicitly
  DoubleCRT(const Context& _context, const IndexSet& indexSet);

  // Assignment operator, the following two lines are equivalent:
  //    DoubleCRT dCRT(poly, context, indexSet);
  // or
  //    DoubleCRT dCRT(context, indexSet); dCRT = poly;

  DoubleCRT& operator=(const DoubleCRT& other);

  // Copy only the primes in s \intersect other.getIndexSet()
  //  void partialCopy(const DoubleCRT& other, const IndexSet& s);

  DoubleCRT& operator=(const zzX& poly);
  DoubleCRT& operator=(const NTL::ZZX& poly);
  DoubleCRT& operator=(const NTL::ZZ& num);
  DoubleCRT& operator=(const long num)
  {
    *this = NTL::to_ZZ(num);
    return *this;
  }

  //! Get one row of a polynomial
  long getOneRow(NTL::Vec<long>& row, long idx, bool positive = false) const;
  long getOneRow(NTL::zz_pX& row, long idx) const; // This affects NTL's modulus

  //! @brief Recovering the polynomial in coefficient representation.
  //! This yields an integer polynomial with coefficients in [-P/2,P/2],
  //! unless the positive flag is set to true, in which case we get
  //! coefficients in [0,P-1] (P is the product of all moduli used).
  //! Using the optional IndexSet param we compute the polynomial
  //! reduced modulo the product of only the primes in that set.
  void toPoly(NTL::ZZX& p, const IndexSet& s, bool positive = false) const;
  void toPoly(NTL::ZZX& p, bool positive = false) const;

  // The variant toPolyMod has another argument, which is a modulus Q, and it
  // computes toPoly() mod Q. This is offered as a separate function in the
  // hope that one day we will figure out a more efficient method of computing
  // this. Right now it is not implemented
  //
  // void toPolyMod(ZZX& p, const ZZ &Q, const IndexSet& s) const;

  bool operator==(const DoubleCRT& other) const
  {
    assertEq(&context,
             &other.context,
             "Cannot compare DoubleCRTs with different context");
    return map == other.map;
  }

  bool operator!=(const DoubleCRT& other) const { return !(*this == other); }

  // @brief Set to zero
  DoubleCRT& SetZero()
  {
    *this = NTL::ZZ::zero();
    return *this;
  }

  // @brief Set to one
  DoubleCRT& SetOne()
  {
    *this = 1;
    return *this;
  }

  //! @brief Break into n digits,according to the primeSets in context.digits.
  //! See Section 3.1.6 of the design document (re-linearization)
  //! Returns the sum of the canonical embedding of the digits
  NTL::xdouble breakIntoDigits(std::vector<DoubleCRT>& dgts) const;

  //! @brief Expand the index set by s1.
  //! It is assumed that s1 is disjoint from the current index set.
  //! If poly_p != 0, then *poly_p will first be set to the result of applying
  //! toPoly.
  void addPrimes(const IndexSet& s1, NTL::ZZX* poly_p = 0);

  //! @brief Expand index set by s1, and multiply by Prod_{q in s1}.
  //! s1 is disjoint from the current index set, returns log(product).
  double addPrimesAndScale(const IndexSet& s1);

  //! @brief Remove s1 from the index set
  void removePrimes(const IndexSet& s1) { map.remove(s1); }

  //! @ brief make prime set equal to s1
  void setPrimes(const IndexSet& s1)
  {
    addPrimes(s1 / getIndexSet());
    removePrimes(getIndexSet() / s1);
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

  DoubleCRT& operator+=(const DoubleCRT& other) { return Op(other, AddFun()); }

  DoubleCRT& operator+=(const NTL::ZZX& poly) { return Op(poly, AddFun()); }

  DoubleCRT& operator+=(const NTL::ZZ& num) { return Op(num, AddFun()); }

  DoubleCRT& operator+=(long num) { return Op(NTL::to_ZZ(num), AddFun()); }

  DoubleCRT& operator-=(const DoubleCRT& other) { return Op(other, SubFun()); }

  DoubleCRT& operator-=(const NTL::ZZX& poly) { return Op(poly, SubFun()); }

  DoubleCRT& operator-=(const NTL::ZZ& num) { return Op(num, SubFun()); }

  DoubleCRT& operator-=(long num) { return Op(NTL::to_ZZ(num), SubFun()); }

  // These are the prefix versions, ++dcrt and --dcrt.
  DoubleCRT& operator++() { return (*this += 1); };
  DoubleCRT& operator--() { return (*this -= 1); };

  // These are the postfix versions -- return type is void,
  // so it is offered just for style...
  void operator++(int) { *this += 1; };
  void operator--(int) { *this -= 1; };

  // Multiplication
  DoubleCRT& operator*=(const DoubleCRT& other)
  {
    // return Op(other,MulFun());
    return do_mul(other);
  }

  DoubleCRT& operator*=(const NTL::ZZX& poly) { return Op(poly, MulFun()); }

  DoubleCRT& operator*=(const NTL::ZZ& num) { return Op(num, MulFun()); }

  DoubleCRT& operator*=(long num) { return Op(NTL::to_ZZ(num), MulFun()); }

  // NOTE: the matchIndexSets business in the following routines
  // is DEPRECATED.  The requirement is that the prime set of other
  // must contain the prime set of *this.  If not, an exception is raised.
  // Also, if matchIndexSets == true and the prime set of *this does
  // not contain the prime set of other, an exception is also raised.

  // Procedural equivalents, supporting also the matchIndexSets flag
  void Add(const DoubleCRT& other, bool matchIndexSets = true)
  {
    Op(other, AddFun(), matchIndexSets);
  }

  void Sub(const DoubleCRT& other, bool matchIndexSets = true)
  {
    Op(other, SubFun(), matchIndexSets);
  }

  void Mul(const DoubleCRT& other, bool matchIndexSets = true)
  {
    // Op(other, MulFun(), matchIndexSets);
    do_mul(other, matchIndexSets);
  }

  // Division by constant
  DoubleCRT& operator/=(const NTL::ZZ& num);
  DoubleCRT& operator/=(long num) { return (*this /= NTL::to_ZZ(num)); }

  //! @brief Small-exponent polynomial exponentiation
  void Exp(long k);

  //! Apply the automorphism F(X) --> F(X^k)  (with gcd(k,m)=1)
  void automorph(long k);
  DoubleCRT& operator>>=(long k)
  {
    automorph(k);
    return *this;
  }

  //! Compute the complex conjugate, the same as automorph(m-1)
  void complexConj();
  ///@}

  // Utilities

  const Context& getContext() const { return context; }
  const IndexMap<NTL::vec_long>& getMap() const { return map; }
  const IndexSet& getIndexSet() const { return map.getIndexSet(); }

  // Choose random DoubleCRT's, either at random or with small/Gaussian
  // coefficients.

  //! @brief Fills each row i with random ints mod pi, uses NTL's PRG
  void randomize(const NTL::ZZ* seed = nullptr);

  //! Sampling routines:
  //! Each of these return a high probability bound on L-infty norm
  //! of canonical embedding

  //! @brief Coefficients are -1/0/1, Prob[0]=1/2
  double sampleSmall();
  double sampleSmallBounded();

  //! @brief Coefficients are -1/0/1 with pre-specified number of nonzeros
  double sampleHWt(long Hwt);
  double sampleHWtBounded(long Hwt);

  //! @brief Coefficients are Gaussians
  //! Return a high probability bound on L-infty norm of canonical embedding
  double sampleGaussian(double stdev = 0.0);
  double sampleGaussianBounded(double stdev = 0.0);

  //! @brief Coefficients are uniform in [-B..B]
  double sampleUniform(long B);
  NTL::xdouble sampleUniform(const NTL::ZZ& B);

  // used to implement modulus switching
  void scaleDownToSet(const IndexSet& s, long ptxtSpace, NTL::ZZX& delta);

  void FFT(const NTL::ZZX& poly, const IndexSet& s);
  void FFT(const zzX& poly, const IndexSet& s);
  // for internal use

  void reduce() const {} // place-holder for consistent with AltCRT

  // Raw I/O
  void read(std::istream& str);
  void write(std::ostream& str) const;

  // I/O: ONLY the matrix is outputted/recovered, not the moduli chain!! An
  // error is raised on input if this is not consistent with the current chain

  friend std::ostream& operator<<(std::ostream& s, const DoubleCRT& d);
  friend std::istream& operator>>(std::istream& s, DoubleCRT& d);
};

inline void conv(DoubleCRT& d, const NTL::ZZX& p) { d = p; }

// FIXME-IndexSet
#if 0
inline DoubleCRT to_DoubleCRT(const NTL::ZZX& p) {
  return DoubleCRT(p);
}
#endif

inline void conv(NTL::ZZX& p, const DoubleCRT& d) { d.toPoly(p); }

inline NTL::ZZX to_ZZX(const DoubleCRT& d)
{
  NTL::ZZX p;
  d.toPoly(p);
  return p;
}

typedef std::shared_ptr<DoubleCRT> DCRTptr;
typedef std::shared_ptr<NTL::ZZX> ZZXptr;

} // namespace helib

#endif // #ifndef HELIB_DOUBLECRT_H
