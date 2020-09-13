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
#ifndef HELIB_CONTEXT_H
#define HELIB_CONTEXT_H
/**
 * @file Context.h
 * @brief Keeps the parameters of an instance of the cryptosystem
 **/
#include <helib/PAlgebra.h>
#include <helib/CModulus.h>
#include <helib/IndexSet.h>
#include <helib/recryption.h>
#include <helib/primeChain.h>
#include <helib/powerful.h>
#include <helib/apiAttributes.h>
#include <helib/range.h>

#include <NTL/Lazy.h>

namespace helib {

/**
 * @brief Returns smallest parameter m satisfying various constraints:
 * @param k security parameter
 * @param L number of levels
 * @param c number of columns in key switching matrices
 * @param p characteristic of plaintext space
 * @param d embedding degree (d ==0 or d==1 => no constraint)
 * @param s at least that many plaintext slots
 * @param chosen_m preselected value of m (0 => not preselected)
 * Fails with an error message if no suitable m is found
 * prints an informative message if verbose == true
 **/
long FindM(long k,
           long nBits,
           long c,
           long p,
           long d,
           long s,
           long chosen_m,
           bool verbose = false);

class EncryptedArray;
struct PolyModRing;
/**
 * @class Context
 * @brief Maintaining the parameters
 **/
class Context
{
  std::vector<Cmodulus> moduli; // Cmodulus objects for the different primes
  // This is private since the implementation assumes that the list of
  // primes only grows and no prime is ever modified or removed.

public:
  // Context is meant for convenience, not encapsulation: Most data
  // members are public and can be initialized by the application program.

  //! @brief The structure of Zm*
  PAlgebra zMStar;

  //! @brief The structure of Z[X]/(Phi_m(X),p^r)
  PAlgebraMod alMod;

  //! @brief A default EncryptedArray
  std::shared_ptr<const EncryptedArray> ea;

  std::shared_ptr<const PowerfulDCRT> pwfl_converter;

  /** @brief The structure of a single slot of the plaintext space.
   *
   * This will be Z[X]/(G(x),p^r) for some irreducible factor G of Phi_m(X).
   **/
  std::shared_ptr<PolyModRing> slotRing;

  //! @brief sqrt(variance) of the LWE error (default=3.2)
  NTL::xdouble stdev;

  //======================= high probability bounds ================
  double scale; // default = 10

  //! erfc(scale/sqrt(2)) * phi(m) should be less than some negligible
  //! parameter epsilon.
  //! The default value of 10 should be good enough for most applications.
  //! NOTE: -log(erfc(8/sqrt(2)))/log(2)  = 49.5
  //!       -log(erfc(10/sqrt(2)))/log(2) = 75.8
  //!       -log(erfc(11/sqrt(2)))/log(2) = 91.1
  //!       -log(erfc(12/sqrt(2)))/log(2) =107.8

  //! The way this is used is as follows. If we have a normal random
  //! variable X with variance sigma^2, then the probability that
  //! that X lies outside the interval [-scale*sigma, scale*sigma] is
  //! delta=erfc(scale/sqrt(2)). We will usually apply the union bound
  //! to a vector of phi(m) such random variables (one for each primitive
  //! m-th root of unity), so that the probability that that the L-infty
  //! norm exceeds scale*sigma is at most epsilon=phim*delta. Thus,
  //! scale*sigma will be used as a high-probability bound on the
  //! L-infty norm of such vectors.

  //=======================================

  //! Assume the polynomial f(x) = sum_{i < k} f_i x^i is chosen so
  //! that each f_i is chosen uniformly and independently from the
  //! interval [-magBound, magBound], and that k = degBound.
  //! This returns a bound B such that the L-infty norm
  //! of the canonical embedding exceeds B with probability at most
  //! epsilon.

  // NOTE: this is a bit heuristic: we assume that if we evaluate
  // f at a primitive root of unity, then we get something that well
  // approximates a normal random variable with the same variance,
  // which is equal to the sum of the variances of the individual
  // f_i's, which is (2*magBound)^2/12 = magBound^2/3.
  // We then multiply the sqrt of the variance by scale to get
  // the high probability bound.

  double noiseBoundForUniform(double magBound, long degBound) const
  {
    return scale * std::sqrt(double(degBound) / 3.0) * magBound;
  }

  NTL::xdouble noiseBoundForUniform(NTL::xdouble magBound, long degBound) const
  {
    return scale * std::sqrt(double(degBound) / 3.0) * magBound;
  }

  //=======================================

  //! Assume the polynomial f(x) = sum_{i < k} f_i x^i is chosen so
  //! that each f_i is chosen uniformly and independently from the
  //! from the set of balanced residues modulo the given modulus.
  //! This returns a bound B such that the L-infty norm
  //! of the canonical embedding exceeds B with probability at most
  //! epsilon.

  // NOTE: for odd modulus, this means each f_i is uniformly distributed
  // over { -floor(modulus/2), ..., floor(modulus/2) }.
  // For even modulus, this means each f_i is uniformly distributed
  // over { modulus/2, ..., modulus/2 }, except that the two endpoints
  // (which represent the same residue class) occur with half the
  // probability of the others.

  // NOTE: this is a bit heuristic: we assume that if we evaluate
  // f at a primitive root of unity, then we get something that well
  // approximates a normal random variable with the same variance,
  // which is equal to the sum of the variances of the individual
  // f_i's, which is (modulus)^2/12 + 1/6 for even modulus,
  // and is at most (modulus^2)/12 for odd modulus.
  // We then multiply the sqrt of the variance by scale to get
  // the high probability bound.

  // NOTE: this is slightly more accurate that just calling
  // noiseBoundForUniform with magBound=modulus/2.

  double noiseBoundForMod(long modulus, long degBound) const
  {
    double var = fsquare(modulus) / 12.0;
    if (modulus % 2 == 0)
      var += 1.0 / 6.0;

    return scale * std::sqrt(degBound * var);
  }

  //=======================================

  //! Assume the polynomial f(x) = sum_{i < k} f_i x^i is chosen
  //! so that each f_i is chosen uniformly and independently from
  //! N(0, sigma^2), and that k = degBound.
  //! This returns a bound B such that the L-infty norm
  //! of the canonical embedding exceeds B with probability at most
  //! epsilon.

  // NOTE: if we evaluate f at a primitive root of unity,
  // then we get a normal random variable variance degBound * sigma^2.
  // We then multiply the sqrt of the variance by scale to get
  // the high probability bound.

  double noiseBoundForGaussian(double sigma, long degBound) const
  {
    return scale * std::sqrt(double(degBound)) * sigma;
  }

  //=======================================

  //! Assume the polynomial f(x) = sum_{i < k} f_i x^i is chosen
  //! so that each f_i is zero with probability 1-prob, 1 with probability
  //! prob/2, and -1 with probability prob/2.
  //! This returns a bound B such that the L-infty norm
  //! of the canonical embedding exceeds B with probability at most
  //! epsilon.

  // NOTE: this is a bit heuristic: we assume that if we evaluate
  // f at a primitive root of unity, then we get something that
  // well approximates a normal random variable with the same variance,
  // which is equal to the sum of the individual variances,
  // which is degBound*prob.
  // We then multiply the sqrt of the variance by scale to get
  // the high probability bound.

  double noiseBoundForSmall(double prob, long degBound) const
  {
    return scale * std::sqrt(double(degBound)) * std::sqrt(prob);
  }

  //=======================================

  //! Assume the polynomial f(x) = sum_{i < k} f_i x^i is chosen
  //! hwt coefficients are chosen to \pm 1, and the remainder zero.
  //! This returns a bound B such that the L-infty norm
  //! of the canonical embedding exceeds B with probability at most
  //! epsilon.

  // NOTE: this is a bit heuristic: we assume that if we evaluate
  // f at a primitive root of unity, then we get something that
  // well approximates a normal random variable with the same variance,
  // which is hwt.
  // We then multiply the sqrt of the variance by scale to get
  // the high probability bound.

  // NOTE: degBound is not used here, but I include it
  // for consistency with the other noiseBound routines

  double noiseBoundForHWt(long hwt, UNUSED long degBound) const
  {
    return scale * std::sqrt(double(hwt));
  }

  //=======================================

  //! This computes a high probability bound on the L-infty norm
  //! of x0+s*x1 in the pwrfl basis, assuming is chosen with coeffs
  //! in the pwrfl basis uniformly and independently dist'd over [-1/2,1/2],
  //! x0 has arbitrary coeffs over [-1/2,1/2] in the pwrfl basis,
  //! and assuming s is chosen with skHwt nonzero coeffs mod X^m-1
  //! in the power basis (uniformly and independently over {-1,1}).
  //! The bound should be satisfied with probability epsilon.

  //! NOTE: this is a bit heuristic. See design document for details.

  //! NOTE: this is still valid even when m is a power of 2

  double stdDevForRecryption(long skHwt = 0) const
  {
    if (!skHwt)
      skHwt = rcData.skHwt;
    // the default reverts to rcData.skHwt, *not* rcData.defSkHwt

    long k = zMStar.getNFactors();
    // number of prime factors of m

    long m = zMStar.getM();
    long phim = zMStar.getPhiM();

    double mrat = double(phim) / double(m);

    return std::sqrt(mrat * double(skHwt) * double(1L << k) / 3.0) * 0.5;
  }

  double boundForRecryption(long skHwt = 0) const
  {
    double c_m = zMStar.get_cM();
    // multiply by this fudge factor

    return 0.5 + c_m * scale * stdDevForRecryption(skHwt);
  }

  /**
   * The "ciphertext primes" are the "normal" primes that are used to
   * represent the public encryption key and ciphertexts. These are all
   * "large" single=precision primes, or bit-size roughly NTL_SP_SIZE bits.
   **/
  IndexSet ctxtPrimes;

  //! A disjoint set of primes, used for key switching. See section 3.1.6
  //! in the design document (key-switching). These too are "large"
  //! single=precision primes, or bit-size close to NTL_SP_SIZE bits.
  IndexSet specialPrimes;

  //! Yet a third set of primes, aimed at allowing modulus-switching with
  //! higher resolution. These are somewhat smaller single-precision
  //! primes, of size from NTL_SP_SIZE-20 to NTL_SP_SIZE-1.
  IndexSet smallPrimes;

  //! A helper table to map required modulo-sizes to primeSets
  ModuliSizes modSizes;
  void setModSizeTable() { modSizes.init(*this); }

  /**
   * @brief The set of primes for the digits.
   *
   * The different columns in any key-switching matrix contain encryptions
   * of multiplies of the secret key, sk, B1*sk, B2*B1*sk, B3*B2*B1*sk,...
   * with each Bi a product of a few "non-special" primes in the chain. The
   * digits data member indicate which primes correspond to each of the Bi's.
   * These are all IndexSet objects, whose union is the subset ctxtPrimes.
   *
   * The number of Bi's is one less than the number of columns in the key
   * switching matrices (since the 1st column encrypts sk, without any Bi's),
   * but we keep in the digits std::vector also an entry for the primes that do
   * not participate in any Bi (so digits.size() is the same as the number
   * of columns in the key switching matrices).
   * See section 3.1.6 in the design document (key-switching).
   **/
  // Digits of ctxt/columns of key-switching matrix
  std::vector<IndexSet> digits;

  //! Bootstrapping-related data in the context
  // includes both thin and thick
  ThinRecryptData rcData;

  /******************************************************************/
  // constructor
  Context(unsigned long m,
          unsigned long p,
          unsigned long r,
          const std::vector<long>& gens = std::vector<long>(),
          const std::vector<long>& ords = std::vector<long>());

  void makeBootstrappable(const NTL::Vec<long>& mvec,
                          long skWht = 0,
                          bool build_cache = false,
                          bool alsoThick = true)
  {
    rcData.init(*this, mvec, alsoThick, skWht, build_cache);
  }

  bool isBootstrappable() const { return rcData.alMod != nullptr; }

  IndexSet fullPrimes() const { return ctxtPrimes | specialPrimes; }

  IndexSet allPrimes() const
  {
    return smallPrimes | ctxtPrimes | specialPrimes;
  }

  // returns first nprimes ctxtPrimes
  IndexSet getCtxtPrimes(long nprimes) const
  {
    long first = ctxtPrimes.first();
    long last = std::min(ctxtPrimes.last(), first + nprimes - 1);
    return IndexSet(first, last);
  }

  // FIXME: replacement for bitsPerLevel...placeholder for now
  long BPL() const { return 30; }

  bool operator==(const Context& other) const;
  bool operator!=(const Context& other) const { return !(*this == other); }

  //! @brief The ith small prime in the modulus chain
  long ithPrime(unsigned long i) const
  {
    return (i < moduli.size()) ? moduli[i].getQ() : 0;
  }

  //! @brief Cmodulus object corresponding to ith small prime in the chain
  const Cmodulus& ithModulus(unsigned long i) const { return moduli[i]; }

  //! @brief Total number of small prime in the chain
  long numPrimes() const { return moduli.size(); }

  //! @brief Is num divisible by any of the primes in the chain?
  bool isZeroDivisor(const NTL::ZZ& num) const
  {
    for (long i : range(moduli.size()))
      if (divide(num, moduli[i].getQ()))
        return true;
    return false;
  }

  //! @brief Is p already in the chain?
  bool inChain(long p) const
  {
    for (long i : range(moduli.size()))
      if (p == moduli[i].getQ())
        return true;
    return false;
  }

  ///@{
  //! @brief The product of all the primes in the given set
  void productOfPrimes(NTL::ZZ& p, const IndexSet& s) const;
  NTL::ZZ productOfPrimes(const IndexSet& s) const
  {
    NTL::ZZ p;
    productOfPrimes(p, s);
    return p;
  }
  ///@}

  // FIXME: run-time error when ithPrime(i) returns 0
  //! @brief Returns the natural logarithm of the ith prime
  double logOfPrime(unsigned long i) const { return log(ithPrime(i)); }

  //! @brief Returns the natural logarithm of productOfPrimes(s)
  double logOfProduct(const IndexSet& s) const
  {
    if (s.last() >= numPrimes())
      throw RuntimeError("Context::logOfProduct: IndexSet has too many rows");

    double ans = 0.0;
    for (long i : s)
      ans += logOfPrime(i);
    return ans;
  }

  //! @brief Size in bits of Q.
  long bitSizeOfQ() const
  {
    IndexSet primes = ctxtPrimes | specialPrimes;
    return std::ceil(logOfProduct(primes) / log(2.0));
  }

  //! @brief An estimate for the security-level. This has a lower bound of 0.
  double securityLevel() const
  {
    long phim = zMStar.getPhiM();
    IndexSet primes = ctxtPrimes | specialPrimes;

    if (primes.card() == 0) {
      throw LogicError(
          "Security level cannot be determined as modulus chain is empty.");
    }

    double bitsize = logOfProduct(primes) / log(2.0);
    double ret = (7.2 * phim / bitsize - 110);
    return ret < 0.0 ? 0.0 : ret; // If ret is negative then return 0.0
  }

  //! @brief print out algebra and other important info
  void printout(std::ostream& out = std::cout) const;

  //! @brief Just add the given prime to the chain
  void AddSmallPrime(long q);
  void AddCtxtPrime(long q);
  void AddSpecialPrime(long q);

  ///@{
  /**
     @name I/O routines

  To write out all the data associated with a context, do the following:

  \code
    writeContextBase(str, context);
    str << context;
  \endcode

  The first function call writes out just [m p r gens ords], which is the
  data needed to invoke the context constructor.

  The second call writes out all other information, including the
  stdev field, the prime sequence (including which primes are "special"),
  and the digits info.

  To read in all the data associated with a context, do the following:

  \code
    unsigned long m, p, r;
    std::vector<long> gens, ords;

    readContextBase(str, m, p, r, gens, ords);

    Context context(m, p, r, gens, ords);

    str >> context;
  \endcode

  The call to readContextBase just reads the values m, p, r and the set
  of generators in Zm* /(p) and their order. Then, after constructing the
  context, the >> operator reads in and attaches all other information.
  **/

  //! @brief write [m p r] data
  friend void writeContextBase(std::ostream& str, const Context& context);

  //! @brief Write all other data
  friend std::ostream& operator<<(std::ostream& str, const Context& context);

  //! @brief read [m p r] data, needed to construct context
  friend void readContextBase(std::istream& str,
                              unsigned long& m,
                              unsigned long& p,
                              unsigned long& r,
                              std::vector<long>& gens,
                              std::vector<long>& ords);

  //! @brief read all other data associated with context
  friend std::istream& operator>>(std::istream& str, Context& context);
  ///@}

  friend void writeContextBinary(std::ostream& str, const Context& context);
  friend void readContextBinary(std::istream& str, Context& context);
};

//! @brief write [m p r gens ords] data
void writeContextBase(std::ostream& s, const Context& context);
//! @brief read [m p r gens ords] data, needed to construct context
void readContextBase(std::istream& s,
                     unsigned long& m,
                     unsigned long& p,
                     unsigned long& r,
                     std::vector<long>& gens,
                     std::vector<long>& ords);
std::unique_ptr<Context> buildContextFromAscii(std::istream& str);

//! @brief write [m p r gens ords] data
void writeContextBaseBinary(std::ostream& str, const Context& context);
void writeContextBinary(std::ostream& str, const Context& context);

//! @brief read [m p r gens ords] data, needed to construct context
void readContextBaseBinary(std::istream& s,
                           unsigned long& m,
                           unsigned long& p,
                           unsigned long& r,
                           std::vector<long>& gens,
                           std::vector<long>& ords);
std::unique_ptr<Context> buildContextFromBinary(std::istream& str);
void readContextBinary(std::istream& str, Context& context);

// Build modulus chain with nBits worth of ctxt primes,
// using nDgts digits in key-switching.
// The willBeBootstrappable and skHwt parameters are needed to get around some
// some circularity when making the context boostrappable.
// If you later call context.makeBootstrappable with a given value
// of skHwt, you should first buildModChain with willBeBootstrappable
// set to true and the given value of skHwt.
// FIXME: We should really have a simpler way to do this.
// resolution ... FIXME

void buildModChain(Context& context,
                   long nBits,
                   long nDgts = 3,
                   bool willBeBootstrappable = false,
                   long skHwt = 0,
                   long resolution = 3,
                   long bitsInSpecialPrimes = 0);
// should be called if after you build the mod chain in some way
// *other* than calling buildModChain.
void endBuildModChain(Context& context);

///@}
// Should point to the "current" context
extern Context* activeContext;

} // namespace helib

#endif // ifndef HELIB_CONTEXT_H
