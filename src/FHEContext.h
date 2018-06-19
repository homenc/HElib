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
#ifndef _FHEcontext_H_
#define _FHEcontext_H_
/**
 * @file FHEContext.h
 * @brief Keeps the parameters of an instance of the cryptosystem
 **/

#include "PAlgebra.h"
#include "CModulus.h"
#include "IndexSet.h"
#include "recryption.h"

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
long FindM(long k, long L, long c, long p, long d, long s, long chosen_m, bool verbose=false);

// FIXME: The size of primes in the chain should be computed at run-time
#if (NTL_SP_NBITS<44)
#define FHE_p2Size NTL_SP_NBITS
#else
#define FHE_p2Size 44
#endif
#define FHE_p2Bound (1L<<FHE_p2Size)
#define FHE_pSize (FHE_p2Size/2) /* The size of levels in the chain */

class EncryptedArray;
/**
 * @class FHEcontext
 * @brief Maintaining the parameters
 **/
class FHEcontext {
  vector<Cmodulus> moduli;    // Cmodulus objects for the different primes
  // This is private since the implementation assumes that the list of
  // primes only grows and no prime is ever modified or removed.

public:
  // FHEContext is meant for convenience, not encapsulation: Most data
  // members are public and can be initialized by the application program.

  //! @brief The structure of Zm*
  PAlgebra zMStar;

  //! @brief The structure of Z[X]/(Phi_m(X),p^r)
  PAlgebraMod alMod;

  //! @breif A default EncryptedArray
  const EncryptedArray* ea;

  //! @brief sqrt(variance) of the LWE error (default=3.2)
  xdouble stdev;

  //! @brief number of bits per level
  long bitsPerLevel;
  /**
   * @brief The "ciphertext primes", used for fresh ciphertexts.
   *
   * The public encryption key and "fresh" ciphertexts are encrypted relative
   * to only a subset of the primes, to allow for mod-UP during key-switching.
   * See section 3.1.6 in the design document (key-switching). 
   * In ctxtPrimes we keep the indexes of this subset. Namely, for a ciphertext
   * part p in a fresh ciphertext we have p.getMap().getIndexSet()==ctxtPrimes.
   * It is assumed that all the "ciphertext primes" are roughly the same size,
   * except perhaps the first one (with index 0), which could be smaller.
   **/
  IndexSet ctxtPrimes;

  //! @brief All the other primes in the chain.
  //!
  //! For convenience, we also keep in specialPrimes the complemeting subset,
  //! i.e., specialPrimes = [0,numPrimes()-1] setminus ctxtPrimes.
  IndexSet specialPrimes;

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
   * but we keep in the digits vector also an entry for the primes that do
   * not participate in any Bi (so digits.size() is the same as the number
   * of columns in the key switching matrices).
   * See section 3.1.6 in the design document (key-switching).
  **/
  vector<IndexSet> digits; // digits of ctxt/columns of key-switching matrix

  long fftPrimeCount;

  //! Bootstrapping-related data in the context
  RecryptData rcData;
  ThinRecryptData trcData;

  /******************************************************************/
  ~FHEcontext(); // destructor
  FHEcontext(unsigned long m, unsigned long p, unsigned long r,
             const vector<long>& gens = vector<long>(), 
             const vector<long>& ords = vector<long>() );  // constructor

  void makeBootstrappable(const Vec<long>& mvec, long skWht=0,
			  bool conservative=false, bool build_cache=false)
  { 
    rcData.init(*this, mvec, skWht, conservative, build_cache); 
    trcData.init(*this, mvec, skWht, conservative, build_cache); 
  }

  bool isBootstrappable() const { return (rcData.alMod != NULL); }

  bool operator==(const FHEcontext& other) const;
  bool operator!=(const FHEcontext& other) const { return !(*this==other); }

  //! @brief The ith small prime in the modulus chain
  long ithPrime(unsigned long i) const 
  { return (i<moduli.size())? moduli[i].getQ() :0; }

  //! @brief Cmodulus object corresponding to ith small prime in the chain
  const Cmodulus& ithModulus(unsigned long i) const { return moduli[i]; }

  //! @brief Total number of small prime in the chain
  long numPrimes() const { return moduli.size(); }

  //! @brief Is num divisible by any of the primes in the chain?
  bool isZeroDivisor(const ZZ& num) const {
    for (unsigned long i=0; i<moduli.size(); i++) 
      if (divide(num,moduli[i].getQ())) return true;
    return false;
  }

  //! @brief Is p already in the chain?
  bool inChain(long p) const {
    for (unsigned long i=0; i<moduli.size(); i++) 
      if (p==moduli[i].getQ()) return true;
    return false;
  }

  ///@{
  //! @brief The product of all the primes in the given set
  void productOfPrimes(ZZ& p, const IndexSet& s) const;
  ZZ productOfPrimes(const IndexSet& s) const {
    ZZ p;
    productOfPrimes(p,s);
    return p;
  }
  ///@}

  // FIXME: run-time error when ithPrime(i) returns 0
  //! @brief Returns the natural logarithm of the ith prime
  double logOfPrime(unsigned long i) const { return log(ithPrime(i)); }

  //! @brief Returns the natural logarithm of productOfPrimes(s)
  double logOfProduct(const IndexSet& s) const {
    if (s.last() >= numPrimes())
      Error("FHEContext::logOfProduct: IndexSet has too many rows");

    double ans = 0.0;
    for (long i = s.first(); i <= s.last(); i = s.next(i))
      ans += logOfPrime(i);
    return ans;
  }

  //! @brief An estimate for the security-level
  double securityLevel() const {
    long phim = zMStar.getPhiM();
    IndexSet allPrimes(0,numPrimes()-1);
    double bitsize = logOfProduct(allPrimes)/log(2.0);
    return (7.2*phim/bitsize -110);
  }

  //! @brief Find the next prime and add it to the chain
  long AddPrime(long startFrom, long delta, bool special);

  //! @brief Add an FFT prime to the chain, if it's not already there
  //! returns the value of the prime
  long AddFFTPrime(bool special); 

  //! @brief Test if the chain contains a "half-size" ciphertext prime
  // If it exists, the half-size prime must be the first cipehrtext prime.
  // All other primes are assumed to have roughly the same size.
  bool containsSmallPrime() const {
    if (card(ctxtPrimes)<2) return false;
    long fst = ctxtPrimes.first();
    long scnd= ctxtPrimes.next(fst);
    return (logOfPrime(fst) < (0.75 * logOfPrime(scnd)));
  }
  
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
    vector<long> gens, ords;
    
    readContextBase(str, m, p, r, gens, ords);

    FHEcontext context(m, p, r, gens, ords);

    str >> context;
  \endcode

  The call to readContextBase just reads the values m, p, r and the set
  of generators in Zm* /(p) and their order. Then, after constructing the
  context, the >> operator reads in and attaches all other information.
  **/

  //! @brief write [m p r] data
  friend void writeContextBase(ostream& str, const FHEcontext& context);

  //! @brief Write all other data
  friend ostream& operator<< (ostream &str, const FHEcontext& context);

  //! @brief read [m p r] data, needed to construct context
  friend void readContextBase(istream& str, unsigned long& m, unsigned long& p, unsigned long& r,
			      vector<long>& gens, vector<long>& ords);

  //! @brief read all other data associated with context
  friend istream& operator>> (istream &str, FHEcontext& context);
  ///@}

  friend void writeContextBinary(ostream& str, const FHEcontext& context);
  friend void readContextBinary(istream& str, FHEcontext& context);

};

//! @brief write [m p r gens ords] data
void writeContextBase(ostream& s, const FHEcontext& context);
//! @brief read [m p r gens ords] data, needed to construct context
void readContextBase(istream& s, unsigned long& m,
                     unsigned long& p, unsigned long& r,
		     vector<long>& gens, vector<long>& ords);
std::unique_ptr<FHEcontext> buildContextFromAscii(istream& str);

//! @brief write [m p r gens ords] data
void writeContextBaseBinary(ostream& str, const FHEcontext& context);
void writeContextBinary(ostream& str, const FHEcontext& context);

//! @brief read [m p r gens ords] data, needed to construct context
void readContextBaseBinary(istream& s, unsigned long& m,
                           unsigned long& p, unsigned long& r,
                           vector<long>& gens, vector<long>& ords);
std::unique_ptr<FHEcontext> buildContextFromBinary(istream& str);
void readContextBinary(istream& str, FHEcontext& context);


// VJS: compiler seems to need these declarations out here...wtf...

//@{
//! @name Convenience routines for generating the modulus chain

//! @brief Adds several primes to the chain. If byNumber=true then totalSize
//! specifies the number of primes to add. If byNumber=false then totalSize
//! specifies the target naturals log all the added primes.
//! Returns natural log of the product of all added primes.
double AddManyPrimes(FHEcontext& context, double totalSize, 
		     bool byNumber, bool special=false);

//! @brief Adds to the chain primes whose product is at least e^totalSize, 
//! Returns natural log of the product of all added primes.
inline double AddPrimesBySize(FHEcontext& context, double totalSize,
			      bool special=false)
{
  return AddManyPrimes(context, totalSize, false, special);
}

//! @brief Adds n primes to the chain
//! Returns natural log of the product of all added primes.
inline double AddPrimesByNumber(FHEcontext& context, long nPrimes, 
				bool special=false) 
{
  return AddManyPrimes(context, (double)nPrimes, true, special);
}

//! @brief Build modulus chain with nLevels levels, using c digits in key-switching
void buildModChain(FHEcontext &context, long nLevels, long c=3,
                   long extraBits=0);
//FIXME: The extraBits params is a hack, used to get around some
//       circularity when making the context boostrappable

///@}
extern FHEcontext* activeContext; // Should point to the "current" context
#endif // ifndef _FHEcontext_H_
