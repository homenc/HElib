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

#ifndef _FHEcontext_H_
#define _FHEcontext_H_
/**
 * @file FHEcontext.h
 * @brief Keeps the parameters of an instance of the cryptosystem
 **/
#include <NTL/xdouble.h>
#include "PAlgebra.h"
#include "CModulus.h"
#include "IndexSet.h"

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

  //! @brief The structure of Z[X]/(Phi_m(X),2)
  PAlgebraMod alMod;

  //! @brief sqrt(variance) of the LWE error (default=3.2)
  xdouble stdev;

  /**
   * @brief The "ciphertext primes", used for fresh ciphertexts.
   *
   * The public encryption key and "fresh" ciphertexts are encrypted relative
   * to only a subset of the primes, to allow for mod-UP during key-switching.
   * See section 3.1.6 in the design document (key-switching).
   * In ctxtPrimes we keep the indexes of this subset. Namely, for a ciphertext
   * part p in a fresh ciphertext we have p.getMap().getIndexSet()==ctxtPrimes.
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

  // Constructors must ensure that alMod points to zMStar

  // constructor
  FHEcontext(unsigned long m, unsigned long p, unsigned long r): zMStar(m, p), alMod(zMStar, r)
  { stdev=3.2; }

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

  //! @brief Add p to the chain, if it's not already there
  void AddPrime(long p, bool special); 


  ///@{
  /**
     @name I/O routines

  To write out all the data associated with a context, do the following:

  \code
    writeContextBase(str, context);
    str << context;
  \endcode

  The first function call writes out just [m p r], which is the data needed
  to invoke the context constructor.

  The second call writes out all other information, including the
  stdev field, the prime sequence (including which primes are "special"),
  and the digits info.

  To read in all the data associated with a context, do the following:

  \code
    unsigned long m, p, r;
    readContextBase(str, m, p, r);

    FHEcontext context(m, p, r);

    str >> context;
  \endcode

  The call to readContextBase just reads the values m, p, r. Then, after
  constructing the context, the >> operator reads in and attaches all other
  information.
  **/

  //! @brief write [m p r] data
  friend void writeContextBase(ostream& str, const FHEcontext& context);

  //! @brief Write all other data
  friend ostream& operator<< (ostream &str, const FHEcontext& context);

  //! @brief read [m p r] data, needed to construct context
  friend void readContextBase(istream& str, unsigned long& m, unsigned long& p, unsigned long& r);

  //! @brief read all other data associated with context
  friend istream& operator>> (istream &str, FHEcontext& context);
  ///@}
};

//! @brief write [m p r] data
void writeContextBase(ostream& s, const FHEcontext& context);
//! @brief read [m p r] data, needed to construct context
void readContextBase(istream& s, unsigned long& m, unsigned long& p, unsigned long& r);

// VJS: compiler seems to need these declarations out here...wtf...

//@{
//! @name Convenience routines for generating the modulus chain

//! @brief Adds to the chain primes whose product is at least e^totalSize, 
//! returns natural log of the product of all added primes
double AddPrimesBySize(FHEcontext& context, double totalSize,
		       bool special=false);

//! @brief Adds n primes to the chain
//! returns natural log of the product of all added primes
double AddPrimesByNumber(FHEcontext& context, long nPrimes, 
			 long startAt=1,
			 bool special=false);

//! @brief Build modulus chain for nLvls levels, using c digits in key-switching
void buildModChain(FHEcontext &context, long nLvls, long c=3);
///@}
extern FHEcontext* activeContext; // Points to the "current" context
#endif // ifndef _FHEcontext_H_
