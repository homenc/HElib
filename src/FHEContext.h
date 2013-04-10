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

#include <NTL/xdouble.h>
#include "PAlgebra.h"
#include "CModulus.h"
#include "IndexSet.h"



long FindM(long k, long L, long c, long p, long d, long s, long chosen_m, bool verbose=false);
// FindM returns smallest modulus m satisfying various constraints:
//    k = security parameter
//    L = # of levels 
//    c = # of columns in key switching matrix
//    p = characteristic of plaintext space
//    d = embedding degree (d == 0 or d== 1 => no constraint)
//    s = min slot size
//    chosen_m = preselected m (0 => not preselected)
// Fails with an error message if no suitable m is found
// prints an informative message if verbose == true

class FHEcontext {
  vector<Cmodulus> moduli;    // Cmodulus objects for the different primes
  // This is private since the implementation assumes that the list of
  // primes only grows and no prime is ever modified or removed.

public:
  // FHEContext is meant for convenience, not encapsulation: Most data
  // members are public and can be initialized by the application program.

  PAlgebra zMStar;       // The structure of Zm*
  PAlgebraMod alMod; // The structure of Z[X]/(Phi_m(X),2)

  // The public encryption key and "fresh" ciphertexts are encrypted relative
  // to only a subset of the prime, to allow for mod-UP during key-switching.
  // In ctxtPrimes we keep the indexes of this subset. Namely, for a ciphertext
  // part p in a fresh ciphertext we have p.getMap().getIndexSet()==ctxtPrimes.
  // For convenience, we also keep in specialPrimes the complemeting subset,
  // i.e., specialPrimes = [0,numPrimes()-1] \setminus ctxtPrimes
  IndexSet ctxtPrimes;
  IndexSet specialPrimes;

  // The different columns in any key-switching matrix contain encryptions
  // of multiplies of the secret key, sk, B1*sk, B2*B1*sk, B3*B2*B1*sk,...
  // with each Bi a product of a few "non-special" primes in the chain. The
  // digits data member indicate which primes correspond to each of the Bi's.
  // These are all IndexSet objects, whose union is the subset ctxtPrimes.
  // The number of Bi's is one less than the number of columns in the key
  // switching matrices (since the 1st column encrypts sk, without any Bi's),
  // but we keep in the digits vector also an entry for the primes that do
  // not participate in any Bi (so digits.size() is the same as the number
  // of columns in the key switching matrices).

  vector<IndexSet> digits; // digits of ctxt/columns of key-switching matrix

  xdouble stdev; // stdev is sqrt(variance) of the LWE error (default=3.2)

  /* Removed from here -- Shai
  //  double VarKeySw; // added noise varaince due to key-switching
  //  double VarModSw; // added noise due to modulus-switching

  The VarModSw and VarKeySw quantities depend on "synamic" quantities such
  as the Hamming weight of the secret key in question, the plaintext space
  the size of digits into which we break the cipehrtext, etc. Hance these
  are not system-wide parameters, instead they are computed in run-time when
  we run the corresponding key-switch/mod-switch operation.
  In principle the sizes of primes in the modulus-chain may depend on this
  added-noise, but we leave it to the application to handle this. One way
  to set these sizes is start with primes as small as possible (for the
  given value of m), then experiment to test these sizes.
  */

  // Constructors must ensure that alMod points to zMStar

  // constructor
  FHEcontext(unsigned m, unsigned p, unsigned r): zMStar(m, p), alMod(zMStar, r)
  { stdev=3.2; }

  bool operator==(const FHEcontext& other) const;
  bool operator!=(const FHEcontext& other) const { return !(*this==other); }

  long ithPrime(unsigned i) const 
  { return (i<moduli.size())? moduli[i].getQ() :0; }

  const Cmodulus& ithModulus(unsigned i) const { return moduli[i]; }

  long numPrimes() const { return moduli.size(); }

  // is num divisible by primes in the chain?
  bool isZeroDivisor(const ZZ& num) const {
    for (unsigned i=0; i<moduli.size(); i++) 
      if (divide(num,moduli[i].getQ())) return true;
    return false;
  }

  bool inChain(long p) const {
    for (unsigned i=0; i<moduli.size(); i++) 
      if (p==moduli[i].getQ()) return true;
    return false;
  }

  // The product of all the primes in the given set
  void productOfPrimes(ZZ& p, const IndexSet& s) const;

  ZZ productOfPrimes(const IndexSet& s) const {
    ZZ p;
    productOfPrimes(p,s);
    return p;
  }

  // FIXME: run-time error when ithPrime(i) returns 0
  double logOfPrime(unsigned i) const { return log(ithPrime(i)); }

  // returns the natural log of productOfPrimes(s)
  double logOfProduct(const IndexSet& s) const {
    if (s.last() >= numPrimes())
      Error("FHEContext::logOfProduct: IndexSet has too many rows");

    double ans = 0.0;
    for (long i = s.first(); i <= s.last(); i = s.next(i))
      ans += logOfPrime(i);
    return ans;
  }

  void AddPrime(long p, bool special); 


  /*************************************************

  I/O routines

  To write out all the data associated with a context,
  do the following:

    writeContextBase(str, context);
    str << context;

  The first function call writes out just [m p r], which is the
  data needed to invoke the context constructor.

  The second call writes out all other information, including the
  stdev field, the prime sequence (including which primes are "special"),
  and the digits info.

  To read in all the data associated with a context, 
  do the following:

    unsigned m, p, r;
    readContextBase(str, m, p, r);

    FHEcontext context(m, p, r);

    str >> context;

  The call to readContextBase just reads the values m, p, r.
  Then, after constructing the context, the >> operator reads
  in and attaches all other information.

  **************************************************/
  


  friend void writeContextBase(ostream& str, const FHEcontext& context);
  // write [m p r] data

  friend ostream& operator<< (ostream &str, const FHEcontext& context);
  // write all other data

  friend void readContextBase(istream& str, unsigned& m, unsigned& p, unsigned& r);
  // read [m p r] data, needed to construct context

  friend istream& operator>> (istream &str, FHEcontext& context);
  // read all other data associated with context
  
};

void writeContextBase(ostream& s, const FHEcontext& context);
void readContextBase(istream& s, unsigned& m, unsigned& p, unsigned& r);
// VJS: compiler seems to need these declarations out here...wtf...


// Convenience routines

// Adds to the chain primes whose product is at least e^totalSize, 
// returns natural log of the product of all added primes
double AddPrimesBySize(FHEcontext& context, double totalSize,
		       bool special=false);

// Adds nPrimes primes to the chain
// returns natural log of the product of all added primes
double AddPrimesByNumber(FHEcontext& context, long nPrimes, 
			 long startAt=1,
			 bool special=false);

// Build modulus chain for nLvls levels, using nDgts digits in key-switching
void buildModChain(FHEcontext &context, long nLvls, long nDgts=3);

extern FHEcontext* activeContext; // Points to the "current" context

#if 0
#include <map>
// Wrapper functions for the context store

bool ContextExists(unsigned m); // Does the store include a context for m?

// Get context for m, raise error if no such context exists
const FHEcontext& ContextByParam(unsigned m);

// Insert a new context for m, does nothing if a context for m already exists
void NewContext(unsigned m, const FHEcontext& contxt);
#endif

#endif // ifndef _FHEcontext_H_
