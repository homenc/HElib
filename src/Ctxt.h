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
#ifndef _Ctxt_H_
#define _Ctxt_H_
/**
 * @file Ctxt.h
 * @brief Declerations of a BGV-type cipehrtext and key-switching matrices
 *
 * A ciphertext is a vector of "ciphertext parts", each part consists of
 * a polynomial (element of polynomial ring R_Q) and a "handle" describing
 * the secret-key polynomial that this part multiplies during decryption.
 * For example:
 * + A "canonical" ciphertext has two parts, the first part multiplies 1
 *   and the second multiplies the "base" secret key s.
 * + When you multiply two canonical ciphertexts you get a 3-part
 *   ciphertext, with parts corresponding to 1, s, and s^2.
 * + When you apply automorphism X->X^t to a generic ciphertext, then
 *   - the part corresponding to 1 still remains wrt 1
 *   - every other part corresponding to some s' will now be corresponding
 *     to the polynomial s'(X^t) mod Phi_m(X)
 *
 * This type of representation lets you in principle add ciphertexts that
 * are defined with respect to different keys:
 * + For parts of the two cipehrtexts that point to the same secret-key
 *   polynomial, you just add the two Double-CRT polynomials
 * + Parts in one ciphertext that do not have counter-part in the other
 *   cipehrtext will just be included in the result intact.
 * For example, you have the ciphertexts
 *    C1 = (a relative to 1, b relative to s)
 *    C2 = (u relative to 1, v relative to s(X^3))
 * Then their sum will be
 *    C1+C2 = (a+u relative to 1, b relative to s, v relative to s(X^3))
 *
 * Similarly, in principle you can also multiply arbitrary cipehrtexts, even
 * ones that are defined with respect to different keys, and the result will
 * be defined with respect to the tensor product of the two keys. 
 *
 * The current implementation is more restrictive, however. It requires that
 * a ciphertext has one part wrt 1, that for every r>=1 there is at most one
 * part wrt to s^r(X^t) (for some t), and that the r's are consecutive. For
 * example you cannot have parts wrt (1,s,s^3) without having a part wrt s^2.
 *
 * It follows that you can only add/multiply ciphertexts if one of the two
 * lists of handles is a prefix of the other. For example, one can add a
 * ciphertext wrt (1,s(X^2)) to another wrt (1,s(X^2),s^2(X^2)), but not
 * to another ciphertext wrt (1,s).
 **/ 
#include <vector>
#include <NTL/xdouble.h>
#include "DoubleCRT.h"

class KeySwitch;
class FHEPubKey;
class FHESecKey;

/**
 * @class SKHandle
 * @brief A handle, describing the secret-key element that "matches" a part, of the form s^r(X^t).
 **/
class SKHandle {
  long powerOfS, powerOfX, secretKeyID;

public:
  friend class Ctxt;

  SKHandle(long newPowerOfS=0, long newPowerOfX=1, long newSecretKeyID=0)  
  {
    powerOfS = newPowerOfS;
    powerOfX = newPowerOfX;
    secretKeyID = newSecretKeyID;
  }

  //! @brief Set powerOfS=powerOfX=1
  void setBase(long newSecretKeyID=-1) 
  {
    powerOfS = 1;
    powerOfX = 1;
    if (newSecretKeyID >= 0) secretKeyID = newSecretKeyID;
  }
 
  //! @brief Is powerOfS==powerOfX==1?
  bool isBase(long ofKeyID=0) const
  {
    // If ofKeyID<0, only check that this is base of some key,
    // otherwise check that this is base of the given key
    return powerOfS == 1 && powerOfX == 1 && 
      (ofKeyID<0 || secretKeyID == ofKeyID);
  }

  //! @brief Set powerOfS=0, powerOfX=1
  void setOne(long newSecretKeyID=-1) 
  {
    powerOfS = 0;
    powerOfX = 1;
    if (newSecretKeyID >= 0) secretKeyID = newSecretKeyID;
  }

  //! @brief Is powerOfS==0?
  bool isOne() const
  {
    return powerOfS == 0;
  }

  bool operator==(const SKHandle& other) const
  {
    if (powerOfS==0 && other.powerOfS==0) return true;
    return powerOfS==other.powerOfS &&
	   powerOfX==other.powerOfX &&
           secretKeyID==other.secretKeyID;
  }

  bool operator!=(const SKHandle& other) const {return !(*this==other);}

  /* access functions */

  long getPowerOfS() const { return powerOfS; }
  long getPowerOfX() const { return powerOfX; }
  long getSecretKeyID() const { return secretKeyID; }

  /**
   * @brief Computes the "product" of two handles.
   *
   * The key-ID's and powers of X must match, else an error state arises,
   * which is represented using a key-ID of -1 and returning false. Also,
   * note that inputs may alias outputs.
   * 
   * To detremine if the resulting handle canbe re-liearized using
   * some key-switchingmatrices from the public key, use the method
   * pubKey.haveKeySWmatrix(handle,handle.secretKeyID), from the class
   * FHEPubKey in FHE.h
  */
  bool mul(const SKHandle& a, const SKHandle& b) 
  {
    // If either of the inputs is one, the output equals to the other input
    if (a.isOne()) {
      *this = b;
      return (b.secretKeyID >= 0);
    }
    if (b.isOne()) {
      *this = a;
      return (a.secretKeyID >= 0);
    }

    if (a.secretKeyID == -1 || b.secretKeyID == -1) {
      secretKeyID = -1; // -1 will be used to indicate an "error state"
      return false;
    }

    if (a.secretKeyID != b.secretKeyID) {
      secretKeyID = -1; 
      return false;
    }

    if (a.powerOfX != b.powerOfX) {
      secretKeyID = -1;
      return false;
    }

    secretKeyID = a.secretKeyID;
    powerOfX = a.powerOfX;
    powerOfS = a.powerOfS + b.powerOfS;
    return true;
  }

  friend istream& operator>>(istream& s, SKHandle& handle);
};
inline ostream& operator<<(ostream& s, const SKHandle& handle)
{
  return s << "[" << handle.getPowerOfS()    << " " << handle.getPowerOfX()
	   << " " << handle.getSecretKeyID() << "]";
}

/**
 * @class CtxtPart
 * @brief One entry in a ciphertext vector
 * 
 * A cipehrtext part consists of a polynomial (element of the ring R_Q)
 * and a handle to the corresponding secret-key polynomial.
 **/
class CtxtPart: public DoubleCRT {
public:
  //! @brief The handle is a public data member
  SKHandle skHandle; // The secret-key polynomial corresponding to this part

  bool operator==(const CtxtPart& other) const;
  bool operator!=(const CtxtPart& other) const {return !(*this==other);}

  // Copy constructor from the base class

  explicit
  CtxtPart(const FHEcontext& _context): DoubleCRT(_context)
  { skHandle.setOne(); }

  CtxtPart(const FHEcontext& _context, const IndexSet& s): 
    DoubleCRT(_context,s) { skHandle.setOne(); }

  CtxtPart(const FHEcontext& _context, const IndexSet& s, 
	   const SKHandle& otherHandle): 
    DoubleCRT(_context,s), skHandle(otherHandle) {}

  explicit
  CtxtPart(const DoubleCRT& other): DoubleCRT(other) { skHandle.setOne(); }

  CtxtPart(const DoubleCRT& other, const SKHandle& otherHandle): 
    DoubleCRT(other), skHandle(otherHandle) {}
};
istream& operator>>(istream& s, CtxtPart& p);
ostream& operator<<(ostream& s, const CtxtPart& p);

/**
 * @class Ctxt
 * @brief A Ctxt object holds a single cipehrtext
 *
 * The class Ctxt includes a vector<CtxtPart>: For a Ctxt c, c[i] is the i'th
 * ciphertext part, which can be used also as a DoubleCRT object (since
 * CtxtPart is derived from DoubleCRT). By convention, c[0], the first
 * CtxtPart object in the vector, has skHndl that points to 1 (i.e., it
 * is just added in upon decryption, without being multiplied by anything).
 * We maintain the invariance that all the parts of a ciphertext are defined
 * relative to the same set of primes.
 *
 * A ciphertext contains also pointers to the general parameters of this
 * FHE instance and the public key, and an estimate of the noise variance.
 * The noise variance is determined by the norm of the canonical embedding
 * of the noise polynomials, namely their evaluations in roots of the ring
 * polynomial (which are the complex primitive roots of unity). We consider
 * each such evaluation point as a random variable, and estimate the variances
 * of these variables. This estimate is heuristic, assuming that various
 * quantities "behave like independent random variables".
 * The variance is added on addition, multiplied on multiplications, remains
 * unchanged for automorphism, and is roughly scaled down by mod-switching
 * with some added factor, and similarly scaled up by key-switching with some
 * added factor. The noiseVar data member of the class keeps the esitmated
 * variance.
 **/
class Ctxt {
  friend class FHEPubKey;
  friend class FHESecKey;

  const FHEcontext& context; // points to the parameters of this FHE instance
  const FHEPubKey& pubKey;   // points to the public encryption key;
  vector<CtxtPart> parts;    // the ciphertexe parts
  IndexSet primeSet; // the primes relative to which the parts are defined
  long ptxtSpace;    // plaintext space for this ciphertext (either p or p^r)
  xdouble noiseVar;  // estimating the noise variance in this ciphertext

  // Create a tensor product of c1,c2. It is assumed that *this,c1,c2
  // are defined relative to the same set of primes and plaintext space,
  // and that *this DOES NOT point to the same object as c1,c2
  void tensorProduct(const Ctxt& c1, const Ctxt& c2);

  // Add/subtract a ciphertext part to/from a ciphertext. These are private
  // methods, they cannot update the noiseVar estimate so they must be called
  // from a procedure that will eventually update that estimate.
  Ctxt& operator-=(const CtxtPart& part) { subPart(part); return *this; }
  Ctxt& operator+=(const CtxtPart& part) { addPart(part); return *this; }

  // Procedureal versions with additional parameter
  void subPart(const CtxtPart& part, bool matchPrimeSet=false)
  { subPart(part, part.skHandle, matchPrimeSet); }
  void addPart(const CtxtPart& part, bool matchPrimeSet=false)
  { addPart(part, part.skHandle, matchPrimeSet); }

  void subPart(const DoubleCRT& part, 
	       const SKHandle& handle, bool matchPrimeSet=false) 
  { addPart(part, handle, matchPrimeSet, true); }
  void addPart(const DoubleCRT& part, const SKHandle& handle, 
	       bool matchPrimeSet=false, bool negative=false);

  // Takes as arguments a ciphertext-part p relative to s' and a key-switching
  // matrix W = W[s'->s], use W to switch p relative to (1,s), and add the
  // result to *this.
  void keySwitchPart(const CtxtPart& p, const KeySwitch& W);

  long getPartIndexByHandle(const SKHandle& hanle) const {
    for (size_t i=0; i<parts.size(); i++) 
      if (parts[i].skHandle==hanle) return i;
    return -1;
  }

public:
  Ctxt(const FHEPubKey& newPubKey, long newPtxtSpace=2); // constructor
  Ctxt& operator=(const Ctxt& other);  // assignment

  bool operator==(const Ctxt& other) const { return equalsTo(other); }
  bool operator!=(const Ctxt& other) const { return !equalsTo(other); }

  // a procedural variant with an additional parameter
  bool equalsTo(const Ctxt& other, bool comparePkeys=true) const;

  void clear() { // set as an empty ciphertext
    primeSet.clear();
    parts.clear();
    noiseVar = to_xdouble(0.0);
  }

  // Encryption and decryption are done by the friends FHE[Pub|Sec]Key

  //! @name Ciphertext arithmetic
  ///@{
  void negate();

 // Add/subtract aonther ciphertext
  Ctxt& operator+=(const Ctxt& other) { addCtxt(other); return *this; }
  Ctxt& operator-=(const Ctxt& other) { addCtxt(other,true); return *this; }
  void addCtxt(const Ctxt& other, bool negative=false);

  Ctxt& operator*=(const Ctxt& other); // Multiply by aonther ciphertext
  void automorph(long k); // Apply automorphism F(X) -> F(X^k) (gcd(k,m)=1)
  Ctxt& operator>>=(long k) { automorph(k); return *this; }

  //! @brief automorphism with re-lienarization
  void smartAutomorph(long k);
  // Apply F(X)->F(X^k) followed by re-liearization. The automorphism is
  // possibly evaluated via a sequence of steps, to ensure that we can
  // re-linearize the result of every step.

  // Add a constant polynomial. If the size is not given, we use
  // phi(m)*ptxtSpace^2 as the default value.
  void addConstant(const DoubleCRT& dcrt, double size=0.0);
  void addConstant(const ZZX& poly, double size=0.0)
  { addConstant(DoubleCRT(poly,context,primeSet),size); }

  // Multiply-by-constant. If the size is not given, we use
  // phi(m)*ptxtSpace^2 as the default value.
  void multByConstant(const DoubleCRT& dcrt, double size=0.0);
  void multByConstant(const ZZX& poly, double size=0.0)
  { multByConstant(DoubleCRT(poly,context,primeSet),size); }

  // Higher-level multiply routines
  void multiplyBy(const Ctxt& other);
  void multiplyBy2(const Ctxt& other1, const Ctxt& other2);
  void square() { multiplyBy(*this); }
  void cube() { multiplyBy2(*this, *this); }
  ///@}

  //! @name Ciphertext maintenance
  ///@{
  void reLinearize(long keyIdx=0);
          // key-switch to (1,s_i), s_i is the base key with index keyIdx

  //! @brief Estimate the added noise variance
  xdouble modSwitchAddedNoiseVar() const;

  //! @brief Modulus-switching up (to a larger modulus).
  //! Must have primeSet <= s, and s must contain
  //! either all the special primes or none of them.
  void modUpToSet(const IndexSet &s);

  //! @brief Modulus-switching down (to a smaller modulus).
  //! mod-switch down to primeSet \intersect s, after this call we have
  //! primeSet<=s. s must contain either all special primes or none of them
  void modDownToSet(const IndexSet &s);

  //! @brief Fidn the "natural level" of a cipehrtext.
  //! Find the highest IndexSet so that mod-switching down to that set results
  //! in the dominant noise term being the additive term due to rounding
  void findBaseSet(IndexSet& s) const;
  ///@}

  //! @name Utility methods
  ///@{

  //! @brief A canonical ciphertext has handles pointing to (1,s)
  bool inCanonicalForm(long keyID=0) const {
    return (parts.size()==2 && 
	    parts[0].skHandle.isOne() && parts[1].skHandle.isBase(keyID));
  }

  //! @brief Would this ciphertext be decrypted without errors?
  bool isCorrect() const {
    ZZ q = context.productOfPrimes(primeSet);
    return (to_xdouble(q) > sqrt(noiseVar)*2);
  }
  const FHEcontext& getContext() const { return context; }
  const FHEPubKey& getPubKey() const   { return pubKey; }
  const IndexSet& getPrimeSet() const  { return primeSet; }
  const xdouble& getNoiseVar() const   { return noiseVar; }
  const long getPtxtSpace() const      { return ptxtSpace;}
  const long getKeyID() const;

  //! @brief How many primes in the "base-set" for that ciphertext
  const long getLevel() const { 
    IndexSet s;
    findBaseSet(s);
    return card(s);
  }

  //! @brief Returns log(noise-variance)/2 - log(q)
  double log_of_ratio() const
  {return (log(getNoiseVar())/2 - context.logOfProduct(getPrimeSet()));}
  ///@}
  friend istream& operator>>(istream& str, Ctxt& ctxt);
  friend ostream& operator<<(ostream& str, const Ctxt& ctxt);
};

inline IndexSet baseSetOf(const Ctxt& c) { 
  IndexSet s; c.findBaseSet(s); return s; 
}


// print to cerr some info about ciphertext
void CheckCtxt(const Ctxt& c, const char* label);

#endif // ifndef _Ctxt_H_
