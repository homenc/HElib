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
#ifndef _EncryptedArray_H_
#define _EncryptedArray_H_
/**
 * @file EncryptedArray.h
 * @brief Data-movement operations on encrypted arrays of slots
 */

#include <NTL/Lazy.h>
#include <NTL/pair.h>
#include <NTL/SmartPtr.h>
#include "FHE.h"
#include "timing.h"


// DIRT: we're using undocumented NTL interfaces here
//   also...this probably should be defined in NTL, anyway....
#define FHE_MORE_UNWRAPARGS(n) NTL_SEPARATOR_##n NTL_REPEATER_##n(NTL_UNWRAPARG)




// these are used to implement NewPlaintextArray stuff routines

#define PA_BOILER \
    const PAlgebraModDerived<type>& tab = ea.getTab(); \
    const RX& G = ea.getG(); \
    long n = ea.size(); \
    long d = ea.getDegree(); \
    vector<RX>& data = pa.getData<type>(); \
    RBak bak; bak.save(); tab.restoreContext(); \


#define CPA_BOILER \
    const PAlgebraModDerived<type>& tab = ea.getTab(); \
    const RX& G = ea.getG(); \
    long n = ea.size(); \
    long d = ea.getDegree(); \
    const vector<RX>& data = pa.getData<type>(); \
    RBak bak; bak.save(); tab.restoreContext(); \






class NewPlaintextArray; // forward reference
class EncryptedArray; // forward reference

/**
 * @class EncryptedArrayBase
 * @brief virtual class for data-movement operations on arrays of slots
 *
 * An object ea of type EncryptedArray stores information about an
 * FHEcontext context, and a monic polynomial G.  If context defines
 * parameters m, p, and r, then ea is a helper abject
 * that supports encoding/decoding and encryption/decryption
 * of vectors of plaintext slots over the ring (Z/(p^r)[X])/(G). 
 *
 * The polynomial G should be irreducble over Z/(p^r) (this is not checked).
 * The degree of G should divide the multiplicative order of p modulo m
 * (this is checked). Currently, the following restriction is imposed:
 *
 * either r == 1 or deg(G) == 1 or G == factors[0].
 * 
 * ea stores objects in the polynomial ring Z/(p^r)[X].
 * 
 * Just as for the class PAlegebraMod, if p == 2 and r == 1, then these
 * polynomials are represented as GF2X's, and otherwise as zz_pX's.
 * Thus, the types of these objects are not determined until run time.
 * As such, we need to use a class heirarchy, which mirrors that of
 * PAlgebraMod, as follows.
 * 
 * EncryptedArrayBase is a virtual class
 * 
 * EncryptedArrayDerived<type> is a derived template class, where
 * type is either PA_GF2 or PA_zz_p.
 *
 * The class EncryptedArray is a simple wrapper around a smart pointer to
 * an EncryptedArrayBase object: copying an EncryptedArray object results
 * is a "deep copy" of the underlying object of the derived class.
****************************************************************/

class EncryptedArrayBase {  // purely abstract interface 
public:
  virtual ~EncryptedArrayBase() {}

  virtual EncryptedArrayBase* clone() const = 0;
  // makes this usable with cloned_ptr

  virtual PA_tag getTag() const = 0;

  virtual const FHEcontext& getContext() const = 0;
  virtual const PAlgebra& getPAlgebra() const = 0;
  virtual const long getDegree() const = 0;

  //! @brief Right rotation as a linear array.
  //! E.g., rotating ctxt=Enc(1 2 3 ... n) by k=1 gives Enc(n 1 2 ... n-1)
  virtual void rotate(Ctxt& ctxt, long k) const = 0; 

  //! @brief Non-cyclic right shift with zero fill
  //! E.g., shifting ctxt=Enc(1 2 3 ... n) by k=1 gives Enc(0 1  2... n-1)
  virtual void shift(Ctxt& ctxt, long k) const = 0;

  //! @brief right-rotate k positions along the i'th dimension
  //! @param dc means "don't care", which means that the caller guarantees
  //! that only zero elements rotate off the end -- this allows for some
  //! optimizations that would not otherwise be possible
  virtual void rotate1D(Ctxt& ctxt, long i, long k, bool dc=false) const = 0; 

  //! @brief Right shift k positions along the i'th dimension with zero fill
  virtual void shift1D(Ctxt& ctxt, long i, long k) const = 0; 




  ///@{
  //! @name Encoding/decoding methods
  // encode/decode arrays into plaintext polynomials
  virtual void encode(zzX& ptxt, const vector< long >& array) const = 0;
  virtual void encode(zzX& ptxt, const vector< zzX >& array) const = 0;
  virtual void encode(zzX& ptxt, const NewPlaintextArray& array) const = 0;

  virtual void encode(ZZX& ptxt, const vector< long >& array) const = 0;
  virtual void encode(ZZX& ptxt, const vector< ZZX >& array) const = 0;
  virtual void encode(ZZX& ptxt, const NewPlaintextArray& array) const = 0;

  virtual void decode(vector< long  >& array, const ZZX& ptxt) const = 0;
  virtual void decode(vector< ZZX  >& array, const ZZX& ptxt) const = 0;
  virtual void decode(NewPlaintextArray& array, const ZZX& ptxt) const = 0;

  virtual void random(vector< long >& array) const = 0;
  virtual void random(vector< ZZX >& array) const = 0;

  // FIXME: Inefficient implementation, calls usual decode and returns one slot
  long decode1Slot(const ZZX& ptxt, long i) const
  { vector< long > v; decode(v, ptxt); return v.at(i); }
  void decode1Slot(ZZX& slot, const ZZX& ptxt, long i) const
  { vector< ZZX > v; decode(v, ptxt); slot=v.at(i); }

  //! @brief Encodes a vector with 1 at position i and 0 everywhere else
  virtual void encodeUnitSelector(ZZX& ptxt, long i) const = 0;
  ///@}

  ///@{
  //! @name Encoding+encryption/decryption+decoding
  virtual void encrypt(Ctxt& ctxt, const FHEPubKey& pKey, const vector< long >& ptxt) const = 0;
  virtual void encrypt(Ctxt& ctxt, const FHEPubKey& pKey, const vector< ZZX >& ptxt) const = 0;
  virtual void encrypt(Ctxt& ctxt, const FHEPubKey& pKey, const NewPlaintextArray& ptxt) const = 0;
  virtual void decrypt(const Ctxt& ctxt, const FHESecKey& sKey, vector< long >& ptxt) const = 0;
  virtual void decrypt(const Ctxt& ctxt, const FHESecKey& sKey, vector< ZZX >& ptxt) const = 0;
  virtual void decrypt(const Ctxt& ctxt, const FHESecKey& sKey, NewPlaintextArray& ptxt) const = 0;

  // Also secret-key encryption, for convenience
  virtual void skEncrypt(Ctxt& ctxt, const FHESecKey& sKey, const vector< long >& ptxt, long skIdx=0) const = 0;
  virtual void skEncrypt(Ctxt& ctxt, const FHESecKey& sKey, const vector< ZZX >& ptxt, long skIdx=0) const = 0;
  virtual void skEncrypt(Ctxt& ctxt, const FHESecKey& sKey, const NewPlaintextArray& ptxt, long skIdx=0) const = 0;

  // FIXME: Inefficient implementation, calls usual decrypt and returns one slot
  long decrypt1Slot(const Ctxt& ctxt, const FHESecKey& sKey, long i) const
  { vector< long > v; decrypt(ctxt, sKey, v); return v.at(i); }
  void decrypt1Slot(ZZX& slot, const Ctxt& ctxt, const FHESecKey& sKey, long i) const
  { vector< ZZX > v; decrypt(ctxt, sKey, v); slot = v.at(i); }
  ///@}

  //@{
  //! MUX: ctxt1 = ctxt1*selector + ctxt2*(1-selector)
  virtual void select(Ctxt& ctxt1, const Ctxt& ctxt2, const vector< long >& selector) const = 0;
  virtual void select(Ctxt& ctxt1, const Ctxt& ctxt2, const vector< ZZX >& selector) const = 0;
  virtual void select(Ctxt& ctxt1, const Ctxt& ctxt2, const NewPlaintextArray& selector) const = 0;
  //@}

  //! @brief Linearized polynomials.
  //! L describes a linear map M by describing its action on the standard
  //! power basis: M(x^j mod G) = (L[j] mod G), for j = 0..d-1.  
  //! The result is a coefficient vector C for the linearized polynomial
  //! representing M: a polynoamial h in Z/(p^r)[X] of degree < d is sent to
  //! \f[
  //!  M(h(X) \bmod G)= \sum_{i=0}^{d-1}(C[j] \cdot h(X^{p^j}))\bmod G).
  //! \f]
  virtual void buildLinPolyCoeffs(vector<ZZX>& C, const vector<ZZX>& L) const=0;

  // restore contexts mod p and mod G
  virtual void restoreContext() const = 0;
  virtual void restoreContextForG() const = 0;

  /* some non-virtual convenience functions */

  //! @brief Total size (# of slots) of hypercube
  long size() const { 
    return getPAlgebra().getNSlots(); 
  } 

  //! @brief Number of dimensions of hypercube
  long dimension() const { 
    return getPAlgebra().numOfGens(); 
  }

  //! @brief Size of given dimension
  long sizeOfDimension(long i) const {
    return getPAlgebra().OrderOf(i);
  }

  //! @brief Is rotations in given dimension a "native" operation?
  bool nativeDimension(long i) const {
    return getPAlgebra().SameOrd(i);
  }

  //! @brief returns coordinate of index k along the i'th dimension
  long coordinate(long i, long k) const {
    return getPAlgebra().coordinate(i, k); 
  }

  //! @brief adds offset to index k in the i'th dimension
  long addCoord(long i, long k, long offset) const {
    return getPAlgebra().addCoord(i, k, offset);
  }


  //! @brief rotate an array by offset in the i'th dimension
  //! (output should not alias input)
  template<class U> void rotate1D(vector<U>& out, const vector<U>& in,
                                  long i, long offset) const {
    assert(lsize(in) == size());
    out.resize(in.size());
    for (long j = 0; j < size(); j++)
      out[addCoord(i, j, offset)] = in[j]; 
  }
};

/**
 * @class EncryptedArrayDerived
 * @brief Derived concrete implementation of EncryptedArrayBase
 */
template<class type> class EncryptedArrayDerived : public EncryptedArrayBase {
public:
  PA_INJECT(type)

private:
  const FHEcontext& context;
  const PAlgebraModDerived<type>& tab;

  // FIXME: all of these types should be copyable
  // out of context, but NTL 8.0 still does not copy
  // matrix copies out of context correctly, as it
  // relies on plain SetLength...I need to fix this 
  // in NTL
  //
  MappingData<type> mappingData; // MappingData is defined in PAlgebra.h

  Lazy< Mat<RE> > linPolyMatrix;

  Lazy< Pair< Mat<R>, Mat<R> > > normalBasisMatrices;
  // a is the matrix, b is its inverse

public:
  explicit
  EncryptedArrayDerived(const FHEcontext& _context, const RX& _G,
			const PAlgebraMod& _tab);

  EncryptedArrayDerived(const EncryptedArrayDerived& other) // copy constructor
    : context(other.context), tab(other.tab)
  {
    RBak bak; bak.save(); tab.restoreContext();
    REBak ebak; ebak.save(); other.mappingData.restoreContextForG();
    mappingData = other.mappingData;
    linPolyMatrix = other.linPolyMatrix;
    normalBasisMatrices = other.normalBasisMatrices;
  }

  EncryptedArrayDerived& operator=(const EncryptedArrayDerived& other) // assignment
  {
    if (this == &other) return *this;
    assert(&context == &other.context);
    assert(&tab == &other.tab);

    RBak bak; bak.save(); tab.restoreContext();
    mappingData = other.mappingData;
    linPolyMatrix = other.linPolyMatrix;
    normalBasisMatrices = other.normalBasisMatrices;
    return *this;
  }

  virtual EncryptedArrayBase* clone() const { return new EncryptedArrayDerived(*this); }

  virtual PA_tag getTag() const { return tag; }
  // tag is defined in PA_INJECT, see PAlgebra.h

  template<template <class> class T, class... Args>
  void dispatch(Args&&... args) const
  {
    T<type>::apply(*this, std::forward<Args>(args)...);
  }


  const RX& getG() const { return mappingData.getG(); }

  const Mat<R>& getNormalBasisMatrix() const {
    if (!normalBasisMatrices.built()) initNormalBasisMatrix(); 
    return normalBasisMatrices->a;
  }

  const Mat<R>& getNormalBasisMatrixInverse() const {
    if (!normalBasisMatrices.built()) initNormalBasisMatrix(); 
    return normalBasisMatrices->b;
  }

  void initNormalBasisMatrix() const;

  virtual void restoreContext() const { tab.restoreContext(); }
  virtual void restoreContextForG() const { mappingData.restoreContextForG(); }


  virtual const FHEcontext& getContext() const override { return context; }
  virtual const PAlgebra& getPAlgebra() const override { return tab.getZMStar(); }
  virtual const long getDegree() const { return mappingData.getDegG(); }
  const PAlgebraModDerived<type>& getTab() const { return tab; }

  virtual void rotate(Ctxt& ctxt, long k) const;
  virtual void shift(Ctxt& ctxt, long k) const;
  virtual void rotate1D(Ctxt& ctxt, long i, long k, bool dc=false) const;
  template<class U> void // avoid this being "hidden" by other rotate1D's
    rotate1D(vector<U>& out, const vector<U>& in, long i, long offset) const
    { EncryptedArrayBase::rotate1D(out, in, i, offset); }
  virtual void shift1D(Ctxt& ctxt, long i, long k) const;

  virtual void encode(ZZX& ptxt, const vector< long >& array) const
    { genericEncode(ptxt, array);  }

  virtual void encode(zzX& ptxt, const vector< long >& array) const
    { genericEncode(ptxt, array);  }

  virtual void encode(ZZX& ptxt, const vector< ZZX >& array) const
    {  genericEncode(ptxt, array); }

  virtual void encode(zzX& ptxt, const vector< zzX >& array) const
    {  genericEncode(ptxt, array); }

  virtual void encode(ZZX& ptxt, const NewPlaintextArray& array) const;
  virtual void encode(zzX& ptxt, const NewPlaintextArray& array) const;

  virtual void encodeUnitSelector(zzX& ptxt, long i) const;
  virtual void encodeUnitSelector(ZZX& ptxt, long i) const;

  virtual void decode(vector< long  >& array, const ZZX& ptxt) const
    { genericDecode(array, ptxt); }

  virtual void decode(vector< ZZX  >& array, const ZZX& ptxt) const
    { genericDecode(array, ptxt); }

  virtual void decode(NewPlaintextArray& array, const ZZX& ptxt) const;
  virtual void decode(NewPlaintextArray& array, const zzX& ptxt) const;

  virtual void random(vector< long  >& array) const
    { genericRandom(array); } // choose at random and convert to vector<long>

  virtual void random(vector< ZZX  >& array) const
    { genericRandom(array); } // choose at random and convert to vector<ZZX>

  virtual void encrypt(Ctxt& ctxt, const FHEPubKey& pKey, const vector< long >& ptxt) const
    { genericEncrypt(ctxt, pKey, ptxt); }

  virtual void encrypt(Ctxt& ctxt, const FHEPubKey& pKey, const vector< ZZX >& ptxt) const
    { genericEncrypt(ctxt, pKey, ptxt); }

  virtual void encrypt(Ctxt& ctxt, const FHEPubKey& pKey, const NewPlaintextArray& ptxt) const
    { genericEncrypt(ctxt, pKey, ptxt); }

  virtual void decrypt(const Ctxt& ctxt, const FHESecKey& sKey, vector< long >& ptxt) const
    { genericDecrypt(ctxt, sKey, ptxt);
      if (ctxt.getPtxtSpace()<tab.getPPowR()) {
	for (long i=0; i<(long)ptxt.size(); i++)
	  ptxt[i] %= ctxt.getPtxtSpace();
      }
    }

  virtual void decrypt(const Ctxt& ctxt, const FHESecKey& sKey, vector< ZZX >& ptxt) const
    { genericDecrypt(ctxt, sKey, ptxt);
      if (ctxt.getPtxtSpace()<tab.getPPowR()) {
	for (long i=0; i<(long)ptxt.size(); i++)
	  PolyRed(ptxt[i], ctxt.getPtxtSpace(),/*abs=*/true);
      }
    }


  virtual void decrypt(const Ctxt& ctxt, const FHESecKey& sKey, NewPlaintextArray& ptxt) const
  { genericDecrypt(ctxt, sKey, ptxt); 
    // FIXME: Redudc mod the ciphertext plaintext space as above
    }

  virtual void skEncrypt(Ctxt& ctxt, const FHESecKey& sKey, const vector< long >& ptxt, long skIdx=0) const
    { genericSkEncrypt(ctxt, sKey, ptxt, skIdx); }

  virtual void skEncrypt(Ctxt& ctxt, const FHESecKey& sKey, const vector< ZZX >& ptxt, long skIdx=0) const
    { genericSkEncrypt(ctxt, sKey, ptxt, skIdx); }


  virtual void skEncrypt(Ctxt& ctxt, const FHESecKey& sKey, const NewPlaintextArray& ptxt, long skIdx=0) const
    { genericSkEncrypt(ctxt, sKey, ptxt, skIdx); }


  virtual void select(Ctxt& ctxt1, const Ctxt& ctxt2, const vector< long >& selector) const
    { genericSelect(ctxt1, ctxt2, selector); }

  virtual void select(Ctxt& ctxt1, const Ctxt& ctxt2, const vector< ZZX >& selector) const
    { genericSelect(ctxt1, ctxt2, selector); }


  virtual void select(Ctxt& ctxt1, const Ctxt& ctxt2, const NewPlaintextArray& selector) const
    { genericSelect(ctxt1, ctxt2, selector); }

  virtual void buildLinPolyCoeffs(vector<ZZX>& C, const vector<ZZX>& L) const;

  /* the following are specialized methods, used to work over extension fields...they assume 
     the modulus context is already set
   */

  void encode(zzX& ptxt, const vector< RX >& array) const;
  void decode(vector< RX  >& array, const zzX& ptxt) const;

  void encode(ZZX& ptxt, const vector< RX >& array) const;
  void decode(vector< RX  >& array, const ZZX& ptxt) const;

  void encode(RX& ptxt, const vector< RX >& array) const;
  void decode(vector< RX  >& array, const RX& ptxt) const;

  // Choose random polynomial of the right degree, coeffs in GF2 or zz_p
  void random(vector< RX  >& array) const
  { 
    array.resize(size()); 
    for (long i=0; i<size(); i++) NTL::random(array[i], getDegree());
  }

  void encrypt(Ctxt& ctxt, const FHEPubKey& pKey, const vector< RX >& ptxt) const
    { genericEncrypt(ctxt, pKey, ptxt); }

  void decrypt(const Ctxt& ctxt, const FHESecKey& sKey, vector< RX >& ptxt) const
    { genericDecrypt(ctxt, sKey, ptxt); }

  void skEncrypt(Ctxt& ctxt, const FHESecKey& sKey, const vector< RX >& ptxt, long skIdx=0) const
    { genericSkEncrypt(ctxt, sKey, ptxt, skIdx); }

  virtual void buildLinPolyCoeffs(vector<RX>& C, const vector<RX>& L) const;


private:

  /* helper template methods, to avoid repetitive code */

  template<class T> 
  void genericEncode(ZZX& ptxt, const T& array) const
  {
    RBak bak; bak.save(); tab.restoreContext();

    vector< RX > array1;
    convert(array1, array);
    encode(ptxt, array1);
  }

  template<class T> 
  void genericEncode(zzX& ptxt, const T& array) const
  {
    RBak bak; bak.save(); tab.restoreContext();

    vector< RX > array1;
    convert(array1, array);
    encode(ptxt, array1);
  }

  template<class T>
  void genericDecode(T& array, const ZZX& ptxt) const
  {
    RBak bak; bak.save(); tab.restoreContext();

    vector< RX > array1;
    decode(array1, ptxt);
    convert(array, array1);
  }

  template<class T>
  void genericRandom(T& array) const // T is vector<long> or vector<ZZX>
  {
    RBak bak; bak.save(); tab.restoreContext(); // backup NTL modulus

    vector< RX > array1;    // RX is GF2X or zz_pX
    random(array1);         // choose random coefficients from GF2/zz_p
    convert(array, array1); // convert to type T (see NumbTh.h)
  }

  template<class T>
  void genericEncrypt(Ctxt& ctxt, const FHEPubKey& pKey, 
                      const T& array) const
  {
    assert(&context == &ctxt.getContext());
    ZZX pp;
    encode(pp, array); // Convert the array of slots into a plaintext polynomial
    pKey.Encrypt(ctxt, pp, tab.getPPowR()); // encrypt the plaintext polynomial
  }

  template<class T>
  void genericDecrypt(const Ctxt& ctxt, const FHESecKey& sKey, 
                      T& array) const
  {
    assert(&context == &ctxt.getContext());
    ZZX pp;
    sKey.Decrypt(pp, ctxt);
    decode(array, pp);
  }

  template<class T>
  void genericSkEncrypt(Ctxt& ctxt, const FHESecKey& sKey, 
                      const T& array, long skIdx=0) const
  {
    assert(&context == &ctxt.getContext());
    ZZX pp;
    encode(pp, array); // Convert the array of slots into a plaintext polynomial
    sKey.Encrypt(ctxt, pp, tab.getPPowR(), skIdx); // encrypt the plaintext polynomial
  }


  template<class T>
  void genericSelect(Ctxt& ctxt1, const Ctxt& ctxt2,
			    const T& selector) const
  {
    if (&ctxt1 == &ctxt2) return; // nothing to do

    assert(&context == &ctxt1.getContext() && &context == &ctxt2.getContext());
    ZZX poly; 
    encode(poly,selector);                             // encode as polynomial
    DoubleCRT dcrt(poly, context, ctxt1.getPrimeSet());// convert to DoubleCRT

    ctxt1.multByConstant(dcrt); // keep only the slots with 1's

    dcrt -= 1;      // convert 1 => 0, 0 => -1
    Ctxt tmp=ctxt2; // a copy of ctxt2
    tmp.multByConstant(dcrt);// keep (but negate) only the slots with 0's
    ctxt1 -= tmp;
  }
};


// plaintextAutomorph: Compute b(X) = a(X^k) mod Phi_m(X).
template <class RX, class RXModulus>
void plaintextAutomorph(RX& bb, const RX& a, long k, long m, const RXModulus& PhimX)
{
  // compute b(X) = a(X^k) mod (X^m-1)
  if (k == 1 || deg(a) <= 0) {
    bb = a;
    return;
  }

  RX b;
  b.SetLength(m);
  mulmod_precon_t precon = PrepMulModPrecon(k, m);
  for (long j = 0; j <= deg(a); j++) 
    b[MulModPrecon(j, k, m, precon)] = a[j]; // b[j*k mod m] = a[j]
  b.normalize();

  rem(bb, b, PhimX); // reduce modulo the m'th cyclotomic
}

// same as above, but k = g_i^j mod m.
// also works with i == ea.getPalgebra().numOfGens(),
// which means Frobenius

template<class RX, class type>
void plaintextAutomorph(RX& b, const RX& a, long i, long j, 
                        const EncryptedArrayDerived<type>& ea)
{
  const PAlgebra& zMStar = ea.getPAlgebra();
  const auto& F = ea.getTab().getPhimXMod();
  long k = zMStar.genToPow(i, j);
  long m = zMStar.getM();
  plaintextAutomorph(b, a, k, m, F);
}


//! @brief A "factory" for building EncryptedArrays
EncryptedArrayBase* buildEncryptedArray(const FHEcontext& context,
					const ZZX& G, const PAlgebraMod& alMod);



//! @class EncryptedArray
//! @brief A simple wrapper for a smart pointer to an EncryptedArrayBase.
//! This is the interface that higher-level code should use
class EncryptedArray {
private:
  const PAlgebraMod& alMod;
  cloned_ptr<EncryptedArrayBase> rep;

public:

  //! constructor: G defaults to the monomial X, PAlgebraMod from context
  EncryptedArray(const FHEcontext& context, const ZZX& G = ZZX(1, 1))
    : alMod(context.alMod), rep(buildEncryptedArray(context,G,context.alMod))
  { }
  //! constructor: G defaults to F0, PAlgebraMod explicitly given
  EncryptedArray(const FHEcontext& context, const PAlgebraMod& _alMod)
    : alMod(_alMod), 
      rep(buildEncryptedArray(context, _alMod.getFactorsOverZZ()[0], _alMod))
  { }

  // copy constructor: 

  EncryptedArray& operator=(const EncryptedArray& other) {
    if (this == &other) return *this;
    assert(&alMod== &other.alMod);
    rep = other.rep;
    return *this;
  }

  //! @brief downcast operator
  //! example: const EncryptedArrayDerived<PA_GF2>& rep = ea.getDerived(PA_GF2());
  template<class type> 
  const EncryptedArrayDerived<type>& getDerived(type) const
  { return dynamic_cast< const EncryptedArrayDerived<type>& >( *rep ); }


  ///@{
  //! @name Direct access to EncryptedArrayBase methods

  PA_tag getTag() const { return rep->getTag(); }

  template<template <class> class T, class... Args>
  void dispatch(Args&&... args) const
  {
    switch (getTag()) {
      case PA_GF2_tag: {
        const EncryptedArrayDerived<PA_GF2> *p = 
          static_cast< const EncryptedArrayDerived<PA_GF2> *>(rep.get_ptr());
        p->dispatch<T>(std::forward<Args>(args)...);
        break;
      }
      case PA_zz_p_tag: {
        const EncryptedArrayDerived<PA_zz_p> *p = 
          static_cast< const EncryptedArrayDerived<PA_zz_p> *>(rep.get_ptr());
        p->dispatch<T>(std::forward<Args>(args)...);
        break;
      }
      default: TerminalError("bad tag"); 
    }
  }



  const FHEcontext& getContext() const { return rep->getContext(); }
  const PAlgebraMod& getAlMod() const { return alMod; }
  const PAlgebra& getPAlgebra() const { return rep->getPAlgebra(); }
  const long getDegree() const { return rep->getDegree(); }
  void rotate(Ctxt& ctxt, long k) const { rep->rotate(ctxt, k); }
  void shift(Ctxt& ctxt, long k) const { rep->shift(ctxt, k); }
  void rotate1D(Ctxt& ctxt, long i, long k, bool dc=false) const { rep->rotate1D(ctxt, i, k, dc); }
  void shift1D(Ctxt& ctxt, long i, long k) const { rep->shift1D(ctxt, i, k); }

  void encode(ZZX& ptxt, const vector< long >& array) const 
    { rep->encode(ptxt, array); }
  void encode(ZZX& ptxt, const vector< ZZX >& array) const 
    { rep->encode(ptxt, array); }
  void encode(ZZX& ptxt, const NewPlaintextArray& array) const 
    { rep->encode(ptxt, array); }

  void encode(zzX& ptxt, const vector< long >& array) const 
    { rep->encode(ptxt, array); }
  void encode(zzX& ptxt, const vector< zzX >& array) const 
    { rep->encode(ptxt, array); }
  void encode(zzX& ptxt, const NewPlaintextArray& array) const 
    { rep->encode(ptxt, array); }

  void encodeUnitSelector(ZZX& ptxt, long i) const
    { rep->encodeUnitSelector(ptxt, i); }

  void decode(vector< long  >& array, const ZZX& ptxt) const 
    { rep->decode(array, ptxt); }
  void decode(vector< ZZX  >& array, const ZZX& ptxt) const 
    { rep->decode(array, ptxt); }
  void decode(NewPlaintextArray& array, const ZZX& ptxt) const 
    { rep->decode(array, ptxt); }

  void random(vector< long  >& array) const
    { rep->random(array); }
  void random(vector< ZZX  >& array) const
    { rep->random(array); }

  void encrypt(Ctxt& ctxt, const FHEPubKey& pKey, const vector< long >& ptxt) const 
    { rep->encrypt(ctxt, pKey, ptxt); }
  void encrypt(Ctxt& ctxt, const FHEPubKey& pKey, const vector< ZZX >& ptxt) const 
    { rep->encrypt(ctxt, pKey, ptxt); }
  void encrypt(Ctxt& ctxt, const FHEPubKey& pKey, const NewPlaintextArray& ptxt) const 
    { rep->encrypt(ctxt, pKey, ptxt); }


  void decrypt(const Ctxt& ctxt, const FHESecKey& sKey, vector< long >& ptxt) const 
    { rep->decrypt(ctxt, sKey, ptxt); }
  void decrypt(const Ctxt& ctxt, const FHESecKey& sKey, vector< ZZX >& ptxt) const 
    { rep->decrypt(ctxt, sKey, ptxt); }
  void decrypt(const Ctxt& ctxt, const FHESecKey& sKey, NewPlaintextArray& ptxt) const
    { rep->decrypt(ctxt, sKey, ptxt); }


  void skEncrypt(Ctxt& ctxt, const FHESecKey& sKey, const vector< long >& ptxt, long skIdx=0) const 
    { rep->skEncrypt(ctxt, sKey, ptxt, skIdx); }
  void skEncrypt(Ctxt& ctxt, const FHESecKey& sKey, const vector< ZZX >& ptxt, long skIdx=0) const 
    { rep->skEncrypt(ctxt, sKey, ptxt, skIdx); }
  void skEncrypt(Ctxt& ctxt, const FHESecKey& sKey, const NewPlaintextArray& ptxt, long skIdx=0) const 
    { rep->skEncrypt(ctxt, sKey, ptxt, skIdx); }


  void select(Ctxt& ctxt1, const Ctxt& ctxt2, const vector< long >& selector) const 
    { rep->select(ctxt1, ctxt2, selector); }
  void select(Ctxt& ctxt1, const Ctxt& ctxt2, const vector< ZZX >& selector) const 
    { rep->select(ctxt1, ctxt2, selector); }
  void select(Ctxt& ctxt1, const Ctxt& ctxt2, const NewPlaintextArray& selector) const
    { rep->select(ctxt1, ctxt2, selector); }

  void buildLinPolyCoeffs(vector<ZZX>& C, const vector<ZZX>& L) const
    { rep->buildLinPolyCoeffs(C, L); }

  void restoreContext() const { rep->restoreContext(); }
  void restoreContextForG() const { rep->restoreContextForG(); }

  long size() const { return rep->size(); } 
  long dimension() const { return rep->dimension(); }
  long sizeOfDimension(long i) const { return rep->sizeOfDimension(i); }
  long nativeDimension(long i) const {return rep->nativeDimension(i); }
  long coordinate(long i, long k) const { return rep->coordinate(i, k); }
  long addCoord(long i, long k, long offset) const { return rep->addCoord(i, k, offset); }


  //! @brief rotate an array by offset in the i'th dimension
  //! (output should not alias input)
  template<class U> void rotate1D(vector<U>& out, const vector<U>& in,
                                  long i, long offset) const {
    rep->rotate1D(out, in, i, offset);
  }
  ///@}
};



// NewPlaintaxtArray


class NewPlaintextArrayBase { // purely abstract interface
public:
  virtual ~NewPlaintextArrayBase() {}
  virtual void print(ostream& s) const = 0;

};


template<class type> class NewPlaintextArrayDerived : public NewPlaintextArrayBase {
public:
  PA_INJECT(type)

  vector< RX > data;

  virtual void print(ostream& s) const { s << data; }

};


class NewPlaintextArray {
private:

  CloneablePtr<NewPlaintextArrayBase> rep;

  template<class type>
  class ConstructorImpl {
  public:
    PA_INJECT(type)

    static void apply(const EncryptedArrayDerived<type>& ea, NewPlaintextArray& pa)
    {
      CloneablePtr< NewPlaintextArrayDerived<type> > ptr =
         MakeCloneable< NewPlaintextArrayDerived<type> >();
      ptr->data.resize(ea.size());
      pa.rep = ptr;
    }
  };

public:
  
  NewPlaintextArray(const EncryptedArray& ea)  
    { ea.dispatch<ConstructorImpl>(*this); }

  NewPlaintextArray(const NewPlaintextArray& other) : rep(other.rep.clone()) { }
  NewPlaintextArray& operator=(const NewPlaintextArray& other) 
    { rep = other.rep.clone(); return *this; }

  template<class type>
    vector<typename type::RX>& getData() 
    { return (dynamic_cast< NewPlaintextArrayDerived<type>& >(*rep)).data; }


  template<class type>
    const vector<typename type::RX>& getData() const
    { return (dynamic_cast< NewPlaintextArrayDerived<type>& >(*rep)).data; }


  void print(ostream& s) const { rep->print(s); }

};

inline 
ostream& operator<<(ostream& s, const NewPlaintextArray& pa)
{  pa.print(s); return s; }


void rotate(const EncryptedArray& ea, NewPlaintextArray& pa, long k);
void shift(const EncryptedArray& ea, NewPlaintextArray& pa, long k);

void encode(const EncryptedArray& ea, NewPlaintextArray& pa, const vector<long>& array);
void encode(const EncryptedArray& ea, NewPlaintextArray& pa, const vector<ZZX>& array);
void encode(const EncryptedArray& ea, NewPlaintextArray& pa, long val);
void encode(const EncryptedArray& ea, NewPlaintextArray& pa, const ZZX& val);

void random(const EncryptedArray& ea, NewPlaintextArray& pa);

void decode(const EncryptedArray& ea, vector<long>& array, const NewPlaintextArray& pa);
void decode(const EncryptedArray& ea, vector<ZZX>& array, const NewPlaintextArray& pa);

bool equals(const EncryptedArray& ea, const NewPlaintextArray& pa, const NewPlaintextArray& other);
bool equals(const EncryptedArray& ea, const NewPlaintextArray& pa, const vector<long>& other);
bool equals(const EncryptedArray& ea, const NewPlaintextArray& pa, const vector<ZZX>& other);

void add(const EncryptedArray& ea, NewPlaintextArray& pa, const NewPlaintextArray& other);
void sub(const EncryptedArray& ea, NewPlaintextArray& pa, const NewPlaintextArray& other);
void mul(const EncryptedArray& ea, NewPlaintextArray& pa, const NewPlaintextArray& other);
void negate(const EncryptedArray& ea, NewPlaintextArray& pa);


void frobeniusAutomorph(const EncryptedArray& ea, NewPlaintextArray& pa, long j);
void frobeniusAutomorph(const EncryptedArray& ea, NewPlaintextArray& pa, const Vec<long>& vec);

void applyPerm(const EncryptedArray& ea, NewPlaintextArray& pa, const Vec<long>& pi);

void power(const EncryptedArray& ea, NewPlaintextArray& pa, long e);




// Following are functions for performing "higher level"
// operations on "encrypted arrays".  There is really no
// reason for these to be members of the EncryptedArray class,
// so they are just declared as separate functions.

//! @brief A ctxt that encrypts \f$(x_1, ..., x_n)\f$ is replaced by an
//! encryption of \f$(y_1, ..., y_n)\f$, where \f$y_i = sum_{j\le i} x_j\f$.
void runningSums(const EncryptedArray& ea, Ctxt& ctxt);
// The implementation uses O(log n) shift operations.


//! @brief A ctxt that encrypts \f$(x_1, ..., x_n)\f$ is replaced by an
//! encryption of \f$(y, ..., y)\$, where \f$y = sum_{j=1}^n x_j.\f$
void totalSums(const EncryptedArray& ea, Ctxt& ctxt);


//! @brief Map all non-zero slots to 1, leaving zero slots as zero.
//! Assumes that r=1, and that all the slots contain elements from GF(p^d).
void mapTo01(const EncryptedArray& ea, Ctxt& ctxt);
// Implemented in eqtesting.cpp. We compute
//             x^{p^d-1} = x^{(1+p+...+p^{d-1})*(p-1)}
// by setting y=x^{p-1} and then outputting y * y^p * ... * y^{p^{d-1}},
// with exponentiation to powers of p done via Frobenius.


//! @brief (only for p=2, r=1), test if prefixes of bits in slots are all zero.
//! Set slot j of res[i] to 0 if bits 0..i of j'th slot in ctxt are all zero,
//! else sets slot j of res[i] to 1.
//! It is assumed that res and the res[i]'s are initialized by the caller.
void incrementalZeroTest(Ctxt* res[], const EncryptedArray& ea,
			 const Ctxt& ctxt, long n);
// Complexity: O(d + n log d) smart automorphisms
//             O(n d) 

/*************** End linear transformation functions ****************/
/********************************************************************/

///@{
/**
 * @name Apply linearized polynomials to a ciphertext.
 *
 * Example usage: The map L selects just the even coefficients
 * \code
 *     long d = ea.getDegree();
 *     vector<ZZX> L(d);
 *     for (long j = 0; j < d; j++)
 *       if (j % 2 == 0) L[j] = ZZX(j, 1);
 *
 *     vector<ZZX> C;
 *     ea.buildLinPolyCoeffs(C, L); 
 *     applyLinPoly1(ea, ctxt, C);
 * \endcode
 */

//! @brief Apply the same linear transformation to all the slots
//! @param C is the output of ea.buildLinPolyCoeffs;
void applyLinPoly1(const EncryptedArray& ea, Ctxt& ctxt, const vector<ZZX>& C);

//! @brief Apply different transformations to different slots
//! @param Cvec is a vector of length ea.size(), each entry of which
//!        is the output of ea.buildLinPolyCoeffs; 
void applyLinPolyMany(const EncryptedArray& ea, Ctxt& ctxt, 
                      const vector< vector<ZZX> >& Cvec);

//! @brief a low-level variant:
//! @param encodedCoeffs has all the linPoly coeffs encoded  in slots;
//!        different transformations can be encoded in different slots
template<class P>  // P can be ZZX or DoubleCRT
void applyLinPolyLL(Ctxt& ctxt, const vector<P>& encodedC, long d);
///@}

#endif /* ifdef _EncryptedArray_H_ */
