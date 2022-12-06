/* Copyright (C) 2012-2021 IBM Corp.
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
#ifndef HELIB_ENCRYPTEDARRAY_H
#define HELIB_ENCRYPTEDARRAY_H
/**
 * @file EncryptedArray.h
 * @brief Data-movement operations on encrypted arrays of slots
 */
#include <exception>
#include <cmath>
#include <complex>
#include <NTL/Lazy.h>
#include <NTL/pair.h>
#include <NTL/SmartPtr.h>

#include <helib/DoubleCRT.h>
#include <helib/Context.h>
#include <helib/Ctxt.h>
#include <helib/keys.h>
#include <helib/exceptions.h>
#include <helib/log.h>

namespace helib {

typedef std::complex<double> cx_double;

// DIRT: we're using undocumented NTL interfaces here
//   also...this probably should be defined in NTL, anyway....
#define HELIB_MORE_UNWRAPARGS(n)                                               \
  NTL_SEPARATOR_##n NTL_REPEATER_##n(NTL_UNWRAPARG)

// these are used to implement PlaintextArray stuff routines

// NOTE: Variables have been marked as UNUSED to silence the unused variable
// warnings due to unused injected variables. It has not been possible to mark
// the whole section with:
// _Pragma("GCC diagnostic ignored \"-Wunused-variable\"") because GCC 5.4 does
// not honor such _Pragma.
// NOTE: Consider marking the section with _Pragma removing the UNUSED tag when
// GCC 5.4 won't be supported anymore.
#define PA_BOILER(type)                                                        \
  const PAlgebraModDerived<type>& tab = ea.getTab();                           \
  UNUSED const RX& G = ea.getG();                                              \
  UNUSED long n = ea.size();                                                   \
  UNUSED long d = ea.getDegree();                                              \
  std::vector<RX>& data = pa.getData<type>();                                  \
  RBak bak;                                                                    \
  bak.save();                                                                  \
  tab.restoreContext();

// NOTE: Variables have been marked as UNUSED to silence the unused variable
// warnings due to unused injected variables. It has not been possible to mark
// the whole section with:
// _Pragma("GCC diagnostic ignored \"-Wunused-variable\"") because GCC 5.4 does
// not honor such _Pragma.
// NOTE: Consider marking the section with _Pragma removing the UNUSED tag when
// GCC 5.4 won't be supported anymore.
#define CPA_BOILER(type)                                                       \
  const PAlgebraModDerived<type>& tab = ea.getTab();                           \
  UNUSED const RX& G = ea.getG();                                              \
  UNUSED long n = ea.size();                                                   \
  UNUSED long d = ea.getDegree();                                              \
  const std::vector<RX>& data = pa.getData<type>();                            \
  RBak bak;                                                                    \
  bak.save();                                                                  \
  tab.restoreContext();

class PlaintextArray; // forward reference
class PtxtArray;      // forward reference
class EncryptedArray; // forward reference

typedef EncryptedArray View;
// New and improved name for EncryptedArray.
// Documentation should use this name.

/**
 * @class EncryptedArrayBase
 * @brief virtual class for data-movement operations on arrays of slots
 *
 * An object ea of type EncryptedArray stores information about an
 * Context context, and a monic polynomial G.  If context defines
 * parameters m, p, and r, then ea is a helper abject
 * that supports encoding/decoding and encryption/decryption
 * of std::vectors of plaintext slots over the ring (Z/(p^r)[X])/(G).
 *
 * The polynomial G should be irreducible over Z/(p^r) (this is not checked).
 * The degree of G should divide the multiplicative order of p modulo m
 * (this is checked). Currently, the following restriction is imposed:
 *
 * either r == 1 or deg(G) == 1 or G == factors[0].
 *
 * ea stores objects in the polynomial ring Z/(p^r)[X].
 *
 * Just as for the class PAlgebraMod, if p == 2 and r == 1, then these
 * polynomials are represented as GF2X's, and otherwise as zz_pX's.
 * Thus, the types of these objects are not determined until run time.
 * As such, we need to use a class hierarchy, which mirrors that of
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
 **/

class EncryptedArrayBase
{ // purely abstract interface
public:
  virtual ~EncryptedArrayBase() {}

  virtual EncryptedArrayBase* clone() const = 0;
  // makes this usable with ClonedPtr

  virtual PA_tag getTag() const = 0;

  virtual const Context& getContext() const = 0;
  virtual const PAlgebra& getPAlgebra() const = 0;
  virtual long getDegree() const = 0;
  virtual long getP2R() const = 0;

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
  virtual void rotate1D(Ctxt& ctxt, long i, long k, bool dc = false) const = 0;

  //! @brief Right shift k positions along the i'th dimension with zero fill
  virtual void shift1D(Ctxt& ctxt, long i, long k) const = 0;

  ///@{
  //! @name Encoding/decoding methods
  // encode/decode arrays into plaintext polynomials
  // These methods are working for some of the derived classes (throwing
  // otherwise)
  virtual void encode(zzX& ptxt, const std::vector<long>& array) const = 0;
  virtual void encode(NTL::ZZX& ptxt, const std::vector<long>& array) const = 0;

  virtual void encode(zzX& ptxt, const std::vector<NTL::ZZX>& array) const = 0;
  virtual void encode(NTL::ZZX& ptxt,
                      const std::vector<NTL::ZZX>& array) const = 0;

  virtual void encode(zzX& ptxt, const PlaintextArray& array) const = 0;
  virtual void encode(NTL::ZZX& ptxt, const PlaintextArray& array) const = 0;

  virtual void encode(zzX& ptxt, const std::vector<zzX>& array) const = 0;

  //=============== new EncodedPtxt interfaces

  // BGV only
  virtual void encode(EncodedPtxt& eptxt,
                      const std::vector<NTL::ZZX>& array) const = 0;
  virtual void encode(EncodedPtxt& eptxt,
                      const std::vector<long>& array) const = 0;

  // CKKS only
  // mag: defaults to Norm(array).
  // prec: defaults to r=getAlMod().getR(), which
  // is usually the same as context.getPrecision().

  // mag should be an upper bound on Norm(array).
  // If an encoding will be encrypted, the user may wish
  // to hide Norm(array) by setting mag to some data-independent
  // upper bound. A warning is issued if Norm(array) > mag.

  // The encoding will normally have an accuracy of 2^{-prec}, meaning that
  // Norm(array - decode(encode(array))) <= 2^{-prec}.
  // Note that prec may be positive, negative, or zero.
  // The exact logic is a bit heuristic, and a warning is
  // issued if the the accuracy exceeds 2^{-prec}.

  // NOTE: Norm above is the infinity (i.e., max) norm.

  virtual void encode(EncodedPtxt& eptxt,
                      const std::vector<cx_double>& array,
                      double mag = -1,
                      OptLong prec = OptLong()) const = 0;

  virtual void encode(EncodedPtxt& eptxt,
                      const std::vector<double>& array,
                      double mag = -1,
                      OptLong prec = OptLong()) const = 0;

  // BGV and CKKS
  virtual void encode(EncodedPtxt& eptxt,
                      const PlaintextArray& array,
                      double mag = -1,
                      OptLong prec = OptLong()) const = 0;
  // NOTE: for BGV, mag,prec must be defaulted

  virtual void encode(EncodedPtxt& eptxt,
                      const std::vector<bool>& array) const = 0;
  // NOTE: for CKKS, mag,prec are default

  virtual void encodeUnitSelector(EncodedPtxt& eptxt, long i) const = 0;
  // NOTE: for CKKS, mag,prec are default

  virtual double defaultScale(UNUSED double err,
                              UNUSED OptLong prec = OptLong()) const
  {
    throw LogicError("function not implemented");
  }

  virtual double defaultErr() const
  {
    throw LogicError("function not implemented");
  }

  //================================

  // These methods are working for some of the derived classes (throwing
  // otherwise)
  virtual void decode(std::vector<long>& array, const NTL::ZZX& ptxt) const = 0;
  virtual void decode(std::vector<NTL::ZZX>& array,
                      const NTL::ZZX& ptxt) const = 0;
  virtual void decode(PlaintextArray& array, const NTL::ZZX& ptxt) const = 0;

  virtual void random(std::vector<long>& array) const = 0; // must be defined

  // These methods are working for some of the derived classes (throwing
  // otherwise)
  virtual void random(std::vector<NTL::ZZX>& array) const = 0;

  // FIXME: Inefficient implementation, calls usual decode and returns one slot
  long decode1Slot(const NTL::ZZX& ptxt, long i) const
  {
    std::vector<long> v;
    decode(v, ptxt);
    return v.at(i);
  }
  void decode1Slot(NTL::ZZX& slot, const NTL::ZZX& ptxt, long i) const
  {
    std::vector<NTL::ZZX> v;
    decode(v, ptxt);
    slot = v.at(i);
  }

  //! @brief Encodes a std::vector with 1 at position i and 0 everywhere else
  virtual void encodeUnitSelector(zzX& ptxt, long i) const = 0;
  ///@}

  ///@{

  // These methods are working for some of the derived calsses (throwing
  // otherwise)
  virtual void decrypt(const Ctxt& ctxt,
                       const SecKey& sKey,
                       std::vector<long>& ptxt) const = 0;

  virtual void decrypt(const Ctxt& ctxt,
                       const SecKey& sKey,
                       std::vector<NTL::ZZX>& ptxt) const = 0;

  virtual void decrypt(const Ctxt& ctxt,
                       const SecKey& sKey,
                       std::vector<double>& ptxt,
                       OptLong prec = OptLong()) const = 0;

  virtual void decrypt(const Ctxt& ctxt,
                       const SecKey& sKey,
                       std::vector<cx_double>& ptxt,
                       OptLong prec = OptLong()) const = 0;

  virtual void rawDecrypt(const Ctxt& ctxt,
                          const SecKey& sKey,
                          std::vector<cx_double>& ptxt) const = 0;

  virtual void rawDecrypt(const Ctxt& ctxt,
                          const SecKey& sKey,
                          std::vector<double>& ptxt) const = 0;

  virtual void decrypt(const Ctxt& ctxt,
                       const SecKey& sKey,
                       PlaintextArray& ptxt,
                       OptLong prec = OptLong()) const = 0;

  virtual void decryptComplex(const Ctxt& ctxt,
                              const SecKey& sKey,
                              PlaintextArray& ptxt,
                              OptLong prec = OptLong()) const = 0;

  virtual void decryptReal(const Ctxt& ctxt,
                           const SecKey& sKey,
                           PlaintextArray& ptxt,
                           OptLong prec = OptLong()) const = 0;

  virtual void rawDecrypt(const Ctxt& ctxt,
                          const SecKey& sKey,
                          PlaintextArray& ptxt) const = 0;

  virtual void rawDecryptComplex(const Ctxt& ctxt,
                                 const SecKey& sKey,
                                 PlaintextArray& ptxt) const = 0;

  virtual void rawDecryptReal(const Ctxt& ctxt,
                              const SecKey& sKey,
                              PlaintextArray& ptxt) const = 0;

  // FIXME: Inefficient implementation, calls usual decrypt and returns one slot
  long decrypt1Slot(const Ctxt& ctxt, const SecKey& sKey, long i) const
  {
    std::vector<long> v;
    decrypt(ctxt, sKey, v);
    return v.at(i);
  }
  void decrypt1Slot(NTL::ZZX& slot,
                    const Ctxt& ctxt,
                    const SecKey& sKey,
                    long i) const
  {
    std::vector<NTL::ZZX> v;
    decrypt(ctxt, sKey, v);
    slot = v.at(i);
  }
  ///@}

  //! @brief Linearized polynomials.
  //! L describes a linear map M by describing its action on the standard
  //! power basis: M(x^j mod G) = (L[j] mod G), for j = 0..d-1.
  //! The result is a coefficient vector C for the linearized polynomial
  //! representing M: a polynomial h in Z/(p^r)[X] of degree < d is sent to
  //! \f[
  //!  M(h(X) \bmod G)= \sum_{j=0}^{d-1}(C[j] \cdot h(X^{p^j}))\bmod G).
  //! \f]
  // These methods are working for some of the derived calsses (throwing
  // otherwise)
  virtual void buildLinPolyCoeffs(std::vector<NTL::ZZX>& C,
                                  const std::vector<NTL::ZZX>& L) const = 0;

  // restore contexts mod p and mod G
  virtual void restoreContext() const {}
  virtual void restoreContextForG() const {}

  /* some non-virtual convenience functions */

  //! @brief Total size (# of slots) of hypercube
  long size() const { return getPAlgebra().getNSlots(); }

  //! @brief Number of dimensions of hypercube
  long dimension() const { return getPAlgebra().numOfGens(); }

  //! @brief Size of given dimension
  long sizeOfDimension(long i) const { return getPAlgebra().OrderOf(i); }

  //! @brief Is rotations in given dimension a "native" operation?
  bool nativeDimension(long i) const { return getPAlgebra().SameOrd(i); }

  //! @brief returns coordinate of index k along the i'th dimension
  long coordinate(long i, long k) const
  {
    return getPAlgebra().coordinate(i, k);
  }

  //! @brief adds offset to index k in the i'th dimension
  long addCoord(long i, long k, long offset) const
  {
    return getPAlgebra().addCoord(i, k, offset);
  }

  //! @brief rotate an array by offset in the i'th dimension
  //! (output should not alias input)
  template <typename U>
  void rotate1D(std::vector<U>& out,
                const std::vector<U>& in,
                long i,
                long offset) const
  {
    assertEq(lsize(in),
             size(),
             "Input vector has wrong size (must equal EncryptedArray::size())");
    out.resize(in.size());
    for (long j : range(size()))
      out[addCoord(i, j, offset)] = in[j];
  }
};

/**
 * @class EncryptedArrayDerived
 * @brief Derived concrete implementation of EncryptedArrayBase
 */
template <typename type>
class EncryptedArrayDerived : public EncryptedArrayBase
{
public:
  PA_INJECT(type)

private:
  const Context& context;
  const PAlgebraModDerived<type>& tab;

  // FIXME: all of these types should be copyable
  // out of context, but NTL 8.0 still does not copy
  // matrix copies out of context correctly, as it
  // relies on plain SetLength...I need to fix this
  // in NTL
  //
  MappingData<type> mappingData; // MappingData is defined in PAlgebra.h

  NTL::Lazy<NTL::Mat<RE>> linPolyMatrix;

  NTL::Lazy<NTL::Pair<NTL::Mat<R>, NTL::Mat<R>>> normalBasisMatrices;
  // a is the matrix, b is its inverse

public:
  explicit EncryptedArrayDerived(const Context& _context,
                                 const RX& _G,
                                 const PAlgebraMod& _tab);

  EncryptedArrayDerived(const EncryptedArrayDerived& other) // copy constructor
      :
      context(other.context), tab(other.tab)
  {
    RBak bak;
    bak.save();
    tab.restoreContext();
    REBak ebak;
    ebak.save();
    other.mappingData.restoreContextForG();
    mappingData = other.mappingData;
    linPolyMatrix = other.linPolyMatrix;
    normalBasisMatrices = other.normalBasisMatrices;
  }

  EncryptedArrayDerived& operator=(const EncryptedArrayDerived&) = delete;

  virtual EncryptedArrayBase* clone() const override
  {
    return new EncryptedArrayDerived(*this);
  }

  virtual PA_tag getTag() const override { return tag; }
  // tag is defined in PA_INJECT, see PAlgebra.h

  template <template <typename> class T, typename... Args>
  void dispatch(Args&&... args) const
  {
    T<type>::apply(*this, std::forward<Args>(args)...);
  }

  const RX& getG() const { return mappingData.getG(); }

  const NTL::Mat<R>& getNormalBasisMatrix() const
  {
    if (!normalBasisMatrices.built())
      initNormalBasisMatrix();
    return normalBasisMatrices->a;
  }

  const NTL::Mat<R>& getNormalBasisMatrixInverse() const
  {
    if (!normalBasisMatrices.built())
      initNormalBasisMatrix();
    return normalBasisMatrices->b;
  }

  void initNormalBasisMatrix() const;

  virtual void restoreContext() const override { tab.restoreContext(); }
  virtual void restoreContextForG() const override
  {
    mappingData.restoreContextForG();
  }

  virtual const Context& getContext() const override { return context; }
  virtual const PAlgebra& getPAlgebra() const override
  {
    return tab.getZMStar();
  }
  virtual long getDegree() const override { return mappingData.getDegG(); }
  const PAlgebraModDerived<type>& getTab() const { return tab; }

  virtual void rotate(Ctxt& ctxt, long k) const override;
  virtual void shift(Ctxt& ctxt, long k) const override;
  virtual void rotate1D(Ctxt& ctxt,
                        long i,
                        long k,
                        bool dc = false) const override;

  long getP2R() const override { return getTab().getPPowR(); }

  // avoid this being "hidden" by other rotate1D's
  template <typename U>
  void rotate1D(std::vector<U>& out,
                const std::vector<U>& in,
                long i,
                long offset) const
  {
    EncryptedArrayBase::rotate1D(out, in, i, offset);
  }
  virtual void shift1D(Ctxt& ctxt, long i, long k) const override;

  /* Begin CKKS functions. They will simply throw here. */
  /**
   * @brief Unimplemented decrypt function for CKKS. It will always
   * throw helib::LogicError
   * @param ctxt Unused.
   * @param sKey Unused.
   * @param ptxt Unused.
   * @param prec Unused.
   */
  // TODO: document this better (especially the prec parameter)
  void decrypt(UNUSED const Ctxt& ctxt,
               UNUSED const SecKey& sKey,
               UNUSED std::vector<double>& ptxt,
               UNUSED OptLong prec = OptLong()) const override
  {
    throw LogicError("Unimplemented: "
                     "EncryptedArrayDerived::decrypt for CKKS type");
  }
  /**
   * @brief Unimplemented decrypt function for CKKS. It will always
   * throw helib::LogicError
   * @param ctxt Unused.
   * @param sKey Unused.
   * @param ptxt Unused.
   * @param prec Unused.
   */
  // TODO: document this better (especially the prec parameter)
  void decrypt(UNUSED const Ctxt& ctxt,
               UNUSED const SecKey& sKey,
               UNUSED std::vector<cx_double>& ptxt,
               UNUSED OptLong prec = OptLong()) const override
  {
    throw LogicError("Unimplemented: "
                     "EncryptedArrayDerived::decrypt for CKKS type");
  }

  void rawDecrypt(UNUSED const Ctxt& ctxt,
                  UNUSED const SecKey& sKey,
                  UNUSED std::vector<cx_double>& ptxt) const override
  {
    throw LogicError("Unimplemented: function only available for CKKS");
  }

  void rawDecrypt(UNUSED const Ctxt& ctxt,
                  UNUSED const SecKey& sKey,
                  UNUSED std::vector<double>& ptxt) const override
  {
    throw LogicError("Unimplemented: function only available for CKKS");
  }

  void rawDecrypt(UNUSED const Ctxt& ctxt,
                  UNUSED const SecKey& sKey,
                  UNUSED PlaintextArray& ptxt) const override
  {
    throw LogicError("function not implemented");
  }

  void decryptComplex(UNUSED const Ctxt& ctxt,
                      UNUSED const SecKey& sKey,
                      UNUSED PlaintextArray& ptxt,
                      UNUSED OptLong prec = OptLong()) const override
  {
    throw LogicError("function not implemented");
  }

  void rawDecryptComplex(UNUSED const Ctxt& ctxt,
                         UNUSED const SecKey& sKey,
                         UNUSED PlaintextArray& ptxt) const override
  {
    throw LogicError("function not implemented");
  }

  void decryptReal(UNUSED const Ctxt& ctxt,
                   UNUSED const SecKey& sKey,
                   UNUSED PlaintextArray& ptxt,
                   UNUSED OptLong prec = OptLong()) const override
  {
    throw LogicError("function not implemented");
  }

  void rawDecryptReal(UNUSED const Ctxt& ctxt,
                      UNUSED const SecKey& sKey,
                      UNUSED PlaintextArray& ptxt) const override
  {
    throw LogicError("function not implemented");
  }
  /* End CKKS functions. */

  virtual void encode(NTL::ZZX& ptxt,
                      const std::vector<long>& array) const override
  {
    genericEncode(ptxt, array);
  }

  virtual void encode(zzX& ptxt, const std::vector<long>& array) const override
  {
    genericEncode(ptxt, array);
  }

  virtual void encode(NTL::ZZX& ptxt,
                      const std::vector<NTL::ZZX>& array) const override
  {
    genericEncode(ptxt, array);
  }

  virtual void encode(zzX& ptxt, const std::vector<zzX>& array) const override
  {
    genericEncode(ptxt, array);
  }

  virtual void encode(NTL::ZZX& ptxt,
                      const PlaintextArray& array) const override;
  virtual void encode(zzX& ptxt, const PlaintextArray& array) const override;

  virtual void encodeUnitSelector(zzX& ptxt, long i) const override;

  virtual void encode(zzX& ptxt,
                      const std::vector<NTL::ZZX>& array) const override
  {
    NTL::ZZX tmp;
    encode(tmp, array);
    convert(ptxt, tmp);
  }

  //================ new EncodedPtxt interfaces

  virtual void encode(EncodedPtxt& eptxt,
                      const std::vector<NTL::ZZX>& array) const override
  {
    zzX poly;
    encode(poly, array);
    eptxt.resetBGV(poly, getP2R(), getContext());
  }

  virtual void encode(EncodedPtxt& eptxt,
                      const std::vector<long>& array) const override
  {
    zzX poly;
    encode(poly, array);
    eptxt.resetBGV(poly, getP2R(), getContext());
  }

  virtual void encode(UNUSED EncodedPtxt& eptxt,
                      UNUSED const std::vector<cx_double>& array,
                      UNUSED double mag = -1,
                      UNUSED OptLong prec = OptLong()) const override
  {
    throw LogicError("function not implemented for BGV");
  }

  virtual void encode(UNUSED EncodedPtxt& eptxt,
                      UNUSED const std::vector<double>& array,
                      UNUSED double mag = -1,
                      UNUSED OptLong prec = OptLong()) const override
  {
    throw LogicError("function not implemented for BGV");
  }

  virtual void encode(EncodedPtxt& eptxt,
                      const PlaintextArray& array,
                      double mag = -1,
                      OptLong prec = OptLong()) const override
  {
    assertTrue(mag < 0 && !prec.isDefined(),
               "BGV encoding: mag,prec set must be defaulted");
    zzX poly;
    encode(poly, array);
    eptxt.resetBGV(poly, getP2R(), getContext());
  }

  virtual void encode(EncodedPtxt& eptxt,
                      const std::vector<bool>& array) const override
  {
    std::vector<long> array1;
    convert(array1, array);
    encode(eptxt, array1);
  }

  virtual void encodeUnitSelector(EncodedPtxt& eptxt, long i) const override
  {
    zzX poly;
    encodeUnitSelector(poly, i);
    eptxt.resetBGV(poly, getP2R(), getContext());
  }

  //===========================================

  virtual void decode(std::vector<long>& array,
                      const NTL::ZZX& ptxt) const override
  {
    genericDecode(array, ptxt);
  }

  virtual void decode(std::vector<NTL::ZZX>& array,
                      const NTL::ZZX& ptxt) const override
  {
    genericDecode(array, ptxt);
  }

  virtual void decode(PlaintextArray& array,
                      const NTL::ZZX& ptxt) const override;
  virtual void decode(PlaintextArray& array, const zzX& ptxt) const;

  // choose at random and convert to std::vector<long>
  virtual void random(std::vector<long>& array) const override
  {
    genericRandom(array);
  }

  // choose at random and convert to std::vector<ZZX>
  virtual void random(std::vector<NTL::ZZX>& array) const override
  {
    genericRandom(array);
  }

  virtual void decrypt(const Ctxt& ctxt,
                       const SecKey& sKey,
                       std::vector<long>& ptxt) const override
  {
    genericDecrypt(ctxt, sKey, ptxt);
    if (ctxt.getPtxtSpace() < getP2R()) {
      helib::Warning("EncryptedArray::decrypt: reducing plaintext modulus");
      for (long i = 0; i < (long)ptxt.size(); i++)
        ptxt[i] %= ctxt.getPtxtSpace();
    }
  }

  virtual void decrypt(const Ctxt& ctxt,
                       const SecKey& sKey,
                       std::vector<NTL::ZZX>& ptxt) const override
  {
    genericDecrypt(ctxt, sKey, ptxt);
    if (ctxt.getPtxtSpace() < getP2R()) {
      helib::Warning("EncryptedArray::decrypt: reducing plaintext modulus");
      for (long i = 0; i < (long)ptxt.size(); i++)
        PolyRed(ptxt[i], ctxt.getPtxtSpace(), /*abs=*/true);
    }
  }

  virtual void decrypt(const Ctxt& ctxt,
                       const SecKey& sKey,
                       PlaintextArray& ptxt,
                       OptLong prec = OptLong()) const override
  {
    if (prec.isDefined())
      throw LogicError("EncryptedArray::decrypt: the precision parameter "
                       "(prec) must be defaulted");

    genericDecrypt(ctxt, sKey, ptxt);
    if (ctxt.getPtxtSpace() < getP2R()) {
      throw LogicError("EncryptedArray::decrypt: bad plaintext modulus");
      // FIXME: Reduce mod the ciphertext plaintext space as above
      // What we would have to do is:
      //   1. convert the PlaintextArray to a vector of ZZX's
      //   2. call PolyRed as above to reduce coefficients
      //   3. convert the vector of ZZX's back to a PlaintextArray
      // But do we *really* want to do this?
    }
  }

  virtual void buildLinPolyCoeffs(
      std::vector<NTL::ZZX>& C,
      const std::vector<NTL::ZZX>& L) const override;

  /* the following are specialized methods, used to work over extension
     fields... they assume the modulus context is already set
   */

  void encode(zzX& ptxt, const std::vector<RX>& array) const;
  void decode(std::vector<RX>& array, const zzX& ptxt) const;

  void encode(NTL::ZZX& ptxt, const std::vector<RX>& array) const;
  void decode(std::vector<RX>& array, const NTL::ZZX& ptxt) const;

  void encode(RX& ptxt, const std::vector<RX>& array) const;
  void decode(std::vector<RX>& array, const RX& ptxt) const;

  // Choose random polynomial of the right degree, coeffs in GF2 or zz_p
  void random(std::vector<RX>& array) const
  {
    array.resize(size());
    for (long i = 0; i < size(); i++)
      NTL::random(array[i], getDegree());
  }

  void decrypt(const Ctxt& ctxt,
               const SecKey& sKey,
               std::vector<RX>& ptxt) const
  {
    genericDecrypt(ctxt, sKey, ptxt);
  }

  virtual void buildLinPolyCoeffs(std::vector<RX>& C,
                                  const std::vector<RX>& L) const;

private:
  /* helper template methods, to avoid repetitive code */

  template <typename T>
  void genericEncode(NTL::ZZX& ptxt, const T& array) const
  {
    RBak bak;
    bak.save();
    tab.restoreContext();

    std::vector<RX> array1;
    convert(array1, array);
    encode(ptxt, array1);
  }

  template <typename T>
  void genericEncode(zzX& ptxt, const T& array) const
  {
    RBak bak;
    bak.save();
    tab.restoreContext();

    std::vector<RX> array1;
    convert(array1, array);
    encode(ptxt, array1);
  }

  template <typename T>
  void genericDecode(T& array, const NTL::ZZX& ptxt) const
  {
    RBak bak;
    bak.save();
    tab.restoreContext();

    std::vector<RX> array1;
    decode(array1, ptxt);
    convert(array, array1);
  }

  // T is std::vector<long> or std::vector<ZZX>
  template <typename T>
  void genericRandom(T& array) const
  {
    static_assert(std::is_same<T, std::vector<long>>::value ||
                      std::is_same<T, std::vector<NTL::ZZX>>::value,
                  "Template type T must be either std::vector<long> or "
                  "std::vector<NTL::ZZX>.");
    RBak bak;
    // backup NTL modulus
    bak.save();
    tab.restoreContext();

    std::vector<RX> array1; // RX is GF2X or zz_pX
    random(array1);         // choose random coefficients from GF2/zz_p
    convert(array, array1); // convert to type T (see NumbTh.h)
  }

  template <typename T>
  void genericDecrypt(const Ctxt& ctxt, const SecKey& sKey, T& array) const
  {
    assertEq(&context,
             &ctxt.getContext(),
             "Cannot decrypt when ciphertext has different context than "
             "EncryptedArray");
    NTL::ZZX pp;
    sKey.Decrypt(pp, ctxt);
    decode(array, pp);
  }
};

//! A different derived class to be used for the approximate-numbers scheme
template <>
class EncryptedArrayDerived<PA_cx> : public EncryptedArrayBase
{
  const Context& context;
  const PAlgebraModCx& alMod;

  zzX iEncoded; // an encoded plaintext with i in all the slots
  // VJS-FIXME: this is a bad idea...see commenets in EaCx.cpp

public:
  static double roundedSize(double x)
  {
    long rounded = std::ceil(std::fabs(x));
    if (rounded < 1)
      rounded = 1;
    return double(1L << NTL::NumBits(rounded - 1));
  }

  double encodei(zzX& ptxt, long precision = -1) const; // encode i in all slots

  explicit EncryptedArrayDerived(const Context& _context) :
      context(_context), alMod(context.getAlMod().getCx())
  {
    clear(iEncoded);
  }
  EncryptedArrayDerived(const Context& _context, const PAlgebraModCx& _alMod) :
      context(_context), alMod(_alMod)
  {
    clear(iEncoded);
  }

  EncryptedArrayBase* clone() const override
  {
    return new EncryptedArrayDerived(*this);
  }

  const zzX& getiEncoded() const;
  PA_tag getTag() const override { return PA_cx_tag; }
  const Context& getContext() const override { return context; }
  const PAlgebra& getPAlgebra() const override { return alMod.getZMStar(); }
  long getDegree() const override { return 2L; }

  void rotate(Ctxt& ctxt, long k) const override;
  void shift(Ctxt& ctxt, long k) const override;
  void rotate1D(Ctxt& ctxt, long i, long k, bool dc = false) const override;
  void shift1D(Ctxt& ctxt, long i, long k) const override;

  long getP2R() const override { return alMod.getPPowR(); }

  // the following help with some template code
  cx_double getG() const { return 0.0; }
  const PAlgebraModCx getTab() const { return alMod; }

  template <template <typename> class T, typename... Args>
  void dispatch(Args&&... args) const
  {
    T<PA_cx>::apply(*this, std::forward<Args>(args)...);
  }

  /* Begin BGV functions. They will simply throw here. */
  // encode
  /**
   * @brief Unimplemented encode function for BGV. It will always throw
   * helib::LogicError.
   * @param ptxt Unused.
   * @param array Unused.
   */
  void encode(UNUSED zzX& ptxt,
              UNUSED const std::vector<long>& array) const override
  {
    throw LogicError("Unimplemented: EncryptedArrayCx::encode for BGV type");
  }
  /**
   * @brief Unimplemented encode function for BGV. It will always throw
   * helib::LogicError.
   * @param ptxt Unused.
   * @param array Unused.
   */
  void encode(UNUSED NTL::ZZX& ptxt,
              UNUSED const std::vector<long>& array) const override
  {
    throw LogicError("Unimplemented: EncryptedArrayCx::encode for BGV type");
  }
  /**
   * @brief Unimplemented encode function for BGV. It will always throw
   * helib::LogicError.
   * @param ptxt Unused.
   * @param array Unused.
   */
  void encode(UNUSED zzX& ptxt,
              UNUSED const std::vector<zzX>& array) const override
  {
    throw LogicError("Unimplemented: EncryptedArrayCx::encode for BGV type");
  }
  /**
   * @brief Unimplemented encode function for BGV. It will always throw
   * helib::LogicError.
   * @param ptxt Unused.
   * @param array Unused.
   */
  void encode(UNUSED zzX& ptxt,
              UNUSED const PlaintextArray& array) const override
  {
    throw LogicError("Unimplemented: EncryptedArrayCx::encode for BGV type");
  }
  /**
   * @brief Unimplemented encode function for BGV. It will always throw
   * helib::LogicError.
   * @param ptxt Unused.
   * @param array Unused.
   */
  void encode(UNUSED NTL::ZZX& ptxt,
              UNUSED const std::vector<NTL::ZZX>& array) const override
  {
    throw LogicError("Unimplemented: EncryptedArrayCx::encode for BGV type");
  }
  /**
   * @brief Unimplemented encode function for BGV. It will always throw
   * helib::LogicError.
   * @param ptxt Unused.
   * @param array Unused.
   */
  void encode(UNUSED NTL::ZZX& ptxt,
              UNUSED const PlaintextArray& array) const override
  {
    throw LogicError("Unimplemented: EncryptedArrayCx::encode for BGV type");
  }

  void encode(UNUSED zzX& ptxt,
              UNUSED const std::vector<NTL::ZZX>& array) const override
  {
    throw LogicError("Unimplemented: EncryptedArrayCx::encode for BGV type");
  }

  // decode
  /**
   * @brief Unimplemented decode function for BGV. It will always throw
   * helib::LogicError.
   * @param array Unused.
   * @param ptxt Unused.
   */
  void decode(UNUSED std::vector<long>& array,
              UNUSED const NTL::ZZX& ptxt) const override
  {
    throw LogicError("Unimplemented: EncryptedArrayCx::decode for BGV type");
  }
  /**
   * @brief Unimplemented decode function for BGV. It will always throw
   * helib::LogicError.
   * @param array Unused.
   * @param ptxt Unused.
   */
  void decode(UNUSED std::vector<NTL::ZZX>& array,
              UNUSED const NTL::ZZX& ptxt) const override
  {
    throw LogicError("Unimplemented: EncryptedArrayCx::decode for BGV type");
  }
  /**
   * @brief Unimplemented decode function for BGV. It will always throw
   * helib::LogicError.
   * @param array Unused.
   * @param ptxt Unused.
   */
  void decode(UNUSED PlaintextArray& array,
              UNUSED const NTL::ZZX& ptxt) const override
  {
    throw LogicError("Unimplemented: EncryptedArrayCx::decode for BGV type");
  }

  // random
  /**
   * @brief Unimplemented random function for BGV. It will always throw
   * helib::LogicError.
   * @param array Unused.
   */
  void random(UNUSED std::vector<NTL::ZZX>& array) const override
  {
    throw LogicError("Unimplemented: EncryptedArrayCx::decode for BGV type");
  }

  // decrypt
  /**
   * @brief Unimplemented decrypt function for BGV. It will always throw
   * helib::LogicError.
   * @param ctxt Unused.
   * @param sKey Unused.
   * @param ptxt Unused.
   */
  void decrypt(UNUSED const Ctxt& ctxt,
               UNUSED const SecKey& sKey,
               UNUSED std::vector<long>& ptxt) const override
  {
    throw LogicError("Unimplemented: EncryptedArrayCx::decrypt for BGV type");
  }
  /**
   * @brief Unimplemented decrypt function for BGV. It will always throw
   * helib::LogicError.
   * @param ctxt Unused.
   * @param sKey Unused.
   * @param ptxt Unused.
   */
  void decrypt(UNUSED const Ctxt& ctxt,
               UNUSED const SecKey& sKey,
               UNUSED std::vector<NTL::ZZX>& ptxt) const override
  {
    throw LogicError("Unimplemented: EncryptedArrayCx::decrypt for BGV type");
  }

  // buildLinPolyCoeffs
  /**
   * @brief Unimplemented buildLinPolyCoeffs function for BGV. It will always
   * throw helib::LogicError.
   * @param C Unused.
   * @param L Unused.
   */
  void buildLinPolyCoeffs(UNUSED std::vector<NTL::ZZX>& C,
                          UNUSED const std::vector<NTL::ZZX>& L) const override
  {
    throw LogicError("Unimplemented: "
                     "EncryptedArrayCx::buildLinPolyCoeffs for BGV type");
  }
  /* End BGV functions. */

  // These EaCx-specific encoding routines return the
  // scaling factor that was used in the encoding routine
  double encode(zzX& ptxt,
                const std::vector<cx_double>& array,
                double useThisSize,
                long precision = -1) const;
  double encode(zzX& ptxt,
                const std::vector<double>& array,
                double useThisSize,
                long precision = -1) const
  {
    std::vector<cx_double> tmp;
    convert(tmp, array);
    return encode(ptxt, tmp, useThisSize, precision);
  }
  double encode(zzX& ptxt,
                const std::vector<long>& array,
                double useThisSize,
                long precision = -1) const
  {
    std::vector<cx_double> tmp;
    convert(tmp, array);
    return encode(ptxt, tmp, useThisSize, precision);
  }

  /**
   * @brief Encode a `Ptxt` object into a `zzX`.
   * @tparam Scheme Encryption scheme to be used (either `BGV` or `CKKS`).
   * @param out Polynomial to encode into.
   * @param ptxt Plaintext `Ptxt` object to encode.
   * @param useThisSize Size to use.
   * @param precision Precision to use.
   * @return The scaling factor used in the encoding routine.
   **/
  template <typename Scheme>
  double encode(zzX& out,
                const Ptxt<Scheme>& ptxt,
                double useThisSize,
                long precision = -1) const
  {
    return encode(out, ptxt.getSlotRepr(), useThisSize, precision);
  }

  // VJS-FIXME: why do some encode functions have a
  // default value for useThisSize and others do not???
  // It's very confusing

  // VJS-FIXME: these routine have a number of issues and should
  // be deprecated in favor of the new EncodedPtxt-based routines

  /**
   * @brief Function for encoding a `double` into a `zzX`.
   * @deprecated This routine has a number of issues and is deprecated in favor
   * of the new `EncodedPtxt`-based routine.\n
   * Please use `PtxtArray::encode()` instead.
   **/
  double encode(zzX& ptxt,
                double aSingleNumber,
                double useThisSize = -1,
                long precision = -1) const;

  /**
   * @brief Function for encoding a `double` into a `zzX`.
   * @deprecated This routine has a number of issues and is deprecated in favor
   * of the new `EncodedPtxt`-based routine.\n
   * Please use `PtxtArray::encode()` instead.
   **/
  template <typename PTXT>
  double encode(NTL::ZZX& ptxt,
                const PTXT& pt,
                double useThisSize = -1,
                long precision = -1) const
  {
    zzX tmp;
    double f = encode(tmp, pt, useThisSize, precision);
    convert(ptxt, tmp);
    return f;
  }

  //========= new EncodedPtxt interface

  virtual void encode(UNUSED EncodedPtxt& eptxt,
                      UNUSED const std::vector<NTL::ZZX>& array) const override
  {
    throw LogicError("function not implemented for CKKS");
  }

  virtual void encode(UNUSED EncodedPtxt& eptxt,
                      UNUSED const std::vector<long>& array) const override
  {
    throw LogicError("function not implemented for CKKS");
  }

  virtual void encode(EncodedPtxt& eptxt,
                      const std::vector<cx_double>& array,
                      double mag = -1,
                      OptLong prec = OptLong()) const override;
  // implemented in EaCx.cpp

  virtual void encode(EncodedPtxt& eptxt,
                      const std::vector<double>& array,
                      double mag = -1,
                      OptLong prec = OptLong()) const override
  {
    std::vector<cx_double> array1;
    convert(array1, array);
    encode(eptxt, array1, mag, prec);
  }

  virtual void encode(EncodedPtxt& eptxt,
                      const PlaintextArray& array,
                      double mag = -1,
                      OptLong prec = OptLong()) const override;
  // implemented in EaCx.cpp

  virtual void encode(EncodedPtxt& eptxt,
                      const std::vector<bool>& array) const override
  {
    std::vector<cx_double> array1;
    convert(array1, array);
    encode(eptxt, array1);
  }

  virtual void encodeUnitSelector(EncodedPtxt& eptxt, long i) const override
  {
    std::vector<cx_double> array(this->size(), cx_double(0.0));
    array.at(i) = 1.0;
    encode(eptxt, array);
    // VJS-FIXME: this could be much more efficient
  }

  //==========================================

  void encryptOneNum(Ctxt& ctxt,
                     const PubKey& key,
                     double num,
                     double useThisSize = -1,
                     long precision = -1) const
  {
    assertEq(&getContext(),
             &ctxt.getContext(),
             "Cannot encrypt when ciphertext has different context than "
             "EncryptedArray");
    if (useThisSize <= 0.0)
      useThisSize = roundedSize(num); // rounded to power of two
    zzX pp;                           // Convert num into a plaintext polynomial
    double f = encode(pp, num, useThisSize, precision);

    key.CKKSencrypt(ctxt, pp, useThisSize, f); // encrypt resulting polynomial
  }

  template <typename PTXT>
  void encrypt(Ctxt& ctxt,
               const PubKey& key,
               const PTXT& ptxt,
               double useThisSize,
               long precision = -1) const
  {
    assertEq(&getContext(),
             &ctxt.getContext(),
             "Cannot encrypt when ciphertext has different context than "
             "EncryptedArray");
    zzX pp;
    double f = encode(pp, ptxt, useThisSize, precision);
    // Convert into a polynomial
    key.CKKSencrypt(ctxt, pp, useThisSize, f); // encrypt the polynomial
  }

  template <typename PTXT>
  void encrypt(Ctxt& ctxt, const PubKey& key, const PTXT& ptxt) const
  {
    encrypt(ctxt, key, ptxt, -1.0, -1);
  }

  // The methods below override EncryptedArrayBase, they use
  // the default size=0 and precision=0, which yield size=1
  // and precision=2^{-alMod.getR()-1}
  void encodeUnitSelector(zzX& ptxt, long i) const override
  {
    std::vector<cx_double> v(this->size(), cx_double(0.0));
    v.at(i) = cx_double(1.0, 0.0);
    encode(ptxt, v, /*size=*/1.0, /*default precision*/ -1);
    // VJS-FIXME: this could be much more efficient
  } // The implicit scaling factor is encodeScalingFactor() below

  // A bound on the rounding error for encoding
  double encodeRoundingError() const
  {
    const Context& context = getContext();
    long phim = context.getPhiM();

    // VJS-NOTE: I changed m to phi(m).
    // VJS-FIXME: for the power of two case, noiseBoundForUniform
    // is a bit too pessimistic, as this is the circularly symmetric
    // case.
    return context.noiseBoundForUniform(0.5, phim);
  }
  // The scaling factor to use when encoding/decoding plaintext elements
  long encodeScalingFactor(long precision = -1, double roundErr = -1.0) const
  {
    assertTrue<InvalidArgument>(precision < NTL_SP_BOUND,
                                "Precision exceeds max single precision bound");
    if (precision <= 0)
      precision = (1L << alMod.getR());
    if (roundErr < 0)
      roundErr = encodeRoundingError();

    // VJS-FIXME: the computation of f and/or return value could overflow

    long f = std::ceil(precision * roundErr);
    // We round the factor up to the next power of two
    return (1L << NTL::NextPowerOfTwo(f));
  }

  virtual double defaultErr() const override
  {
    const Context& context = getContext();
    long phim = context.getPhiM();

    // VJS-FIXME: For the power of two case, noiseBoundForUniform
    // is a bit too pessimistic, as this is the circularly symmetric
    // case.

    // VJS-NOTE: this is extremey heuristic: we are assuming
    // that the coefficients mod 1 are modeled as uniform
    // on [0,1].

    return context.noiseBoundForUniform(0.5, phim);
  }

  virtual double defaultScale(double err,
                              OptLong prec = OptLong()) const override
  {
    if (err < 1.0)
      err = 1.0;
    long r = alMod.getR(); // default r-value
    if (prec.isDefined())
      r = prec; // override if necessary

    // we want to compute
    //   2^(ceil(log2(err*2^r))) = 2^(ceil(log2(err) + r))
    //                           = 2^(r + ceil(log2(err)))
    int e;
    std::frexp(1 / err, &e);
    // we have 2^{e-1} <= 1/err < 2^e, i.e.,
    //         2^{-e} < err <= 2^{-e+1}
    // so ceil(log2(err)) = -e+1
    return std::ldexp(1.0, r - e + 1);
  }

  void decode(std::vector<cx_double>& array,
              const zzX& ptxt,
              double scaling) const;

  void decode(std::vector<cx_double>& array,
              const NTL::ZZX& ptxt,
              double scaling) const
  {
    zzX tmp;
    convert(tmp, ptxt);
    decode(array, tmp, scaling);
  }

  void decode(std::vector<double>& array, const zzX& ptxt, double scaling) const
  {
    std::vector<cx_double> v;
    decode(v, ptxt, scaling);
    project(array, v);
  }

  void decode(std::vector<double>& array,
              const NTL::ZZX& ptxt,
              double scaling) const
  {
    std::vector<cx_double> v;
    decode(v, ptxt, scaling);
    project(array, v);
  }

  /**
   * @deprecated This routine has a number of issues and is deprecated in favor
   * of either `RandomComplex()` or `PtxtArray::random()`.
   **/
  void random(std::vector<cx_double>& array, double rad = 1.0) const;
  void random(std::vector<double>& array, double rad = 1.0) const
  {
    std::vector<cx_double> v;
    random(v, rad);
    project(array, v);
  }
  void random(std::vector<long>& array) const override
  {
    std::vector<cx_double> v;
    random(v, 1.0);
    project_and_round(array, v);
  }

  void decrypt(const Ctxt& ctxt,
               const SecKey& sKey,
               std::vector<cx_double>& ptxt,
               OptLong prec = OptLong()) const override;

  void decrypt(const Ctxt& ctxt,
               const SecKey& sKey,
               std::vector<double>& ptxt,
               OptLong prec = OptLong()) const override;

  void rawDecrypt(const Ctxt& ctxt,
                  const SecKey& sKey,
                  std::vector<cx_double>& ptxt) const override;

  void rawDecrypt(const Ctxt& ctxt,
                  const SecKey& sKey,
                  std::vector<double>& ptxt) const override;

  void decrypt(const Ctxt& ctxt,
               const SecKey& sKey,
               PlaintextArray& ptxt,
               OptLong prec = OptLong()) const override
  {
    decryptReal(ctxt, sKey, ptxt, prec);
  }

  void rawDecrypt(const Ctxt& ctxt,
                  const SecKey& sKey,
                  PlaintextArray& ptxt) const override
  {
    rawDecryptReal(ctxt, sKey, ptxt);
  }

  void decryptComplex(const Ctxt& ctxt,
                      const SecKey& sKey,
                      PlaintextArray& ptxt,
                      OptLong prec = OptLong()) const override;

  void rawDecryptComplex(const Ctxt& ctxt,
                         const SecKey& sKey,
                         PlaintextArray& ptxt) const override;

  void decryptReal(const Ctxt& ctxt,
                   const SecKey& sKey,
                   PlaintextArray& ptxt,
                   OptLong prec = OptLong()) const override;

  void rawDecryptReal(const Ctxt& ctxt,
                      const SecKey& sKey,
                      PlaintextArray& ptxt) const override;

  /**
   * @brief Decrypt ciphertext to a plaintext relative to a specific scheme.
   * @tparam Scheme Encryption scheme to be used (either `BGV` or `CKKS`).
   * @param ctxt Ciphertext to decrypt.
   * @param sKey Secret key to be used for decryption.
   * @param ptxt Plaintext into which to decrypt.
   * @param prec `CKKS` precision to be used (must be defaulted if Scheme is
   *`BGV`). Decrypt a `Ctxt` ciphertext object to a `Ptxt` plaintext one
   *relative to a specific scheme.
   **/

  // VJS-FIXME: something seems odd here. This code clearly only
  // works for CKKS because ptxtArray has type vector<cx_double>.
  // However, Scheme could be CKKS or BGV.  Is there a specialized
  // version of this somewhere?
  // TODO: document this better (especially the prec parameter)
  template <typename Scheme>
  void decrypt(const Ctxt& ctxt,
               const SecKey& sKey,
               Ptxt<Scheme>& ptxt,
               OptLong prec = OptLong()) const
  {
    std::vector<cx_double> ptxtArray;
    decrypt(ctxt, sKey, ptxtArray, prec);
    ptxt.setData(std::move(ptxtArray));
  }

  void extractRealPart(Ctxt& c) const;

  /**
   * @brief Extract the real part of a `CKKS` plaintext.
   * @tparam Scheme Encryption scheme to be used (must be `CKKS`).
   * @param p Plaintext on which to operate.
   **/
  template <typename Scheme>
  void extractRealPart(Ptxt<Scheme>& p) const
  {
    p = p.real();
  }

  /**
   * @brief Extract the imaginary part of a `CKKS` plaintext.
   * @tparam Scheme Encryption scheme to be used (must be `CKKS`).
   * @param p Plaintext on which to operate.
   **/
  template <typename Scheme>
  void extractImPart(Ptxt<Scheme>& p) const
  {
    p = p.imag();
  }

  //! Note: If called with dcrt==nullptr, extractImPart will perform FFT's
  //! when encoding i as a DoubleCRT object. If called with dcrt!=nullptr,
  //! it assumes that dcrt points to an object that encodes i.
  void extractImPart(Ctxt& c, DoubleCRT* dcrt = nullptr) const;

  //! @name Linearized polynomials for EncryptedArrayCx
  ///@{
  //! buildLinPolyCoeffs returns in C two encoded constants such that the
  //! linear transformation(s) defined as L(1) = oneImage and L(i)=iImage
  //! can be computed as:      L(x) = C[0]*x + C[1]*conjugate(x).
  //! Once C is computed, we can apply this L to a ciphertext by calling
  //! applyLinPolyLL(ctxt, C, 2).
  //! Alternatively, we can convert C to a vector of two DoubleCRT objects,
  //! then call applyLinPolyLL(ctxt, dcrtVec, 2). This lets us compute the
  //! DoubleCRT object just once, then use them many times.

  //! First variant: same linear transformation in all the slots
  void buildLinPolyCoeffs(std::vector<zzX>& C,
                          const cx_double& oneImage,
                          const cx_double& iImage,
                          long precision = 0) const;

  //! Second variant: different linear transformation in each slots
  void buildLinPolyCoeffs(std::vector<zzX>& C,
                          const std::vector<cx_double>& oneImages,
                          const std::vector<cx_double>& iImages,
                          long precision = 0) const;
  ///@}
};

typedef EncryptedArrayDerived<PA_cx> EncryptedArrayCx;

// plaintextAutomorph: Compute b(X) = a(X^k) mod Phi_m(X).
template <typename RX, typename RXModulus>
void plaintextAutomorph(RX& bb,
                        const RX& a,
                        long k,
                        long m,
                        const RXModulus& PhimX)
{
  // compute b(X) = a(X^k) mod (X^m-1)
  if (k == 1 || deg(a) <= 0) {
    bb = a;
    return;
  }

  RX b;
  b.SetLength(m);
  NTL::mulmod_precon_t precon = NTL::PrepMulModPrecon(k, m);
  for (long j = 0; j <= deg(a); j++)
    b[NTL::MulModPrecon(j, k, m, precon)] = a[j]; // b[j*k mod m] = a[j]
  b.normalize();

  rem(bb, b, PhimX); // reduce modulo the m'th cyclotomic
}

// same as above, but k = g_i^j mod m.
// also works with i == ea.getPalgebra().numOfGens(),
// which means Frobenius

template <typename RX, typename type>
void plaintextAutomorph(RX& b,
                        const RX& a,
                        long i,
                        long j,
                        const EncryptedArrayDerived<type>& ea)
{
  const PAlgebra& zMStar = ea.getPAlgebra();
  const auto& F = ea.getTab().getPhimXMod();
  long k = zMStar.genToPow(i, j);
  long m = zMStar.getM();
  plaintextAutomorph(b, a, k, m, F);
}

//! @brief A "factory" for building EncryptedArrays
EncryptedArrayBase* buildEncryptedArray(const Context& context,
                                        const PAlgebraMod& alMod,
                                        const NTL::ZZX& G = NTL::ZZX::zero());

//! @class EncryptedArray
//! @brief A simple wrapper for a smart pointer to an EncryptedArrayBase.
//! This is the interface that higher-level code should use
class EncryptedArray
{
private:
  const PAlgebraMod& alMod;
  ClonedPtr<EncryptedArrayBase> rep;

public:
  //! constructor: G defaults to the monomial X, PAlgebraMod from context
  EncryptedArray(const Context& context, const NTL::ZZX& G = NTL::ZZX(1, 1)) :
      alMod(context.getAlMod()),
      rep(buildEncryptedArray(context, context.getAlMod(), G))
  {}
  //! constructor: G defaults to F0, PAlgebraMod explicitly given
  EncryptedArray(const Context& context, const PAlgebraMod& _alMod) :
      alMod(_alMod), rep(buildEncryptedArray(context, _alMod))
  {}

  // NOTES:
  //  (1) the second constructor is provided mainly for BGV bootstrapping
  //  (2) we do not currently provide a constructor that allows
  //      the user to select both G and alMod, but this could be added

  // copy constructor:

#if 1
  EncryptedArray& operator=(const EncryptedArray& other) = delete;
#else
  EncryptedArray& operator=(const EncryptedArray& other)
  {
    if (this == &other)
      return *this;
    assertEq(&alMod,
             &other.alMod,
             "Cannot assign EncryptedArrays with different algebras");
    rep = other.rep;
    return *this;
  }
#endif

  //! @brief downcast operator
  //! example: const EncryptedArrayDerived<PA_GF2>& rep =
  //! ea.getDerived(PA_GF2());
  template <typename type>
  const EncryptedArrayDerived<type>& getDerived(type) const
  {
    return dynamic_cast<const EncryptedArrayDerived<type>&>(*rep);
  }

  const EncryptedArrayCx& getCx() const
  {
    return dynamic_cast<const EncryptedArrayCx&>(*rep);
  }

  ///@{
  //! @name Direct access to EncryptedArrayBase methods

  PA_tag getTag() const { return rep->getTag(); }
  bool isCKKS() const { return getTag() == PA_cx_tag; }

  template <template <typename> class T, typename... Args>
  void dispatch(Args&&... args) const
  {
    switch (getTag()) {
    case PA_GF2_tag: {
      const EncryptedArrayDerived<PA_GF2>* p =
          static_cast<const EncryptedArrayDerived<PA_GF2>*>(rep.get());
      p->dispatch<T>(std::forward<Args>(args)...);
      break;
    }
    case PA_zz_p_tag: {
      const EncryptedArrayDerived<PA_zz_p>* p =
          static_cast<const EncryptedArrayDerived<PA_zz_p>*>(rep.get());
      p->dispatch<T>(std::forward<Args>(args)...);
      break;
    }
    case PA_cx_tag: {
      const EncryptedArrayDerived<PA_cx>* p =
          static_cast<const EncryptedArrayDerived<PA_cx>*>(rep.get());
      p->dispatch<T>(std::forward<Args>(args)...);
      break;
    }
    default:
      throw RuntimeError("EncryptedArray: bad tag");
    }
  }

  const Context& getContext() const { return rep->getContext(); }
  const PAlgebraMod& getAlMod() const { return alMod; }
  const PAlgebra& getPAlgebra() const { return rep->getPAlgebra(); }
  long getDegree() const { return rep->getDegree(); }
  void rotate(Ctxt& ctxt, long k) const { rep->rotate(ctxt, k); }
  void shift(Ctxt& ctxt, long k) const { rep->shift(ctxt, k); }
  void rotate1D(Ctxt& ctxt, long i, long k, bool dc = false) const
  {
    rep->rotate1D(ctxt, i, k, dc);
  }
  void shift1D(Ctxt& ctxt, long i, long k) const { rep->shift1D(ctxt, i, k); }

  void encode(zzX& ptxt, const std::vector<long>& array) const
  {
    rep->encode(ptxt, array);
  }

  void encode(NTL::ZZX& ptxt, const std::vector<long>& array) const
  {
    rep->encode(ptxt, array);
  }

  void encode(zzX& ptxt, const std::vector<zzX>& array) const
  {
    rep->encode(ptxt, array);
  }

  void encode(zzX& ptxt, const PlaintextArray& array) const
  {
    rep->encode(ptxt, array);
  }

  void encode(NTL::ZZX& ptxt, const std::vector<NTL::ZZX>& array) const
  {
    rep->encode(ptxt, array);
  }

  void encode(NTL::ZZX& ptxt, const PlaintextArray& array) const
  {
    rep->encode(ptxt, array);
  }

  void encode(zzX& ptxt, const std::vector<NTL::ZZX>& array) const
  {
    rep->encode(ptxt, array);
  }

  //=============== new EncodedPtxt interfaces

  void encode(EncodedPtxt& eptxt, const std::vector<NTL::ZZX>& array) const
  {
    rep->encode(eptxt, array);
  }

  void encode(EncodedPtxt& eptxt, const std::vector<long>& array) const
  {
    rep->encode(eptxt, array);
  }

  void encode(EncodedPtxt& eptxt,
              const std::vector<cx_double>& array,
              double mag = -1,
              OptLong prec = OptLong()) const
  {
    rep->encode(eptxt, array, mag, prec);
  }

  void encode(EncodedPtxt& eptxt,
              const std::vector<double>& array,
              double mag = -1,
              OptLong prec = OptLong()) const
  {
    rep->encode(eptxt, array, mag, prec);
  }

  void encode(EncodedPtxt& eptxt,
              const PlaintextArray& array,
              double mag = -1,
              OptLong prec = OptLong()) const
  {
    rep->encode(eptxt, array, mag, prec);
  }

  void encode(EncodedPtxt& eptxt, const std::vector<bool>& array) const
  {
    rep->encode(eptxt, array);
  }

  void encodeUnitSelector(EncodedPtxt& eptxt, long i) const
  {
    rep->encodeUnitSelector(eptxt, i);
  }

  //================================

  void encodeUnitSelector(zzX& ptxt, long i) const
  {
    rep->encodeUnitSelector(ptxt, i);
  }

  template <typename PTXT, typename ARRAY>
  void decode(ARRAY& array, const PTXT& ptxt) const
  {
    rep->decode(array, ptxt);
  }

  template <typename T>
  void random(std::vector<T>& array) const
  {
    rep->random(array);
  }

  //========= new encryption interfaces
  // NOTE: I provide both a "legacy" interface that includes
  // Ctxt and PubKey, and a cleaner interface that includes
  // only Ctxt.

  // BGV only
  void encrypt(Ctxt& ctxt,
               const PubKey& key,
               const std::vector<NTL::ZZX>& array) const
  {
    EncodedPtxt eptxt;
    encode(eptxt, array);
    key.Encrypt(ctxt, eptxt);
  }

  void encrypt(Ctxt& ctxt, const std::vector<NTL::ZZX>& array) const
  {
    encrypt(ctxt, ctxt.getPubKey(), array);
  }

  void encrypt(Ctxt& ctxt,
               const PubKey& key,
               const std::vector<long>& array) const
  {
    EncodedPtxt eptxt;
    encode(eptxt, array);
    key.Encrypt(ctxt, eptxt);
  }

  void encrypt(Ctxt& ctxt, const std::vector<long>& array) const
  {
    encrypt(ctxt, ctxt.getPubKey(), array);
  }

  // CKKS only
  void encrypt(Ctxt& ctxt,
               const PubKey& key,
               const std::vector<cx_double>& array,
               double mag,
               OptLong prec = OptLong()) const
  {
    if (mag < 0)
      throw LogicError("CKKS encryption: mag must be set to non-default");
    EncodedPtxt eptxt;
    encode(eptxt, array, mag, prec);
    key.Encrypt(ctxt, eptxt);
  }

  void encrypt(Ctxt& ctxt,
               const std::vector<cx_double>& array,
               UNUSED double mag,
               OptLong prec = OptLong()) const
  {
    encrypt(ctxt, ctxt.getPubKey(), array, prec);
  }

  void encrypt(Ctxt& ctxt,
               const PubKey& key,
               const std::vector<double>& array,
               double mag,
               OptLong prec = OptLong()) const
  {
    if (mag < 0)
      throw LogicError("CKKS encryption: mag must be set to non-default");
    EncodedPtxt eptxt;
    encode(eptxt, array, mag, prec);
    key.Encrypt(ctxt, eptxt);
  }

  void encrypt(Ctxt& ctxt,
               const std::vector<double>& array,
               double mag,
               OptLong prec = OptLong()) const
  {
    encrypt(ctxt, ctxt.getPubKey(), array, mag, prec);
  }

  // BGV and CKKS
  void encrypt(Ctxt& ctxt,
               const PubKey& key,
               const PlaintextArray& array,
               double mag = -1,
               OptLong prec = OptLong()) const
  // NOTES: (1) for BGV, mag,prec must be defaulted;
  // (2) for CKKS, mag must be set to non-default value
  {
    if (getTag() == PA_cx_tag && mag < 0)
      throw LogicError("CKKS encryption: mag must be set to non-default");
    EncodedPtxt eptxt;
    encode(eptxt, array, mag, prec);
    key.Encrypt(ctxt, eptxt);
  }

  void encrypt(Ctxt& ctxt,
               const PlaintextArray& array,
               double mag = -1,
               OptLong prec = OptLong()) const
  {
    encrypt(ctxt, ctxt.getPubKey(), array, mag, prec);
  }

  //=========================

  template <typename T>
  void decrypt(const Ctxt& ctxt, const SecKey& sKey, T& ptxt) const
  {
    rep->decrypt(ctxt, sKey, ptxt);
  }

  template <typename T>
  void decrypt(const Ctxt& ctxt,
               const SecKey& sKey,
               T& ptxt,
               OptLong prec) const
  {
    rep->decrypt(ctxt, sKey, ptxt, prec);
  }

  template <typename T>
  void rawDecrypt(const Ctxt& ctxt, const SecKey& sKey, T& ptxt) const
  {
    rep->rawDecrypt(ctxt, sKey, ptxt);
  }

  void decryptComplex(const Ctxt& ctxt,
                      const SecKey& sKey,
                      PlaintextArray& ptxt,
                      OptLong prec = OptLong()) const
  {
    rep->decryptComplex(ctxt, sKey, ptxt, prec);
  }

  void rawDecryptComplex(const Ctxt& ctxt,
                         const SecKey& sKey,
                         PlaintextArray& ptxt) const
  {
    rep->rawDecryptComplex(ctxt, sKey, ptxt);
  }

  void decryptReal(const Ctxt& ctxt,
                   const SecKey& sKey,
                   PlaintextArray& ptxt,
                   OptLong prec = OptLong()) const
  {
    rep->decryptReal(ctxt, sKey, ptxt, prec);
  }

  void rawDecryptReal(const Ctxt& ctxt,
                      const SecKey& sKey,
                      PlaintextArray& ptxt) const
  {
    rep->rawDecryptReal(ctxt, sKey, ptxt);
  }

  void buildLinPolyCoeffs(std::vector<NTL::ZZX>& C,
                          const std::vector<NTL::ZZX>& L) const
  {
    rep->buildLinPolyCoeffs(C, L);
  }

  void restoreContext() const { rep->restoreContext(); }
  void restoreContextForG() const { rep->restoreContextForG(); }

  long size() const { return rep->size(); }
  long dimension() const { return rep->dimension(); }
  long sizeOfDimension(long i) const { return rep->sizeOfDimension(i); }
  long nativeDimension(long i) const { return rep->nativeDimension(i); }
  long coordinate(long i, long k) const { return rep->coordinate(i, k); }
  long addCoord(long i, long k, long offset) const
  {
    return rep->addCoord(i, k, offset);
  }

  //! @brief rotate an array by offset in the i'th dimension
  //! (output should not alias input)
  template <typename U>
  void rotate1D(std::vector<U>& out,
                const std::vector<U>& in,
                long i,
                long offset) const
  {
    rep->rotate1D(out, in, i, offset);
  }
  ///@}
};

// Convenience routines to avoid EncryptedArray's when
// using the "default" one in Context

inline void rotate(Ctxt& ctxt, long k)
{
  ctxt.getContext().getView().rotate(ctxt, k);
}

inline void shift(Ctxt& ctxt, long k)
{
  ctxt.getContext().getView().shift(ctxt, k);
}

inline void rotate1D(Ctxt& ctxt, long i, long k, bool dc = false)
{
  ctxt.getContext().getView().rotate1D(ctxt, i, k, dc);
}

inline void shift1D(Ctxt& ctxt, long i, long k)
{
  ctxt.getContext().getView().shift1D(ctxt, i, k);
}

// PlaintextArray

class PlaintextArrayBase
{ // purely abstract interface
public:
  virtual ~PlaintextArrayBase() {}
  virtual void print(std::ostream& s) const = 0;
};

template <typename type>
class PlaintextArrayDerived : public PlaintextArrayBase
{
public:
  PA_INJECT(type)

  std::vector<RX> data;

  virtual void print(std::ostream& s) const { s << data; }
};

/**
 * @deprecated There is a somewhat "friendlier" interface available in
 * `PtxtArray` as it carries with it a reference to an `EncryptedArray`.
 * It is recommended that `PlaintextArray` be deprecated in favor of
 * `PtxtArray`.
 **/
class PlaintextArray
{
private:
  NTL::CloneablePtr<PlaintextArrayBase> rep;

  template <typename type>
  class ConstructorImpl
  {
  public:
    PA_INJECT(type)

    static void apply(const EncryptedArrayDerived<type>& ea, PlaintextArray& pa)
    {
      NTL::CloneablePtr<PlaintextArrayDerived<type>> ptr =
          NTL::MakeCloneable<PlaintextArrayDerived<type>>();
      ptr->data.resize(ea.size());
      pa.rep = ptr;
    }
  };

public:
  PlaintextArray(const EncryptedArray& ea)
  {
    ea.dispatch<ConstructorImpl>(*this);
  }

  PlaintextArray(const PlaintextArray& other) : rep(other.rep.clone()) {}
  PlaintextArray& operator=(const PlaintextArray& other)
  {
    rep = other.rep.clone();
    return *this;
  }

  template <typename type>
  std::vector<typename type::RX>& getData()
  {
    return (dynamic_cast<PlaintextArrayDerived<type>&>(*rep)).data;
  }

  template <typename type>
  const std::vector<typename type::RX>& getData() const
  {
    return (dynamic_cast<PlaintextArrayDerived<type>&>(*rep)).data;
  }

  void print(std::ostream& s) const { rep->print(s); }
};

inline std::ostream& operator<<(std::ostream& s, const PlaintextArray& pa)
{
  pa.print(s);
  return s;
}

void rotate(const EncryptedArray& ea, PlaintextArray& pa, long k);
void shift(const EncryptedArray& ea, PlaintextArray& pa, long k);

void rotate1D(const EncryptedArray& ea, PlaintextArray& pa, long i, long k);
void shift1D(const EncryptedArray& ea, PlaintextArray& pa, long i, long k);

void encode(const EncryptedArray& ea,
            PlaintextArray& pa,
            const std::vector<long>& array);
void encode(const EncryptedArray& ea,
            PlaintextArray& pa,
            const std::vector<NTL::ZZX>& array);
void encode(const EncryptedArray& ea,
            PlaintextArray& pa,
            const std::vector<cx_double>& array);
void encode(const EncryptedArray& ea,
            PlaintextArray& pa,
            const std::vector<double>& array);

void encode(const EncryptedArray& ea, PlaintextArray& pa, long val);
void encode(const EncryptedArray& ea, PlaintextArray& pa, const NTL::ZZX& val);
void encode(const EncryptedArray& ea, PlaintextArray& pa, double val);
void encode(const EncryptedArray& ea, PlaintextArray& pa, cx_double val);

void decode(const EncryptedArray& ea,
            std::vector<long>& array,
            const PlaintextArray& pa);
void decode(const EncryptedArray& ea,
            std::vector<NTL::ZZX>& array,
            const PlaintextArray& pa);
void decode(const EncryptedArray& ea,
            std::vector<cx_double>& array,
            const PlaintextArray& pa);
void decode(const EncryptedArray& ea,
            std::vector<double>& array,
            const PlaintextArray& pa);

void random(const EncryptedArray& ea, PlaintextArray& pa);
void randomReal(const EncryptedArray& ea, PlaintextArray& pa);
void randomComplex(const EncryptedArray& ea, PlaintextArray& pa);

bool equals(const EncryptedArray& ea,
            const PlaintextArray& pa,
            const PlaintextArray& other);
bool equals(const EncryptedArray& ea,
            const PlaintextArray& pa,
            const std::vector<long>& other);
bool equals(const EncryptedArray& ea,
            const PlaintextArray& pa,
            const std::vector<NTL::ZZX>& other);

bool equals(const EncryptedArray& ea,
            const PlaintextArray& pa,
            const PlaintextArray& other,
            double tolerance,
            double floor);

void add(const EncryptedArray& ea,
         PlaintextArray& pa,
         const PlaintextArray& other);
void sub(const EncryptedArray& ea,
         PlaintextArray& pa,
         const PlaintextArray& other);
void mul(const EncryptedArray& ea,
         PlaintextArray& pa,
         const PlaintextArray& other);
void negate(const EncryptedArray& ea, PlaintextArray& pa);

void frobeniusAutomorph(const EncryptedArray& ea, PlaintextArray& pa, long j);
void frobeniusAutomorph(const EncryptedArray& ea,
                        PlaintextArray& pa,
                        const NTL::Vec<long>& vec);

void extractRealPart(const EncryptedArray& ea, PlaintextArray& pa);
void extractImPart(const EncryptedArray& ea, PlaintextArray& pa);

void applyPerm(const EncryptedArray& ea,
               PlaintextArray& pa,
               const NTL::Vec<long>& pi);

void power(const EncryptedArray& ea, PlaintextArray& pa, long e);

double Norm(const EncryptedArray& ea, const PlaintextArray& pa);
double Distance(const EncryptedArray& ea,
                const PlaintextArray& pa,
                const PlaintextArray& other);

void totalSums(const EncryptedArray& ea, PlaintextArray& pa);
void runningSums(const EncryptedArray& ea, PlaintextArray& pa);

//=====================================

// PtxtArray is a somewhat "friendlier" interface than
// PlaintextArray, as it carries with it a reference to
// an EncryptedArray.  It is recommended that PlaintextArray
// is deprecated in favor of PtxtArray.

class PtxtArray
{
public:
  /**
   * @brief Class label to be added to JSON serialization as object type
   * information.
   */
  static constexpr std::string_view typeName = "PtxtArray";

  // These two data fields should really be private, but there are
  // a lot of internal functions that need to access them
  const EncryptedArray& ea;
  PlaintextArray pa;

  explicit PtxtArray(const EncryptedArray& ea_) : ea(ea_), pa(ea) {}

  explicit PtxtArray(const Context& context) : ea(context.getView()), pa(ea) {}

  // copy constructor: default
  PtxtArray(const PtxtArray&) = default;

  // templates that allows construction via convert:
  // T can be any type supported by convert(PtxtArray,T)
  template <class T>
  PtxtArray(const EncryptedArray& ea, const T& t) : PtxtArray(ea)
  {
    load(t);
  }

  template <class T>
  PtxtArray(const Context& context, const T& t) : PtxtArray(context)
  {
    load(t);
  }

  PtxtArray& operator=(const PtxtArray& other)
  {
    assertTrue(&ea == &other.ea, "PtxtArray: inconsistent assignment");
    pa = other.pa;
    return *this;
  }

  // template that allows assignment via convert:
  // T can be any type supported by PxtArray::load(T)
  template <class T>
  PtxtArray& operator=(const T& t)
  {
    load(t);
    return *this;
  }

  const EncryptedArray& getView() const { return ea; }
  const EncryptedArray& getEA() const { return ea; }

  long size() const { return ea.size(); }

  // direct encode, encrypt, and decrypt methods
  void encode(EncodedPtxt& eptxt,
              double mag = -1,
              OptLong prec = OptLong()) const
  {
    if (ea.isCKKS())
      ea.encode(eptxt, pa, mag, prec);
    else
      ea.encode(eptxt, pa); // ignore mag,prec for BGV
  }

  void encrypt(Ctxt& ctxt, double mag = -1, OptLong prec = OptLong()) const
  {
    if (ea.isCKKS()) {
      if (mag < 0)
        mag = NextPow2(Norm(pa.getData<PA_cx>()));
      // if mag is defaulted, set it to 2^(ceil(log2(max(Norm(pa),1))))
      ea.encrypt(ctxt, pa, mag, prec);
    } else {
      ea.encrypt(ctxt, pa); // ignore mag,prec for BGV
    }
  }

  void decrypt(const Ctxt& ctxt, const SecKey& sKey, OptLong prec = OptLong())
  {
    if (ea.isCKKS())
      ea.decrypt(ctxt, sKey, pa, prec);
    else
      ea.decrypt(ctxt, sKey, pa);
  }

  void rawDecrypt(const Ctxt& ctxt, const SecKey& sKey)
  {
    ea.rawDecrypt(ctxt, sKey, pa);
  }

  void decryptComplex(const Ctxt& ctxt,
                      const SecKey& sKey,
                      OptLong prec = OptLong())
  {
    ea.decryptComplex(ctxt, sKey, pa, prec);
  }

  void rawDecryptComplex(const Ctxt& ctxt, const SecKey& sKey)
  {
    ea.rawDecryptComplex(ctxt, sKey, pa);
  }

  void decryptReal(const Ctxt& ctxt,
                   const SecKey& sKey,
                   OptLong prec = OptLong())
  {
    ea.decryptReal(ctxt, sKey, pa, prec);
  }

  void rawDecryptReal(const Ctxt& ctxt, const SecKey& sKey)
  {
    ea.rawDecryptReal(ctxt, sKey, pa);
  }

  void randomReal() { helib::randomReal(ea, pa); }

  void randomComplex() { helib::randomComplex(ea, pa); }

  void random() { helib::random(ea, pa); }

  //======== load ========
  // Puts vector or scalar data into a PtxtArray

  void load(const std::vector<int>& array)
  {
    std::vector<long> array1;
    convert(array1, array);
    helib::encode(ea, pa, array1);
  }

  void load(const std::vector<long>& array) { helib::encode(ea, pa, array); }

  void load(const std::vector<NTL::ZZX>& array)
  {
    helib::encode(ea, pa, array);
  }

  void load(const std::vector<cx_double>& array)
  {
    helib::encode(ea, pa, array);
  }

  void load(const std::vector<double>& array) { helib::encode(ea, pa, array); }

  void load(int val) { helib::encode(ea, pa, long(val)); }

  void load(long val) { helib::encode(ea, pa, val); }

  void load(const NTL::ZZX& val) { helib::encode(ea, pa, val); }

  void load(double val) { helib::encode(ea, pa, val); }

  void load(cx_double val) { helib::encode(ea, pa, val); }

  // additional convenience conversions from NTL types (BGV only)
  // NTL vectors of NTL ring types
  // Implementation note: these all go via conversion to
  // std::vector<NTL::ZZX> to enforce BGV.  This is not
  // the most efficient way to do this, and if it
  // becomes a bottleneck, we can revisit this implementation
  // (without changing semantics).
  void load(const NTL::Vec<NTL::GF2>& vec)
  {
    std::vector<NTL::ZZX> v;
    convert(v, vec);
    load(v);
  }
  void load(const NTL::Vec<NTL::GF2X>& vec)
  {
    std::vector<NTL::ZZX> v;
    convert(v, vec);
    load(v);
  }
  void load(const NTL::Vec<NTL::zz_p>& vec)
  {
    std::vector<NTL::ZZX> v;
    convert(v, vec);
    load(v);
  }
  void load(const NTL::Vec<NTL::zz_pX>& vec)
  {
    std::vector<NTL::ZZX> v;
    convert(v, vec);
    load(v);
  }

  // NTL scalar ring types
  // Implementation note: these all go via conversion to
  // NTL::ZZX to enforce BGV
  void load(NTL::GF2 scalar)
  {
    NTL::ZZX s;
    convert(s, scalar);
    load(s);
  }
  void load(const NTL::GF2X& scalar)
  {
    NTL::ZZX s;
    convert(s, scalar);
    load(s);
  }
  void load(NTL::zz_p scalar)
  {
    NTL::ZZX s;
    convert(s, scalar);
    load(s);
  }
  void load(const NTL::zz_pX& scalar)
  {
    NTL::ZZX s;
    convert(s, scalar);
    load(s);
  }

  //============== store ============
  // Puts data into a std::vector aka `unload`

  void store(std::vector<long>& array) const { decode(ea, array, pa); }

  void store(std::vector<NTL::ZZX>& array) const { decode(ea, array, pa); }

  void store(std::vector<cx_double>& array) const { decode(ea, array, pa); }

  void store(std::vector<double>& array) const { decode(ea, array, pa); }

  //===============================

  // this is here for consistency with Ctxt class
  void negate() { helib::negate(ea, pa); }

  /**
   * @brief Function to serialize `this` `PtxtArray`.
   * @param os Output `std::ostream`.
   * @note `PtxtArray` `context` is not serialized, see note of `readJSON`.
   *
   * The output stream will be a JSON where the `PtxtArray` content will be
   * serialized in the `slots` field.\n
   * Each slot of `PtxtArray` will be serialized in an element of such list by
   * the JSON serializer function determined by the scheme.\n
   * For example if we have a plaintext `pa` such that `pa[0]=slot0`,
   * `pa[1]=slot1`, `pa[2]=slot2`, and `pa[i]=0` for `i>2`, it will be
   * serialized as '['slot0', 'slot1', 'slot2', `0`, `0` ...]'.
   **/
  void writeToJSON(std::ostream& os) const;

  /**
   * @brief Function to serialize `this` `PtxtArray`.
   * @return The `JsonWrapper` containing the serialization.
   * @note `PtxtArray` `context` is not serialized, see note of `readJSON`.
   *
   * The output JsonWrapper will be a JSON where the `PtxtArray` content will
   * be serialized in the `slots` field.\n
   * Each slot of `PtxtArray` will be serialized in an element of such list by
   * the JSON serializer function determined by the scheme.\n
   * For example if we have a plaintext `pa` such that `pa[0]=slot0`,
   * `pa[1]=slot1`, `pa[2]=slot2`, and `pa[i]=0` for `i>2`, it will be
   * serialized as '['slot0', 'slot1', 'slot2', `0`, `0`, ...]'.
   **/
  JsonWrapper writeToJSON() const;

  /**
   * @brief Function to deserialize and return a `PtxtArray` from a JSON
   * stream.
   * @param is Input `std::istream`.
   * @throws IOError if the stream is badly formatted (i.e. it does not contain
   * a valid JSON).
   * @code
   * PtxtArray my_pa = PtxtArray::readFromJSON(std::cin, context);
   * @endcode
   *
   * The input stream has to contain a valid typed JSON value.\n
   * Each element of the content list will be deserialized as a slot of the type
   * determined by the scheme.\n
   * If the number of tokens in the slot list is less than the number of slots,
   * the remaining slots will be padded by 0.\n
   * For example a slot list '['slot0', 'slot1', 'slot2']' will be deserialized
   * as a plaintext `pa` where `pa[0]=slot0`, `pa[1]=slot1`, `pa[2]=slot2`, and
   * `pa[i]=0` for `i>2`.
   **/
  static PtxtArray readFromJSON(std::istream& is, const Context& context);

  /**
   * @brief Function to deserialize and return a `PtxtArray` from a
   * `JsonWrapper` object.
   * @param jw `JsonWrapper` containing the serialized object.
   * @throws IOError if the `JsonWrapper` object does not contains a valid
   * serialization of a `PtxtArray`.
   * @code
   * PtxtArray my_pa = PtxtArray::readFromJSON(..., context);
   * @endcode
   *
   * The `JsonWrapper` must contain a valid `PtxtArray` serialization.\n
   * Each element of the content list will be deserialized as a slot of the
   * type determined by the scheme.\n
   * If the number of tokens in the slot list is less than the number of slots,
   * the remaining slots will be padded by 0.\n
   * For example a slot list '['slot0', 'slot1', 'slot2']' will be deserialized
   * as a plaintext `pa` where `pa[0]=slot0`, `pa[1]=slot1`, `pa[2]=slot2` and
   * `pa[i]=0` for `i>2`.
   **/
  static PtxtArray readFromJSON(const JsonWrapper& jw, const Context& context);

  /**
   * @brief In-place function to deserialize a `PtxtArray` from a JSON stream.
   * @param is Input `std::istream`.
   * @throws IOError if the stream is badly formatted (i.e. it does not contain
   * a valid JSON).
   * @note `this` must be constructed with an appropriate context @b BEFORE
   * calling this function. For example,
   * @code
   * PtxtArray my_pa(context);
   * my_pa.readJSON(std::cin);
   * @endcode
   *
   * The input stream has to contain a valid typed JSON value.\n
   * Each element of the content list will be deserialized as a slot of the type
   * determined by the scheme.\n
   * If the number of tokens in the slot list is less than the number of slots,
   * the remaining slots will be padded by 0.\n
   * For example a slot list '['slot0', 'slot1', 'slot2']' will be deserialized
   * as a plaintext `pa` where `pa[0]=slot0`, `pa[1]=slot1`, `pa[2]=slot2`, and
   * `pa[i]=0` for `i>2`.
   **/
  void readJSON(std::istream& is);

  /**
   * @brief In-place function to deserialize a `PtxtArray` from a `JsonWrapper`
   * object.
   * @param jw `JsonWrapper` containing the serialized object.
   * @throws IOError if the `JsonWrapper` object does not contain a valid
   * serialization of a `PtxtArray`.
   * @note `this` must be constructed with an appropriate context @b BEFORE
   * calling this function. For example,
   * @code
   * PtxtArray my_pa(context);
   * my_pa.readJSON(...);
   * @endcode
   *
   * The `JsonWrapper` must contain a valid `PtxtArray` serialization.\n
   * Each element of the content list will be deserialized as a slot of the
   * type determined by the scheme.\n
   * If the number of tokens in the slot list is less than the number of slots,
   * the remaining slots will be padded by 0.\n
   * For example a slot list '['slot0', 'slot1', 'slot2']' will be deserialized
   * as a plaintext `pa` where `pa[0]=slot0`, `pa[1]=slot1`, `pa[2]=slot2` and
   * `pa[i]=0` for `i>2`.
   **/
  void readJSON(const JsonWrapper& jw);

  friend std::istream& operator>>(std::istream& is, PtxtArray& pa);

  friend std::ostream& operator<<(std::ostream& os, const PtxtArray& pa);
};

inline void rotate(PtxtArray& a, long k) { rotate(a.ea, a.pa, k); }

inline void shift(PtxtArray& a, long k) { shift(a.ea, a.pa, k); }

inline void rotate1D(PtxtArray& a, long i, long k)
{
  rotate1D(a.ea, a.pa, i, k);
}

inline void shift1D(PtxtArray& a, long i, long k) { shift1D(a.ea, a.pa, i, k); }

inline bool operator==(const PtxtArray& a, const PtxtArray& b)
{
  assertTrue(&a.ea == &b.ea, "PtxtArray: inconsistent operation");
  return equals(a.ea, a.pa, b.pa);
}

inline bool operator!=(const PtxtArray& a, const PtxtArray& b)
{
  assertTrue(&a.ea == &b.ea, "PtxtArray: inconsistent operation");
  return !equals(a.ea, a.pa, b.pa);
}

inline PtxtArray& operator+=(PtxtArray& a, const PtxtArray& b)
{
  assertTrue(&a.ea == &b.ea, "PtxtArray: inconsistent operation");
  add(a.ea, a.pa, b.pa);
  return a;
}

template <class T>
auto operator+=(PtxtArray& a, const T& b) -> decltype(a.load(b), a)
// SFINAE: this allows operator+= to be more easily overloaded
// PtxtArray& operator+=(PtxtArray& a, const T& b)
{
  return a += PtxtArray(a.ea, b);
}

inline PtxtArray& operator-=(PtxtArray& a, const PtxtArray& b)
{
  assertTrue(&a.ea == &b.ea, "PtxtArray: inconsistent operation");
  sub(a.ea, a.pa, b.pa);
  return a;
}

template <class T>
auto operator-=(PtxtArray& a, const T& b) -> decltype(a.load(b), a)
// SFINAE: this allows operator-= to be more easily overloaded
// PtxtArray& operator-=(PtxtArray& a, const T& b)
{
  return a -= PtxtArray(a.ea, b);
}

inline PtxtArray& operator*=(PtxtArray& a, const PtxtArray& b)
{
  assertTrue(&a.ea == &b.ea, "PtxtArray: inconsistent operation");
  mul(a.ea, a.pa, b.pa);
  return a;
}

template <class T>
auto operator*=(PtxtArray& a, const T& b) -> decltype(a.load(b), a)
// SFINAE: this allows operator*= to be more easily overloaded
// PtxtArray& operator*=(PtxtArray& a, const T& b)
{
  return a *= PtxtArray(a.ea, b);
}

inline void frobeniusAutomorph(PtxtArray& a, long j)
{
  frobeniusAutomorph(a.ea, a.pa, j);
}

inline void frobeniusAutomorph(PtxtArray& a, const NTL::Vec<long>& vec)
{
  frobeniusAutomorph(a.ea, a.pa, vec);
}

inline void conjugate(PtxtArray& a) { frobeniusAutomorph(a, 1); }

inline void extractRealPart(PtxtArray& a) { extractRealPart(a.ea, a.pa); }

inline void extractImPart(PtxtArray& a) { extractImPart(a.ea, a.pa); }

inline void applyPerm(PtxtArray& a, const NTL::Vec<long>& pi)
{
  applyPerm(a.ea, a.pa, pi);
}

inline void power(PtxtArray& a, long e) { power(a.ea, a.pa, e); }

// For CKKS, returns max norm of slots, for BGV the trivial norm
// (i.e., 0 if zero, 1 otherwise)
inline double Norm(const PtxtArray& a) { return Norm(a.ea, a.pa); }

inline double Distance(const PtxtArray& a, const PtxtArray& b)
{
  assertTrue(&a.ea == &b.ea, "PtxtArray: inconsistent operation");
  return Distance(a.ea, a.pa, b.pa);
}

inline void totalSums(PtxtArray& a) { totalSums(a.ea, a.pa); }
inline void runningSums(PtxtArray& a) { runningSums(a.ea, a.pa); }

//=====================================

// Following are functions for performing "higher level"
// operations on "encrypted arrays".  There is really no
// reason for these to be members of the EncryptedArray class,
// so they are just declared as separate functions.

//! @brief A ctxt that encrypts \f$(x_1, ..., x_n)\f$ is replaced by an
//! encryption of \f$(y_1, ..., y_n)\f$, where \f$y_i = sum_{j\le i} x_j\f$.
void runningSums(const EncryptedArray& ea, Ctxt& ctxt);
// The implementation uses O(log n) shift operations.

inline void runningSums(Ctxt& ctxt)
{
  runningSums(ctxt.getContext().getView(), ctxt);
}

//! @brief A ctxt that encrypts \f$(x_1, ..., x_n)\f$ is replaced by an
//! encryption of \f$(y, ..., y)\$, where \f$y = sum_{j=1}^n x_j.\f$
void totalSums(const EncryptedArray& ea, Ctxt& ctxt);

inline void totalSums(Ctxt& ctxt)
{
  totalSums(ctxt.getContext().getView(), ctxt);
}

//! @brief Map all non-zero slots to 1, leaving zero slots as zero.
//! Assumes that r=1, and that all the slots contain elements from GF(p^d).
void mapTo01(const EncryptedArray& ea, Ctxt& ctxt, bool multithread = true);
// Implemented in eqtesting.cpp. We compute
//             x^{p^d-1} = x^{(1+p+...+p^{d-1})*(p-1)}
// by setting y=x^{p-1} and then outputting y * y^p * ... * y^{p^{d-1}},
// with exponentiation to powers of p done via Frobenius.

//! @brief Map all non-zero slots to 1, leaving zero slots as zero.
template <typename Scheme>
void mapTo01(const EncryptedArray&, Ptxt<Scheme>& ptxt);

//! @brief (only for p=2, r=1), test if prefixes of bits in slots are all zero.
//! Set slot j of res[i] to 0 if bits 0..i of j'th slot in ctxt are all zero,
//! else sets slot j of res[i] to 1.
//! It is assumed that res and the res[i]'s are initialized by the caller.
void incrementalZeroTest(Ctxt* res[],
                         const EncryptedArray& ea,
                         const Ctxt& ctxt,
                         long n);
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
 *     std::vector<ZZX> L(d);
 *     for (long j = 0; j < d; j++)
 *       if (j % 2 == 0) L[j] = ZZX(j, 1);
 *
 *     std::vector<ZZX> C;
 *     ea.buildLinPolyCoeffs(C, L);
 *     applyLinPoly1(ea, ctxt, C);
 * \endcode
 */

// VJS-FIXME: we plan on re-implementing this logic in matmul.{h,cpp}.
// The reasons are two-fold:
//   (1) the interfaces should look a lot more like the matmul interfaces,
//       instead of a completely different interface
//   (2) we want to exploit better algorithms, like BS/GS and hoisting,
//       which we don't do now.

//! @brief Apply the same linear transformation to all the slots
//! @param C is the output of ea.buildLinPolyCoeffs;
void applyLinPoly1(const EncryptedArray& ea,
                   Ctxt& ctxt,
                   const std::vector<NTL::ZZX>& C);

//! @brief Apply different transformations to different slots
//! @param Cvec is a std::vector of length ea.size(), each entry of which
//!        is the output of ea.buildLinPolyCoeffs;
void applyLinPolyMany(const EncryptedArray& ea,
                      Ctxt& ctxt,
                      const std::vector<std::vector<NTL::ZZX>>& Cvec);

//! @brief a low-level variant:
//! @param encodedCoeffs has all the linPoly coeffs encoded  in slots;
//!        different transformations can be encoded in different slots
template <typename P> // P can be ZZX or DoubleCRT
void applyLinPolyLL(Ctxt& ctxt, const std::vector<P>& encodedC, long d);
///@}

// Helper class for unimplemented pa_impl classes

template <typename type>
struct pa_no_impl
{
  template <typename... Args>
  static void apply(UNUSED Args&&... args)
  {
    throw helib::LogicError("function not implemented");
  }
};

#define HELIB_NO_CKKS_IMPL(impl)                                               \
  template <>                                                                  \
  class impl<PA_cx> : public pa_no_impl<PA_cx>                                 \
  {};

} // namespace helib

#endif // ifndef HELIB_ENCRYPTEDARRAY_H
