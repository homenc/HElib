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
/**
 * @file powerful.h
 * @brief The "powerful basis" of a cyclotomic ring
 **/

#include "NumbTh.h"
#include "cloned_ptr.h"
#include "bluestein.h"
#include "hypercube.h"
#include "DoubleCRT.h"
#include "FHEContext.h"

//! @class PowerfulTranslationIndexes
//! @brief Holds index tables for translation between powerful and zz_pX
class PowerfulTranslationIndexes {
 public:
  long m;     // product of mi's
  long phim;  // phi(m) = product of phi(mi)'s
  Vec<long> mvec;   // mvec[i] = mi
  Vec<long> phivec; // phivec[i] = phi(mi)
  Vec<long> divvec; // divvec[i] = m/mi
  Vec<long> invvec; // invvec[i] = (m/mi)^{-1} mod mi

  CubeSignature longSig;  // hypercube of dimensions mi
  CubeSignature shortSig; // hypercube of dimensions phi(mi)

  Vec<long> polyToCubeMap; // index translation tables
  Vec<long> cubeToPolyMap;
  Vec<long> shortToLongMap;

  Vec<ZZX> cycVec; // cycvec[i] = Phi_mi(X)
  ZZX phimX;

  PowerfulTranslationIndexes(const Vec<long>& mv);
};

/**
 * @class PowerfulConversion
 * @brief Conversion between powerful representation in R_m/(q) and zz_pX
 *
 * Usage pattern is as follows:
 *
 *    // compute tables for index translation
 *    PowerfulTranslationIndexes ind(mvec); // mvec is some factorization of m
 *
 *    // ... set the current zz_p::modulus to something before initializing
 *    PowerfulConversion pConv(ind);
 *
 *    // Alternatively use
 *    //   PowerfulConversion pConv(); pConv.initPConv(ind);
 *    // Only the latter call needs zz_p::modulus to be defined
 *
 *    // A powerful basis is defined wrt same modulus and cube signature
 *    HyperCube<zz_p> powerful(pConv.getShortSig());
 *
 *    // ... some code here to initialize powerful
 *    // code can also do other stuff, perhaps changing zz_p::modulus
 *
 *    pConv.restoreModulus(); // restore  zz_p::modulus
 *    zz_pX poly;             // defined relative to the same modulus
 *    pConv.powerfulToPoly(poly, powerful);
 *
 *    // ... some more code here, perhaps changing zz_p::modulus again
 *
 *    pConv.restoreModulus(); // restore  zz_p::modulus
 *    pConv.polyToPowerful(powerful, poly);
 **/
class PowerfulConversion {
  const PowerfulTranslationIndexes* indexes;
  zz_pContext zzpContext;
  Vec<zz_pXModulus> cycVec_p; // cycvec[i] = Phi_mi(X)
  zz_pXModulus phimX_p;

public:

  PowerfulConversion(): indexes(NULL) {}

  explicit PowerfulConversion(const PowerfulTranslationIndexes& ind):
  indexes(NULL) { initPConv(ind); }

  void initPConv(const PowerfulTranslationIndexes& ind)
  {
    if (indexes!=NULL) return; // cannot re-initialize a non-NULL object
    indexes = &ind;

    cycVec_p.SetLength(ind.cycVec.length());
    zzpContext.save(); // store the current modulus
    for (long i=0; i<ind.cycVec.length(); i++) {
      cycVec_p[i] = conv<zz_pX>(ind.cycVec[i]); // convert to zz_pXModulus
    }
    phimX_p = conv<zz_pX>(ind.phimX); // convert to zz_pXModulus
  }

  void restoreModulus() const { zzpContext.restore(); }
  const CubeSignature& getLongSig() const { return indexes->longSig; }
  const CubeSignature& getShortSig() const { return indexes->shortSig; }

  //! The conversion routines return the value of the modulus q.
  //! It is assumed that the modulus is already set before calling them
  long powerfulToPoly(zz_pX& poly, const HyperCube<zz_p>& powerful) const;
  long polyToPowerful(HyperCube<zz_p>& powerful, const zz_pX& poly) const;
};

/**
 * @class PowerfulDCRT
 * @brief Conversion between powerful representation, DoubleCRT, and ZZX
 **/
class PowerfulDCRT {
  const FHEcontext& context; // points to the context for the DoubleCRT's

  PowerfulTranslationIndexes indexes; // modulus-independent tables

  // a vector of PowerfulConversion tables, one for each modulus in the chain
  Vec<PowerfulConversion> pConvVec;

public:
  PowerfulDCRT(const FHEcontext& _context, const Vec<long>& mvec);

  const PowerfulTranslationIndexes& getIndexTranslation() const
  { return indexes; }
  const PowerfulConversion& getPConv(long i) const
  { return pConvVec.at(i); }

  void dcrtToPowerful(Vec<ZZ>& powerful, const DoubleCRT& dcrt) const;
  void powerfulToDCRT(DoubleCRT& dcrt, const Vec<ZZ>& powerful) const;

  // If the IndexSet is omitted, default to all the primes in the chain
  void ZZXtoPowerful(Vec<ZZ>& powerful, const ZZX& poly,
		     IndexSet s=IndexSet::emptySet()    ) const;
  void powerfulToZZX(ZZX& poly, const Vec<ZZ>& powerful,
		     IndexSet s=IndexSet::emptySet()    ) const;
};







/********************************************************************/
/****************    UNUSED CODE - COMMENTED OUT   ******************/
/********************************************************************/

#if 0
//! For vec[d] = (a_d , b_d), returns \prod_d a_d^{b_d}
long computeProd(const Vec< Pair<long, long> >& vec);

//! For x=p^e, returns phi(p^e) = (p-1) p^{e-1}
inline long computePhi(const Pair<long, long>& x)
{
  long p = x.a;
  long e = x.b;
  return power_long(p, e - 1) * (p-1);
}

//! For factors[d] = (p_d, e_d), computes
//! phiVec[d] = phi(p_d^{e_d}) = (p_d-1) p_i{e_d-1}
void computePhiVec(Vec<long>& phiVec,
                   const Vec< Pair<long, long> >& factors);

//! For powVec[d] = p_d^{e_d}, m = \prod_d p_d^{e_d}, computes
//! divVec[d] = m/p_d^{e_d}
void computeDivVec(Vec<long>& divVec, long m,
                   const Vec<long>& powVec);

//! For divVec[d] = m/p_d^{e_d}, powVec[d] = p^{e_d}, computes
//! invVec[d] = divVec[d]^{-1} mod powVec[d]
void computeInvVec(Vec<long>& invVec,
                   const Vec<long>& divVec, const Vec<long>& powVec);


//! shortSig is a CubeSignature for (phi(m_1), .., phi(m_k)),
//! longSig is a CubeSignature for (m_1, ..., m_k).
//! computes shortToLongMap[i] that maps an index i 
//! with respect to shortSig to the corresponding index
//! with respect to longSig.
void computeShortToLongMap(Vec<long>& shortToLongMap, 
                           const CubeSignature& shortSig, 
                           const CubeSignature& longSig);

//! Computes the inverse of the shortToLongMap, computed above.
//! "undefined" entries are initialzed to -1.
void computeLongToShortMap(Vec<long>& longToShortMap,
                           long m,
                           const Vec<long>& shortToLongMap);

//! This routine recursively reduces each hypercolumn
//! in dimension d (viewed as a coeff vector) by Phi_{m_d}(X)
//! If one starts with a cube of dimension (m_1, ..., m_k),
//! one ends up with a cube that is effectively of dimension
//! phi(m_1, ..., m_k). Viewed as an element of the ring
//! F_p[X_1,...,X_k]/(Phi_{m_1}(X_1), ..., Phi_{m_k}(X_k)),
//! the cube remains unchanged.
void recursiveReduce(const CubeSlice<zz_p>& s, 
                     const Vec<zz_pXModulus>& cycVec, 
                     long d,
                     zz_pX& tmp1,
                     zz_pX& tmp2);

//! This routine implements the isomorphism from
//! F_p[X]/(Phi_m(X)) to F_p[X_1, ..., X_k]/(Phi_{m_1}(X_1),...,Phi_{m_k}(X_k))
//! The input is poly, which must be of degree < m, and the
//! output is cube, which is a HyperCube of dimension (phi(m_1), ..., phi(m_k)).
//! The caller is responsible to supply "scratch space" in the
//! form of a HyperCube tmpCube of dimension (m_1, ..., m_k).
void convertPolyToPowerful(HyperCube<zz_p>& cube, 
                           HyperCube<zz_p>& tmpCube, 
                           const zz_pX& poly,
                           const Vec<zz_pXModulus>& cycVec,
                           const Vec<long>& polyToCubeMap,
                           const Vec<long>& shortToLongMap);

//! This implements the inverse of the above isomorphism.
void convertPowerfulToPoly(zz_pX& poly,
                           const HyperCube<zz_p>& cube,
                           long m,
                           const Vec<long>& shortToLongMap,
                           const Vec<long>& cubeToPolyMap,
                           const zz_pXModulus& phimX);
#endif
/********************************************************************/
#if 0
//! For x = (p,e), returns p^e
inline long computePow(const Pair<long, long>& x)
{
  long p = x.a;
  long e = x.b;
  return power_long(p, e);
}

//! For factors[d] = (p_d, e_d), computes powVec[d] = p_d^{e_d}
void computePowVec(Vec<long>& powVec, 
                   const Vec< Pair<long, long> >& factors);

//! this maps an index j in [phi(m)] to a vector
//! representing the powerful basis coordinates
void mapIndexToPowerful(Vec<long>& pow, long j, const Vec<long>& phiVec);

//! @deprecated For powVec[d]=p_d^{e_d}, cycVec[d] = Phi_{p_d^{e_d}}(X) mod p
void computeCycVec(Vec<zz_pXModulus>& cycVec, const Vec<long>& powVec);

//! m = m_1 ... m_k, m_d = p_d^{e_d} 
//! powVec[d] = m_d
//! invVec[d] = (m/m_d)^{-1} mod m_d
//! computes polyToCubeMap[i] and cubeToPolyMap[i] 
//!   where polyToCubeMap[i] is the index of (i_1, ..., i_k)
//!   in the cube with CubeSignature longSig = (m_1, ..., m_k),
//!   and (i_1, ..., i_k) is the unique tuple satistifying
//!   i = i_1 (m/m_1) + ... + i_k (m/m_k) mod m
//! and
//!   cubeToPolyMap is the inverse map.
void computePowerToCubeMap(Vec<long>& polyToCubeMap,
                           Vec<long>& cubeToPolyMap,
                           long m,
                           const Vec<long>& powVec,
                           const Vec<long>& invVec,
                           const CubeSignature& longSig);

//! powVec[d] = m_d = p_d^{e_d}
//! computes multiEvalPoints[d] as a vector of length phi(m_d)
//!   containing base^{(m/m_d) j} for j in Z_{m_d}^*
void computeMultiEvalPoints(Vec< Vec<zz_p> >& multiEvalPoints,
                            const zz_p& base,
                            long m,
                            const Vec<long>& powVec,
                            const Vec<long>& phiVec);

//! computes linearEvalPoints[i] = base^i, i in Z_m^*
void computeLinearEvalPoints(Vec<zz_p>& linearEvalPoints,
                             const zz_p& base,
                             long m, long phim);

//! powVec[d] = m_d = p_d^{e_d}
//! computes compressedIndex[d] as a vector of length m_d,
//!   where compressedIndex[d][j] = -1 if GCD(j, m_d) != 1,
//!   and otherwise is set to the relative numerical position
//!   of j among Z_{m_d}^*
void computeCompressedIndex(Vec< Vec<long> >& compressedIndex,
                            const Vec<long>& powVec);

//! computes powToCompressedIndexMap[i] as -1 if GCD(i, m) != 1, 
//!   and otherwise as the index of the point (j_1, ..., j_k)
//!   relative to a a cube of dimension (phi(m_1), ..., phi(m_k)),
//!   where each j_d is the compressed index of i_d = i mod m_d.
void computePowToCompressedIndexMap(Vec<long>& powToCompressedIndexMap,
                                    long m,
                                    const Vec<long>& powVec,
                                    const Vec< Vec<long> >& compressedIndex,
                                    const CubeSignature& shortSig);

void recursiveEval(const CubeSlice<zz_p>& s,
                   const Vec< Vec<zz_p> >& multiEvalPoints,
                   long d,
                   zz_pX& tmp1,
                   Vec<zz_p>& tmp2);

inline void eval(HyperCube<zz_p>& cube,
		 const Vec< Vec<zz_p> >& multiEvalPoints)
{
   zz_pX tmp1;
   Vec<zz_p> tmp2;

   recursiveEval(CubeSlice<zz_p>(cube), multiEvalPoints, 0, tmp1, tmp2);
} 

void mapPowerfulToPoly(ZZX& poly, 
                       const Vec<long>& pow, 
                       const Vec<long>& divVec,
                       long m,
                       const ZZX& phimX);
#endif
/********************************************************************/
#if 0
//! @class FFTHelper
//! @brief Helper class for FFT over NTL::zz_p
//!
//! The class FFTHelper is used to help perform FFT over zz_p.
//! The constructor supplies an element x in zz_p and an integer m,
//! where x has order m and x has a square root in zz-p.
//! Auxilliary data structures are constricted to support
//! evaluation and interpolation at the points x^i for i in Z_m^*
//!
//! It is assumed that the zz_p-context is set prior to all
//! constructor and method invocations.
class FFTHelper {

private:
  long m;
  zz_p m_inv;
  zz_p root, iroot;
  zz_pXModulus phimx;
  Vec<bool> coprime;
  long phim;

  mutable zz_pX powers, ipowers;
  mutable Vec<mulmod_precon_t> powers_aux, ipowers_aux;
  mutable fftRep Rb, iRb;
  mutable fftrep_aux Rb_aux, iRb_aux;
  mutable fftRep Ra; 
  mutable zz_pX tmp;

public:
  FFTHelper(long _m, zz_p x);

  void FFT(const zz_pX& f, Vec<zz_p>& v) const;
    // compute v = { f(x^i) }_{i in Z_m^*}
  void iFFT(zz_pX& f, const Vec<zz_p>& v, bool normalize = true) const;
    // computes inverse transform. If !normalize, result is scaled by m

  const zz_p& get_m_inv() const { return m_inv; }
  // useful for aggregate normalization

};


//! powVec[d] = m_d = p_d^{e_d}
//! computes multiEvalPoints[d] as an FFTHelper for base^{m/m_d}
void computeMultiEvalPoints(Vec< copied_ptr<FFTHelper> >& multiEvalPoints,
                            const zz_p& base,
                            long m,
                            const Vec<long>& powVec,
                            const Vec<long>& phiVec);

void recursiveEval(const CubeSlice<zz_p>& s,
                   const Vec< copied_ptr<FFTHelper> >& multiEvalPoints,
                   long d,
                   zz_pX& tmp1,
                   Vec<zz_p>& tmp2);

inline void eval(HyperCube<zz_p>& cube,
          const Vec< copied_ptr<FFTHelper> >& multiEvalPoints)
{
   zz_pX tmp1;
   Vec<zz_p> tmp2;

   recursiveEval(CubeSlice<zz_p>(cube), multiEvalPoints, 0, tmp1, tmp2);
} 

void recursiveInterp(const CubeSlice<zz_p>& s,
                     const Vec< copied_ptr<FFTHelper> >& multiEvalPoints,
                     long d,
                     zz_pX& tmp1,
                     Vec<zz_p>& tmp2);

void interp(HyperCube<zz_p>& cube,
	    const Vec< copied_ptr<FFTHelper> >& multiEvalPoints);
#endif
