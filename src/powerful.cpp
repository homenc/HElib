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
/* a preliminary test program for playing around with Peikert's "powerful"
 * basis. EXPERIMENTAL CODE, not usable yet
 */

#include <helib/powerful.h>

namespace helib {

// powVec[d] = p_d^{e_d}, m = \prod_d p_d^{e_d}
// computes divVec[d] = m/p_d^{e_d}
inline void computeDivVec(NTL::Vec<long>& divVec,
                          long m,
                          const NTL::Vec<long>& powVec)
{
  long k = powVec.length();
  divVec.SetLength(k);

  for (long d = 0; d < k; d++)
    divVec[d] = m / powVec[d];
}

// divVec[d] = m/p_d^{e_d}, powVec[d] = p^{e_d}
// computes invVec[d] = divVec[d]^{-1} mod powVec[d]
inline void computeInvVec(NTL::Vec<long>& invVec,
                          const NTL::Vec<long>& divVec,
                          const NTL::Vec<long>& powVec)
{
  long k = divVec.length();
  invVec.SetLength(k);

  for (long d = 0; d < k; d++) {
    long t1 = divVec[d] % powVec[d];
    long t2 = NTL::InvMod(t1, powVec[d]);
    invVec[d] = t2;
  }
}

// m = m_1 ... m_k, m_d = p_d^{e_d}
// powVec[d] = m_d
// invVec[d] = (m/m_d)^{-1} mod m_d
// computes polyToCubeMap[i] and cubeToPolyMap[i]
//   where polyToCubeMap[i] is the index of (i_1, ..., i_k)
//   in the cube with CubeSignature longSig = (m_1, ..., m_k),
//   and (i_1, ..., i_k) is the unique tuple satisfying
//   i = i_1 (m/m_1) + ... + i_k (m/m_k) mod m
// and
//   cubeToPolyMap is the inverse map.
static void computePowerToCubeMap(NTL::Vec<long>& polyToCubeMap,
                                  NTL::Vec<long>& cubeToPolyMap,
                                  long m,
                                  const NTL::Vec<long>& powVec,
                                  const NTL::Vec<long>& invVec,
                                  const CubeSignature& longSig)
{
  long k = powVec.length();

  polyToCubeMap.SetLength(m);
  cubeToPolyMap.SetLength(m);

  for (long i = 0; i < m; i++) {
    long j = 0;
    for (long d = 0; d < k; d++) {
      long i_d = NTL::MulMod((i % powVec[d]), invVec[d], powVec[d]);
      j += i_d * longSig.getProd(d + 1);
    }
    polyToCubeMap[i] = j;
    cubeToPolyMap[j] = i;
  }
}

// shortSig is a CubeSignature for (phi(m_1), .., phi(m_k)),
// longSig is a CubeSignature for (m_1, ..., m_k).
// computes shortToLongMap[i] that maps an index i
// with respect to shortSig to the corresponding index
// with respect to longSig.
static void computeShortToLongMap(NTL::Vec<long>& shortToLongMap,
                                  const CubeSignature& shortSig,
                                  const CubeSignature& longSig)
{
  long phim = shortSig.getSize();
  long k = shortSig.getNumDims();

  shortToLongMap.SetLength(phim);

  for (long i = 0; i < phim; i++) {
    long j = 0;
    for (long d = 0; d < k; d++) {
      long i_d = shortSig.getCoord(i, d);
      j += i_d * longSig.getProd(d + 1);
    }

    shortToLongMap[i] = j;
  }
}

// This routine recursively reduces each hypercolumn
// in dimension d (viewed as a coeff vector) by Phi_{m_d}(X)
// If one starts with a cube of dimension (m_1, ..., m_k),
// one ends up with a cube that is effectively of dimension
// phi(m_1, ..., m_k). Viewed as an element of the ring
// F_p[X_1,...,X_k]/(Phi_{m_1}(X_1), ..., Phi_{m_k}(X_k)),
// the cube remains unchanged.
static void recursiveReduce(const CubeSlice<NTL::zz_p>& s,
                            const NTL::Vec<NTL::zz_pXModulus>& cycVec,
                            long d,
                            NTL::zz_pX& tmp1,
                            NTL::zz_pX& tmp2)
{
  long numDims = s.getNumDims();
  assertTrue(numDims > 0l, "CubeSlice s has negative number of dimensions");

  long deg0 = deg(cycVec[d]);

  long posBnd = s.getProd(1);
  for (long pos = 0; pos < posBnd; pos++) {
    getHyperColumn(tmp1.rep, s, pos);
    tmp1.normalize();

    // tmp2 may not be normalized, so clear it first
    clear(tmp2);

    rem(tmp2, tmp1, cycVec[d]);

    // now pad tmp2.rep with zeros to length deg0...
    // tmp2 may not be normalized
    long len = tmp2.rep.length();
    tmp2.rep.SetLength(deg0);
    for (long i = len; i < deg0; i++)
      tmp2.rep[i] = 0;

    setHyperColumn(tmp2.rep, s, pos);
  }

  if (numDims == 1)
    return;

  for (long i = 0; i < deg0; i++)
    recursiveReduce(CubeSlice<NTL::zz_p>(s, i), cycVec, d + 1, tmp1, tmp2);
}

PowerfulTranslationIndexes::PowerfulTranslationIndexes(
    const NTL::Vec<long>& mv) :
    mvec(mv) // copy the vector of factors
{
  // mvec contains the prime-power factorization of m = \prod_{i=1}^k mi
  long nfactors = mvec.length(); // = k
  m = computeProd(mvec);         // compute m itself

  // phivec holds phi(mi) for all factors mi
  phivec.SetLength(nfactors);
  for (long i = 0; i < nfactors; i++)
    phivec[i] = phi_N(mvec[i]);
  phim = computeProd(phivec); // phi(m) = prod_i phi(mi)

  computeDivVec(divvec, m, mvec); // divvec[i] = m/mi

  computeInvVec(invvec, divvec, mvec); // invvec[i] = (m/mi)^{-1} mod mi

  // Let (i_1,...,i_k) be the representation of i in base
  // (m/m1,...,m/mk), namely i = i_1 (m/m_1)+...+i_k (m/m_k) mod m.
  // Then polyToCubeMap[i] is the lexicographic index of the tuple
  // (i_1,...,i_k) in the cube with dimensions (m_1, ..., m_k).
  // cubeToPolyMap is the inverse map, polyToCubeMap[cubeToPolyMap[j]]=j.
  longSig.initSignature(mvec);
  shortSig.initSignature(phivec);
  computePowerToCubeMap(polyToCubeMap, cubeToPolyMap, m, mvec, invvec, longSig);

  // shortSig is a CubeSignature for (phi(m_1),..., phi(m_k)), and longSig
  // is a CubeSignature for (m_1, ..., m_k). shortToLongMap[i] maps an
  // index i wrt shortSig to an index i' wrt longSig so that both indexes
  // correspond to the same tuple (i_1,...,i_k).
  computeShortToLongMap(shortToLongMap, shortSig, longSig);

  cycVec.SetLength(nfactors);
  for (long d = 0; d < nfactors; d++)
    cycVec[d] = Cyclotomic(mvec[d]);

  phimX = Cyclotomic(m);
}

// This routine implements the isomorphism from F_p[X]/(Phi_m(X)) to
// F_p[X_1, ..., X_k]/(Phi_{m_1}(X_1), ..., Phi_{m_k}(X_k)). The input
// is a polynomial mod q, which must be of degree < m. The output is a
// HyperCube of dimension (phi(m_1), ..., phi(m_k)).
//
// It is assumed that the current modulus is already set.
// For convenience, this method returns the value of the modulus q.
long PowerfulConversion::polyToPowerful(HyperCube<NTL::zz_p>& powerful,
                                        const NTL::zz_pX& poly) const
{
  HyperCube<NTL::zz_p> tmpCube(getLongSig());

  long n = deg(poly);
  assertTrue(n < indexes->m,
             "Degree of polynomial poly is greater or equal than indexes->m");

  for (long i = 0; i <= n; i++)
    tmpCube[indexes->polyToCubeMap[i]] = poly[i];

  for (long i = n + 1; i < indexes->m; i++)
    tmpCube[indexes->polyToCubeMap[i]] = 0;

  NTL::zz_pX tmp1, tmp2;
  recursiveReduce(CubeSlice<NTL::zz_p>(tmpCube), cycVec_p, 0, tmp1, tmp2);

  for (long i = 0; i < indexes->phim; i++)
    powerful[i] = tmpCube[indexes->shortToLongMap[i]];

  return NTL::zz_p::modulus();
}

long PowerfulConversion::powerfulToPoly(
    NTL::zz_pX& poly,
    const HyperCube<NTL::zz_p>& powerful) const
{
  //  convertPowerfulToPoly(poly, powerful, indexes->m, indexes.shortToLongMap,
  //			indexes->cubeToPolyMap, phimX_p);
  // a temporary degree-(m-1) polynomial, initialized to all-zero
  NTL::zz_pX tmp;
  tmp.SetLength(indexes->m);
  for (long i = 0; i < indexes->m; i++)
    tmp[i] = 0;

  // copy the coefficient from hypercube in the right order
  for (long i = 0; i < indexes->phim; i++)
    tmp[indexes->cubeToPolyMap[indexes->shortToLongMap[i]]] = powerful[i];
  // FIXME: these two maps could be composed into a single map

  tmp.normalize();
  rem(poly, tmp, phimX_p); // reduce modulo Phi_m(X)

  return NTL::zz_p::modulus();
}

PowerfulDCRT::PowerfulDCRT(const Context& _context,
                           const NTL::Vec<long>& mvec) :
    context(_context), indexes(mvec)
{
  triv = (mvec.length() == 1);
  // triv means poly and pwfl bases are identical

  if (triv)
    return;

  double polyNormBnd = context.zMStar.getPolyNormBnd();

  // double normBnd = context.zMStar.getNormBnd();
  // we compute this from the given factorization, just in case
  // it is not a prime-power factorization
  double normBnd = 1.0;
  for (long fac : mvec)
    normBnd *= calcPolyNormBnd(fac);

  long phim = context.zMStar.getPhiM();

  to_pwfl_excess_bits = long((log(normBnd) + log(double(phim))) / log(2.0)) + 5;
  to_poly_excess_bits =
      long((log(polyNormBnd) + log(double(phim))) / log(2.0)) + 5;
  // NOTE: we are using the following:
  //   |poly(a)| <= polyNormBnd*|canon(a)| <= polyNormBnd*phim*|pwfl(a)|
  //   |pwfl(a)| <= normBnd*|canon(a)| <= normBnd*phim*|poly(a)|
  // where |.| is the infty norm.

  // Why +5?  +1 to get the ceiling, +1 for rounding errors, +1 for the extra
  // factor of 2 needed to CRT, +1 because NumBits(product) is not quite right,
  // and another +1 for good measure

  NTL::ZZ prod = context.productOfPrimes(context.allPrimes());
  long max_bits = (prod.size() + 1) * NTL_ZZ_NBITS; // +1 to be safe
  if (max_bits < 256)
    max_bits = 256; // just to make sure we can handle
                    // various smallish values

  long target_bits = max_bits + to_pwfl_excess_bits + to_poly_excess_bits;
  // we can afford back and forth conversions

  NTL::zz_pPush push; // backup NTL's current modulus

  long n = 0; // number of primes we need
  NTL::ZZ prod1{1};

  do {
    NTL::zz_p::FFTInit(n);
    long p = NTL::zz_p::modulus();
    prod1 *= p;
    product_bits.append(NumBits(prod1));
    n++;
  } while (product_bits[n - 1] < target_bits);

  // initialize the modulus-dependent tables for all the moduli
  pConvVec.SetLength(n);         // allocate space
  for (long i = 0; i < n; i++) { // initialize
    NTL::zz_p::FFTInit(i);
    pConvVec[i].initPConv(indexes); // initialize tables
  }
}

void PowerfulDCRT::ZZXtoPowerful(NTL::Vec<NTL::ZZ>& out,
                                 const NTL::ZZX& poly) const
{
  long phim = context.zMStar.getPhiM();

  if (triv) {
    NTL::VectorCopy(out, poly, phim);
    return;
  }

  long max_sz = 0;
  for (const NTL::ZZ& coeff : poly.rep) {
    long sz = coeff.size();
    if (max_sz < sz)
      max_sz = sz;
  }

  long target_bits = max_sz * NTL_ZZ_NBITS + to_pwfl_excess_bits;

  long n = product_bits.length();
  long m = 1;
  while (m <= n && product_bits[m - 1] < target_bits)
    m++;

  if (m > n)
    throw LogicError("ZZXtoPowerful: not enough primes");

  NTL::zz_pPush push; // backup NTL's current modulus
  NTL::ZZ product{1};

  NTL::Vec<NTL::ZZ> res;
  res.SetLength(phim);

  for (long i : range(m)) {
    pConvVec[i].restoreModulus();
    NTL::zz_pX oneRowPoly;
    NTL::conv(oneRowPoly, poly);
    HyperCube<NTL::zz_p> oneRowPwrfl(indexes.shortSig);
    pConvVec[i].polyToPowerful(oneRowPwrfl, oneRowPoly);
    NTL::CRT(res, product, oneRowPwrfl.getData());
  }

  out = res;
}

void PowerfulDCRT::powerfulToZZX(NTL::ZZX& poly,
                                 const NTL::Vec<NTL::ZZ>& powerful) const
{
  if (triv) {
    NTL::conv(poly, powerful);
    return;
  }

  long max_sz = 0;
  for (const NTL::ZZ& coeff : powerful) {
    long sz = coeff.size();
    if (max_sz < sz)
      max_sz = sz;
  }

  long target_bits = max_sz * NTL_ZZ_NBITS + to_poly_excess_bits;

  long n = product_bits.length();
  long m = 1;
  while (m <= n && product_bits[m - 1] < target_bits)
    m++;

  if (m > n)
    throw LogicError("powerfulToZZX: not enough primes");

  NTL::zz_pPush push; // backup NTL's current modulus
  clear(poly);
  NTL::ZZ product{1};

  for (long i : range(m)) {
    pConvVec[i].restoreModulus();
    HyperCube<NTL::zz_p> oneRowPwrfl(indexes.shortSig);
    NTL::conv(oneRowPwrfl.getData(), powerful);
    NTL::zz_pX oneRowPoly;
    pConvVec[i].powerfulToPoly(oneRowPoly, oneRowPwrfl);
    NTL::CRT(poly, product, oneRowPoly);
  }
}

void PowerfulDCRT::dcrtToPowerful(NTL::Vec<NTL::ZZ>& powerful,
                                  const DoubleCRT& dcrt) const
{
  NTL::ZZX poly;
  dcrt.toPoly(poly);

  if (triv) {
    long phim = context.zMStar.getPhiM();
    NTL::VectorCopy(powerful, poly, phim);
    return;
  }

  NTL::Vec<NTL::ZZ> pwfl;
  this->ZZXtoPowerful(pwfl, poly);
  NTL::ZZ Q = context.productOfPrimes(dcrt.getIndexSet());
  vecRed(powerful, pwfl, Q, /*abs=*/false);
  // reduce to interval [-Q/2,+Q/2]
}

/********************************************************************/
/****************    UNUSED CODE - COMMENTED OUT   ******************/
/********************************************************************/
#if 0

static void convertPolyToPowerful(HyperCube<NTL::zz_p>& cube,
                           HyperCube<NTL::zz_p>& tmpCube,
                           const NTL::zz_pX& poly,
                           const NTL::Vec<NTL::zz_pXModulus>& cycVec,
                           const NTL::Vec<long>& polyToCubeMap,
                           const NTL::Vec<long>& shortToLongMap)
{
   long m = tmpCube.getSize();
   long phim = cube.getSize();
   long n = deg(poly);

   assertTrue(n < m, "Degree of polynomial poly must be less than size of the hypercube tmpCube");

   for (long i = 0; i <= n; i++)
      tmpCube[polyToCubeMap[i]] = poly[i];

   for (long i = n+1; i < m; i++)
      tmpCube[polyToCubeMap[i]] = 0;

   NTL::zz_pX tmp1, tmp2;
   recursiveReduce(CubeSlice<NTL::zz_p>(tmpCube), cycVec, 0, tmp1, tmp2);

   for (long i = 0; i < phim; i++)
      cube[i] = tmpCube[shortToLongMap[i]];
}

// This implements the inverse of the above isomorphism.
static void convertPowerfulToPoly(NTL::zz_pX& poly,
                           const HyperCube<NTL::zz_p>& cube,
                           long m,
                           const NTL::Vec<long>& shortToLongMap,
                           const NTL::Vec<long>& cubeToPolyMap,
                           const NTL::zz_pXModulus& phimX)
{
   long phim = cube.getSize();

   NTL::zz_pX tmp;

   tmp.SetLength(m);

   for (long i = 0; i < m; i++)
      tmp[i] = 0;

   for (long i = 0; i < phim; i++)
      tmp[cubeToPolyMap[shortToLongMap[i]]] = cube[i];
      // FIXME: these two maps could be composed into a single map

   tmp.normalize();

   rem(poly, tmp, phimX);
}

// powVec[d] = p_d^{e_d}
// cycVec[d] = Phi_{p_d^{e_d}}(X) mod p
void computeCycVec(NTL::Vec<NTL::zz_pXModulus>& cycVec, const NTL::Vec<long>& powVec)
{
  long k = powVec.length();
  cycVec.SetLength(k);

  for (long d = 0; d < k; d++) {
    NTL::ZZX PhimX = Cyclotomic(powVec[d]);
    cycVec[d] = NTL::conv<NTL::zz_pX>(PhimX);
  }
}

// factors[d] = (p_d, e_d)
// computes phiVec[d] = phi(p_d^{e_d}) = (p_d-1) p_i{e_d-1}
static void computePhiVec(NTL::Vec<long>& phiVec,
			  const NTL::Vec< Pair<long, long> >& factors)
{
  long k = factors.length();
  phiVec.SetLength(k);

  for (long d = 0; d < k; d++)
    phiVec[d] = computePhi(factors[d]);
}

void mapIndexToPowerful(NTL::Vec<long>& pow, long j, const NTL::Vec<long>& phiVec)
// this maps an index j in [phi(m)] to a vector
// representing the powerful basis coordinates

{
  long k = phiVec.length();
  long phim = computeProd(phiVec);
  assertInRange(j, 0l, phim, "Index j is not in [0, computeProd(phiVec))");

  pow.SetLength(k);

  for (long i = k-1; i >= 0; i--) {
    pow[i] = j % phiVec[i];
    j = (j - pow[i])/phiVec[i];
  }
}


void mapPowerfulToPoly(NTL::ZZX& poly,
                       const NTL::Vec<long>& pow,
                       const NTL::Vec<long>& divVec,
                       long m,
                       const NTL::ZZX& phimX)
{
  long k = pow.length();
  assertEq(divVec.length(), k, "pow and divVec have different sizes");

  long j = 0;
  for (long i = 0; i < k; i++)
    j += pow[i] * divVec[i];

  j %= m;

  NTL::ZZX f = NTL::ZZX(j, 1);

  poly = f % phimX;
}

// powVec[d] = p_d^{e_d}
// cycVec[d] = Phi_{p_d^{e_d}}(X) mod p
void computeCycVec(NTL::Vec<NTL::zz_pXModulus>& cycVec, const NTL::Vec<long>& powVec)
{
  long k = powVec.length();
  cycVec.SetLength(k);

  for (long d = 0; d < k; d++) {
    NTL::ZZX PhimX = Cyclotomic(powVec[d]);
    cycVec[d] = NTL::conv<NTL::zz_pX>(PhimX);
  }
}

// vec[d] = (a_d , b_d)
// returns \prod_d a_d^{b_d}
long computeProd(const NTL::Vec< Pair<long, long> >& vec)
{
  long prod = 1;
  long k = vec.length();
  for (long d = 0; d < k; d++) {
    prod = prod * computePow(vec[d]);
  }
  return prod;
}


// factors[d] = (p_d, e_d)
// computes powVec[d] = p_d^{e_d}
void computePowVec(NTL::Vec<long>& powVec,
                   const NTL::Vec< Pair<long, long> >& factors)
{
  long k = factors.length();
  powVec.SetLength(k);
  for (long d = 0; d < k; d++)
    powVec[d] = computePow(factors[d]);
}

// Computes the inverse of the shortToLongMap, computed above.
// "undefined" entries are initialzed to -1.
void computeLongToShortMap(NTL::Vec<long>& longToShortMap,
                           long m,
                           const NTL::Vec<long>& shortToLongMap)
{
   long n = shortToLongMap.length();

   longToShortMap.SetLength(m);

   for (long i = 0; i < m; i++) longToShortMap[i] = -1;

   for (long j = 0; j < n; j++) {
      long i = shortToLongMap[j];
      longToShortMap[i] = j;
   }
}


// powVec[d] = m_d = p_d^{e_d}
// computes multiEvalPoints[d] as a vector of length phi(m_d)
//   containing base^{(m/m_d) j} for j in Z_{m_d}^*
void computeMultiEvalPoints(NTL::Vec< NTL::Vec<NTL::zz_p> >& multiEvalPoints,
                            const NTL::zz_p& base,
                            long m,
                            const NTL::Vec<long>& powVec,
                            const NTL::Vec<long>& phiVec)
{
   long k = powVec.length();

   multiEvalPoints.SetLength(k);

   for (long d = 0; d < k; d++) {
      long m_d = powVec[d];
      long phi_d = phiVec[d];
      long count = 0;

      NTL::zz_p pow = NTL::conv<NTL::zz_p>(1);
      NTL::zz_p mult = power(base, m/m_d);

      multiEvalPoints[d].SetLength(phi_d);

      for (long j = 0; j < m_d; j++) {
         if (GCD(j, m_d) == 1) {
            multiEvalPoints[d][count] = pow;
            count++;
         }
         pow = pow * mult;
      }
   }
}

// computes linearEvalPoints[i] = base^i, i in Z_m^*
void computeLinearEvalPoints(NTL::Vec<NTL::zz_p>& linearEvalPoints,
                             const NTL::zz_p& base,
                             long m, long phim)
{
   linearEvalPoints.SetLength(phim);
   NTL::zz_p pow = NTL::conv<NTL::zz_p>(1);

   for (long i = 0, j = 0; i < m; i++) {
      if (GCD(i, m) == 1) linearEvalPoints[j++] = pow;
      pow = pow * base;
   }
}

// powVec[d] = m_d = p_d^{e_d}
// computes compressedIndex[d] as a vector of length m_d,
//   where compressedIndex[d][j] = -1 if GCD(j, m_d) != 1,
//   and otherwise is set to the relative numerical position
//   of j among Z_{m_d}^*
void computeCompressedIndex(NTL::Vec< NTL::Vec<long> >& compressedIndex,
                            const NTL::Vec<long>& powVec)
{
   long k = powVec.length();

   compressedIndex.SetLength(k);

   for (long d = 0; d < k; d++) {
      long m_d = powVec[d];
      long count = 0;
      compressedIndex[d].SetLength(m_d);
      for (long j = 0; j < m_d; j++) {
         if (GCD(j, m_d) == 1) {
            compressedIndex[d][j] = count;
            count++;
         }
         else {
            compressedIndex[d][j] = -1;
         }
      }
   }
}


// computes powToCompressedIndexMap[i] as
//   -1 if GCD(i, m) != 1,
//   and otherwise as the index of the point (j_1, ..., j_k)
//   relative to a a cube of dimension (phi(m_1), ..., phi(m_k)),
//   where each j_d is the compressed index of i_d = i mod m_d.

void computePowToCompressedIndexMap(NTL::Vec<long>& powToCompressedIndexMap,
                                    long m,
                                    const NTL::Vec<long>& powVec,
                                    const NTL::Vec< NTL::Vec<long> >& compressedIndex,
                                    const CubeSignature& shortSig
                                   )
{
   long k = powVec.length();
   powToCompressedIndexMap.SetLength(m);

   for (long i = 0; i < m; i++) {
      if (GCD(i, m) != 1)
         powToCompressedIndexMap[i] = -1;
      else {
         long j = 0;
         for (long d = 0; d < k; d++) {
            long i_d = i % powVec[d];
            long j_d = compressedIndex[d][i_d];
            j += j_d * shortSig.getProd(d+1);
         }
         powToCompressedIndexMap[i] = j;
      }
   }
}

void recursiveEval(const CubeSlice<NTL::zz_p>& s,
                   const NTL::Vec< NTL::Vec<NTL::zz_p> >& multiEvalPoints,
                   long d,
                   NTL::zz_pX& tmp1,
                   NTL::Vec<NTL::zz_p>& tmp2)
{
   long numDims = s.getNumDims();
   assertTrue(numDims > 0, "CubeSlice s has negative dimension number");

   if (numDims > 1) {
      long dim0 = s.getDim(0);
      for (long i = 0; i < dim0; i++)
         recursiveEval(CubeSlice<NTL::zz_p>(s, i), multiEvalPoints, d+1, tmp1, tmp2);
   }

   long posBnd = s.getProd(1);
   for (long pos = 0; pos < posBnd; pos++) {
      getHyperColumn(tmp1.rep, s, pos);
      tmp1.normalize();
      eval(tmp2, tmp1, multiEvalPoints[d]);
      setHyperColumn(tmp2, s, pos);
   }

}
#endif
/********************************************************************/
#if 0
// Implementation of FFTHelper
#include <helib/bluestein.h>
#include <helib/clonedPtr.h>

FFTHelper::FFTHelper(long _m, NTL::zz_p x)
{
  m = _m;
  m_inv = 1/NTL::conv<NTL::zz_p>(m);
  root = NTL::conv<NTL::zz_p>( SqrRootMod( NTL::conv<NTL::ZZ>(x), NTL::conv<NTL::ZZ>(NTL::zz_p::modulus())) );
    // NOTE: the previous line is a pain because NTL does not have
    // a single-precision variant of SqrRootMod...
  iroot = 1/root;

  phim = 0;
  coprime.SetLength(m);
  for (long i = 0; i < m; i++) {
    coprime[i] = (GCD(i, m) == 1);
    if (coprime[i]) phim++;
  }

  build(phimx, NTL::conv<NTL::zz_pX>( Cyclotomic(m) ));
}


void FFTHelper::FFT(const NTL::zz_pX& f, NTL::Vec<NTL::zz_p>& v) const
{
  tmp = f;
  BluesteinFFT(tmp, m, root, powers, powers_aux, Rb, Rb_aux, Ra);
  v.SetLength(phim);

  for (long i = 0, j = 0; i < m; i++)
    if (coprime[i]) v[j++] = coeff(tmp, i);
}

void FFTHelper::iFFT(NTL::zz_pX& f, const NTL::Vec<NTL::zz_p>& v, bool normalize) const
{
  tmp.rep.SetLength(m);
  for (long i = 0, j = 0; i < m; i++) {
    if (coprime[i]) tmp.rep[i] = v[j++];
  }
  tmp.normalize();

  BluesteinFFT(tmp, m, iroot, ipowers, ipowers_aux, iRb, iRb_aux, Ra);

  rem(f, tmp, phimx);

  if (normalize) f *= m_inv;
}

// powVec[d] = m_d = p_d^{e_d}
// computes multiEvalPoints[d] as an FFTHelper for base^{m/m_d}
void computeMultiEvalPoints(NTL::Vec< copied_ptr<FFTHelper> >& multiEvalPoints,
                            const NTL::zz_p& base,
                            long m,
                            const NTL::Vec<long>& powVec,
                            const NTL::Vec<long>& phiVec)
{
   long k = powVec.length();

   multiEvalPoints.SetLength(k);

   for (long d = 0; d < k; d++) {
      long m_d = powVec[d];
      multiEvalPoints[d].set_ptr(new FFTHelper(m_d, power(base, m/m_d)));
   }

}


void recursiveEval(const CubeSlice<NTL::zz_p>& s,
                   const NTL::Vec< copied_ptr<FFTHelper> >& multiEvalPoints,
                   long d,
                   NTL::zz_pX& tmp1,
                   NTL::Vec<NTL::zz_p>& tmp2)
{
   long numDims = s.getNumDims();
   assertTrue(numDims > 0, "CubeSlice s has negative dimension number");

   if (numDims > 1) {
      long dim0 = s.getDim(0);
      for (long i = 0; i < dim0; i++)
         recursiveEval(CubeSlice<NTL::zz_p>(s, i), multiEvalPoints, d+1, tmp1, tmp2);
   }

   long posBnd = s.getProd(1);
   for (long pos = 0; pos < posBnd; pos++) {
      getHyperColumn(tmp1.rep, s, pos);
      tmp1.normalize();
      multiEvalPoints[d]->FFT(tmp1, tmp2);
      setHyperColumn(tmp2, s, pos);
   }

}

void recursiveInterp(const CubeSlice<NTL::zz_p>& s,
                     const NTL::Vec< copied_ptr<FFTHelper> >& multiEvalPoints,
                     long d,
                     NTL::zz_pX& tmp1,
                     NTL::Vec<NTL::zz_p>& tmp2)
{
  long numDims = s.getNumDims();
  assertTrue(numDims > 0, "CubeSlice s has negative dimension number");

  long posBnd = s.getProd(1);
  for (long pos = 0; pos < posBnd; pos++) {
    getHyperColumn(tmp2, s, pos);
    multiEvalPoints[d]->iFFT(tmp1, tmp2, false); // do not normalize
    setHyperColumn(tmp1.rep, s, pos, NTL::zz_p::zero());
  }

  if (numDims > 1) {
    long dim0 = s.getDim(0);
    for (long i = 0; i < dim0; i++)
      recursiveInterp(CubeSlice<NTL::zz_p>(s, i), multiEvalPoints, d+1, tmp1, tmp2);
  }
}

void interp(HyperCube<NTL::zz_p>& cube,
          const NTL::Vec< copied_ptr<FFTHelper> >& multiEvalPoints)
{
   NTL::zz_pX tmp1;
   NTL::Vec<NTL::zz_p> tmp2;

   recursiveInterp(CubeSlice<NTL::zz_p>(cube), multiEvalPoints, 0, tmp1, tmp2);

   // result is not normalized, so we fix that now...

   long k = multiEvalPoints.length();

   NTL::zz_p m_inv = NTL::conv<NTL::zz_p>(1);
   for (long d = 0; d < k; d++) m_inv *= multiEvalPoints[d]->get_m_inv();

   cube.getData() *= m_inv;
}
#endif

} // namespace helib
