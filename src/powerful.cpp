/* a prelimary test program for playing around with Peikert's "powerful" basis.
 * EXPERIMENTAL CODE, not usable yet
 */
#include "NumbTh.h"
#include "bluestein.h"
#include "cloned_ptr.h"
#include "powerful.h"
#include "hypercube.h"

// Implementation of FFTHelper

FFTHelper::FFTHelper(long _m, zz_p x)
{
  m = _m;
  m_inv = 1/conv<zz_p>(m);
  root = conv<zz_p>( SqrRootMod( conv<ZZ>(x), conv<ZZ>(zz_p::modulus())) );
    // NOTE: the previous line is a pain because NTL does not have
    // a single-precision variant of SqrRootMod...
  iroot = 1/root;

  phim = 0;
  coprime.SetLength(m);
  for (long i = 0; i < m; i++) {
    coprime[i] = (GCD(i, m) == 1); 
    if (coprime[i]) phim++;
  }

  build(phimx, conv<zz_pX>( Cyclotomic(m) ));
}


void FFTHelper::FFT(const zz_pX& f, Vec<zz_p>& v) const
{
  tmp = f;
  BluesteinFFT(tmp, m, root, powers, powers_aux, Rb, Rb_aux, Ra);
  v.SetLength(phim);

  for (long i = 0, j = 0; i < m; i++)
    if (coprime[i]) v[j++] = coeff(tmp, i);
}

void FFTHelper::iFFT(zz_pX& f, const Vec<zz_p>& v, bool normalize) const
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

// returns \prod_d vec[d]
long computeProd(const Vec<long>& vec)
{
  long prod = 1;
  long k = vec.length();
  for (long d = 0; d < k; d++)
    prod = prod * vec[d];
  return prod;
}

// vec[d] = (a_d , b_d)
// returns \prod_d a_d^{b_d}
long computeProd(const Vec< Pair<long, long> >& vec)
{
  long prod = 1;
  long k = vec.length();
  for (long d = 0; d < k; d++) {
    prod = prod * computePow(vec[d]);
  }
  return prod;
}

// factors[d] = (p_d, e_d)
// computes phiVec[d] = phi(p_d^{e_d}) = (p_d-1) p_i{e_d-1}
void computePhiVec(Vec<long>& phiVec, 
                   const Vec< Pair<long, long> >& factors)
{
  long k = factors.length();
  phiVec.SetLength(k);

  for (long d = 0; d < k; d++) 
    phiVec[d] = computePhi(factors[d]);
}

// factors[d] = (p_d, e_d)
// computes powVec[d] = p_d^{e_d}
void computePowVec(Vec<long>& powVec, 
                   const Vec< Pair<long, long> >& factors)
{
  long k = factors.length();
  powVec.SetLength(k);
  for (long d = 0; d < k; d++)
    powVec[d] = computePow(factors[d]);
}


// powVec[d] = p_d^{e_d}, m = \prod_d p_d^{e_d}
// computes divVec[d] = m/p_d^{e_d}
void computeDivVec(Vec<long>& divVec, long m,
                   const Vec<long>& powVec)
{
  long k = powVec.length();
  divVec.SetLength(k);

  for (long d = 0; d < k; d++)
    divVec[d] = m/powVec[d];
}


// divVec[d] = m/p_d^{e_d}, powVec[d] = p^{e_d}
// computes invVec[d] = divVec[d]^{-1} mod powVec[d]
void computeInvVec(Vec<long>& invVec,
                   const Vec<long>& divVec, const Vec<long>& powVec)
{
  long k = divVec.length();
  invVec.SetLength(k);

  for (long d = 0; d < k; d++) {
    long t1 = divVec[d] % powVec[d];
    long t2 = InvMod(t1, powVec[d]);
    invVec[d] = t2;
  }
}

// powVec[d] = p_d^{e_d}
// cycVec[d] = Phi_{p_d^{e_d}}(X) mod p
void computeCycVec(Vec<zz_pX>& cycVec, const Vec<long>& powVec)
{
  long k = powVec.length();
  cycVec.SetLength(k);

  for (long d = 0; d < k; d++) {
    ZZX PhimX = Cyclotomic(powVec[d]);
    cycVec[d] = conv<zz_pX>(PhimX);
  }
}


// m = m_1 ... m_k, m_d = p_d^{e_d} 
// powVec[d] = m_d
// invVec[d] = (m/m_d)^{-1} mod m_d
// computes polyToCubeMap[i] and cubeToPolyMap[i] 
//   where polyToCubeMap[i] is the index of (i_1, ..., i_k)
//   in the cube with CubeSignature longSig = (m_1, ..., m_k),
//   and (i_1, ..., i_k) is the unique tuple satistifying
//   i = i_1 (m/m_1) + ... + i_k (m/m_k) mod m
// and
//   cubeToPolyMap is the inverse map.
void computePowerToCubeMap(Vec<long>& polyToCubeMap,
                           Vec<long>& cubeToPolyMap,
                           long m,
                           const Vec<long>& powVec,
                           const Vec<long>& invVec,
                           const CubeSignature& longSig)
{
   long k = powVec.length();

   polyToCubeMap.SetLength(m);
   cubeToPolyMap.SetLength(m);

   for (long i = 0; i < m; i++) {
      long j = 0;
      for (long d = 0; d < k; d++) {
         long i_d = MulMod((i % powVec[d]), invVec[d], powVec[d]);
         j += i_d * longSig.getProd(d+1);
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
void computeShortToLongMap(Vec<long>& shortToLongMap, 
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
         j += i_d * longSig.getProd(d+1);
      }

      shortToLongMap[i] = j;
   }
}


// Computes the inverse of the shortToLongMap, computed above.
// "undefined" entries are initialzed to -1.
void computeLongToShortMap(Vec<long>& longToShortMap,
                           long m,
                           const Vec<long>& shortToLongMap)
{
   long n = shortToLongMap.length();

   longToShortMap.SetLength(m);

   for (long i = 0; i < m; i++) longToShortMap[i] = -1;

   for (long j = 0; j < n; j++) {
      long i = shortToLongMap[j];
      longToShortMap[i] = j;
   }
}



// This routine recursively reduces each hypercolumn
// in dimension d (viewed as a coeff vector) by Phi_{m_d}(X)
// If one starts with a cube of dimension (m_1, ..., m_k),
// one ends up with a cube that is effectively of dimension
// phi(m_1, ..., m_k). Viewed as an element of the ring
// F_p[X_1,...,X_k]/(Phi_{m_1}(X_1), ..., Phi_{m_k}(X_k)),
// the cube remains unchanged.
void recursiveReduce(const CubeSlice<zz_p>& s, 
                     const Vec<zz_pX>& cycVec, 
                     long d,
                     zz_pX& tmp1,
                     zz_pX& tmp2)
{
   long numDims = s.getNumDims();
   assert(numDims > 0);

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
      for (long i = len; i < deg0; i++) tmp2.rep[i] = 0;

      setHyperColumn(tmp2.rep, s, pos);
   }

   if (numDims == 1) return;

   for (long i = 0; i < deg0; i++) 
      recursiveReduce(CubeSlice<zz_p>(s, i), cycVec, d+1, tmp1, tmp2);

}

// This routine implements the isomorphism from
// F_p[X]/(Phi_m(X)) to F_p[X_1, ..., X_k]/(Phi_{m_1}(X_1), ..., Phi_{m_k}(X_k))
// The input is poly, which must be of degree < m, and the
// output is cube, which is a HyperCube of dimension (phi(m_1), ..., phi(m_k)).
// The caller is responsible to supply "scratch space" in the
// form of a HyperCube tmpCube of dimension (m_1, ..., m_k).
void convertPolyToPowerful(HyperCube<zz_p>& cube, 
                           HyperCube<zz_p>& tmpCube, 
                           const zz_pX& poly,
                           const Vec<zz_pX>& cycVec,
                           const Vec<long>& polyToCubeMap,
                           const Vec<long>& shortToLongMap)
{
   long m = tmpCube.getSize();
   long phim = cube.getSize();
   long n = deg(poly);
 
   assert(n < m);

   for (long i = 0; i <= n; i++)
      tmpCube[polyToCubeMap[i]] = poly[i];

   for (long i = n+1; i < m; i++)
      tmpCube[polyToCubeMap[i]] = 0;

   zz_pX tmp1, tmp2;
   recursiveReduce(CubeSlice<zz_p>(tmpCube), cycVec, 0, tmp1, tmp2);

   for (long i = 0; i < phim; i++)
      cube[i] = tmpCube[shortToLongMap[i]];
}


// This implements the inverse of the above isomorphism.
void convertPowerfulToPoly(zz_pX& poly,
                           const HyperCube<zz_p>& cube,
                           long m,
                           const Vec<long>& shortToLongMap,
                           const Vec<long>& cubeToPolyMap,
                           const zz_pX& phimX)
{
   long phim = cube.getSize();

   zz_pX tmp;

   tmp.SetLength(m);
  
   for (long i = 0; i < m; i++)
      tmp[i] = 0;

   for (long i = 0; i < phim; i++) 
      tmp[cubeToPolyMap[shortToLongMap[i]]] = cube[i];
      // FIXME: these two maps could be composed into a single map

   tmp.normalize();

   rem(poly, tmp, phimX);
}



// powVec[d] = m_d = p_d^{e_d}
// computes multiEvalPoints[d] as a vector of length phi(m_d)
//   containing base^{(m/m_d) j} for j in Z_{m_d}^*
void computeMultiEvalPoints(Vec< Vec<zz_p> >& multiEvalPoints,
                            const zz_p& base,
                            long m,
                            const Vec<long>& powVec,
                            const Vec<long>& phiVec)
{
   long k = powVec.length();

   multiEvalPoints.SetLength(k);

   for (long d = 0; d < k; d++) {
      long m_d = powVec[d];
      long phi_d = phiVec[d];
      long count = 0;

      zz_p pow = conv<zz_p>(1);
      zz_p mult = power(base, m/m_d);

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

// powVec[d] = m_d = p_d^{e_d}
// computes multiEvalPoints[d] as an FFTHelper for base^{m/m_d}
void computeMultiEvalPoints(Vec< copied_ptr<FFTHelper> >& multiEvalPoints,
                            const zz_p& base,
                            long m,
                            const Vec<long>& powVec,
                            const Vec<long>& phiVec)
{
   long k = powVec.length();

   multiEvalPoints.SetLength(k);

   for (long d = 0; d < k; d++) {
      long m_d = powVec[d];
      multiEvalPoints[d].set_ptr(new FFTHelper(m_d, power(base, m/m_d))); 
   }
   
}


// computes linearEvalPoints[i] = base^i, i in Z_m^*
void computeLinearEvalPoints(Vec<zz_p>& linearEvalPoints,
                             const zz_p& base,
                             long m, long phim)
{
   linearEvalPoints.SetLength(phim);
   zz_p pow = conv<zz_p>(1);

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
void computeCompressedIndex(Vec< Vec<long> >& compressedIndex,
                            const Vec<long>& powVec)
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

void computePowToCompressedIndexMap(Vec<long>& powToCompressedIndexMap,
                                    long m,
                                    const Vec<long>& powVec,
                                    const Vec< Vec<long> >& compressedIndex,
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


void recursiveEval(const CubeSlice<zz_p>& s,
                   const Vec< Vec<zz_p> >& multiEvalPoints,
                   long d,
                   zz_pX& tmp1,
                   Vec<zz_p>& tmp2)
{
   long numDims = s.getNumDims();
   assert(numDims > 0);

   if (numDims > 1) {
      long dim0 = s.getDim(0);
      for (long i = 0; i < dim0; i++)
         recursiveEval(CubeSlice<zz_p>(s, i), multiEvalPoints, d+1, tmp1, tmp2);
   }

   long posBnd = s.getProd(1);
   for (long pos = 0; pos < posBnd; pos++) {
      getHyperColumn(tmp1.rep, s, pos);
      tmp1.normalize();
      eval(tmp2, tmp1, multiEvalPoints[d]);
      setHyperColumn(tmp2, s, pos);
   }

}


void recursiveEval(const CubeSlice<zz_p>& s,
                   const Vec< copied_ptr<FFTHelper> >& multiEvalPoints,
                   long d,
                   zz_pX& tmp1,
                   Vec<zz_p>& tmp2)
{
   long numDims = s.getNumDims();
   assert(numDims > 0);

   if (numDims > 1) {
      long dim0 = s.getDim(0);
      for (long i = 0; i < dim0; i++)
         recursiveEval(CubeSlice<zz_p>(s, i), multiEvalPoints, d+1, tmp1, tmp2);
   }

   long posBnd = s.getProd(1);
   for (long pos = 0; pos < posBnd; pos++) {
      getHyperColumn(tmp1.rep, s, pos);
      tmp1.normalize();
      multiEvalPoints[d]->FFT(tmp1, tmp2);
      setHyperColumn(tmp2, s, pos);
   }

}


void recursiveInterp(const CubeSlice<zz_p>& s,
                     const Vec< copied_ptr<FFTHelper> >& multiEvalPoints,
                     long d,
                     zz_pX& tmp1,
                     Vec<zz_p>& tmp2)
{
   long numDims = s.getNumDims();
   assert(numDims > 0);

   long posBnd = s.getProd(1);
   for (long pos = 0; pos < posBnd; pos++) {
      getHyperColumn(tmp2, s, pos);
      multiEvalPoints[d]->iFFT(tmp1, tmp2, false); // do not normalize
      setHyperColumn(tmp1.rep, s, pos, zz_p::zero());
   }

   if (numDims > 1) {
      long dim0 = s.getDim(0);
      for (long i = 0; i < dim0; i++)
         recursiveInterp(CubeSlice<zz_p>(s, i), multiEvalPoints, d+1, tmp1, tmp2);
   }

}


void interp(HyperCube<zz_p>& cube,
          const Vec< copied_ptr<FFTHelper> >& multiEvalPoints)
{
   zz_pX tmp1;
   Vec<zz_p> tmp2;

   recursiveInterp(CubeSlice<zz_p>(cube), multiEvalPoints, 0, tmp1, tmp2);

   // result is not normalized, so we fix that now...

   long k = multiEvalPoints.length();

   zz_p m_inv = conv<zz_p>(1);
   for (long d = 0; d < k; d++) m_inv *= multiEvalPoints[d]->get_m_inv();

   cube.getData() *= m_inv;

} 


// -------------------------------------

        
void mapIndexToPowerful(Vec<long>& pow, long j, const Vec<long>& phiVec)
// this maps an index j in [phi(m)] to a vector
// representing the powerful basis coordinates

{
  long k = phiVec.length();
  long phim = computeProd(phiVec);
  assert(j >= 0 && j < phim);

  pow.SetLength(k);

  for (long i = k-1; i >= 0; i--) {
    pow[i] = j % phiVec[i];
    j = (j - pow[i])/phiVec[i];
  }
}


void mapPowerfulToPoly(ZZX& poly, 
                       const Vec<long>& pow, 
                       const Vec<long>& divVec,
                       long m,
                       const ZZX& phimX)
{
  long k = pow.length();
  assert(divVec.length() == k);

  long j = 0;
  for (long i = 0; i < k; i++)
    j += pow[i] * divVec[i];

  j %= m;

  ZZX f = ZZX(j, 1);

  poly = f % phimX;
}
