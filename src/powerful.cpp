/* a prelimary test program for playing around with Peikert's "powerful" basis.
 * EXPERIMENTAL CODE, not usable yet
 */

#include "NumbTh.h"
#include <NTL/lzz_pX.h>

#include <cassert>
#include <iomanip>

using namespace std;
using namespace NTL;


// class CubSignature: such an object is initialized
// with a vector of dimensions for a hypercube, and
// some auxilliary data is computed

class CubeSignature {
private:
   Vec<long> dims;  // dims[i] is the size along the i'th diemnsion
   Vec<long> prods; // prods[i] = \prod_{j=i}^{n-1} dims[i]
   long ndims;
   long size;

   CubeSignature(); // disabled

public:
   CubeSignature(const Vec<long>& _dims)
   {
      dims = _dims;
      ndims = dims.length();
      assert(ndims > 0);

      
      prods.SetLength(ndims+1);
      prods[ndims] = 1;
      for (long i = ndims-1; i >= 0; i--) {
         assert(dims[i] > 0);
         prods[i] = dims[i]*prods[i+1];
      }

      size = prods[0];
   }

   // total size of cube
   long getSize() const { return size; }

   // number of dimensions
   long getNumDims() const { return ndims; }

   // size of dimension d
   long getDim(long d) const { return dims.at(d); }

   // product of sizes of dimensions d, d+1, ...
   long getProd(long d) const { return prods.at(d);}

   // get coordinate in dimension d of index i
   long getCoord(long i, long d) const {
      assert(i >= 0 && i < size);
   
      return (i % prods.at(d)) / prods.at(d+1); 
   }

   // add offset to coordinate in dimension d of index i
   long addCoord(long i, long d, long offset) const {
      assert(i >= 0 && i < size);
      
      offset = offset % dims.at(d);
      if (offset < 0) offset += dims.at(d);

      long i_d = getCoord(i, d);
      long i_d1 = (i_d + offset) % dims.at(d);

      long i1 = i + (i_d1 - i_d) * prods.at(d+1);

      return i1;
   }

};
   


template<class T>
class HyperCube;  // forward reference


// A ConstCubeSlice acts like a pointer to a lower dimensional
// constant subcube of a hypercube.  It is initialized using a reference
// to a hypercube, which must remain alive during the lifetime
// of the slice, to prevent dangling pointers.

// The subclass CubeSlice works with non-constant cubes and subcubes.

template<class T>
class ConstCubeSlice {
private:
   const HyperCube<T>* cube;
   long dimOffset;
   long sizeOffset;

   ConstCubeSlice(); // disabled

public:

   // initialize the slice to the full cube
   explicit ConstCubeSlice(const HyperCube<T>& _cube);

   // initialize the slice to point to the i-th subcube
   // of the cube pointed to by other
   ConstCubeSlice(const ConstCubeSlice& other, long i);

   // use default copy constructor and assignment operators,
   // which means shallow copy


   // the following mimic the corresponding methods
   // in the HyperCube class, restricted to the slice

   // total size 
   long getSize() const;

   // number of dimensions 
   long getNumDims() const;

   // size of dimension d
   long getDim(long d) const;

   // product of sizes of dimensions d, d+1, ...
   long getProd(long d) const;

   // get coordinate in dimension d of index i
   long getCoord(long i, long d) const;

   // add offset to coordinate in dimension d of index i
   long addCoord(long i, long d, long offset) const;

   // read-only reference to element at position i, with bounds check 
   const T& at(long i) const;

   // read-only reference to element at position i, without bounds check 
   const T& operator[](long i) const;

};

template<class T>
class CubeSlice : public ConstCubeSlice<T> {
private:
   CubeSlice(); // disabled
public:

   // initialize the slice to the full cube
   explicit CubeSlice(HyperCube<T>& _cube);

   // initialize the slice to point to the i-th subcube
   // of the cube pointed to by other
   CubeSlice(const CubeSlice<T>& other, long i);

   // deep copy of a slice: copies other into this 
   void copy(const ConstCubeSlice<T>& other) const;

   // reference to element at position i, with bounds check 
   T& at(long i) const;

   // reference to element at position i, without bounds check 
   T& operator[](long i) const;

};

// The class HyperCube<T> represents a multi-dimensional cube.
// Such an object is initialzied with a CubeSignature: a reference
// to the signature is stored with the cube, and so the signature
// must remain alive during the lifetime of the cube, to
// prevent dangling pointers.



template<class T>
class HyperCube {
private:
   const CubeSignature& sig;
   Vec<T> data;

   HyperCube(); // disable default constructor

public:
   // initialzie a HyperCube with a CubeSignature
   HyperCube(const CubeSignature& _sig) : sig(_sig) {
      data.SetLength(sig.getSize());
   }

   // use default copy constructor 

   // assignment: signatures must be the same
   HyperCube& operator=(const HyperCube<T>& other)
   {
      assert(&this->sig == &other.sig);
      data = other.data;
   }


   // const ref to signature
   const CubeSignature& getSig() const { return sig; }

   // total size of cube
   long getSize() const { return sig.getSize(); }

   // number of dimensions
   long getNumDims() const { return sig.getNumDims(); }

   // size of dimension d
   long getDim(long d) const { return sig.getDim(d); }

   // product of sizes of dimensions d, d+1, ...
   long getProd(long d) const { return sig.getProd(d);}

   // get coordinate in dimension d of index i
   long getCoord(long i, long d) const { return sig.getCoord(i, d); }

   // add offset to coordinate in dimension d of index i
   long addCoord(long i, long d, long offset) const { return sig.addCoord(i, d, offset); }

   // reference to element at position i, with bounds check 
   T& at(long i) { return data.at(i); }

   // reference to element at position i, without bounds check 
   T& operator[](long i) { return data[i]; }

   // read-only reference to element at position i, with bounds check 
   const T& at(long i) const { return data.at(i); }

   // read-only reference to element at position i, without bounds check 
   const T& operator[](long i) const { return data[i]; }

};



// Implementation of ConstCubeSlice

template<class T>
ConstCubeSlice<T>::ConstCubeSlice(const HyperCube<T>& _cube)
{
   cube = &_cube;
   dimOffset = 0;
   sizeOffset = 0;
}


template<class T>
ConstCubeSlice<T>::ConstCubeSlice(const ConstCubeSlice<T>& other, long i) 
{
   cube = other.cube;
   dimOffset = other.dimOffset + 1;

   assert(dimOffset <= cube->getNumDims());
   // allow zero-dimensional slice

   assert(i >= 0 && i < cube->getDim(other.dimOffset));

   sizeOffset = other.sizeOffset + i*cube->getProd(dimOffset);
}


template<class T>
long ConstCubeSlice<T>::getSize() const 
{
   return cube->getProd(dimOffset);
}


template<class T>
long ConstCubeSlice<T>::getNumDims() const
{
   return cube->getNumDims() - dimOffset;
}


template<class T>
long ConstCubeSlice<T>::getDim(long d) const
{
   return cube->getDim(d + dimOffset);
}


template<class T>
long ConstCubeSlice<T>::getProd(long d) const
{
   return cube->getProd(d + dimOffset);
}

template<class T>
long ConstCubeSlice<T>::getCoord(long i, long d) const
{
   assert(i >= 0 && i < getSize());
   return cube->getCoord(i + sizeOffset, d + dimOffset);
}

template<class T>
long ConstCubeSlice<T>::addCoord(long i, long d, long offset) const
{
   assert(i >= 0 && i < getSize());
   return cube->addCoord(i + sizeOffset, d + dimOffset, offset);
}


template<class T>
const T& ConstCubeSlice<T>::at(long i) const
{
   assert(i >= 0 && i < getSize());
   return (*cube)[i + sizeOffset];
}


template<class T>
const T& ConstCubeSlice<T>::operator[](long i) const
{
   return (*cube)[i + sizeOffset];
}


// Implementation of CubeSlice

template<class T>
CubeSlice<T>::CubeSlice(HyperCube<T>& _cube) : ConstCubeSlice<T>(_cube) {}

template<class T>
CubeSlice<T>::CubeSlice(const CubeSlice<T>& other, long i) : ConstCubeSlice<T>(other, i) {}

template<class T>
T& CubeSlice<T>::at(long i) const
{
   return const_cast<T&>(this->ConstCubeSlice<T>::at(i));
}

template<class T>
T& CubeSlice<T>::operator[](long i) const
{
   return const_cast<T&>(this->ConstCubeSlice<T>::operator[](i));
}

template<class T>
void CubeSlice<T>::copy(const ConstCubeSlice<T>& other) const
{
   long n = this->getSize();

   // we only check that the sizes match
   assert(n == other.getSize());

   T *dst = &(*this)[0];
   const T *src = &other[0];

   for (long i = 0; i < n; i++)
      dst[i] = src[i];
}

// getHyperColumn reads out a (multi-dimensional) from
// a slice.  The parameter pos specifies the position of the column,
// which must be in the range 0 <= pos < s.getProd(1).
// The vector v is filled with values whose coordinate in the lower
// dimensional subcube is equal to pos.  The length of v will be
// set to s.getDim(0).

template<class T>
void getHyperColumn(Vec<T>& v, const ConstCubeSlice<T>& s, long pos)
{
   long m = s.getProd(1);
   long n = s.getDim(0);

   assert(pos >= 0 && pos < m);
   v.SetLength(n);

   for (long i = 0; i < n; i++)
      v[i] = s[pos + i*m];
}

// setHyperColumn does the reverse of getHyperColumn, setting the column
// to the given vector

template<class T>
void setHyperColumn(const Vec<T>& v, const CubeSlice<T>& s, long pos)
{
   long m = s.getProd(1);
   long n = s.getDim(0);

   assert(pos >= 0 && pos < m);
   if (v.length() < n) n = v.length();

   for (long i = 0; i < n; i++)
      s[pos + i*m] = v[i];
}



template<class T>
void print3D(const HyperCube<T>& c) 
{
   assert(c.getNumDims() == 3);

   ConstCubeSlice<T> s0(c);

   for (long i = 0; i < s0.getDim(0); i++) {
      
      ConstCubeSlice<T> s1(s0, i);
      for (long j = 0; j < s1.getDim(0); j++) {

         ConstCubeSlice<T> s2(s1, j);
         for (long k = 0; k < s2.getDim(0); k++)
            cout << setw(3) << s2.at(k);

         cout << "\n";
      }

      cout << "\n";
   }
}


// x = p^e
// returns phi(p^e) = (p-1) p^{e-1}
long computePhi(const Pair<long, long>& x)
{
  long p = x.a;
  long e = x.b;
  return power_long(p, e - 1) * (p-1);
}

// x = (p, e)
// returns p^e
long computePow(const Pair<long, long>& x)
{
  long p = x.a;
  long e = x.b;
  return power_long(p, e);
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


// computes linearEvalPoints[i] = base^i, i = 0..m-1
void computeLinearEvalPoints(Vec<zz_p>& linearEvalPoints,
                             const zz_p& base,
                             long m)
{
   linearEvalPoints.SetLength(m);
   zz_p pow = conv<zz_p>(1);

   for (long i = 0; i < m; i++) {
      linearEvalPoints[i] = pow;
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


void eval(HyperCube<zz_p>& cube,
          const Vec< Vec<zz_p> >& multiEvalPoints)
{
   zz_pX tmp1;
   Vec<zz_p> tmp2;

   recursiveEval(CubeSlice<zz_p>(cube), multiEvalPoints, 0, tmp1, tmp2);
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

void usage()
{
  cerr << "bad args\n";
  exit(0);
}

int main(int argc, char *argv[])
{

  argmap_t argmap;

  argmap["q"] = "101";
  argmap["m"] = "100";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage();

  long q = atoi(argmap["q"]);
  long m = atoi(argmap["m"]);

  cout << "q=" << q << "\n";
  cout << "m=" << m << "\n";

  Vec< Pair<long, long> > factors;

  factorize(factors, m);

  cout << factors << "\n";

  Vec<long> phiVec;
  computePhiVec(phiVec, factors);
  cout << phiVec << "\n";

  long phim = computeProd(phiVec);
  cout << phim << "\n";

  Vec<long> powVec;
  computePowVec(powVec, factors);
  cout << powVec << "\n";

  Vec<long> divVec;
  computeDivVec(divVec, m, powVec);
  cout << divVec << "\n";

  Vec<long> invVec;
  computeInvVec(invVec, divVec, powVec);
  cout << invVec << "\n";


  CubeSignature shortSig(phiVec);
  CubeSignature longSig(powVec);

  Vec<long> polyToCubeMap;
  Vec<long> cubeToPolyMap;
  computePowerToCubeMap(polyToCubeMap, cubeToPolyMap, m, powVec, invVec, longSig);
  cout << polyToCubeMap << "\n";
  cout << cubeToPolyMap << "\n";

  Vec<long> shortToLongMap;
  computeShortToLongMap(shortToLongMap, shortSig, longSig);
  cout << shortToLongMap << "\n";

  Vec<long> longToShortMap;
  computeLongToShortMap(longToShortMap, m, shortToLongMap);
  cout << longToShortMap << "\n";
  

  zz_p::init(q);

  Vec<zz_pX> cycVec;
  computeCycVec(cycVec, powVec);
  cout << cycVec << "\n";


  ZZX PhimX = Cyclotomic(m);
  zz_pX phimX = conv<zz_pX>(PhimX);

  cout << phimX << "\n";

  zz_pX poly;
  random(poly, phim);

  HyperCube<zz_p> cube(shortSig);
  HyperCube<zz_p> tmpCube(longSig);

  zz_pX tmp1, tmp2;

  convertPolyToPowerful(cube, tmpCube, poly, cycVec, 
                        polyToCubeMap, shortToLongMap);

  zz_pX poly1;

  convertPowerfulToPoly(poly1, cube, m, shortToLongMap, cubeToPolyMap, phimX);

  if (poly == poly1)
    cout << ":-)\n";
  else {
    cout << ":-(\n";
    cout << poly << "\n";
    cout << poly1 << "\n";
  }

  long lbase = 1;
  while (lbase < q && multOrd(lbase, q) != m) lbase++;

  assert(lbase < q);

  zz_p base = conv<zz_p>(lbase);

  Vec< Vec<zz_p> > multiEvalPoints;
  computeMultiEvalPoints(multiEvalPoints, base, m, powVec, phiVec);

  Vec<zz_p> linearEvalPoints;
  computeLinearEvalPoints(linearEvalPoints, base, m);

  Vec< Vec<long> > compressedIndex;
  computeCompressedIndex(compressedIndex, powVec);

  Vec<long> powToCompressedIndexMap;
  computePowToCompressedIndexMap(powToCompressedIndexMap, m, powVec, 
                                 compressedIndex, shortSig);

  eval(cube, multiEvalPoints);

  Vec<zz_p> res;
  eval(res, poly, linearEvalPoints);

  bool eval_ok = true;

  for (long i = 0; i < m; i++) {
    if (GCD(i, m) == 1) {
      long j = powToCompressedIndexMap[i];
      if (cube[j] != res[i]) eval_ok = false;
    }
  } 


  if (eval_ok)
    cout << "eval ok\n";
  else
    cout << "eval not ok\n";

}

