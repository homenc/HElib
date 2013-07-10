/* a prelimary test program for playing around
 * with Peikert's "powerful" basis.
 */

#include "NumbTh.h"

#include <cassert>
#include <iomanip>

using namespace std;
using namespace NTL;


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

   long getSize() const { return size; }
   long getNumDims() const { return ndims; }
   long getDim(long d) const { return dims.at(d); }
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
// to a hypercube, which must remain alive duriring the lifetime
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

   long getSize() const;
   long getNumDims() const;
   long getDim(long d) const;
   long getProd(long d) const;
   long getCoord(long i, long d) const;
   long addCoord(long i, long d, long offset) const;

   const T& at(long i) const;
   const T& operator[](long i) const;

};

template<class T>
class CubeSlice : public ConstCubeSlice<T> {
private:
   CubeSlice(); // disabled
public:

   // initialize the slice to the full cube
   explicit CubeSlice(HyperCube<T>& _cube);

   CubeSlice(const CubeSlice<T>& other, long i);

   // deep copy of a slice 
   void copy(const ConstCubeSlice<T>& other) const;

   T& at(long i) const;
   T& operator[](long i) const;

};


template<class T>
class HyperCube {
private:
   const CubeSignature& sig;
   Vec<T> data;

   HyperCube(); // disable default constructor

public:
   HyperCube(const CubeSignature& _sig) : sig(_sig) {
      data.SetLength(sig.getSize());
   }

   // use default copy constructor 

   HyperCube& operator=(const HyperCube<T>& other)
   {
      assert(&this->sig == &other.sig);
      data = other.data;
   }


   const CubeSignature& getSig() const { return sig; }

   long getSize() const { return sig.getSize(); }
   long getNumDims() const { return sig.getNumDims(); }
   long getDim(long d) const { return sig.getDim(d); }
   long getProd(long d) const { return sig.getProd(d);}
   long getCoord(long i, long d) const { return sig.getCoord(i, d); }
   long addCoord(long i, long d, long offset) const { return sig.addCoord(i, d, offset); }

   T& at(long i) { return data.at(i); }
   T& operator[](long i) { return data[i]; }

   const T& at(long i) const { return data.at(i); }
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
   assert(v.length() == n);

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






long computePhi(const Pair<long, long>& x)
{
  long p = x.a;
  long e = x.b;
  return power_long(p, e - 1) * (p-1);
}

long computePow(const Pair<long, long>& x)
{
  long p = x.a;
  long e = x.b;
  return power_long(p, e);
}


long computeProd(const Vec<long>& vec)
{
  long prod = 1;
  long k = vec.length();
  for (long i = 0; i < k; i++)
    prod = prod * vec[i];
  return prod;
}

long computeProd(const Vec< Pair<long, long> >& vec)
{
  long prod = 1;
  long k = vec.length();
  for (long i = 0; i < k; i++) {
    prod = prod * computePow(vec[i]);
  }
  return prod;
}

void computePhiVec(Vec<long>& phiVec, 
                   const Vec< Pair<long, long> >& factors)
{
  long k = factors.length();
  phiVec.SetLength(k);

  for (long i = 0; i < k; i++) 
    phiVec[i] = computePhi(factors[i]);
}

void computePowVec(Vec<long>& powVec, 
                   const Vec< Pair<long, long> >& factors)
{
  long k = factors.length();
  powVec.SetLength(k);
  for (long i = 0; i < k; i++)
    powVec[i] = computePow(factors[i]);
}

void computeDivVec(Vec<long>& divVec, long m,
                   const Vec<long>& powVec)
{
  long k = powVec.length();
  divVec.SetLength(k);

  for (long i = 0; i < k; i++)
    divVec[i] = m/powVec[i];
}


void computeInvVec(Vec<long>& invVec,
                   const Vec<long>& divVec, const Vec<long>& powVec)
{
  long k = divVec.length();
  invVec.SetLength(k);

  for (long i = 0; i < k; i++) {
    long t1 = divVec[i] % powVec[i];
    long t2 = InvMod(t1, powVec[i]);
    invVec[i] = t2;
  }
}

        
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

   Vec<long> dims;
   dims.SetLength(3);

   dims[0] = 3;
   dims[1] = 4;
   dims[2] = 5;

   CubeSignature sig(dims);

   HyperCube<double> c(sig);

   for (long i = 0; i < c.getSize(); i++)
      c.at(i) = i;

   print3D(c);

   CubeSlice<double> s0(c);

   CubeSlice<double> s1(s0, 0);
   CubeSlice<double> s2(s0, 2);

   Vec<double> v;
   getHyperColumn(v, s0, 1);
   setHyperColumn(v, s0, 7);

   print3D(c);
   


   




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


  ZZX phimX = Cyclotomic(m);


  for (long j = 0; j < phim; j++) {
    Vec<long> pow;
    mapIndexToPowerful(pow, j, phiVec);
    ZZX poly;
    mapPowerfulToPoly(poly, pow, divVec, m, phimX);
    cout << pow << "  " << poly << "\n";
  }

}

