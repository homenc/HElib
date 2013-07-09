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
   vector<long> dims;  // dims[i] is the size along the i'th diemnsion
   vector<long> prods; // prods[i] = \prod_{j=i}^{n-1} dims[i]
   long ndims;
   long size;

   CubeSignature(); // disabled

public:
   CubeSignature(const vector<long>& _dims)
   {
      dims = _dims;
      ndims = dims.size();
      assert(ndims > 0);

      
      prods.resize(ndims+1);
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
   


class HyperCube;  // forward reference


// A ConstCubeSlice acts like a pointer to a lower dimensional
// constant subcube of a hypercube.  It is initialized using a reference
// to a hypercube, which must remain alive duriring the lifetime
// of the slice, to prevent dangling pointers.

// The subclass CubeSlice works with non-constant cubes and subcubes.

class ConstCubeSlice {
private:
   const HyperCube* cube;
   long dimOffset;
   long sizeOffset;

   ConstCubeSlice(); // disabled

public:

   // initialize the slice to the full cube
   explicit ConstCubeSlice(const HyperCube& _cube);

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

   const double& at(long i) const;
   const double& operator[](long i) const;

};

class CubeSlice : public ConstCubeSlice {
private:
   CubeSlice(); // disabled
public:

   // initialize the slice to the full cube
   explicit CubeSlice(HyperCube& _cube);

   CubeSlice(const CubeSlice& other, long i);

   // deep copy of a slice 
   void copy(const ConstCubeSlice& other) const;

   double& at(long i) const;
   double& operator[](long i) const;

};


class HyperCube {
private:
   const CubeSignature& sig;
   vector<double> data;

   HyperCube(); // disable default constructor

public:
   HyperCube(const CubeSignature& _sig) : sig(_sig) {
      data.resize(sig.getSize());
   }

   // use default copy constructor 

   HyperCube& operator=(const HyperCube& other)
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

   double& at(long i) { return data.at(i); }
   double& operator[](long i) { return data[i]; }

   const double& at(long i) const { return data.at(i); }
   const double& operator[](long i) const { return data[i]; }

};



// Implementation of ConstCubeSlice

ConstCubeSlice::ConstCubeSlice(const HyperCube& _cube)
{
   cube = &_cube;
   dimOffset = 0;
   sizeOffset = 0;
}


ConstCubeSlice::ConstCubeSlice(const ConstCubeSlice& other, long i) 
{
   cube = other.cube;
   dimOffset = other.dimOffset + 1;

   assert(dimOffset <= cube->getNumDims());
   // allow zero-dimensional slice

   assert(i >= 0 && i < cube->getDim(other.dimOffset));

   sizeOffset = other.sizeOffset + i*cube->getProd(dimOffset);
}


long ConstCubeSlice::getSize() const 
{
   return cube->getProd(dimOffset);
}


long ConstCubeSlice::getNumDims() const
{
   return cube->getNumDims() - dimOffset;
}


long ConstCubeSlice::getDim(long d) const
{
   return cube->getDim(d + dimOffset);
}


long ConstCubeSlice::getProd(long d) const
{
   return cube->getProd(d + dimOffset);
}

long ConstCubeSlice::getCoord(long i, long d) const
{
   assert(i >= 0 && i < getSize());
   return cube->getCoord(i + sizeOffset, d + dimOffset);
}

long ConstCubeSlice::addCoord(long i, long d, long offset) const
{
   assert(i >= 0 && i < getSize());
   return cube->addCoord(i + sizeOffset, d + dimOffset, offset);
}


const double& ConstCubeSlice::at(long i) const
{
   assert(i >= 0 && i < getSize());
   return (*cube)[i + sizeOffset];
}

const double& ConstCubeSlice::operator[](long i) const
{
   return (*cube)[i + sizeOffset];
}


// Implementation of CubeSlice

CubeSlice::CubeSlice(HyperCube& _cube) : ConstCubeSlice(_cube) {}
CubeSlice::CubeSlice(const CubeSlice& other, long i) : ConstCubeSlice(other, i) {}

double& CubeSlice::at(long i) const
{
   return const_cast<double&>(this->ConstCubeSlice::at(i));
}

double& CubeSlice::operator[](long i) const
{
   return const_cast<double&>(this->ConstCubeSlice::operator[](i));
}

void CubeSlice::copy(const ConstCubeSlice& other) const
{
   long n = getSize();

   // we only check that the sizes match
   assert(n == other.getSize());

   double *dst = &(*this)[0];
   const double *src = &other[0];

   for (long i = 0; i < n; i++)
      dst[i] = src[i];
}



void print3D(const HyperCube& c) 
{
   assert(c.getNumDims() == 3);

   ConstCubeSlice s0(c);

   for (long i = 0; i < s0.getDim(0); i++) {
      
      ConstCubeSlice s1(s0, i);
      for (long j = 0; j < s1.getDim(0); j++) {

         ConstCubeSlice s2(s1, j);
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

void computeDivVec(Vec<long>& divVec,
                   const Vec< Pair<long, long> >& factors)
{
  long k = factors.length();
  divVec.SetLength(k);

  long m = computeProd(factors);

  for (long i = 0; i < k; i++)
    divVec[i] = m/computePow(factors[i]);
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

   vector<long> dims;
   dims.resize(3);

   dims[0] = 3;
   dims[1] = 4;
   dims[2] = 5;

   CubeSignature sig(dims);

   HyperCube c(sig);

   for (long i = 0; i < c.getSize(); i++)
      c.at(i) = i;

   print3D(c);

   CubeSlice s0(c);

   CubeSlice s1(s0, 0);
   CubeSlice s2(s0, 2);

   s2.copy(s1);

   print3D(c);
   

   exit(0);

   




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

  Vec<long> divVec;
  computeDivVec(divVec, factors);
  cout << divVec << "\n";


  ZZX phimX = Cyclotomic(m);


  for (long j = 0; j < phim; j++) {
    Vec<long> pow;
    mapIndexToPowerful(pow, j, phiVec);
    ZZX poly;
    mapPowerfulToPoly(poly, pow, divVec, m, phimX);
    cout << pow << "  " << poly << "\n";
  }

}

