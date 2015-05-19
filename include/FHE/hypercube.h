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
/**
 * @file hypercube.h
 * @brief Hypercubes and their slices.
 **/
#ifndef _HYPERCUBE_H_
#define _HYPERCUBE_H_

#include "NumbTh.h"

class PAlgebra; // forward decleration

//! @class CubeSignature
//! @brief Holds a vector of dimensions for a hypercube and some additional data
class CubeSignature {
private:
   Vec<long> dims;  // dims[i] is the size along the i'th diemnsion
   Vec<long> prods; // prods[i] = \prod_{j=i}^{n-1} dims[i]
   long ndims;
   long size;

public:
   CubeSignature(): ndims(0) {} // a NULL signature

   CubeSignature(const Vec<long>& _dims): ndims(0) { initSignature(_dims); }

   //! Build a CubeSignature to reflect the hypercube structure of Zm* /(p)
   explicit CubeSignature(const PAlgebra& alg); // in PAlgebra.cpp

   /* When we get C++11 support, we could #include <initializer_list>
    * and then do e.g., CubeSignature s {1,2,3};
   CubeSignature(std::initializer_list<long> _dims)
   {
      dims.SetLength(_dims.size());
      long i;
      std::initializer_list<long>::iterator it;
      for (i=0, it=_dims.begin(); it!=_dims.end(); ++i, ++it)
	dims[i] = *it;

      [...] // continue as above, ndims = dims.length() etc.
   }
   **********************************************************/

   void initSignature(const Vec<long>& _dims)
   {
     assert(ndims == 0); // can only initialize a NULL signature
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

   //! total size of cube
   long getSize() const { return size; }

   //! number of dimensions
   long getNumDims() const { return ndims; }

   //! size of dimension d
   long getDim(long d) const { return dims.at(d); }

   //! product of sizes of dimensions d, d+1, ...
   long getProd(long d) const { return prods.at(d);}

   //! product of sizes of dimensions from, from+1, ..., to-1
   long getProd(long from, long to) const 
   { return prods.at(from)/prods.at(to); }

   //! get coordinate in dimension d of index i
   long getCoord(long i, long d) const {
      assert(i >= 0 && i < size);
   
      return (i % prods.at(d)) / prods.at(d+1); 
   }

   //! add offset to coordinate in dimension d of index i
   long addCoord(long i, long d, long offset) const {
      assert(i >= 0 && i < size);
      
      offset = offset % dims.at(d);
      if (offset < 0) offset += dims.at(d);

      long i_d = getCoord(i, d);
      long i_d1 = (i_d + offset) % dims.at(d);

      long i1 = i + (i_d1 - i_d) * prods.at(d+1);

      return i1;
   }

   //! number of slices
   long numSlices(long d=1) const { return getProd(0, d); }

   //! size of one slice
   long sliceSize(long d=1) const { return getProd(d); }

   //! number of columns
   long numCols() const { return getProd(1); }

   friend ostream& operator<<(ostream &s, const CubeSignature& sig);
};

inline ostream& operator<<(ostream &s, const CubeSignature& sig)
{ return s << sig.dims; }


//! @class HyperCube
//! @brief A multi-dimensional cube.
//!
//! Such an object is initialzied with a CubeSignature: a reference to the
//! signature is stored with the cube, and so the signature must remain alive
//! during the lifetime of the cube, to prevent dangling pointers.
template<class T>
class HyperCube {
private:
   const CubeSignature& sig;
   Vec<T> data;

   HyperCube(); // disable default constructor

public:
   //! initialzie a HyperCube with a CubeSignature
   HyperCube(const CubeSignature& _sig) : sig(_sig) {
      data.FixLength(sig.getSize());
   }

   // use default copy constructor 

   //! assignment: signatures must be the same
   HyperCube& operator=(const HyperCube<T>& other)
   {
      assert(&this->sig == &other.sig);
      data = other.data;
      return *this;
   }


   //! equality testing: signaturees must be the same
   bool operator==(const HyperCube<T>& other) const
   {
      assert(&this->sig == &other.sig);
      return data == other.data;
   }

   bool operator!=(const HyperCube<T>& other) const
   {
      return !(*this == other);
   }


   //! const ref to signature
   const CubeSignature& getSig() const { return sig; }

   //! read/write ref to the data vector.
   //! Note that the length of data is fixed upon construction,
   //! so it cannot be changed through this ref.
   Vec<T>& getData() { return data; }

   //! read-only ref to data vector
   const Vec<T>& getData() const { return data; } 

   //! total size of cube
   long getSize() const { return sig.getSize(); }

   //! number of dimensions
   long getNumDims() const { return sig.getNumDims(); }

   //! size of dimension d
   long getDim(long d) const { return sig.getDim(d); }

   //! product of sizes of dimensions d, d+1, ...
   long getProd(long d) const { return sig.getProd(d);}

   //! product of sizes of dimensions from, from+1, ..., to-1
   long getProd(long from, long to) const { return sig.getProd(from,to);} 

   //! get coordinate in dimension d of index i
   long getCoord(long i, long d) const { return sig.getCoord(i, d); }

   //! add offset to coordinate in dimension d of index i
   long addCoord(long i, long d, long offset) const { return sig.addCoord(i, d, offset); }

   //! number of slices
   long numSlices(long d=1) const { return getProd(0, d); }

   //! size of one slice
   long sliceSize(long d=1) const { return getProd(d); }

   //! number of columns
   long numCols() const { return getProd(1); }

   //! reference to element at position i, with bounds check 
   T& at(long i) { return data.at(i); }

   //! reference to element at position i, without bounds check 
   T& operator[](long i) { return data[i]; }

   //! read-only reference to element at position i, with bounds check 
   const T& at(long i) const { return data.at(i); }

   //! read-only reference to element at position i, without bounds check 
   const T& operator[](long i) const { return data[i]; }

  //! @brief rotate k positions along the i'th dimension
   void rotate1D(long i, long k);

  //! @brief Shift k positions along the i'th dimension with zero fill
   void shift1D(long i, long k);
};

//! @class ConstCubeSlice
//! @brief A constant lower-dimension slice of a hypercube
//!
//! A ConstCubeSlice acts like a pointer to a lower dimensional constant
//! subcube of a hypercube. It is initialized using a reference to a hypercube,
//! which must remain alive during the lifetime of the slice, to prevent
//! dangling pointers.
//! The subclass CubeSlice works also with non-constant cubes and subcubes.
//! 
//! In addition, for greater flexibility, a "slice" may be initialized
//! with a vector and a signature, rather than a cube

template<class T>
class ConstCubeSlice {
private:
   const Vec<T>* data;
   const CubeSignature* sig;
   long dimOffset; // # of "missing dimensions" is this slice vs. the full cube
   long sizeOffset; // index in the cube of the first element in this slice

   ConstCubeSlice(); // disabled

public:

   //! initialize the slice to the full cube
   explicit ConstCubeSlice(const HyperCube<T>& _cube) {
     data = &_cube.getData(); 
     sig = &_cube.getSig(); 
     dimOffset = 0; 
     sizeOffset = 0; 
   }

   ConstCubeSlice(const Vec<T>& _data, const CubeSignature& _sig) { 
     assert(_data.length() == _sig.getSize());
     data = &_data;
     sig = &_sig;
     dimOffset = 0; 
     sizeOffset = 0; 
   }

   //! initialize the slice to point to the i-th subcube (with some
   //! given dimension offset) of the cube pointed to by _cube or bigger.
   ConstCubeSlice(const ConstCubeSlice& bigger, long i, long _dimOffset=1);
   ConstCubeSlice(const HyperCube<T>& _cube, long i, long _dimOffset=1);

   // use default copy constructor and assignment operators,
   // which means shallow copy

   // the following mimic the corresponding methods
   // in the HyperCube class, restricted to the slice

   //! const ref to signature
    const CubeSignature& getSig() const { return *sig; }

   //! total size 
   long getSize() const { return sig->getProd(dimOffset); }

   //! number of dimensions 
   long getNumDims() const { return sig->getNumDims() - dimOffset; }

   //! size of dimension d
   long getDim(long d) const { return sig->getDim(d + dimOffset); }

   //! product of sizes of dimensions d, d+1, ...
   long getProd(long d) const { return sig->getProd(d + dimOffset); }

   //! product of sizes of dimensions from, from+1, ..., to-1
   long getProd(long from, long to) const 
   { return sig->getProd(from + dimOffset, to + dimOffset); } 

   //! get coordinate in dimension d of index i
   long getCoord(long i, long d) const {
     assert(i >= 0 && i < getSize());
     return sig->getCoord(i + sizeOffset, d + dimOffset);
   }

   //! add offset to coordinate in dimension d of index i
   long addCoord(long i, long d, long offset) const {
     assert(i >= 0 && i < getSize());
     return sig->addCoord(i + sizeOffset, d + dimOffset, offset);
   }

   //! number of slices
   long numSlices(long d=1) const { return getProd(0, d); }

   //! size of one slice
   long sliceSize(long d=1) const { return getProd(d); }

   //! number of columns
   long numCols() const { return getProd(1); }

   //! read-only reference to element at position i, with bounds check 
   const T& at(long i) const {
     assert(i >= 0 && i < getSize());
     return (*data)[i + sizeOffset];
   }

   //! read-only reference to element at position i, without bounds check 
   const T& operator[](long i) const { return (*data)[i + sizeOffset]; }
};


//! @class CubeSlice
//! @brief A lower-dimension slice of a hypercube
template<class T>
class CubeSlice : public ConstCubeSlice<T> {
private:
   CubeSlice(); // disabled
public:

   // initialize the slice to the full cube
   explicit CubeSlice(HyperCube<T>& _cube) : ConstCubeSlice<T>(_cube) {}

   CubeSlice(Vec<T>& _data, const CubeSignature& _sig) : ConstCubeSlice<T>(_data, _sig) {}

   // initialize the slice to point to the i-th subcube
   // of the cube pointed to by bigger
   CubeSlice(const CubeSlice<T>& bigger, long i, long _dimOffset=1) :
   ConstCubeSlice<T>(bigger, i, _dimOffset) {}

   CubeSlice(HyperCube<T>& _cube, long i, long _dimOffset=1) :
   ConstCubeSlice<T>(_cube, i, _dimOffset) {}

   // deep copy of a slice: copies other into this 
   void copy(const ConstCubeSlice<T>& other) const;

   // reference to element at position i, with bounds check 
   T& at(long i) const {
     return const_cast<T&>(this->ConstCubeSlice<T>::at(i));
   }

   // reference to element at position i, without bounds check 
   T& operator[](long i) const {
     return const_cast<T&>(this->ConstCubeSlice<T>::operator[](i));
   }
};


//! getHyperColumn reads out a (multi-dimensional) column from a slice. The
//! parameter pos specifies the position of the column, which must be in the
//! range 0 <= pos < s.getProd(1). The vector v is filled with values whose
//! coordinate in the lower dimensional subcube is equal to pos. The length
//! of v will be set to s.getDim(0).
template<class T>
void getHyperColumn(Vec<T>& v, const ConstCubeSlice<T>& s, long pos);

//! setHyperColumn does the reverse of getHyperColumn, setting the column
//! to the given vector
template<class T>
void setHyperColumn(const Vec<T>& v, const CubeSlice<T>& s, long pos);

//! this version of setHyperColumn implicitly pads v with a default value,
//! if v is too short
template<class T>
void setHyperColumn(const Vec<T>& v, const CubeSlice<T>& s, long pos, const T& val);

template<class T>
void print3D(const HyperCube<T>& c);

#endif /* ifndef _HYPERCUBE_H_ */
