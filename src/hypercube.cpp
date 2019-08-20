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
#include "hypercube.h"
#include <iomanip>

NTL_CLIENT

  //! Break an index into the hypercube to index of the
  //! dimension-dim subcube and index inside that subcube.
std::pair<long,long> CubeSignature::breakIndexByDim(long idx, long dim) const
{
  std::pair<long,long> ans;
  ans.first = (idx%prods[dim+1]) + (idx/prods[dim])*prods[dim+1];
  ans.second = (idx % prods[dim]) / prods[dim+1];
  return ans;
}

   //! The inverse of breakIndexByDim
long CubeSignature::assembleIndexByDim(std::pair<long,long> idx, long dim) const
{
  long third = idx.first % prods[dim+1];
  idx.first = (idx.first - third) * dims[dim];
  return idx.first + idx.second*prods[dim+1] + third;
}

// Rotate k positions along the d'th dimension: The content of slot j that has
// coordinate j_d in the d'th dimension is moved to j' that has coordinates the
// same as j in all the dimensions except d, and has coordinate at dimension d
// j'_d = j_d +k mod sz (where sz is the size of the d'th dimension).
template<class T>
void HyperCube<T>::rotate1D(long d, long k)
{
  //OLD: assert(d >=0 && d < getNumDims());
  helib::assertInRange(d, 0l, getNumDims(), "d must be between 0 and number of dimensions", true);

  // Make sure rotation amount is in the range [1,dimSize-1]
  k %= getDim(d);
  if (k == 0) return;
  if (k < 0) k += getDim(d);

  // A simple implementation with a temporary vector
  Vec<T> tmp(INIT_SIZE, getSize());
  for (long i=0; i<getSize(); i++) 
    tmp[addCoord(i, d, k)] = data[i];
  for (long i=0; i<getSize(); i++)
    data[i] = tmp[i];
}

// Shift k positions along the d'th dimension with zero fill
template<class T>
void HyperCube<T>::shift1D(long d, long k)
{
  //OLD: assert(d >=0 && d < getNumDims());
  helib::assertInRange(d, 0l, getNumDims(), "d must be between 0 and number of dimensions");

  // Make sure rotation amount is in the range [1,dimSize-1]
  bool negative = (k<0);
  k %= getDim(d);
  if (k == 0) return;
  if (k < 0) k += getDim(d);

  if (negative) { // left shift, zeros at the end
    for (long i=getSize()-1; i>=0; i--) {
      long iPrime = addCoord(i, d, k);
      if (iPrime < i) data[iPrime] = data[i];
      else            data[iPrime] = T(); // empty
    }
  } else {        // right shift, zeros at the beginning
    for (long i=0; i<getSize(); i++) {
      long iPrime = addCoord(i, d, k);
      if (iPrime > i) data[iPrime] = data[i];
      else            data[iPrime] = T(); // empty
    }
  }
}

// initialize the slice to point to the i-th subcube
// of the cube pointed to by bigger
template<class T>
ConstCubeSlice<T>::ConstCubeSlice(const ConstCubeSlice<T>& bigger, long i,
				  long dOffset)
{
  //OLD: assert(dOffset >= 0 && dOffset <= bigger.getNumDims()); 
  helib::assertInRange(dOffset, 0l, bigger.getNumDims(), "dOffset must be between 0 and bigger.getNumDims()", true); 
  // allow zero-dimensional slice

  //OLD: assert(i >= 0 && i < bigger.getProd(0, dOffset));
  helib::assertInRange(i, 0l, bigger.getProd(0, dOffset), "i must be between 0 and bigger.getProd(0, dOffset)"); 

   data = bigger.data;
   sig = bigger.sig;
   dimOffset = dOffset + bigger.dimOffset;
   sizeOffset = bigger.sizeOffset + i*bigger.getProd(dOffset);
}

template<class T>
ConstCubeSlice<T>::ConstCubeSlice(const HyperCube<T>& _cube, long i,
				  long dOffset)
{
  //OLD: assert(dOffset >= 0 && dOffset <= _cube.getNumDims()); // allow zero-dimensional slice
  helib::assertInRange(dOffset, 0l, _cube.getNumDims(), "dOffset must be non-negative and at most _cube.getNumDims()", true);
  //OLD: assert(i >= 0 && i < _cube.getProd(0,dOffset));
  helib::assertInRange(i, 0l, _cube.getProd(0, dOffset), "i must be non-negative and at most _cube.getProd(0, dOffset)");

   data = &_cube.getData();
   sig = &_cube.getSig();
   dimOffset = dOffset;
   sizeOffset = i*_cube.getProd(dimOffset);
}

// deep copy of a slice: copies other into this 
template<class T>
void CubeSlice<T>::copy(const ConstCubeSlice<T>& other) const
{
   long n = this->getSize();

   // we only check that the sizes match
   //OLD: assert(n == other.getSize());
   helib::assertEq(n, other.getSize(), "Cube sizes do not match");

   T *dst = &(*this)[0];
   const T *src = &other[0];

   for (long i = 0; i < n; i++)
      dst[i] = src[i];
}

// getHyperColumn reads out a (multi-dimensional) column from
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

   //OLD: assert(pos >= 0 && pos < m);
   helib::assertInRange(pos, 0l, m, "pos must be between 0 and s.getProd(1)");
   v.SetLength(n);

   T* vp = &v[0];
   const T* sp = &s[0];

   for (long i = 0; i < n; i++)
      vp[i] = sp[pos + i*m];
}

// setHyperColumn does the reverse of getHyperColumn, setting the column
// to the given vector

template<class T>
void setHyperColumn(const Vec<T>& v, const CubeSlice<T>& s, long pos)
{
   long m = s.getProd(1);
   long n = s.getDim(0);

   //OLD: assert(pos >= 0 && pos < m);
   helib::assertInRange(pos, 0l, m, "pos must be between 0 and s.getProd(1)");
   if (v.length() < n) n = v.length();

   const T* vp = &v[0];
   T* sp = &s[0];

   for (long i = 0; i < n; i++)
      sp[pos + i*m] = vp[i];
}


// this version of setHyperColumn implicitly pads v with a default value,
// if v is too short

template<class T>
void setHyperColumn(const Vec<T>& v, const CubeSlice<T>& s, long pos, const T& val)
{
   long m = s.getProd(1);
   long n = s.getDim(0);
   long n1 = n;

   //OLD: assert(pos >= 0 && pos < m);
   helib::assertInRange(pos, 0l, m, "pos must be between 0 and s.getProd(1)");
   if (v.length() < n) n1 = v.length();

   const T* vp = &v[0];
   T* sp = &s[0];

   for (long i = 0; i < n1; i++)
      sp[pos + i*m] = vp[i];

   for (long i = n1; i < n; i++)
      sp[pos + i*m] = val;
}



template<class T>
void print3D(const HyperCube<T>& c) 
{
   //OLD: assert(c.getNumDims() == 3);
  helib::assertEq(c.getNumDims(), 3l, "Cube must be 3-dimensional for call to print3D");

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

// Explicit instantiations for long and for zz_p
template class HyperCube<long>;
template class ConstCubeSlice<long>;
template class CubeSlice<long>;
template void getHyperColumn(Vec<long>& v, const ConstCubeSlice<long>& s, long pos);
template void setHyperColumn(const Vec<long>& v, const CubeSlice<long>& s, long pos);
template void setHyperColumn(const Vec<long>& v, const CubeSlice<long>& s, long pos, const long& val);
template void print3D(const HyperCube<long>& c);


#include <NTL/lzz_p.h>

template class HyperCube<NTL::zz_p>;
template class ConstCubeSlice<NTL::zz_p>;
template class CubeSlice<NTL::zz_p>;
template void getHyperColumn(Vec<NTL::zz_p>& v, const ConstCubeSlice<NTL::zz_p>& s, long pos);
template void setHyperColumn(const Vec<NTL::zz_p>& v, const CubeSlice<NTL::zz_p>& s, long pos);
template void setHyperColumn(const Vec<NTL::zz_p>& v, const CubeSlice<NTL::zz_p>& s, long pos, const NTL::zz_p& val);
template void print3D(const HyperCube<NTL::zz_p>& c);
