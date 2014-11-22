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
#include "hypercube.h"
#include <iomanip>


// Rotate k positions along the d'th dimension: The content of slot j that has
// coordinate j_d in the d'th dimension is moved to j' that has coordinates the
// same as j in all the dimensions except d, and has coordinate at dimension d
// j'_d = j_d +k mod sz (where sz is the size of the d'th dimension).
template<class T>
void HyperCube<T>::rotate1D(long d, long k)
{
  assert(d >=0 && d < getNumDims());

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
  assert(d >=0 && d < getNumDims());

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
   assert(dOffset >= 0 && dOffset <= bigger.getNumDims()); 
   // allow zero-dimensional slice

   assert(i >= 0 && i < bigger.getProd(0, dOffset));

   data = bigger.data;
   sig = bigger.sig;
   dimOffset = dOffset + bigger.dimOffset;
   sizeOffset = bigger.sizeOffset + i*bigger.getProd(dOffset);
}

template<class T>
ConstCubeSlice<T>::ConstCubeSlice(const HyperCube<T>& _cube, long i,
				  long dOffset)
{
   assert(dOffset >= 0 && dOffset <= _cube.getNumDims()); // allow zero-dimensional slice
   assert(i >= 0 && i < _cube.getProd(0,dOffset));

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
   assert(n == other.getSize());

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

   assert(pos >= 0 && pos < m);
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

   assert(pos >= 0 && pos < m);
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

   assert(pos >= 0 && pos < m);
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
