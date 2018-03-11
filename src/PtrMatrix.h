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
#ifndef _PTRMATRIX_H
#define _PTRMATRIX_H
/**
 * @file PtrMatrix.h
 * @brief Convenience class templates providing a unified interface
 *    for a matrix of objects, returning pointers to these objects.
 **/
#include "PtrVector.h"

//! @brief An abstract class for an array of PtrVectors
template<typename T>
struct PtrMatrix {
  virtual PtrVector<T>& operator[](long) =0;             // returns a row
  virtual const PtrVector<T>& operator[](long) const =0; // returns a row
  virtual long size() const =0;        // How many rows
  virtual void resize(long newSize)    // reset the number of rows
  { throw(std::logic_error("Cannot resize generic PtrMatrix")); }
  virtual ~PtrMatrix(){}

  // Return a pointer to some non-Null T, if it can find one.
  // This is convenient since T may not have an empty constructor
  virtual const T* ptr2nonNull() const
  {
    for (long i=0; i<size(); i++) {
      const T* pt = (*this)[i].ptr2nonNull();
      if (pt!=nullptr) return pt;
    }
    return nullptr;
  }
};

template<typename T> long lsize(const PtrMatrix<T>& v) {return v.size();}
template<typename T> void resize(PtrMatrix<T>& v, long newSize)
{ v.resize(newSize); }
template<typename T> void setLengthZero(PtrMatrix<T>& v){v.resize(0);}
// implementation of resize function below

// This header provides some implementations of these interfaces, but
// users can define their own as needed. The ones defined here are:

//struct PtrMatrix_Vec;    // NTL::Vec<NTL::Vec<T>>
//struct PtrMatrix_vector; // std::vector<std::vector<T>>
//struct PtrMatrix_ptVec;    // NTL::Vec<NTL::Vec<T>*>
//struct PtrMatrix_ptvector; // std::vector<std::vector<T>*>

#include <initializer_list>
template<typename T>
const T* ptr2nonNull(std::initializer_list<const PtrVector<T>*> list)
{
  for (auto elem : list) {
    const T* ptr = elem->ptr2nonNull();
    if (ptr != nullptr) return ptr;
  }
  return nullptr;
}


/*******************************************************************/
/* Implementation details: applications should not care about them */
/*******************************************************************/

//! @brief An implementation of PtrMatrix using Vec< Vec<T> >
template<typename T>
struct PtrMatrix_Vec : PtrMatrix<T> {
  NTL::Vec< NTL::Vec<T> >& buffer;
  std::vector< PtrVector_VecT<T> > rows;
  // rows[i] is a PtrVector_VecT<T> object 'pointing' to buffer[i]
  // the above uses std::vector to be able to use emplace

  PtrMatrix_Vec(NTL::Vec< NTL::Vec<T> >& mat): buffer(mat)
  {
    rows.reserve(lsize(mat));         // allocate memory
    for (int i=0; i<lsize(mat); i++)  // initialize
      rows.emplace_back(buffer[i]);
  }
  PtrVector<T>& operator[](long i) override             // returns a row
  { return rows[i]; }
  const PtrVector<T>& operator[](long i) const override // returns a row
  { return rows[i]; }
  long size() const override { return lsize(rows); }    // How many rows
  void resize(long newSize) override         // reset the number of rows
  {
    long oldSize = size();
    if (oldSize == newSize) return; // nothing to do

    buffer.SetLength(newSize); // resize buffer, then add/delete 'pointers'
    if (newSize > oldSize) {
      rows.reserve(newSize);
      for (int i=oldSize; i<newSize; i++)
        rows.emplace_back(buffer[i]);
    }
    //    else rows.resize(newSize);
    // Can't shrink without operator=
    else { std::cerr << "Attempt to shrink PtrMatrix_Vec failed\n"; }
  }
};

//! @brief An implementation of PtrMatrix using Vec< Vec<T>* >
template<typename T>
struct PtrMatrix_ptVec : PtrMatrix<T> {
  NTL::Vec< NTL::Vec<T>* >& buffer;
  std::vector< PtrVector_VecT<T> > rows;
  // rows[i] is a PtrVector_VecT<T> object 'pointing' to *buffer[i]
  // the above uses std::vector to be able to use emplace

  PtrMatrix_ptVec(NTL::Vec< NTL::Vec<T>* >& mat): buffer(mat)
  {
    rows.reserve(lsize(mat));         // allocate memory
    for (int i=0; i<lsize(mat); i++)  // initialize
      rows.emplace_back(*(buffer[i]));
  }
  PtrVector<T>& operator[](long i) override             // returns a row
  { return rows[i]; }
  const PtrVector<T>& operator[](long i) const override // returns a row
  { return rows[i]; }
  long size() const override { return lsize(rows); }    // How many rows
};

//! @brief An implementation of PtrMatrix using vector< vector<T> >
template<typename T>
struct PtrMatrix_vector : PtrMatrix<T> {
  std::vector< std::vector<T> >& buffer;
  std::vector< PtrVector_vectorT<T> > rows;
  // rows[i] is a PtrVector_vectorT<T> object 'pointing' to buffer[i]

  PtrMatrix_vector(std::vector< std::vector<T> >& mat): buffer(mat)
  {
    rows.reserve(lsize(mat));         // allocate memory
    for (int i=0; i<lsize(mat); i++)  // initialize
      rows.emplace_back(buffer[i]);
  }
  PtrVector<T>& operator[](long i) override             // returns a row
  { return rows[i]; }
  const PtrVector<T>& operator[](long i) const override // returns a row
  { return rows[i]; }
  long size() const override { return lsize(rows); }    // How many rows
  void resize(long newSize) override         // reset the number of rows
  {
    long oldSize = size();
    if (oldSize == newSize) return; // nothing to do

    buffer.resize(newSize); // resize buffer, then add/delete 'pointers'
    if (newSize > oldSize) {
      rows.reserve(newSize);
      for (int i=oldSize; i<newSize; i++)
        rows.emplace_back(buffer[i]);
    }
    //    else rows.resize(newSize);
    // Can't shrink without operator=
    else { std::cerr << "Attempt to shrink PtrMatrix_vector failed\n"; }
  }
};

//! @brief An implementation of PtrMatrix using vector< vector<T>* >
template<typename T>
struct PtrMatrix_ptvector : PtrMatrix<T> {
  std::vector< std::vector<T>* >& buffer;
  std::vector< PtrVector_vectorT<T> > rows;
  // rows[i] is a PtrVector_vectorT<T> object 'pointing' to *buffer[i]

  PtrMatrix_ptvector(std::vector< std::vector<T>* >& mat): buffer(mat)
  {
    rows.reserve(lsize(mat));         // allocate memory
    for (int i=0; i<lsize(mat); i++)  // initialize
      rows.emplace_back(*(buffer[i]));
  }
  PtrVector<T>& operator[](long i) override             // returns a row
  { return rows[i]; }
  const PtrVector<T>& operator[](long i) const override // returns a row
  { return rows[i]; }
  long size() const override { return lsize(rows); }    // How many rows
};
#endif // _PTRMATRIX_H
