#ifndef _PTRMATRIX_H
#define _PTRMATRIX_H
/** PtrMatrix.h: convenience class templates providing a unified
 * interface for a matrix of objects, returning pointers to these objects.
 */
#include "PtrVector.h"

template<typename T>
struct PtrMatrix {
  virtual PtrVector<T>& operator[](long) =0;             // returns a row
  virtual const PtrVector<T>& operator[](long) const =0; // returns a row
  virtual long size() const =0;        // How many rows
  virtual void resize(long newSize)=0; // reset the number of rows
  virtual ~PtrMatrix(){}

  // Return a pointer to some non-Null T, if it can find one.
  // This is convenient since T may not have an empty constructor
  virtual const T* ptr2nonNull() const
  {
    for (long i=0; i<size(); i++) {
      const PtrVector<T>& row = (*this)[i];
      for (long j=0; j<row.size(); j++) {
        T* pt = row[j];
        if (pt!=nullptr) return pt;
      }
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


/*******************************************************************/
/* Implementation details: applications should not care about them */
/*******************************************************************/

// An implementation using Vec< Vec<T> >
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

// An implementation using Vec< Vec<T>* >
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
  void resize(long newSize) override         // reset the number of rows
  { throw(std::logic_error("Cannot resize PtrMatrix_ptVec")); }
};

// An implementation using vector<vector<T>>
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

// An implementation using vector<vector<T>*>
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
  void resize(long newSize) override         // reset the number of rows
  { throw(std::logic_error("Cannot resize PtrMatrix_ptvector")); }
};
#endif // _PTRMATRIX_H
