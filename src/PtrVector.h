#ifndef _PTRVECTOR_H
#define _PTRVECTOR_H
/** PtrVector.h: convenience class templates providing a unified
 * interface for a collection of objects, returning pointers to these
 * object. ciphertexts.
 */
#include <stdexcept>
#include <climits>
#include <vector>
#include <NTL/vector.h>

template<typename T>
struct PtrVector {
  virtual T* operator[](long) const =0;
    // NOTE: the PtrVector is const, but the pointer T* is not
  virtual long size() const =0;
  virtual void resize(long newSize, const PtrVector* another=nullptr)=0;
  virtual ~PtrVector(){}

  // How many non-null entries there are (beginning at startFrom)
  virtual long numNonNull(long startFrom=0, long last=LONG_MAX) const=0;

  // Return a pointer to some non-Null T, if it can find one.
  // This is convenient since T doesn't have an empty constructor
  virtual const T* ptr2nonNull() const =0;
};

// This header provides five implementations of these interfaces,
// users can define their own as needed. The ones defined here are:

//struct PtrVector_VecT;    // constrcuted as PtrVector_VecT(NTL::Vec<T>)
//struct PtrVector_VecPt;    // constrcuted as PtrVector_VecPt(NTL::Vec<T*>)
//struct PtrVector_vectorT; // constrcuted as PtrVector_vectorT(std::vector<T>)
//struct PtrVector_vectorPt; // constrcuted as PtrVector_vectorPt(std::vector<T*>)

//struct PtrVector_slice;    // A slice, PtrVector_slice(PtrVector, start, length)


/*******************************************************************/
/* Implementation details: applications should not care about them */
/*******************************************************************/

// Convenient functions for use in the pointer-based derived classes
template<typename T>
long countNonNull(const T** ptrs, long first, long last)
{
  long count = 0;
  for (long i=first; i<last; i++) if (ptrs[i]!=nullptr) count++;
  return count;
}
template<typename T>
const T* findNonNull(const T** ptrs, long nn)
{
  for (long i=0; i<nn; i++) if (ptrs[i]!=nullptr) return ptrs[i];
  return nullptr;
}

// An implementation using Vec<T*>
template<typename T>
struct PtrVector_VecPt : PtrVector<T> {
  NTL::Vec<T*>& v;
  PtrVector_VecPt(NTL::Vec<T*>& _v) : v(_v){}
  T* operator[](long i) const override {return v[i];}
  long size() const override     { return v.length(); }

  void resize(long newSize, const PtrVector<T>* another=nullptr) override
  { v.SetLength(newSize, nullptr); }

  long numNonNull(long first=0, long last=LONG_MAX) const override
  {
    if (first<0) first = 0;
    if (last>v.length()) last = v.length();
    return countNonNull((const T**)v.elts(), first, last);
  }
  const T* ptr2nonNull() const override
  { return findNonNull((const T**)v.elts(), v.length()); }
};

// An implementation using vector<T*>
template<typename T>
struct PtrVector_vectorPt : PtrVector<T> {
  std::vector<T*>& v;
  PtrVector_vectorPt(std::vector<T*>& _v) : v(_v){}
  T* operator[](long i) const override {return v[i];}
  long size() const     override { return long(v.size()); }

  void resize(long newSize, const PtrVector<T>* another=nullptr) override
  { v.resize(newSize, nullptr); }

  long numNonNull(long first=0, long last=LONG_MAX) const override
  {
    if (first<0) first = 0;
    if (last>long(v.size())) last = long(v.size());
    return countNonNull((const T**)v.data(), first, last);
  }
  const T* ptr2nonNull() const override
  { return findNonNull((const T**)v.data(), v.size()); }
};

// An implementation using Vec<T>
template<typename T>
struct PtrVector_VecT : PtrVector<T> {
  NTL::Vec<T>& v;
  PtrVector_VecT(NTL::Vec<T>& _v) : v(_v){}
  T* operator[](long i) const override {return &(v[i]);}
  long size() const     override { return v.length(); }

  void resize(long newSize, const PtrVector<T>* another) override
  {
    if (newSize==0) setLengthZero(v);
    else {
      if (another==nullptr) another = this;
      const T* pt = another->ptr2nonNull();
      assert(pt!=nullptr);
      v.SetLength(newSize, *pt);
    }
  }

  long numNonNull(long first=0, long last=LONG_MAX) const override
  {
    if (first<0) first = 0;
    if (last>v.length()) last = v.length();
    return last - first;
  }
  const T* ptr2nonNull() const override
  { return ((v.length()>0)? &(v[0]) : nullptr); }
};

// An implementation using vector<T>
template<typename T>
struct PtrVector_vectorT : PtrVector<T> {
  std::vector<T>& v;
  PtrVector_vectorT(std::vector<T>& _v) : v(_v){}
  T* operator[](long i) const override {return &(v[i]);}
  long size() const     override { return long(v.size()); }

  void resize(long newSize, const PtrVector<T>* another) override
  {
    if (newSize==0) setLengthZero(v);
    else {
      if (another==nullptr) another = this;
      const T* pt = another->ptr2nonNull();
      assert(pt!=nullptr);
      v.resize(newSize, *pt);
    }
  }
  long numNonNull(long first=0, long last=LONG_MAX) const override
  {
    if (first<0) first = 0;
    if (last>long(v.size())) last = long(v.size());
    return last - first;
  }
  const T* ptr2nonNull() const override
  { return ((v.size()>0)? &(v[0]) : nullptr); }
};

template<typename T>
struct PtrVector_slice : PtrVector<T> {
  const PtrVector<T>& orig;
  long start, sz;
  // Special case for a slice of a slice
  PtrVector_slice(const PtrVector_slice<T>& _orig, long from, long _sz=-1)
    : orig(_orig.orig), start(_orig.start +from), sz(_sz)
  {
    if (start<0) start=0;
    else if (start>_orig.orig.size()) start=_orig.orig.size();
    if (sz<0 || sz>_orig.orig.size()-start) sz = _orig.orig.size()-start;
  }
  // The general case: slice of a PtrVector
  PtrVector_slice(const PtrVector<T>& _orig, long from, long _sz=-1)
    : orig(_orig), start(from), sz(_sz)
  {
    if (start<0) start=0;
    else if (start>_orig.size()) start=_orig.size();
    if (sz<0 || sz>_orig.size()-start) sz = _orig.size()-start;
  }
  T* operator[](long i) const override { return orig[i+start]; }
  long size() const override { return sz; }

  void resize(long newSize, const PtrVector<T>* another=nullptr) override
  { throw(std::logic_error("Cannot resize a slice")); }

  long numNonNull(long first=0, long last=LONG_MAX) const override
  { return orig.numNonNull(start+first, start+std::min(sz,last)); }
  const T* ptr2nonNull() const override { return orig.ptr2nonNull(); }
};

#endif // _PTRVECTOR_H
