#ifndef _PTRVECTOR_H
#define _PTRVECTOR_H
/** PtrVector.h: convenience class templates providing a unified
 * interface for a collection of objects, returning pointers to these objects.
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
  virtual void resize(long newSize, const PtrVector* another=nullptr)
  { throw(std::logic_error("Cannot resize a generic PtrVector")); }
  virtual ~PtrVector(){}

  bool isSet(long i) const
  {
    if (i<0 || i>=size()) return false;
    return ((*this)[i] != nullptr);
  }

  // How many non-null entries there are (beginning at startFrom)
  virtual long numNonNull(long first=0, long last=LONG_MAX) const
  {
    if (first<0) first = 0;
    if (last>size()) last = size();
    long count = 0;
    for (long i=first; i<last; i++) if ((*this)[i]!=nullptr) count++;
    return count;
  }

  // Return a pointer to some non-Null T, if it can find one.
  // This is convenient since T may not have an empty constructor
  virtual const T* ptr2nonNull() const
  {
    for (long i=0; i<size(); i++) {
      T* pt = (*this)[i];
      if (pt!=nullptr) return pt;
    }
    return nullptr;
  }
};

// This header provides five implementations of these interfaces, but
// users can define their own as needed. The ones defined here are:

//struct PtrVector_VecT;    // constrcuted as PtrVector_VecT(NTL::Vec<T>)
//struct PtrVector_VecPt;   // constrcuted as PtrVector_VecPt(NTL::Vec<T*>)
//struct PtrVector_vectorT;// constrcuted as PtrVector_vectorT(std::vector<T>)
//struct PtrVector_vectorPt;// constrcuted PtrVector_vectorPt(std::vector<T*>)

//struct PtrVector_slice;// A slice, PtrVector_slice(PtrVector, start, length)


template<typename T> long lsize(const PtrVector<T>& v) {return v.size();}
template<typename T> void setLengthZero(PtrVector<T>& v){v.resize(0);}
template<typename T> void resize(PtrVector<T>& v, long newSize, const T& val);
  // implementation of resize function below
template<typename T> void resize(PtrVector<T>& v, long newSize, const T* val)
{ resize<T>(v, newSize, *val); }


// Templates for element-wies vector-copy

// Generic version for std::vector, NTL::Vec
template<typename V1, typename V2>
void vecCopy(V1& v1, const V2& v2)
{
  int n = lsize(v2);
  if (n==0)
    setLengthZero(v1);
  else {
    resize(v1, n, v2[0]);
    for (int i=0; i<n; i++) v1[i] = v2[i];
  }
}

// Specializations for PtrVector
template<typename V, typename T> // V is either Vec<T> or vector<T>
void vecCopy(V& v1, const PtrVector<T>& v2)
{
  int n = lsize(v2);
  if (n==0)
    setLengthZero(v1);
  else {
    resize(v1, n, *(v2[0]));
    for (int i=0; i<n; i++) v1[i] = *(v2[i]);
  }
}
template<typename V, typename T> // V is either Vec<T> or vector<T>
void vecCopy(PtrVector<T>& v1, const V& v2)
{
  int n = lsize(v2);
  if (n==0)
    setLengthZero(v1);
  else {
    resize(v1, n, v2[0]);
    for (int i=0; i<n; i++) *(v1[i]) = v2[i];
  }
}
template<typename T> // V is either Vec<T> or vector<T>
void vecCopy(PtrVector<T>& v1, const PtrVector<T>& v2)
{
  int n = lsize(v2);
  if (n==0)
    setLengthZero(v1);
  else {
    resize(v1, n, *(v2[0]));
    for (int i=0; i<n; i++) *(v1[i]) = *(v2[i]);
  }
}



/*******************************************************************/
/* Implementation details: applications should not care about them */
/*******************************************************************/

// An implementation using Vec<T*>
template<typename T>
struct PtrVector_VecPt : PtrVector<T> {
  NTL::Vec<T*>& v;
  PtrVector_VecPt(NTL::Vec<T*>& _v) : v(_v){}
  T* operator[](long i) const override {return v[i];}
  long size() const override     { return v.length(); }

  void resize(long newSize, const PtrVector<T>* another=nullptr) override
  { v.SetLength(newSize, nullptr); }
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
      // Try to find a non-null pointer to T that you can give to resize
      if (another==nullptr) another = this;
      const T* pt = another->ptr2nonNull();
      assert(pt!=nullptr);
      v.SetLength(newSize, *pt); // Do the actual resize
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
      // Try to find a non-null pointer to T that you can give to resize
      if (another==nullptr) another = this;
      const T* pt = another->ptr2nonNull();
      assert(pt!=nullptr);
      v.resize(newSize, *pt); // do the actual resize
    }
  }
  long numNonNull(long first=0, long last=LONG_MAX) const override
  {
    if (first<0) first = 0;
    if (last>long(v.size())) last = long(v.size());
    return last - first;
  }
};

// Implementing a slice of a vector
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

  long numNonNull(long first=0, long last=LONG_MAX) const override
  { return orig.numNonNull(start+first, start+std::min(sz,last)); }
  const T* ptr2nonNull() const override { return orig.ptr2nonNull(); }
};

// An implementation using a single T*
template<typename T>
struct PtrVector_Singleton : PtrVector<T> {
  const T* v;
  PtrVector_Singleton(const T* _v) : v(_v){}
  T* operator[](long i) const override {return (i==0)? ((T*)v) : nullptr;}
  long size() const override { return 1L; }
};

template<typename T> void resize(PtrVector<T>& v, long newSize, const T& val)
{
  PtrVector_Singleton<T> t(&val);
  v.resize(newSize, &t);
}

#endif // _PTRVECTOR_H
