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
#ifndef _NumbTh
#define _NumbTh
/**
 * @file NumbTh.h
 * @brief Miscellaneous utility functions.
 **/
#include <vector>
#include <set>
#include <cmath>
#include <cassert>
#include <string>
#include <climits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <ctime>

#include <NTL/version.h>
#include <NTL/ZZ.h>
#include <NTL/ZZX.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ_pX.h>
#include <NTL/xdouble.h>

#include <NTL/mat_GF2.h>
#include <NTL/mat_GF2E.h>
#include <NTL/GF2XFactoring.h>

#include <NTL/mat_lzz_p.h>
#include <NTL/mat_lzz_pE.h>
#include <NTL/lzz_pXFactoring.h>

#include <NTL/GF2EX.h>
#include <NTL/lzz_pEX.h>

#include <NTL/FFT.h>

// Test for the "right version" of NTL (currently 10.0.0)
#if (NTL_MAJOR_VERSION<10)
#error "This version of HElib requires NTL version 10.0.0 or above"
#endif

#include "range.h"

using namespace std;
using namespace NTL;

namespace FHEglobals
{
  //! @brief A dry-run flag
  //! The dry-run option disables most operations, to save time. This lets
  //! us quickly go over the evaluation of a circuit and estimate the
  //! resulting noise magnitude, without having to actually compute anything. 
  extern bool dryRun;

  //! @brief A list of required automorphisms
  //! When non-NULL, causes Ctxt::smartAutomorphism to just record the
  //! requested automorphism rather than actualy performing it. This can
  //! be used to get a list of needed automorphisms for certain operations
  //! and then generate all these key-switching matrices. Should only be
  //! used in conjunction with dryRun=true
  extern std::set<long>* automorphVals; 
  extern std::set<long>* automorphVals2; 
}
inline bool setDryRun(bool toWhat=true) { return (FHEglobals::dryRun=toWhat); }
inline bool isDryRun() { return FHEglobals::dryRun; }

inline void setAutomorphVals(std::set<long>* aVals)
{ FHEglobals::automorphVals=aVals; }
inline bool isSetAutomorphVals() { return FHEglobals::automorphVals!=NULL; }
inline void recordAutomorphVal(long k) { FHEglobals::automorphVals->insert(k); }

inline void setAutomorphVals2(std::set<long>* aVals)
{ FHEglobals::automorphVals2=aVals; }
inline bool isSetAutomorphVals2() { return FHEglobals::automorphVals2!=NULL; }
inline void recordAutomorphVal2(long k) { FHEglobals::automorphVals2->insert(k); }

#if (__cplusplus>199711L)
#include <memory>
#include <unordered_map>
#else
#include <tr1/memory>
#include <tr1/unordered_map>
using namespace tr1;
#endif


//! @typedef
typedef unordered_map<string, const char *> argmap_t;

typedef long LONG; // using this to identify casts that we should
                   // really get rid of at some point in the future
typedef NTL::Vec<long> zzX;

inline
bool IsZero(const zzX& a) { return a.length() == 0; }

inline 
void clear(zzX& a) { a.SetLength(0); }

//! @brief Code for parsing command line arguments.
/**
 * Tries to parse each argument as arg=val, and returns a correspinding map.
 * It returns false if errors were detected, and true otherwise. 
 **/
bool parseArgs(int argc,  char *argv[], argmap_t& argmap);


//! @brief Easier arg parsing
/**
 * Example use:
 *   ArgMapping amap;
 *
 *   long p = 2;
 *   amap.arg("p", p, "doc for p");
 *   long m = 0;
 *   amap.arg("m", m, "doc for m", "undefined"); // special default info
 *   long k = 0;
 *   amap.arg("k", k, "doc for k", NULL); // no default info
 *
 *   amap.parse(argc, argv); // parses and overrides initail values
 *                           // of p and p, returns false on error
 *
 *   amap.documentation(); // returns string with documentation
 *                         // for each parameter, one per line,
 *
 **/


/* doArgProcessing: converts c-string s to value T,
 * returns upon success.  By default, we parse using
 * the istream input operator, except when T = string
 * and just convert without any parsing.
 */

//! \cond FALSE (make doxygen ignore these classes)
template<class T>
bool doArgProcessing(T *value, const char *s)
{
  string ss(s);
  stringstream sss(ss);
  return bool(sss >> *value);
}

bool doArgProcessing(string *value, const char *s);

/* ArgProcessor: virtual base class */

class ArgProcessor {
public:
virtual bool process(const char *s) = 0;
};

/* ArgProcessorDerived: templated subclasses */

template<class T>
class ArgProcessorDerived : public ArgProcessor   {
public:
  T *value;

  virtual bool process(const char *s)
  {
    return doArgProcessing(value, s);
  }

  ArgProcessorDerived(T* _value) : value(_value) {}
};

class ArgMapping {
public:
  unordered_map< string, shared_ptr<ArgProcessor> > map;
  stringstream doc;

  // no documentation
  template<class T>
  void arg(const char *name, T& value) 
  { 
    shared_ptr<ArgProcessor> ap = 
      shared_ptr<ArgProcessor>(new ArgProcessorDerived<T>(&value));

    assert(!map[name]);
    map[name] = ap;
  }

  // documentation + standard default info
  template<class T>
  void arg(const char *name, T& value, const char *doc1) 
  {
    arg(name, value);
    doc << "\t" << name << " \t" << doc1 << "  [ default=" << value << " ]" << "\n";
  }

  // documentation + standard non-standard default info: 
  // NULL => no default info
  template<class T>
  void arg(const char *name, T& value, const char *doc1, const char *info) 
  {
    arg(name, value);
    doc << "\t" << name << " \t" << doc1; 
    if (info) 
      doc << "  [ default=" << info << " ]"  << "\n";
    else
      doc << "\n";
  }

  void note(const char *s)
  {
    doc << "\t\t   " << s << "\n";
  }

  void usage(const char *prog) 
  {
    cerr << "Usage: " << prog << " [ name=value ]...\n";
    cerr << documentation();
    exit(0);
  }

  void parse(int argc, char **argv)
  {
    for (long i = 1; i < argc; i++) {
      const char *x = argv[i];
      long j = 0;
      while (x[j] != '=' && x[j] != '\0') j++; 
      if (x[j] == '\0') usage(argv[0]);
      string name(x, j);
      const char *s = x+j+1;

      shared_ptr<ArgProcessor> ap = map[name];
      if (!ap) return usage(argv[0]);
      if (!ap->process(s)) usage(argv[0]);
    }
  }

  string documentation() 
  {
    return doc.str();
  }
};
//! \endcond


//! @brief Routines for computing mathematically correct mod and div.
//! 
//! mcDiv(a, b) = floor(a / b), mcMod(a, b) = a - b*mcDiv(a, b);
//! in particular, mcMod(a, b) is 0 or has the same sign as b

long mcMod(long a, long b);
long mcDiv(long a, long b);

//! Return multiplicative order of p modulo m, or 0 if GCD(p, m) != 1
long multOrd(long p, long m);


//! @brief Prime power solver.
//!
//! A is an n x n matrix, b is a length n (row) vector, this function finds a
//! solution for the matrix-vector equation x A = b. An error is raised if A
//! is not inverible mod p.
//!
//! NTL's current smallint modulus, zz_p::modulus(), is assumed to be p^r,
//! for p prime, r >= 1 integer.
void ppsolve(vec_zz_pE& x, const mat_zz_pE& A, const vec_zz_pE& b,
             long p, long r); 

//! @brief A version for GF2: must have p == 2 and r == 1
void ppsolve(vec_GF2E& x, const mat_GF2E& A, const vec_GF2E& b,
             long p, long r);

//! @brief Compute the inverse mod p^r of an n x n matrix.
//! 
//! NTL's current smallint modulus zz_p::modulus() is assumed to be p^r for
//! p prime, r >= 1 integer. For the zz_pE variant also zz_pE::modulus() must
//! be initialized. An error is raised if A is not inverible mod p.
void ppInvert(mat_zz_p& X, const mat_zz_p& A, long p, long r);
void ppInvert(mat_zz_pE& X, const mat_zz_pE& A, long p, long r);

// variants for GF2/GF2E to help with template code
inline void ppInvert(mat_GF2& X, const mat_GF2& A, long p, long r)
{ NTL::inv(X, A); }
inline void ppInvert(mat_GF2E& X, const mat_GF2E& A, long p, long r)
{ NTL::inv(X, A); }

void buildLinPolyMatrix(mat_zz_pE& M, long p);
void buildLinPolyMatrix(mat_GF2E& M, long p);

//! @brief Combination of buildLinPolyMatrix and ppsolve.
//!
//! Obtain the linearized polynomial coefficients from a vector L representing 
//! the action of a linear map on the standard basis for zz_pE over zz_p.
//!
//! NTL's current smallint modulus, zz_p::modulus(), is assumed to be p^r,
//! for p prime, r >= 1 integer.
void buildLinPolyCoeffs(vec_zz_pE& C, const vec_zz_pE& L, long p, long r);

//! @brief A version for GF2: must be called with p == 2 and r == 1
void buildLinPolyCoeffs(vec_GF2E& C, const vec_GF2E& L, long p, long r);

//! @brief Apply a linearized polynomial with coefficient vector C.
//!
//! NTL's current smallint modulus, zz_p::modulus(), is assumed to be p^r,
//! for p prime, r >= 1 integer.
void applyLinPoly(zz_pE& beta, const vec_zz_pE& C, const zz_pE& alpha, long p);

//! @brief A version for GF2: must be called with p == 2 and r == 1
void applyLinPoly(GF2E& beta, const vec_GF2E& C, const GF2E& alpha, long p);


//! Base-2 logarithm
inline double log2(const xdouble& x){ return log(x) * 1.442695040889; }
inline double log2(const double x){ return log(x) * 1.442695040889; }

//! @brief Factoring by trial division, only works for N<2^{60}, only the
//! primes are recorded, not their multiplicity.
void factorize(vector<long> &factors, long N);
void factorize(vector<ZZ> &factors, const ZZ& N);


//! @brief Factoring by trial division, only works for N<2^{60}
//! primes and multiplicities are recorded
void factorize(Vec< Pair<long, long> > &factors, long N);

//! @brief Prime-power factorization
void pp_factorize(vector<long>& factors, long N);

//! Compute Phi(N) and also factorize N.
void phiN(long &phiN, vector<long> &facts, long N);
void phiN(ZZ &phiN, vector<ZZ> &facts, const ZZ &N);

//! Compute Phi(N).
long phi_N(long N);

//! Returns in gens a generating set for Zm* /<p>, and in ords the
//! order of these generators. Return value is the order of p in Zm*.
long findGenerators(vector<long>& gens, vector<long>& ords, long m, long p,
                    const vector<long>& candidates=vector<long>());

//! Find e-th root of unity modulo the current modulus.
void FindPrimitiveRoot(zz_p &r, unsigned long e);
void FindPrimitiveRoot(ZZ_p &r, unsigned long e);

//! Compute mobius function (naive method as n is small).
long mobius(long n);

//! Compute cyclotomic polynomial.
ZZX Cyclotomic(long N);

//! Return a degree-d irreducible polynomial mod p
ZZX makeIrredPoly(long p, long d);

//! Find a primitive root modulo N.
long primroot(long N,long phiN);

//! Compute the highest power of p that divides N.
long ord(long N,long p);


// Returns a random mod p polynomial of degree < n
ZZX RandPoly(long n,const ZZ& p);

///@{
/**
 * @brief Reduce all the coefficients of a polynomial modulo q.
 *
 * When abs=false reduce to interval (-q/2,...,q/2), when abs=true reduce
 * to [0,q). When abs=false and q=2, maintains the same sign as the input.
 */
void PolyRed(ZZX& out, const ZZX& in,       long q, bool abs=false);
void PolyRed(ZZX& out, const ZZX& in, const ZZ& q, bool abs=false);
inline void PolyRed(ZZX& F, long q, bool abs=false) { PolyRed(F,F,q,abs); }
inline void PolyRed(ZZX& F, const ZZ& q, bool abs=false)
{ PolyRed(F,F,q,abs); }
void vecRed(Vec<ZZ>& out, const Vec<ZZ>& in, long q, bool abs);
///@}

//! Multiply the polynomial f by the integer a modulo q
void MulMod(ZZX& out, const ZZX& f, long a, long q, bool abs=true);
inline ZZX MulMod(const ZZX& f, long a, long q, bool abs=true) {
  ZZX res;
  MulMod(res, f, a, q, abs);
  return res;
}

///@{
//! @name Some enhanced conversion routines
inline void convert(long& x1, const GF2X& x2)
{
   x1 = rep(ConstTerm(x2));
}
inline void convert(long& x1, const zz_pX& x2)
{
   x1 = rep(ConstTerm(x2));
}
void convert(vec_zz_pE& X, const vector<ZZX>& A);
void convert(mat_zz_pE& X, const vector< vector<ZZX> >& A);
void convert(vector<ZZX>& X, const vec_zz_pE& A);
void convert(vector< vector<ZZX> >& X, const mat_zz_pE& A);
void convert(NTL::Vec<long>& out, const NTL::ZZX& in);
void convert(NTL::Vec<long>& out, const NTL::zz_pX& in);
void convert(NTL::Vec<long>& out, const NTL::GF2X& in);
void convert(NTL::ZZX& out, const NTL::Vec<long>& in);
void convert(NTL::GF2X& out, const NTL::Vec<long>& in);
// right now, this is just a place-holder...it may or may not 
// eventually be further fleshed out

inline void convert(zz_pX& x, const zzX& a)
{
   conv(x.rep, a);
   x.normalize();
}
///@}

//! A generic template that resolves to NTL's conv routine
template<class T1, class T2>
void convert(T1& x1, const T2& x2) 
{
   conv(x1, x2);
}

//! generic vector conversion routines
template<class T1, class T2> 
void convert(vector<T1>& v1, const vector<T2>& v2)
{
   long n = v2.size();
   v1.resize(n);
   for (long i = 0; i < n; i++)
      convert(v1[i], v2[i]);
}

template<class T1, class T2> 
void convert(vector<T1>& v1, const Vec<T2>& v2)
{
   long n = v2.length();
   v1.resize(n);
   for (long i = 0; i < n; i++)
      convert(v1[i], v2[i]);
}

template<class T1, class T2> 
void convert(Vec<T1>& v1, const vector<T2>& v2)
{
   long n = v2.size();
   v1.SetLength(n);
   for (long i = 0; i < n; i++)
      convert(v1[i], v2[i]);
}

template<class T1, class T2>
T1 convert(const T2& v2)
{  
   T1 v1;
   convert(v1, v2);
   return v1;
}


template<class T>
vector<T> vector_replicate(const T& a, long n)
{
   vector<T> res;
   res.resize(n);
   for (long i = 0; i < n; i++) res[i] = a;
   return res;
}

template<class T>
vector<T> Vec_replicate(const T& a, long n)
{
   Vec<T> res;
   res.SetLength(n);
   for (long i = 0; i < n; i++) res[i] = a;
   return res;
}

//! returns \prod_d vec[d]
long computeProd(const Vec<long>& vec);
long computeProd(const vector<long>& vec);

// some useful operations
void mul(vector<ZZX>& x, const vector<ZZX>& a, long b);
void div(vector<ZZX>& x, const vector<ZZX>& a, long b);
void add(vector<ZZX>& x, const vector<ZZX>& a, const vector<ZZX>& b);




//! @brief Finds whether x is an element of the set X of size sz,
//! Returns -1 it not and the location if true
long is_in(long x,int* X,long sz);

//! @brief Returns a CRT coefficient: x = (0 mod p, 1 mod q).
//! If symmetric is set then x \in [-pq/2, pq/2), else x \in [0,pq)
inline long CRTcoeff(long p, long q, bool symmetric=false)
{
  long pInv = InvMod(p,q); // p^-1 mod q \in [0,q)
  if (symmetric && 2*pInv >= q) return p*(pInv-q);
  else                          return p*pInv;
}

/**
 * @brief Incremental integer CRT for vectors.
 * 
 * Expects co-primes p,q with q odd, and such that all the entries in v1 are
 * in [-p/2,p/2). Returns in v1 the CRT of vp mod p and vq mod q, as integers
 * in [-pq/2, pq/2). Uses the formula:
 * \f[               CRT(vp,p,vq,q) = vp + [(vq-vp) * p^{-1}]_q * p, \f]
 * where [...]_q means reduction to the interval [-q/2,q/2). Notice that if
 * q is odd then this is the same as reducing to [-(q-1)/2,(q-1)/2], which
 * means that [...]_q * p is in [-p(q-1)/2, p(q-1)/2], and since vp is in
 * [-p/2,p/2) then the sum is indeed in [-pq/2,pq/2).
 *
 * Return true is both vectors are of the same length, false otherwise
 */
template <class zzvec>     // zzvec can be vec_ZZ, vec_long, or Vec<zz_p>
bool intVecCRT(vec_ZZ& vp, const ZZ& p, const zzvec& vq, long q);

/**
 * @brief Find the index of the (first) largest/smallest element.
 *
 * These procedures are roughly just simpler variants of std::max_element and
 * std::min_element. argmin/argmax are implemented as a template, so the code
 * must be placed in the header file for the comiler to find it. The class T
 * must have an implementation of operator> and operator< for this template to
 * work.
 * @tparam maxFlag A boolean value: true - argmax, false - argmin
 **/
template <class T, bool maxFlag>
long argminmax(vector<T>& v)
{
  if (v.size()<1) return -1; // error: this is an empty array
  unsigned long idx = 0;
  T target = v[0];
  for (unsigned long i=1; i<v.size(); i++)
    if (maxFlag) { if (v[i] > target) { target = v[i]; idx = i;} }
    else         { if (v[i] < target) { target = v[i]; idx = i;} }
  return (long) idx;
}

template <class T> long argmax(vector<T>& v)
{  return argminmax<T,true>(v); }

template <class T> long argmin(vector<T>& v)
{  return argminmax<T,false>(v); }

//! @brief A variant with a specialized comparison function
//! (*moreThan)(a,b) returns the comparison a>b
inline long argmax(vector<long>& v, bool (*moreThan)(long, long))
{
  if (v.size()<1) return -INT_MAX; // error: this is an empty array
  unsigned long idx = 0;
  long target = v[0];
  for (unsigned long i=1; i<v.size(); i++)
    if ((*moreThan)(v[i],target)) { target = v[i]; idx = i; }
  return (long) idx;
}

// Sample polynomials with entries {-1,0,1}. These functions are similar to
// the SampleSmall class from v1, but without a class around it.

// In sampleSmall, 
// sampleHWt, min(Hwt,n) random coefficients are chosen at random in {-1,+1}
// and the others are set to zero. If n=0 then n=poly.deg()+1 is used. 

//! @brief Sample polynomials with entries {-1,0,1}. Each coefficient is 0 with probability 1/2 and +-1 with probability 1/4.
void sampleSmall(ZZX &poly, long n=0);

//! @brief Sample polynomials with entries {-1,0,1} with a given HAming weight.
//!
//! Choose min(Hwt,n) coefficients at random in {-1,+1} and the others are set
//! to zero. If n=0 then n=poly.deg()+1 is used. 
 void sampleHWt(ZZX &poly, long Hwt, long n=0);

//! Sample polynomials with Gaussian coefficients.
void sampleGaussian(ZZX &poly, long n=0, double stdev=1.0);

//! Sample polynomials with coefficients sampled uniformy
//! over [-B..B]
void sampleUniform(ZZX& poly, const ZZ& B, long n=0);


/**
 * @brief Facility for "restoring" the NTL PRG state.
 *
 * NTL's random number generation faciliity is pretty limited, and does not
 * provide a way to save/restore the state of a pseudo-random stream. This
 * class gives us that ability: Constructing a RandomState object uses the PRG
 * to generate 512 bits and stores them. Upon destruction (or an explicit call
 * to restore()), these bits are used to re-set the seed of the PRG. A typical
 * usage of thie class is as follows:
 * \code
 *   {
 *     RandomState r;      // save the random state
 *
 *     SetSeed(something); // set the PRG seed to something
 *     ...                 // more code that uses the new PRG seed
 *
 *   } // The destructor is called implicitly, PRG state is restored
 * \endcode
 **/
class RandomState {
private:
  ZZ state;
  bool restored;

public:
  RandomState() {
    RandomBits(state, 512);
    restored = false;
  }

  //! Restore the PRG state of NTL
  void restore() {
    if (!restored) {
      SetSeed(state);
      restored = true;
    }
  }

  ~RandomState() {
    restore();
  }

private:
  RandomState(const RandomState&); // disable copy constructor
  RandomState& operator=(const RandomState&); // disable assignment
};

//! @brief Advance the input stream beyond white spaces and a single instance of the char cc
void seekPastChar(istream& str, int cc);

//! @brief Reverse a vector in place
template<class T> void reverse(Vec<T>& v, long lo, long hi)
{
  long n = v.length();
  assert(lo >= 0 && lo <= hi && hi < n);

  if (lo >= hi) return;

  for (long i = lo, j = hi; i < j; i++, j--) swap(v[i], v[j]); 
}

//! @brief Rotate a vector in place using swaps
// Example: rotate by 1 means [0 1 2 3] -> [3 0 1 2]
//          rotate by -1 means [0 1 2 3] -> [1 2 3 0]
template<class T> void rotate(Vec<T>& v, long k)
{
  long n = v.length();
  if (n <= 1) return;

  k %= n;
  if (k < 0) k += n;

  if (k == 0) return;

  reverse(v, 0, n-1);
  reverse(v, 0, k-1);
  reverse(v, k, n-1);
}

// An experimental facility...it is annoying that vector::size() is an
// unsigned quantity...this leads to all kinds of annoying warning messages...
//! @brief Size of STL vector as a long (rather than unsigned long)
template <typename T>
inline long lsize(const std::vector<T>& v) {
  return (long) v.size();
}

//! NTL/std compatability

// Utility functions, release memory of std::vector and NTL::Vec
template<typename T> void killVec(std::vector<T>& vec)
{ std::vector<T>().swap(vec); }
template<typename T> void killVec(NTL::Vec<T>& vec)
{ vec.kill(); }

// Set length to zero, but don't necessarily release memory
template<typename T> void setLengthZero(std::vector<T>& vec)
{ if (vec.size()>0) vec.resize(0, vec[0]); }
template<typename T> void setLengthZero(NTL::Vec<T>& vec)
{ if (vec.length()>0) vec.SetLength(0, vec[0]); }

template <typename T>
inline long lsize(const NTL::Vec<T>& v) {
  return v.length();
}


// VJS: I changed the resize functions so that the
// optional arg is handled as in C++11

template <typename T>
void resize(NTL::Vec<T>& v, long sz, const T& val) {
  return v.SetLength(sz, val);
}

template <typename T>
void resize(std::vector<T>& v, long sz, const T& val) {
  return v.resize(sz, val);
}

template <typename T>
void resize(NTL::Vec<T>& v, long sz) {
  return v.SetLength(sz);
}

template <typename T>
void resize(std::vector<T>& v, long sz) {
  return v.resize(sz);
}


//! @brief Testing if two vectors point to the same object
// Believe it or not, this is really the way to do it...
template <typename T1, typename T2>
bool sameObject(const T1* p1, const T2* p2) {
  return dynamic_cast<const void*>(p1) == dynamic_cast<const void*>(p2);
}

//! @brief Modular composition of polynomials: res = g(h) mod f
void ModComp(ZZX& res, const ZZX& g, const ZZX& h, const ZZX& f);

//! @brief Evaluates a modular integer polynomial, returns poly(x) mod p
long polyEvalMod(const ZZX& poly, long x, long p);

//! @brief Interpolate polynomial such that poly(x[i] mod p)=y[i] (mod p^e)
//! It is assumed that the points x[i] are all distinct modulo p
void interpolateMod(ZZX& poly, const vec_long& x, const vec_long& y,
		    long p, long e=1);

//! @brief returns ceiling(a/b); assumes a >=0, b>0, a+b <= MAX_LONG
inline long divc(long a, long b)
{
  return (a + b - 1)/b;
}

///@{
//! @name The size of the coefficient vector of a polynomial.
ZZ sumOfCoeffs(const ZZX& f);  // = f(1)
ZZ largestCoeff(const ZZX& f); // l_infty norm
xdouble coeffsL2Norm(const ZZX& f); // l_2 norm
///@}



//! @class zz_pXModulus1
//! @brief Auxiliary classes to facillitiate faster reduction mod Phi_m(X)
//!        when the input has degree less than m
class zz_pXModulus1 {
public:
   long m;
   zz_pX f;
   long n;

   bool specialLogic;

   long k, k1;
   fftRep R0, R1;

   zz_pXModulus fm; // just in case...

   zz_pXModulus1(long _m, const zz_pX& _f);

   const zz_pXModulus& upcast() const { return fm; } 
};


void rem(zz_pX& r, const zz_pX& a, const zz_pXModulus1& ff);

//! placeholder for pXModulus ...no optimizations
class ZZ_pXModulus1 : public ZZ_pXModulus {
public:
   ZZ_pXModulus1(long _m, const ZZ_pX& _f) : ZZ_pXModulus(_f) { }
   const ZZ_pXModulus& upcast() const { return *this; }
};

template<class T>
ostream& operator<<(ostream& s, vector<T> v)
{
  if (v.size()==0) return (s << "[]");

  s << '[';
  for (long i=0; i<(long) v.size()-1; i++)
    s << v[i] << ' ';
  return (s << v[v.size()-1] << ']');
}

template<class T>
istream& operator>>(istream& s, vector<T>& v)
{
  Vec<T> vv; // read into an NTL vector, then convert
  s >> vv;
  convert(v, vv);
  return s;
}


template<class T>
Vec<T> atoVec(const char *a) 
{
  Vec<T> v;
  string s(a);
  stringstream ss(s);
  ss >> v;
  return v;
}

template<class T>
vector<T> atovector(const char *a)
{
  Vec<T> v1 = atoVec<T>(a);
  vector<T> v2;
  convert(v2, v1);
  return v2;
}




#ifndef NTL_PROVIDES_TRUNC_FFT
// Define truncated FFT routines if not provided by NTL
inline void TofftRep_trunc(fftRep& y, const zz_pX& x, long k, long len,
                    long lo, long hi)
{  TofftRep(y, x, k, lo, hi); }

inline void TofftRep_trunc(fftRep& y, const zz_pX& x, long k, long len)
{ TofftRep_trunc(y, x, k, len, 0, deg(x)); }
#endif


//! Debug printing routines for vectors, ZZX'es, print only a few entries
template<class T> ostream& printVec(ostream& s, const Vec<T>& v,
				    long nCoeffs=40);
ostream& printZZX(ostream& s, const ZZX& poly, long nCoeffs=40);

// NOTE: Maybe NTL should contain conversion routines
// like this for the various polynomial classes?

#if 0
//! @brief stand-in for make_unique, which is C++14, not C++11
template<typename T, typename... Args>
std::unique_ptr<T> build_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}
#endif

//! This should go in NTL some day...
//! Just call as make_lazy(obj, ...) to initialize a lazy object
//! via a call to a constructor T(...)
template<class T, class P, class... Args>
void make_lazy(const Lazy<T,P>& obj, Args&&... args)
{
   typename Lazy<T,P>::Builder builder(obj);
   if (!builder()) return;
   UniquePtr<T,P> ptr;
   ptr.make(std::forward<Args>(args)...);
   builder.move(ptr);
}


#endif
