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
#ifndef _NumbTh
#define _NumbTh
/**
 * @file NumbTh.h
 * @brief Miscellaneous utility functions.
 **/
#include <vector>
#include <cmath>
#include <cassert>
#include <istream>
#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZX.h>
#include <NTL/GF2X.h>
#include <NTL/vec_ZZ.h>
#include <NTL/xdouble.h>
#include <NTL/mat_lzz_pE.h>
#include <NTL/mat_GF2E.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/GF2XFactoring.h>
#include <string>

NTL_CLIENT

#if (__cplusplus>199711L)
#include <unordered_map>
//! @typedef
typedef std::unordered_map<string, const char *> argmap_t;
#else
#include <tr1/unordered_map>
typedef tr1::unordered_map<string, const char *> argmap_t;
#endif


//! @brief Code for parsing command line arguments.
/**
 * Tries to parse each argument as arg=val, and returns a correspinding map.
 * It returns false if errors were detected, and true otherwise. 
 **/
bool parseArgs(int argc,  char *argv[], argmap_t& argmap);


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

//! Compute Phi(N) and also factorize N.
void phiN(long &phiN, vector<long> &facts, long N);
void phiN(ZZ &phiN, vector<ZZ> &facts, const ZZ &N);

//! Compute Phi(N).
long phi_N(long N);

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
///@}

//! A generic template that resolves to NTL's conv routine
template<class T1, class T2>
void convert(T1& x1, const T2& x2) 
{
   conv(x1, x2);
}

//! A generic vector conversion routine
template<class T1, class T2> 
void convert(vector<T1>& v1, const vector<T2>& v2)
{
   long n = v2.size();
   v1.resize(n);
   for (long i = 0; i < n; i++)
      convert(v1[i], v2[i]);
}

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
 * \f[                  CRT(vp,p,vq,q) = vp + [(vq-vp) * p^{-1}]_q * p, \f]
 * where [...]_q means reduction to the interval [-q/2,q/2). Notice that if
 * q is odd then this is the same as reducing to [-(q-1)/2,(q-1)/2], which
 * means that [...]_q * p is in [-p(q-1)/2, p(q-1)/2], and since vp is in
 * [-p/2,p/2) then the sum is indeed in [-pq/2,pq/2).
 *
 * Return true is both vectors are of the same length, false otherwise
 */
template <class zzvec>        // zzvec can be vec_ZZ or vec_long
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
inline long lsize(const vector<T>& v) {
  return (long) v.size();
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



// Auxilliary classes to facillitiate faster reduction mod Phi_m(X)
// when the input has degree less than m



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

// placeholder for ZZ_pX's ...no optimizations

class ZZ_pXModulus1 : public ZZ_pXModulus {
public:
   ZZ_pXModulus1(long _m, const ZZ_pX& _f) : ZZ_pXModulus(_f) { }
   const ZZ_pXModulus& upcast() const { return *this; }
};



#endif
