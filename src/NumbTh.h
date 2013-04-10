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

#include <vector>
#include <cmath>
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



#include <tr1/unordered_map>
#include <string>

NTL_CLIENT



// Code for parsing command line arguments

typedef tr1::unordered_map< string, const char * > argmap_t;


/* parseArgs tries to parse each argument as arg=val,
 * and returns a correspinding map.  It returns false
 * if errors were detected, and true otherwise. 
 */

bool parseArgs(int argc,  char *argv[], argmap_t& argmap);


// return multiplicative order of p modulo m, or 0
// if GCD(p, m) != 1
long multOrd(long p, long m);


// Builds the matrix defining the linearized polynomial transformation
// zz_p::modulus() is assumed to be p^r, for p prime, r >= 1
// Call ppsolve(C, L, M, p, r) to get the coeffecients C for the
// linearized polynomial represented the linear map defined
// by its action on the standard basis for zz_pE over zz_p:
// for i = 0..zz_pE::degree() - 1: x^i -> L[i], where x = (X mod zz_pE::modulus())
void buildLinPolyMatrix(mat_zz_pE& M, long p);

// version for GF2: must be called with p == 2
void buildLinPolyMatrix(mat_GF2E& M, long p);


// prime power solver
// zz_p::modulus() is assumed to be p^r, for p prime, r >= 1
// A is an n x n matrix, b is a length n (row) vector,
// and a solution for the matrix-vector equation x A = b is found.
// If A is not inverible mod p, then error is raised.
void ppsolve(vec_zz_pE& x, const mat_zz_pE& A, const vec_zz_pE& b,
             long p, long r); 

// version for GF2: must have p == 2 and r == 1
void ppsolve(vec_GF2E& x, const mat_GF2E& A, const vec_GF2E& b,
             long p, long r);

// Combines the above two routines to obtain the linearized
// polynomial coefficients from a vector L representing 
// the action of a linear map on the standard basis for zz_pE over zz_p.
// zz_p::modulus() is assumed to be p^r, for p prime, r >= 1
void buildLinPolyCoeffs(vec_zz_pE& C, const vec_zz_pE& L, long p, long r);

// version for GF2: must be called with p == 2 and r == 1
void buildLinPolyCoeffs(vec_GF2E& C, const vec_GF2E& L, long p, long r);

// Apply a linearized polynomial with coefficient vector C.
// zz_p::modulus() is assumed to be p^r, for p prime, r >= 1.
void applyLinPoly(zz_pE& beta, const vec_zz_pE& C, const zz_pE& alpha, long p);

// version for GF2: must be called with p == 2
void applyLinPoly(GF2E& beta, const vec_GF2E& C, const GF2E& alpha, long p);






inline double log2(const xdouble& x){ return log(x) * 1.442695040889; }
inline double log2(const double x){ return log(x) * 1.442695040889; }

// Factoring by trial division, only works for N<2^{60}, only the primes
// are recorded, not their multiplicity. class zz can be long or ZZ
void factorize(vector<long> &factors, long N);
void factorize(vector<ZZ> &factors, const ZZ& N);

/* Compute Phi(N) and also factorize N */
void phiN(long &phiN, vector<long> &facts, long N);
void phiN(ZZ &phiN, vector<ZZ> &facts, const ZZ &N);

/* Compute Phi(N) */
int phi_N(int N);

// Find e-th root of unity modulo the current modulus
void FindPrimitiveRoot(zz_p &r, unsigned e);
void FindPrimitiveRoot(ZZ_p &r, unsigned e);

/* Compute mobius function (naive method as n is small) */
int mobius(int n);

/* Compute cyclotomic polynomial */
ZZX Cyclotomic(int N);

/* Find a primitive root modulo N */
int primroot(int N,int phiN);

int ord(int N,int p);


// Rand mod p poly of degree < n
ZZX RandPoly(int n,const ZZ& p);


/* When abs=false reduce to interval (-q/2,...,q/2), when abs=true reduce
 * to [0,q). When abs=false and q=2, maintains the same sign as the input.
 */
void PolyRed(ZZX& out, const ZZX& in,       int q, bool abs=false);
void PolyRed(ZZX& out, const ZZX& in, const ZZ& q, bool abs=false);
inline void PolyRed(ZZX& F, int q, bool abs=false) { PolyRed(F,F,q,abs); }
inline void PolyRed(ZZX& F, const ZZ& q, bool abs=false)
{ PolyRed(F,F,q,abs); }

// multiply the polynomial f by the integer a modulo q
void MulMod(ZZX& out, const ZZX& f, long a, long q, bool abs=true);
inline ZZX MulMod(const ZZX& f, long a, long q, bool abs=true) {
  ZZX res;
  MulMod(res, f, a, q, abs);
  return res;
}


// some enhanced conversion routines

// specialized instantiations
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


// first a generic template that resolves to NTL's conv routine
template<class T1, class T2>
void convert(T1& x1, const T2& x2) 
{
   conv(x1, x2);
}

//  a generic vector conversion routine
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









/* Finds whether x is an element of the set X
   of size sz, returns -1 it not and the location if true
*/
int is_in(int x,int* X,int sz);

/* Incremental integer CRT for vectors. Expects co-primes p,q with q odd,
 * and such that all the entries in v1 are in [-p/2,p/2). Returns in v1 the
 * CRT of vp mod p and vq mod q, as integers in [-pq/2, pq/2). Uses the
 * formula:
 *                   CRT(vp,p,vq,q) = vp + [(vq-vp)*p^{-1}]_q * p,
 *
 * where [...]_q means reduction to the interval [-q/2,q/2). Notice that if q
 * is odd then this is the same as reducing to [-(q-1)/2,(q-1)/2], which means
 * that [...]_q * p is in [-p(q-1)/2, p(q-1)/2], and since vp is in [-p/2,p/2)
 * then the sum is indeed in [-pq/2,pq/2).
 *
 * return true is both vectors are of the same length, false otherwise
 */
template <class zzvec>        // zzvec can be vec_ZZ or vec_long
bool intVecCRT(vec_ZZ& vp, const ZZ& p, const zzvec& vq, long q);

// argmax(v)/argmin(v) finds the index of the (first) largest/smallest
//   element in the vector v. They are roughly just simpler variants of
//   std::max_element and std::mim_element.
// argmin/argmax are implemented as a template, so the code must be placed
// in the header file for the comiler to find it. The class T must have an
// implementation of operator> and operator< for this template to work.

template <class T, bool maxFlag>
long argminmax(vector<T>& v)
{
  if (v.size()<1) return -1; // error: this is an empty array
  unsigned idx = 0;
  T target = v[0];
  for (unsigned i=1; i<v.size(); i++)
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

// In sampleSmall, each coeff is 0 w/ prob. 1/2 and +-1 w/ prob. 1/4. In
// sampleHWt, min(Hwt,n) random coefficients are chosen at random in {-1,+1}
// and the others are set to zero. If n=0 then n=poly.deg()+1 is used. 

void sampleSmall(ZZX &poly, long n=0);
void sampleHWt(ZZX &poly, long Hwt, long n=0);

// Sample polynomials with Gaussian coefficients. This function is similar
// to the MultivariateGauss class from v1, but without a class around it.

void sampleGaussian(ZZX &poly, long n=0, double stdev=1.0);


// NTL's random number generation faciliity is pretty limited,
// and does not provide a way to save/restore the state of a
// pseudo-random stream.  The following gives us that ability.

class RandomState {
private:
  ZZ state;
  bool restored;

public:
  RandomState() {
    RandomBits(state, 512);
    restored = false;
  }

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

// advance the input stream beyond white spaces and a single instance of cc
void seekPastChar(istream& str, int cc);

// usage: 
//  { RandomState state; ... }
//
// The destructor will restore state at end of scope,
// but this can also be done explicitly via the restore method.


// An experimental facility...it is annoying that vector::size()
// is an unsigned quantity...this leads to all kinds of annoying
// warning messages...

template <typename T>
inline long lsize(const vector<T>& v) {
  return (long) v.size();
}

// Believe it or not, this is the way to test if two pointers 
// really point to the same object

template <typename T1, typename T2>
bool sameObject(const T1* p1, const T2* p2) {
  return dynamic_cast<const void*>(p1) == dynamic_cast<const void*>(p2);
}


// modular composition: res = g(h) mod f
void ModComp(ZZX& res, const ZZX& g, const ZZX& h, const ZZX& f);

// procedures to compute size of the coefficient vector of a polynomial
ZZ sumOfCoeffs(const ZZX& f);  // = f(1)
ZZ largestCoeff(const ZZX& f); // l_infty norm
xdouble coeffsL2Norm(const ZZX& f); // l_2 norm

#endif
