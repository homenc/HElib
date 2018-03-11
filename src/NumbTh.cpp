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
#include "NumbTh.h"

#include "timing.h"

#include <fstream>
#include <cctype>
#include <algorithm>   // defines count(...), min(...)

bool FHEglobals::dryRun = false;

// Code for parsing command line

bool parseArgs(int argc,  char *argv[], argmap_t& argmap)
{
  for (long i = 1; i < argc; i++) {
    char *x = argv[i];
    long j = 0;
    while (x[j] != '=' && x[j] != '\0') j++; 
    if (x[j] == '\0') return false;
    string arg(x, j);
    if (argmap[arg] == NULL) return false;
    argmap[arg] = x+j+1;
  }

  return true;
}

bool doArgProcessing(string *value, const char *s)
{
  *value = string(s);
  return true;
}


// Mathematically correct mod and div, avoids overflow
long mcMod(long a, long b) 
{
   long r = a % b;

   if (r != 0 && (b < 0) != (r < 0))
      return r + b;
   else
      return r;

}

long mcDiv(long a, long b) {

   long r = a % b;
   long q = a / b;

   if (r != 0 && (b < 0) != (r < 0))
      return q + 1;
   else
      return q;
}


// return multiplicative order of p modulo m, or 0 if GCD(p, m) != 1
long multOrd(long p, long m)
{
  if (GCD(p, m) != 1) return 0;

  p = p % m;
  long ord = 1;
  long val = p; 
  while (val != 1) {
    ord++;
    val = MulMod(val, p, m);
  }
  return ord;
}


// returns \prod_d vec[d]
template<class T> static inline long computeProd(const T& vec, long k)
{
  long prod = 1;
  for (long d = 0; d < k; d++)
    prod = prod * vec[d];
  return prod;
}
long computeProd(const Vec<long>& vec) { return computeProd(vec, vec.length());}
long computeProd(const vector<long>& vec) {return computeProd(vec, vec.size());}


// return a degree-d irreducible polynomial mod p
ZZX makeIrredPoly(long p, long d)
{
	assert(d >= 1);
  assert(ProbPrime(p));

  if (d == 1) return ZZX(1, 1); // the monomial X

  zz_pBak bak; bak.save();
  zz_p::init(p);
  return to_ZZX(BuildIrred_zz_pX(d));
}


// Factoring by trial division, only works for N<2^{60}.
// Only the primes are recorded, not their multiplicity
template<class zz> static void factorT(vector<zz> &factors, const zz &N)
{
  factors.resize(0); // reset the factors

  if (N<2) return;   // sanity check

  PrimeSeq s;
  zz n = N;
  while (true) {
    if (ProbPrime(n)) { // we are left with just a single prime
      factors.push_back(n);
      return;
    }
    // if n is a composite, check if the next prime divides it
    long p = s.next();
    if ((n%p)==0) {
      zz pp;
      conv(pp,p);
      factors.push_back(pp);
      do { n /= p; } while ((n%p)==0);
    }
    if (n==1) return;
  }
}
void factorize(vector<long> &factors, long N) { factorT<long>(factors, N);}
void factorize(vector<ZZ> &factors, const ZZ& N) {factorT<ZZ>(factors, N);}

// Returns a list of prime factors and their multiplicity, 
// N = \prod_i factors[i].first^{factors[i].second}
void factorize(Vec< Pair<long, long> > &factors, long N)
{
  factors.SetLength(0);

  if (N < 2) return;

  PrimeSeq s;
  long n = N;
  while (n > 1) {
    if (ProbPrime(n)) { // n itself is a prime, add (n,1) to the list
      append(factors, cons(n, 1L));
      return;
    }

    long p = s.next();
    if ((n % p) == 0) { // p divides n, find its multiplicity
      long e = 1;
      n = n/p;
      while ((n % p) == 0) {
        n = n/p;
        e++;
      }
      append(factors, cons(p, e)); // add (p,e) to the list
    }
  }
}

// Prime-power factorization
void pp_factorize(vector<long>& factors, long N)
{
  Vec< Pair<long, long> > pf;
  factorize(pf,N); // prime factors, N = \prod_i pf[i].first^{pf[i].second}
  factors.resize(pf.length());
  for (long i=0; i<pf.length(); i++)
    factors[i] = power_long(pf[i].a, pf[i].b); // p_i^e_i
}


template<class zz> static void phiNT(zz &phin, vector<zz> &facts, const zz &N)
{
  if (facts.size()==0) factorize(facts,N);

  zz n = N;
  conv(phin,1); // initialize phiN=1
  for (unsigned long i=0; i<facts.size(); i++) {
    zz p = facts[i];
    phin *= (p-1); // first factor of p
    for (n /= p; (n%p)==0; n /= p) phin *= p; // multiple factors of p
  } 
}
// Specific template instantiations for long and ZZ
void phiN(long &pN, vector<long> &fs, long N)  { phiNT<long>(pN,fs,N); }
void phiN(ZZ &pN, vector<ZZ> &fs, const ZZ &N) { phiNT<ZZ>(pN,fs,N);   }

/* Compute Phi(N) */
long phi_N(long N)
{
  long phiN=1,p,e;
  PrimeSeq s;
  while (N!=1)
    { p=s.next();
      e=0;
      while ((N%p)==0) { N=N/p; e++; }
      if (e!=0)
        { phiN=phiN*(p-1)*power_long(p,e-1); }
    }
  return phiN;
}

/* While generating the representation of (Z/mZ)^*, we keep the elements in
 * equivalence classes, and each class has a representative element (called
 * a pivot), which is the smallest element in the class. Initialy each element
 * is in its own class. When we add a new generator g we unify classes if
 * their members are a factor of g from each other, repeating this process
 * until no further unification is possible.
 *
 * We begin by adding p as a generator, thus computing the equivalence
 * classes of (Z/mZ)^* /<p>. Then we repeatedly compute the orders of
 * all elements in the current quotient group, choose the highest-order
 * element and add it as a generator, then recompute the new quotient
 * group and so on, until the remaining quotient group is the trivial
 * one, containing just a a single element. A twist is that when choosing
 * the highest-order generator, we try to find one whose order in the
 * current quotient group is the same as in the original group (Z/mZ)^*.
 **/

// The function conjClasses(classes,g,m) unifies equivalence classes that
// have elements which are a factor of g apart, the pivot of the unified
// class is the smallest element in that class. 
static void conjClasses(vector<long>& classes, long g, long m)
{
  for (long i=0; i<m; i++) {
    if (classes[i]==0) continue; // i \notin (Z/mZ)^*

    if (classes[i]<i) { // i is not a pivot, updated its pivot
      classes[i] = classes[classes[i]];
      continue;
    }

    // If i is a pivot, update other pivots to point to it
    unsigned long j = MulMod(i, g, m);
    while (classes[j] != i) {
      classes[classes[j]]= i; // Merge the equivalence classes of j and i

      // Note: if classes[j]!=j then classes[j] will be updated later,
      //       when we get to i=j and use the code for "i not pivot".

      j = MulMod(j, g, m);
    }
  }
}


// The function compOrder(orders, classes,flag,m) computes the order of elements
// of the quotient group, relative to current equivalent classes. If flag==1
// then also check if the order is the same as in (Z/mZ)^* and store the order
// with negative sign if not.
static void 
compOrder(vector<long>& orders, vector<long>& classes, bool flag, long m)
{
  orders[0] = 0;
  orders[1] = 1;
  for (long i=2; i<m; i++) {
    if (classes[i] <= 1) { // ignore i not in Z_m^* and order-0 elements
      orders[i] = (classes[i]==1)? 1 : 0;
      continue;
    }

    // If not comparing order with (Z/mZ)^*, only compute the order of pivots

    if (!flag && classes[i]<i){          // not a pivot
      orders[i] = orders[classes[i]];
      continue;
    }

    // For an element i>1, the order is at least 2
    long j = MulMod(i, i, m);
    long ord = 2;
    while (classes[j] != 1) {
      j = MulMod(j, i, m); // next element in <i>
      ord++;    // count how many steps until we reach 1
    }

    // When we get here we have classes[j]==1, so if j!=1 it means that the
    // order of i in the quotient group is smaller than its order in the
    // entire group Z_m^*. If the flag is set then we store orders[i] = -ord.
    
    if (flag && j != 1) ord = -ord; // order in Z_m^* is larger than ord
    orders[i] = ord;
  }
}

// Compare numbers based on their absolute value
#if 1

// This version prefers positive numbers over negative
static bool gtAbsVal(long a, long b)
{
  return (abs(a)>abs(b) || (abs(a)==abs(b) && a>b));
}


#else

// This version does not have a preference...
// useful in generating test cases with "bad" dimensions
static bool gtAbsVal(long a, long b)
{
  return (abs(a)>abs(b));
}
#endif

// Returns in gens a generating set for Zm* /<p>, and in ords the
// order of these generators. Return value is the order of p in Zm*.
long findGenerators(vector<long>& gens, vector<long>& ords, long m, long p,
                    const vector<long>& candidates)
{
  gens.clear();
  ords.clear();
  // Compute the generators for (Z/mZ)^*
  vector<long> classes(m);
  vector<long> orders(m);

  for (long i=0; i<m; i++) { // initially each element in its own class
    if (GCD(i,m)!=1) classes[i] = 0; // i is not in (Z/mZ)^*
    else             classes[i] = i;
  }

  // Start building a representation of (Z/mZ)^*, first use the generator p
  conjClasses(classes,p % m,m);  // merge classes that have a factor of p

  // The order of p is the size of the equivalence class of 1
#if 0
  long ordP = std::count(classes.begin(), classes.end(), 1);
    // count(from,to,val) returns # of elements in (from,to) with value=val
#else
  long ordP = 0;
  for (long i = 0; i < lsize(classes); i++)
    if (classes[i] == 1) ordP++;
#endif

  // Compute orders in (Z/mZ)^* /<p> while comparing to (Z/mZ)^*
  long candIdx=0;
  while (true) {
    compOrder(orders,classes,true,m);
    // if the orders of i in Zm* /<p> and Zm* are not the same, then
    // order[i] contains the order in Zm* /<p> with negative sign

    long idx=0;
    if (candIdx<lsize(candidates)) { // try next candidate
      idx = candidates[candIdx++];
      if (abs(orders[idx])<=1)
        idx=0;
    }
    if (idx==0) // no viable candidates supplied externally
      idx = argmax(orders, &gtAbsVal);// find the element with largest order
    long largest = orders[idx];

    if (abs(largest) == 1) break;   // Trivial group, we are done

    // store generator with same order as in (Z/mZ)^*
    gens.push_back(idx);
    ords.push_back(largest);
    conjClasses(classes,idx,m); // merge classes that have a factor of idx
  }
  return ordP;
}

// finding e-th root of unity modulo the current modulus
// VJS: rewritten to be both faster and deterministic,
//  and assumes that current modulus is prime

template<class zp,class zz> void FindPrimRootT(zp &root, unsigned long e)
{
  zz qm1 = zp::modulus()-1;

  assert(qm1 % e == 0);
  
  vector<long> facts;
  factorize(facts,e); // factorization of e

  root = 1;

  for (unsigned long i = 0; i < facts.size(); i++) {
    long p = facts[i];
    long pp = p;
    long ee = e/p;
    while (ee % p == 0) {
      ee = ee/p;
      pp = pp*p;
    }
    // so now we have e = pp * ee, where pp is 
    // the power of p that divides e.
    // Our goal is to find an element of order pp

    PrimeSeq s;
    long q;
    zp qq, qq1;
    long iter = 0;
    do {
      iter++;
      if (iter > 1000000) 
        Error("FindPrimitiveRoot: possible infinite loop?");
      q = s.next();
      conv(qq, q);
      power(qq1, qq, qm1/p);
    } while (qq1 == 1);
    power(qq1, qq, qm1/pp); // qq1 has order pp

    mul(root, root, qq1);
  }

  // independent check that we have an e-th root of unity 
  {
    zp s;

    power(s, root, e);
    if (s != 1) Error("FindPrimitiveRoot: internal error (1)");

    // check that s^{e/p} != 1 for any prime divisor p of e
    for (unsigned long i=0; i<facts.size(); i++) {
      long e2 = e/facts[i];
      power(s, root, e2);   // s = root^{e/p}
      if (s == 1) 
        Error("FindPrimitiveRoot: internal error (2)");
    }
  }
}
// instantiations of the template
void FindPrimitiveRoot(zz_p &r, unsigned long e){FindPrimRootT<zz_p,long>(r,e);}
void FindPrimitiveRoot(ZZ_p &r, unsigned long e){FindPrimRootT<ZZ_p,ZZ>(r,e);}

/* Compute mobius function (naive method as n is small) */
long mobius(long n)
{
  long p,e,arity=0;
  PrimeSeq s;
  while (n!=1)
    { p=s.next();
      e=0;
      while ((n%p==0)) { n=n/p; e++; }
      if (e>1) { return 0; }
      if (e!=0) { arity^=1; }
    }     
  if (arity==0) { return 1; }
  return -1;
}

/* Compute cyclotomic polynomial */
ZZX Cyclotomic(long N)
{
  ZZX Num,Den,G,F;
  NTL::set(Num); NTL::set(Den);
  long m,d;
  for (d=1; d<=N; d++)
    { if ((N%d)==0)
         { clear(G);
           SetCoeff(G,N/d,1); SetCoeff(G,0,-1);
           m=mobius(d);
           if (m==1)       { Num*=G; }
           else if (m==-1) { Den*=G; }
         }
    } 
  F=Num/Den;
  return F;
}

/* Find a primitive root modulo N */
long primroot(long N,long phiN)
{
  long g=2,p;
  PrimeSeq s;
  bool flag=false;

  while (flag==false)
    { flag=true;
      s.reset(1);
      do
        { p=s.next();
          if ((phiN%p)==0)
            { if (PowerMod(g,phiN/p,N)==1)
                { flag=false; }
            }
        }
      while (p<phiN && flag);
      if (flag==false) { g++; }
    }
  return g;
}

long ord(long N,long p)
{
  long o=0;
  while ((N%p)==0)
    { o++;
      N/=p;
    }
  return o;
}

ZZX RandPoly(long n,const ZZ& p)
{ 
  ZZX F; F.SetMaxLength(n);
  ZZ p2;  p2=p>>1;
  for (long i=0; i<n; i++)
    { SetCoeff(F,i,RandomBnd(p)-p2); }
  return F;
}

/* When q=2 maintains the same sign as the input */
void PolyRed(ZZX& out, const ZZX& in, const ZZ& q, bool abs)
{
  // ensure that out has the same degree as in
  out.SetMaxLength(deg(in)+1);               // allocate space if needed
  if (deg(out)>deg(in)) trunc(out,out,deg(in)+1); // remove high degrees

  ZZ q2; q2=q>>1;
  for (long i=0; i<=deg(in); i++)
    { ZZ c=coeff(in,i);
      c %= q;
      if (abs) {
        if (c<0) c += q;
      } 
      else if (q!=2) {
        if (c>q2)  { c=c-q; }
          else if (c<-q2) { c=c+q; }
      }
      else // q=2
        { if (sign(coeff(in,i))!=sign(c))
	    { c=-c; }
        }
      SetCoeff(out,i,c);
    }
}

void PolyRed(ZZX& out, const ZZX& in, long q, bool abs)
{
  // ensure that out has the same degree as in
  out.SetMaxLength(deg(in)+1);               // allocate space if needed
  if (deg(out)>deg(in)) trunc(out,out,deg(in)+1); // remove high degrees

  long q2; q2=q>>1;
  for (long i=0; i<=deg(in); i++)
    { long c=coeff(in,i)%q;
      if (abs)
        { if (c<0) { c=c+q; } }
      else if (q==2)
        { if (coeff(in,i)<0) { c=-c; } }
      else
        { if (c>=q2)  { c=c-q; }
          else if (c<-q2) { c=c+q; }
	}
      SetCoeff(out,i,c);
    }
}

void vecRed(Vec<ZZ>& out, const Vec<ZZ>& in, long q, bool abs)
{
  out.SetLength(in.length());  // allocate space if needed

  for (long i=0; i<in.length(); i++) {
    long c = in[i]%q;
    if (abs)       { if (c<0) c+=q; }
    else if (q==2) { if (in[i]<0) c = -c; }
    else { 
      if (c >= q/2)        c -= q;
      else if (c < -(q/2)) c += q;
    }
    out[i] = c;
  }
}

// multiply the polynomial f by the integer a modulo q
void MulMod(ZZX& out, const ZZX& f, long a, long q, bool abs/*default=true*/)
{
  // ensure that out has the same degree as f
  out.SetMaxLength(deg(f)+1);               // allocate space if needed
  if (deg(out)>deg(f)) trunc(out,out,deg(f)+1); // remove high degrees

  mulmod_precon_t aqinv = PrepMulModPrecon(a, q);
  for (long i=0; i<=deg(f); i++) { 
    long c = rem(coeff(f,i), q);
    c = MulModPrecon(c, a, q, aqinv); // returns c \in [0,q-1]
    if (!abs && c >= q/2)
      c -= q;
    SetCoeff(out,i,c);
  }
}

long is_in(long x,int* X,long sz)
{
  for (long i=0; i<sz; i++)
    { if (x==X[i]) { return i; } }
  return -1;
}

/* Incremental integer CRT for vectors. Expects co-primes p>0,q>0 with q odd,
 * and such that all the entries in vp are in [-p/2,p/2) and all entries in
 * vq are in [0,q-1). Returns in vp the CRT of vp mod p and vq mod q, as
 * integers in [-pq/2, pq/2). Uses the formula:
 *
 *   CRT(vp,p,vq,q) = vp + p*[ (vq-vp)*p^{-1} ]_q
 *
 * where [...]_q means reduction to the interval [-q/2,q/2). As q is odd then
 * this is the same as reducing to [-(q-1)/2,(q-1)/2], hence [...]_q * p is
 * in [-p(q-1)/2, p(q-1)/2], and since vp is in [-p/2,p/2) then the sum is
 * indeed in [-pq/2,pq/2).
 *
 * Returns true if both vectors are of the same length, false otherwise
 */
template <class zzvec>
bool intVecCRT(vec_ZZ& vp, const ZZ& p, const zzvec& vq, long q)
{
  long pInv = InvMod(rem(p,q), q); // p^{-1} mod q
  long n = min(vp.length(),vq.length());
  long q_over_2 = q/2;
  ZZ tmp;
  long vqi;
  mulmod_precon_t pqInv = PrepMulModPrecon(pInv, q);
  for (long i=0; i<n; i++) {
    conv(vqi, vq[i]); // convert to single precision
    long vq_minus_vp_mod_q = SubMod(vqi, rem(vp[i],q), q);

    long delta_times_pInv = MulModPrecon(vq_minus_vp_mod_q, pInv, q, pqInv);
    if (delta_times_pInv > q_over_2) delta_times_pInv -= q;

    mul(tmp, delta_times_pInv, p); // tmp = [(vq_i-vp_i)*p^{-1}]_q * p
    vp[i] += tmp;
  }
  // other entries (if any) are 0 mod q
  for (long i=vq.length(); i<vp.length(); i++) {
    long minus_vp_mod_q = NegateMod(rem(vp[i],q), q);

    long delta_times_pInv = MulModPrecon(minus_vp_mod_q, pInv, q, pqInv);
    if (delta_times_pInv > q_over_2) delta_times_pInv -= q;

    mul(tmp, delta_times_pInv, p); // tmp = [(vq_i-vp_i)*p^{-1}]_q * p
    vp[i] += tmp;
  }
  return (vp.length()==vq.length());
}
// specific instantiations: vq can be vec_long, vec_ZZ, or Vec<zz_p>
template bool intVecCRT(vec_ZZ&, const ZZ&, const vec_ZZ&, long);
template bool intVecCRT(vec_ZZ&, const ZZ&, const vec_long&, long);
template bool intVecCRT(vec_ZZ&, const ZZ&, const Vec<zz_p>&, long);

// MinGW hack
#ifndef lrand48
#if defined(__MINGW32__) || defined(WIN32)
#define drand48() (((double)rand()) / RAND_MAX)
#define lrand48() rand()
#endif
#endif

void sampleHWt(ZZX &poly, long Hwt, long n)
{
  if (n<=0) n=deg(poly)+1; if (n<=0) return;
  clear(poly);          // initialize to zero
  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial

  long b,u,i=0;
  if (Hwt>n) Hwt=n;
  while (i<Hwt) {  // continue until exactly Hwt nonzero coefficients
    u=lrand48()%n; // The next coefficient to choose
    if (IsZero(coeff(poly,u))) { // if we didn't choose it already
      b = lrand48()&2; // b random in {0,2}
      b--;             //   random in {-1,1}
      SetCoeff(poly,u,b);

      i++; // count another nonzero coefficient
    }
  }
  poly.normalize(); // need to call this after we work on the coeffs
}

void sampleSmall(ZZX &poly, long n)
{
  if (n<=0) n=deg(poly)+1; if (n<=0) return;
  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial

  for (long i=0; i<n; i++) {    // Chosse coefficients, one by one
    long u = lrand48();
    if (u&1) {                 // with prob. 1/2 choose between -1 and +1
      u = (u & 2) -1;
      SetCoeff(poly, i, u);
    }
    else SetCoeff(poly, i, 0); // with ptob. 1/2 set to 0
  }
  poly.normalize(); // need to call this after we work on the coeffs
}

void sampleGaussian(ZZX &poly, long n, double stdev)
{
  static double const Pi=4.0*atan(1.0); // Pi=3.1415..
  static long const bignum = 0xfffffff;
  // THREADS: C++11 guarantees these are initialized only once

  if (n<=0) n=deg(poly)+1; if (n<=0) return;
  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial
  for (long i=0; i<n; i++) SetCoeff(poly, i, ZZ::zero());

  // Uses the Box-Muller method to get two Normal(0,stdev^2) variables
  for (long i=0; i<n; i+=2) {
    double r1 = (1+RandomBnd(bignum))/((double)bignum+1);
    double r2 = (1+RandomBnd(bignum))/((double)bignum+1);
    double theta=2*Pi*r1;
    double rr= sqrt(-2.0*log(r2))*stdev;

    assert(rr < 8*stdev); // sanity-check, no more than 8 standard deviations

    // Generate two Gaussians RV's, rounded to integers
    long x = (long) floor(rr*cos(theta) +0.5);
    SetCoeff(poly, i, x);
    if (i+1 < n) {
      x = (long) floor(rr*sin(theta) +0.5);
      SetCoeff(poly, i+1, x);
    }
  }
  poly.normalize(); // need to call this after we work on the coeffs
}

void sampleUniform(ZZX& poly, const ZZ& B, long n)
{
  if (n<=0) n=deg(poly)+1; if (n<=0) return;
  if (B <= 0) {
    clear(poly);
    return;
  }

  poly.SetMaxLength(n); // allocate space for degree-(n-1) polynomial

  ZZ UB, tmp;

  UB =  2*B + 1;
  for (long i = 0; i < n; i++) {
    RandomBnd(tmp, UB);
    tmp -= B; 
    poly.rep[i] = tmp;
  }

  poly.normalize();
}



// ModComp: a pretty lame implementation

void ModComp(ZZX& res, const ZZX& g, const ZZX& h, const ZZX& f)
{
  assert(LeadCoeff(f) == 1);

  ZZX hh = h % f;
  ZZX r = to_ZZX(0);

  for (long i = deg(g); i >= 0; i--) 
    r = (r*hh + coeff(g, i)) % f; 

  res = r;
}

long polyEvalMod(const ZZX& poly, long x, long p)
{
  long ret = 0;
  x %= p; if (x<0) x += p;
  mulmod_precon_t xpinv = PrepMulModPrecon(x, p);
  for (long i=deg(poly); i>=0; i--) {
    long coeff = rem(poly[i], p);
    ret = AddMod(ret, coeff, p);      // Add the coefficient of x^i
    if (i>0) ret = MulModPrecon(ret, x, p, xpinv); // then mult by x
  }
  return ret;
}

static void recursiveInterpolateMod(ZZX& poly, const vec_long& x, vec_long& y,
				    const vec_zz_p& xmod, vec_zz_p& ymod,
				    long p, long p2e)
{
  if (p2e<=1) { // recursion edge condition, mod-1 poly = 0
    clear(poly);
    return;
  }

  // convert y input to zz_p
  for (long j=0; j<y.length(); j++) ymod[j] = to_zz_p(y[j] % p);

  // a polynomial p_i s.t. p_i(x[j]) = i'th p-base digit of poly(x[j])
  zz_pX polyMod;
  interpolate(polyMod, xmod, ymod);    // interpolation modulo p
  ZZX polyTmp; conv(polyTmp, polyMod); // convert to ZZX

  // update ytmp by subtracting the new digit, then dividing by p
  for (long j=0; j<y.length(); j++) {
    y[j] -= polyEvalMod(polyTmp,x[j],p2e); // mod p^e
    if (y[j]<0) y[j] += p2e;
// if (y[j] % p != 0) {
//   cerr << "@@error (p2^e="<<p2e<<"): y["<<j<<"] not divisible by "<<p<< endl;
//   exit(0);
// }
    y[j] /= p;
  } // maybe it's worth optimizing above by using multi-point evaluation

  // recursive call to get the solution of poly'(x)=y mod p^{e-1}
  recursiveInterpolateMod(poly, x, y, xmod, ymod, p, p2e/p);

  // return poly = p*poly' + polyTmp
  poly *= p;
  poly += polyTmp;
}

// Interpolate the integer polynomial such that poly(x[i] mod p)=y[i] (mod p^e)
// It is assumed that the points x[i] are all distinct modulo p
void interpolateMod(ZZX& poly, const vec_long& x, const vec_long& y,
		    long p, long e)
{
  poly = ZZX::zero();       // initialize to zero
  long p2e = power_long(p,e); // p^e

  vec_long ytmp(INIT_SIZE, y.length()); // A temporary writable copy
  for (long j=0; j<y.length(); j++) {
    ytmp[j] = y[j] % p2e;
    if (ytmp[j] < 0) ytmp[j] += p2e;
  }

  zz_pBak bak; bak.save();    // Set the current modulus to p
  zz_p::init(p);

  vec_zz_p xmod(INIT_SIZE, x.length()); // convert to zz_p
  for (long j=0; j<x.length(); j++) xmod[j] = to_zz_p(x[j] % p);

  vec_zz_p ymod(INIT_SIZE, y.length()); // scratch space
  recursiveInterpolateMod(poly, x, ytmp, xmod, ymod, p, p2e);
}

ZZ largestCoeff(const ZZX& f)
{
  ZZ mx = ZZ::zero();
  for (long i=0; i<=deg(f); i++) {
    if (mx < abs(coeff(f,i)))
      mx = abs(coeff(f,i));
  }
  return mx;
}

ZZ sumOfCoeffs(const ZZX& f) // = f(1)
{
  ZZ sum = ZZ::zero();
  for (long i=0; i<=deg(f); i++) sum += coeff(f,i);
  return sum;
}

xdouble coeffsL2Norm(const ZZX& f) // l_2 norm
{
  xdouble s = to_xdouble(0.0);
  for (long i=0; i<=deg(f); i++) {
    xdouble coef = to_xdouble(coeff(f,i));
    s += coef * coef;
  }
  return sqrt(s);
}

// advance the input stream beyond white spaces and a single instance of cc
void seekPastChar(istream& str, int cc)
{
   int c = str.get();
   while (isspace(c)) c = str.get();
   if (c != cc) {
     std::cerr << "Searching for cc='"<<(char)cc<<"' (ascii "<<cc<<")"
	       << ", found c='"<<(char)c<<"' (ascii "<<c<<")\n";
     exit(1);
   }
}

// stuff added relating to linearized polynomials and support routines

// Builds the matrix defining the linearized polynomial transformation.
//
// NTL's current smallint modulus, zz_p::modulus(), is assumed to be p^r,
// for p prime, r >= 1 integer.
//
// After calling this function, one can call ppsolve(C, L, M, p, r) to get
// the coeffecients C for the linearized polynomial represented the linear
// map defined by its action on the standard basis for zz_pE over zz_p:
// for i = 0..zz_pE::degree()-1: x^i -> L[i], where x = (X mod zz_pE::modulus())

void buildLinPolyMatrix(mat_zz_pE& M, long p)
{
   long d = zz_pE::degree();

   M.SetDims(d, d);

   for (long j = 0; j < d; j++) 
      conv(M[0][j], zz_pX(j, 1));

   for (long i = 1; i < d; i++)
      for (long j = 0; j < d; j++)
         M[i][j] = power(M[i-1][j], p);
}

void buildLinPolyMatrix(mat_GF2E& M, long p)
{
   assert(p == 2);

   long d = GF2E::degree();

   M.SetDims(d, d);

   for (long j = 0; j < d; j++) 
      conv(M[0][j], GF2X(j, 1));

   for (long i = 1; i < d; i++)
      for (long j = 0; j < d; j++)
         M[i][j] = power(M[i-1][j], p);
}





// some auxilliary conversion routines

void convert(vec_zz_pE& X, const vector<ZZX>& A)
{
   long n = A.size();
   zz_pX tmp;
   X.SetLength(n);
   for (long i = 0; i < n; i++) {
      conv(tmp, A[i]);
      conv(X[i], tmp); 
   }
} 

void convert(mat_zz_pE& X, const vector< vector<ZZX> >& A)
{
   long n = A.size();

   if (n == 0) {
      long m = X.NumCols();
      X.SetDims(0, m);
      return;
   }

   long m = A[0].size();
   X.SetDims(n, m);

   for (long i = 0; i < n; i++)
      convert(X[i], A[i]);
}

void convert(vector<ZZX>& X, const vec_zz_pE& A)
{
   long n = A.length();
   X.resize(n);
   for (long i = 0; i < n; i++)
      conv(X[i], rep(A[i]));
}

void convert(vector< vector<ZZX> >& X, const mat_zz_pE& A)
{
   long n = A.NumRows();
   X.resize(n);
   for (long i = 0; i < n; i++)
      convert(X[i], A[i]);
}

void convert(NTL::Vec<long>& out, const NTL::ZZX& in)
{
  out.SetLength(in.rep.length());
  for (long i=0; i<out.length(); i++)
    out[i] = conv<long>(in[i]);
}


void convert(NTL::Vec<long>& out, const NTL::zz_pX& in)
{
  out.SetLength(in.rep.length());
  for (long i=0; i<out.length(); i++)
    out[i] = conv<long>(in[i]);
}


void convert(NTL::Vec<long>& out, const NTL::GF2X& in)
{
  out.SetLength(1+deg(in));
  for (long i=0; i<out.length(); i++)
    out[i] = conv<long>(in[i]);
}


void convert(NTL::ZZX& out, const NTL::Vec<long>& in)
{
  out.SetLength(in.length());
  for (long i=0; i<in.length(); i++)
    out[i] = conv<ZZ>(in[i]);
  out.normalize();
}

void convert(NTL::GF2X& out, const NTL::Vec<long>& in)
{
  out.SetLength(in.length());
  for (long i=0; i<in.length(); i++)
    out[i] = conv<GF2>(in[i]);
  out.normalize();
}


void mul(vector<ZZX>& x, const vector<ZZX>& a, long b)
{
   long n = a.size();
   x.resize(n);
   for (long i = 0; i < n; i++) 
      mul(x[i], a[i], b);
}

void div(vector<ZZX>& x, const vector<ZZX>& a, long b)
{
   long n = a.size();
   x.resize(n);
   for (long i = 0; i < n; i++) 
      div(x[i], a[i], b);
}

void add(vector<ZZX>& x, const vector<ZZX>& a, const vector<ZZX>& b)
{
   long n = a.size();
   if (n != (long) b.size()) Error("add: dimension mismatch");
   for (long i = 0; i < n; i++)
      add(x[i], a[i], b[i]);
}

// prime power solver
// zz_p::modulus() is assumed to be p^r, for p prime, r >= 1
// A is an n x n matrix, b is a length n (row) vector,
// and a solution for the matrix-vector equation x A = b is found.
// If A is not inverible mod p, then error is raised.
void ppsolve(vec_zz_pE& x, const mat_zz_pE& A, const vec_zz_pE& b,
             long p, long r) 
{

   if (r == 1) {
      zz_pE det;
      solve(det, x, A, b);
      if (det == 0) Error("ppsolve: matrix not invertible");
      return;
   }

   long n = A.NumRows();
   if (n != A.NumCols()) 
      Error("ppsolve: matrix not square");
   if (n == 0)
      Error("ppsolve: matrix of dimension 0");

   zz_pContext pr_context;
   pr_context.save();

   zz_pEContext prE_context;
   prE_context.save();

   zz_pX G = zz_pE::modulus();

   ZZX GG = to_ZZX(G);

   vector< vector<ZZX> > AA;
   convert(AA, A);

   vector<ZZX> bb;
   convert(bb, b);

   zz_pContext p_context(p);
   p_context.restore();

   zz_pX G1 = to_zz_pX(GG);
   zz_pEContext pE_context(G1);
   pE_context.restore();

   // we are now working mod p...

   // invert A mod p

   mat_zz_pE A1;
   convert(A1, AA);

   mat_zz_pE I1;
   zz_pE det;

   inv(det, I1, A1);
   if (det == 0) {
      Error("ppsolve: matrix not invertible");
   }

   vec_zz_pE b1;
   convert(b1, bb);

   vec_zz_pE y1;
   y1 = b1 * I1;

   vector<ZZX> yy;
   convert(yy, y1);

   // yy is a solution mod p

   for (long k = 1; k < r; k++) {
      // lift solution yy mod p^k to a solution mod p^{k+1}

      pr_context.restore();
      prE_context.restore();
      // we are now working mod p^r

      vec_zz_pE d, y;
      convert(y, yy);

      d = b - y * A;

      vector<ZZX> dd;
      convert(dd, d);

      long pk = power_long(p, k);
      vector<ZZX> ee;
      div(ee, dd, pk);

      p_context.restore();
      pE_context.restore();

      // we are now working mod p

      vec_zz_pE e1;
      convert(e1, ee);
      vec_zz_pE z1;
      z1 = e1 * I1;

      vector<ZZX> zz, ww;
      convert(zz, z1);

      mul(ww, zz, pk);
      add(yy, yy, ww);
   }

   pr_context.restore();
   prE_context.restore();

   convert(x, yy);

   assert(x*A == b);
}

void ppsolve(vec_GF2E& x, const mat_GF2E& A, const vec_GF2E& b,
             long p, long r) 
{
   assert(p == 2 && r == 1);

   GF2E det;
   solve(det, x, A, b);
   if (det == 0) Error("ppsolve: matrix not invertible");
}


// prime power solver
// A is an n x n matrix, we compute its inverse mod p^r. An error is raised
// if A is not inverible mod p. zz_p::modulus() is assumed to be p^r, for
// p prime, r >= 1. Also zz_pE::modulus() is assumed to be initialized.
void ppInvert(mat_zz_pE& X, const mat_zz_pE& A, long p, long r)
{
  if (r == 1) { // use native inversion from NTL
    inv(X, A);    // X = A^{-1}
    return;
  }

  // begin by inverting A modulo p

  // convert to ZZX for a safe transaltion to mod-p objects
  vector< vector<ZZX> > tmp;
  convert(tmp, A);
  { // open a new block for mod-p computation
  ZZX G;
  convert(G, zz_pE::modulus());
  zz_pBak bak_pr; bak_pr.save(); // backup the mod-p^r moduli
  zz_pEBak bak_prE; bak_prE.save();
  zz_p::init(p);   // Set the mod-p moduli
  zz_pE::init(conv<zz_pX>(G));

  mat_zz_pE A1, Inv1;
  convert(A1, tmp);   // Recover A as a mat_zz_pE object modulo p
  inv(Inv1, A1);      // Inv1 = A^{-1} (mod p)
  convert(tmp, Inv1); // convert to ZZX for transaltion to a mod-p^r object
  } // mod-p^r moduli restored on desctuction of bak_pr and bak_prE
  mat_zz_pE XX;
  convert(XX, tmp); // XX = A^{-1} (mod p)

  // Now lift the solution modulo p^r

  // Compute the "correction factor" Z, s.t. XX*A = I - p*Z (mod p^r)
  long n = A.NumRows();
  const mat_zz_pE I = ident_mat_zz_pE(n); // identity matrix
  mat_zz_pE Z = I - XX*A;

  convert(tmp, Z);  // Conver to ZZX to divide by p
  for (long i=0; i<n; i++) for (long j=0; j<n; j++) tmp[i][j] /= p;
  convert(Z, tmp);  // convert back to a mod-p^r object

  // The inverse of A is ( I+(pZ)+(pZ)^2+...+(pZ)^{r-1} )*XX (mod p^r). We use
  // O(log r) products to copmute it as (I+pZ)* (I+(pZ)^2)* (I+(pZ)^4)*...* XX

  long e = NextPowerOfTwo(r); // 2^e is smallest power of two >= r

  Z *= p;                 // = pZ
  mat_zz_pE prod = I + Z; // = I + pZ
  for (long i=1; i<e; i++) {
    sqr(Z, Z);     // = (pZ)^{2^i}
    prod *= (I+Z); // = sum_{j=0}^{2^{i+1}-1} (pZ)^j
  }
  mul(X, prod, XX); // X = A^{-1} mod p^r
  assert(X*A == I);
}

// FIXME: at some point need to make a template for these two functions
// prime power solver
// A is an n x n matrix, we compute its inverse mod p^r. An error is raised
// if A is not inverible mod p. zz_p::modulus() is assumed to be p^r, for
// p prime, r >= 1.
void ppInvert(mat_zz_p& X, const mat_zz_p& A, long p, long r)
{
  if (r == 1) { // use native inversion from NTL
    inv(X, A);    // X = A^{-1}
    return;
  }

  // begin by inverting A modulo p

  // convert to long for a safe transaltion to mod-p objects
  Mat<long> tmp;
  conv(tmp, A);
  { // open a new block for mod-p computation
  zz_pBak bak_pr; bak_pr.save(); // backup the mod-p^r moduli
  zz_p::init(p);   // Set the mod-p moduli

  mat_zz_p A1, Inv1;
  conv(A1, tmp);   // Recover A as a mat_zz_pE object modulo p
  inv(Inv1, A1);      // Inv1 = A^{-1} (mod p)
  conv(tmp, Inv1); // convert to long for transaltion to a mod-p^r object
  } // mod-p^r moduli restored on desctuction of bak_pr and bak_prE
  mat_zz_p XX;
  conv(XX, tmp); // XX = A^{-1} (mod p)

  // Now lift the solution modulo p^r

  // Compute the "correction factor" Z, s.t. XX*A = I - p*Z (mod p^r)
  long n = A.NumRows();
  const mat_zz_p I = ident_mat_zz_p(n); // identity matrix
  mat_zz_p Z = I - XX*A;

  conv(tmp, Z);  // Conver to long to divide by p
  for (long i=0; i<n; i++) for (long j=0; j<n; j++) tmp[i][j] /= p;
  conv(Z, tmp);  // convert back to a mod-p^r object

  // The inverse of A is ( I+(pZ)+(pZ)^2+...+(pZ)^{r-1} )*XX (mod p^r). We use
  // O(log r) products to copmute it as (I+pZ)* (I+(pZ)^2)* (I+(pZ)^4)*...* XX

  long e = NextPowerOfTwo(r); // 2^e is smallest power of two >= r

  Z *= p;                 // = pZ
  mat_zz_p prod = I + Z; // = I + pZ
  for (long i=1; i<e; i++) {
    sqr(Z, Z);     // = (pZ)^{2^i}
    prod *= (I+Z); // = sum_{j=0}^{2^{i+1}-1} (pZ)^j
  }
  mul(X, prod, XX); // X = A^{-1} mod p^r
  assert(X*A == I);
}

void buildLinPolyCoeffs(vec_zz_pE& C_out, const vec_zz_pE& L, long p, long r)
{
   FHE_TIMER_START;
   mat_zz_pE M;
   buildLinPolyMatrix(M, p);

   vec_zz_pE C;
   ppsolve(C, M, L, p, r);

   C_out = C;
   FHE_TIMER_STOP;
}

void buildLinPolyCoeffs(vec_GF2E& C_out, const vec_GF2E& L, long p, long r)
{
   FHE_TIMER_START;
   assert(p == 2 && r == 1);

   mat_GF2E M;
   buildLinPolyMatrix(M, p);

   vec_GF2E C;
   ppsolve(C, M, L, p, r);

   C_out = C;
   FHE_TIMER_STOP;
}

void applyLinPoly(zz_pE& beta, const vec_zz_pE& C, const zz_pE& alpha, long p)
{
   long d = zz_pE::degree();
   assert(d == C.length());

   zz_pE gamma, res;

   gamma = to_zz_pE(zz_pX(1, 1));
   res = C[0]*alpha;
   for (long i = 1; i < d; i++) {
      gamma = power(gamma, p);
      res += C[i]*to_zz_pE(CompMod(rep(alpha), rep(gamma), zz_pE::modulus()));
   }

   beta = res;
}

void applyLinPoly(GF2E& beta, const vec_GF2E& C, const GF2E& alpha, long p)
{
   long d = GF2E::degree();
   assert(d == C.length());

   GF2E gamma, res;

   gamma = to_GF2E(GF2X(1, 1));
   res = C[0]*alpha;
   for (long i = 1; i < d; i++) {
      gamma = power(gamma, p);
      res += C[i]*to_GF2E(CompMod(rep(alpha), rep(gamma), GF2E::modulus()));
   }

   beta = res;
}

// Auxilliary classes to facillitiate faster reduction mod Phi_m(X)
// when the input has degree less than m


static
void LocalCopyReverse(zz_pX& x, const zz_pX& a, long lo, long hi)

   // x[0..hi-lo] = reverse(a[lo..hi]), with zero fill
   // input may not alias output

{
   long i, j, n, m;

   n = hi-lo+1;
   m = a.rep.length();

   x.rep.SetLength(n);

   const zz_p* ap = a.rep.elts();
   zz_p* xp = x.rep.elts();

   for (i = 0; i < n; i++) {
      j = hi-i;
      if (j < 0 || j >= m)
         clear(xp[i]);
      else
         xp[i] = ap[j];
   }

   x.normalize();
} 

static
void LocalCyclicReduce(zz_pX& x, const zz_pX& a, long m)

// computes x = a mod X^m-1

{
   long n = deg(a);
   long i, j;
   zz_p accum;

   if (n < m) {
      x = a;
      return;
   }

   if (&x != &a)
      x.rep.SetLength(m);

   for (i = 0; i < m; i++) {
      accum = a.rep[i];
      for (j = i + m; j <= n; j += m)
         add(accum, accum, a.rep[j]);
      x.rep[i] = accum;
   }

   if (&x == &a)
      x.rep.SetLength(m);

   x.normalize();
}

zz_pXModulus1::zz_pXModulus1(long _m, const zz_pX& _f) 
: m(_m), f(_f), n(deg(f))
{
   assert(m > n);

   specialLogic = (m - n > 10 && m < 2*n);
   build(fm, f);
   
   if (specialLogic) {
      zz_pX P1, P2, P3;

      LocalCopyReverse(P3, f, 0, n);
      InvTrunc(P2, P3, m-n);
      LocalCopyReverse(P1, P2, 0, m-n-1);

      k = NextPowerOfTwo(2*(m-1-n)+1);
      k1 = NextPowerOfTwo(n);

      TofftRep(R0, P1, k); 
      TofftRep(R1, f, k1);
   }
}


void rem(zz_pX& r, const zz_pX& a, const zz_pXModulus1& ff)
{
   if (!ff.specialLogic) {
      rem(r, a, ff.fm);
      return;
   }

   long m = ff.m;
   long n = ff.n;
   long k = ff.k;
   long k1 = ff.k1;
   const fftRep& R0 = ff.R0;
   const fftRep& R1 = ff.R1;

   if (deg(a) < n) {
      r = a;
      return;
   }

   zz_pX P2, P3;

   fftRep R2, R3;

   TofftRep(R2, a, k, n, m-1);
   mul(R2, R2, R0);
   FromfftRep(P3, R2, m-1-n, 2*(m-1-n));
   
   long l = 1L << k1;

   TofftRep(R3, P3, k1);
   mul(R3, R3, R1);
   FromfftRep(P3, R3, 0, n-1);
   LocalCyclicReduce(P2, a, l);
   trunc(P2, P2, n);
   sub(P2, P2, P3);
   r = P2;
}

// Debug printing routines for vectors, ZZX'es, print only a few entries

template<class T> ostream& printVec(ostream& s, const Vec<T>& v,
				    long nCoeffs)
{
  long d = v.length();
  if (d<nCoeffs) return s << v; // just print the whole thing

  // otherwise print only 1st nCoeffs coefficiants
  s << '[';
  for (long i=0; i<nCoeffs-2; i++) s << v[i] << ' ';
  s << "... " << v[d-2] << ' ' << v[d-1] << ']';
  return s;
}
template ostream& printVec(ostream& s, const Vec<zz_p>& v, long nCoeffs);
template ostream& printVec(ostream& s, const Vec<long>& v, long nCoeffs);
template ostream& printVec(ostream& s, const Vec<ZZX>& v, long nCoeffs);

ostream& printZZX(ostream& s, const ZZX& poly, long nCoeffs)
{
  return printVec(s, poly.rep, nCoeffs);
}


