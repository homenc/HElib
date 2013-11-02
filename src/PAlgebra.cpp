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
#include <NTL/ZZ.h>    // defines GCD(...)
#include <NTL/ZZX.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/GF2EXFactoring.h>
#include <NTL/lzz_pEXFactoring.h>


#include <climits>     // defines INT_MAX
#include <cstring>
#include <algorithm>   // defines count(...), min(...)
#include <iostream>
#include <cassert>

#include "NumbTh.h"    // defines argmax(...)
#include "PAlgebra.h"
#include "timing.h"


NTL_CLIENT


// polynomials are sorted lexicographically, with the
// constant term being the "most significant"

template<class RX> bool poly_comp(const RX& a, const RX& b) 
{
  long na = deg(a) + 1;
  long nb = deg(b) + 1;

  long i = 0;
  while (i < na && i < nb && coeff(a, i) == coeff(b, i)) i++;

  if (i < na && i < nb)
    return coeff(a, i) < coeff(b, i);
  else 
    return na < nb;
}

namespace NTL {

// for some weird reason, these need to be in either the std or NTL 
// namespace; otherwise, the compiler won't find them...

bool operator<(GF2 a, GF2 b) { return rep(a) < rep(b); }
bool operator<(zz_p a, zz_p b) { return rep(a) < rep(b); }

bool operator<(const GF2X& a, const GF2X& b) { return poly_comp(a, b); }
bool operator<(const zz_pX& a, const zz_pX& b) { return poly_comp(a, b); }

bool operator<(const GF2E& a, const GF2E& b) { return rep(a) < rep(b); }
bool operator<(const zz_pE& a, const zz_pE& b) { return rep(a) < rep(b); }

bool operator<(const GF2EX& a, const GF2EX& b) { return poly_comp(a, b); }
bool operator<(const zz_pEX& a, const zz_pEX& b) { return poly_comp(a, b); }

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
 * one, containing just a a single element.
 *
 * A twist is that initially we only add an element as a new generator if its
 * order in the current quotient group is the same as in the original group
 * (Z/mZ)^* (these are the gi's). Only after no such elements are available
 * we begin to use generators that do not have the same order as in Z_m^*.
 **/


// The function conjClasses(classes,g,m) unifies equivalence classes that have
// elements which are a factor of g apart, the pivot of the unified class is
// the smallest element in that class. 

static 
void conjClasses(vector<unsigned long>& classes, unsigned long g, unsigned long m)
{
  for (unsigned long i=0; i<m; i++) {
    if (classes[i]==0) continue; // i \notin (Z/mZ)^*

    if (classes[i]<i) { // i is not a pivot, updated its pivot
      classes[i] = classes[classes[i]];
      continue;
    }

    // If i is a pivot, update other pivots to point to it
    unsigned long ii = i;
    unsigned long gg = g;
    unsigned long jj = MulMod(ii, gg, m);
    while (classes[jj] != i) {
      classes[classes[jj]]= i; // Merge the equivalence classes of j and i

      // Note: if classes[j]!=j then classes[j] will be updated later,
      //       when we get to i=j and use the code for "i not pivot".

      jj = MulMod(jj, g, m);
    }
  }
}


// The function compOrder(orders, classes,flag,m) computes the order of elements
// of the quotient group, relative to current equivalent classes. If flag==1
// then also check if the order is the same as in (Z/mZ)^* and store the order
// with negative sign if not.

static 
void compOrder(vector<long>& orders, vector<unsigned long>& classes, bool flag, 
               unsigned long m)
{
  orders[0] = -INT_MAX;
  orders[1] = 0;
  for (unsigned long i=2; i<m; i++) {
    if (classes[i] <= 1) { // ignore i not in Z_m^* and order-0 elements
      orders[i] = (classes[i]==1)? 0 : -INT_MAX;
      continue;
    }

    // If not comparing order with (Z/mZ)^*, only compute the order of pivots

    if (!flag && classes[i]<i){          // not a pivot
      orders[i] = orders[classes[i]];
      continue;
    }

    // For an element i>1, the order is at least 2
    unsigned long ii = i;
    unsigned long jj = MulMod(ii, ii, m);
    long ord = 2;
    while (classes[jj] != 1) {
      jj = MulMod(jj, ii, m); // next element in <i>
      ord++;    // count how many steps until we reach 1
    }

    // When we get here we have classes[j]==1, so if j!=1 it means that the
    // order of i in the quotient group is smaller than its order in the
    // entire group Z_m^*. If the flag is set then we store orders[i] = -ord.
    
    if (flag && jj != 1) ord = -ord; // order in Z_m^* is larger than ord
    orders[i] = ord;
  }
}

bool PAlgebra::operator==(const PAlgebra& other) const
{
  if (m != other.m) return false;
  if (p != other.p) return false;

  return true;
}

bool PAlgebra::nextExpVector(vector<unsigned long>& buffer) const
{
  // increment the vector in lexicographic order
  for (long i=gens.size()-1; i>=0; i--) {
    if (i>=(long)buffer.size()) continue; // sanity check
    // increment current index, set all the ones after it to zero
    if (buffer[i] < OrderOf(i)-1) { 
      buffer[i]++;
      for (unsigned long j=i+1; j<buffer.size(); j++) buffer[j] = 0;
      return true;  // succeeded in incrementing the vector
    }
    // if buffer[i] >= OrderOf(i)-1, mover to previous index i
  }
  return false;     // cannot increment the vector anymore
}

long PAlgebra::coordinate(long i, long k) const
{
  long t = ith_rep(k); // element of Zm^* representing the k'th slot

  // dLog returns the representation of t along the generators, so the
  // i'th entry there is the coordinate relative to i'th geneator
  return dLog(t)[i];
}

long PAlgebra::addCoord(long i, long k, long offset) const
{
  assert(k >= 0 && k < (long) nSlots);
  assert(i >= 0 && i < (long) gens.size());
  
  offset = offset % ((long) OrderOf(i));
  if (offset < 0) offset += OrderOf(i);
  
  long k_i = coordinate(i, k);
  long k_i1 = (k_i + offset) % OrderOf(i);
  
  long k1 = k + (k_i1 - k_i) * prods[i+1];
  
  return k1;
  
}

unsigned long PAlgebra::exponentiate(const vector<unsigned long>& exps,
				bool onlySameOrd) const
{
  unsigned long t = 1;
  unsigned long n = min(exps.size(),gens.size());
  for (unsigned long i=0; i<n; i++) {
    if (onlySameOrd && !SameOrd(i)) continue;
    unsigned long g = PowerMod(gens[i] ,exps[i], m); 
    t = MulMod(t, g, m);
  }
  return t;
}

void PAlgebra::printout() const
{
  unsigned long i;
  cout << "m = " << m << ", p = " << p << ", phi(m) = " << phiM << endl;
  cout << "  ord(p)=" << ordP << endl;
  for (i=0; i<gens.size(); i++) if (gens[i]) {
      cout << "  generator " << gens[i] << " has order ("
           << (SameOrd(i)? "=":"!") << "= Z_m^*) of " 
	   << OrderOf(i) << endl;
  }
  if (qGrpOrd()<100) {
    cout << "  T = [";
    for (i=0; i<T.size(); i++) cout << T[i] << " ";
    cout << "]\n";
  }
}

// Generate the representation of Z_m^* for a given odd integer m
// and plaintext base p
PAlgebra::PAlgebra(unsigned long mm, unsigned long pp)
{
  m = mm;
  p = pp;

  //  assert( (m&1) == 1 );
  assert( ProbPrime(p) );
  //  assert( m > p && (m % p) != 0 );
  assert( (m % p) != 0 );
  assert( m < NTL_SP_BOUND );

  // Compute the generators for (Z/mZ)^*
  vector<unsigned long> classes(m);
  vector<long> orders(m);

  unsigned long i;
  for (i=0; i<m; i++) { // initially each element in its own class
    if (GCD(i,m)!=1) 
      classes[i] = 0; // i is not in (Z/mZ)^*
    else 
      classes[i] = i;
  }

  // Start building a representation of (Z/mZ)^*, first use the generator p
  conjClasses(classes,p,m);  // merge classes that have a factor of p

  // The order of p is the size of the equivalence class of 1
  ordP = (unsigned long) count (classes.begin(), classes.end(), 1);

  // Compute orders in (Z/mZ)^*/<p> while comparing to (Z/mZ)^*
  long idx, largest;
  while (true) {
    compOrder(orders,classes,true,m);
    idx = argmax(orders);      // find the element with largest order
    largest = orders[idx];

    if (largest <= 0) break;   // stop comparing to order in (Z/mZ)^*

    // store generator with same order as in (Z/mZ)^*
    gens.push_back(idx);
    ords.push_back(largest);
    conjClasses(classes,idx,m); // merge classes that have a factor of idx
  }
  // Compute orders in (Z/mZ)^*/<p> without comparing to (Z/mZ)^*
  while (true) {
    compOrder(orders,classes,false,m);
    idx = argmax(orders);      // find the element with largest order
    largest = orders[idx];

    if (largest <= 0) break;   // we have the trivial group, we are done

    // store generator with different order than (Z/mZ)^*
    gens.push_back(idx);
    ords.push_back(-largest);  // store with negative sign
    conjClasses(classes,idx,m);  // merge classes that have a factor of idx
  }

  nSlots = qGrpOrd();
  phiM = ordP * nSlots;

  // Allocate space for the various arrays
  T.resize(nSlots);
  dLogT.resize(nSlots*gens.size());
  Tidx.assign(m,-1);    // allocate m slots, initialize them to -1
  zmsIdx.assign(m,-1);  // allocate m slots, initialize them to -1
  for (i=idx=0; i<m; i++) if (GCD(i,m)==1) zmsIdx[i] = idx++;

  // Now fill the Tidx and dLogT translation tables. We identify an element
  // t\in T with its representation t = \prod_{i=0}^n gi^{ei} mod m (where
  // the gi's are the generators in gens[]) , represent t by the vector of
  // exponents *in reverse order* (en,...,e1,e0), and order these vectors
  // in lexicographic order.

  // buffer is initialized to all-zero, which represents 1=\prod_i gi^0
  vector<unsigned long> buffer(gens.size()); // temporaty holds exponents
  i = idx = 0;
  do {
    unsigned long t = exponentiate(buffer);
    for (unsigned long j=0; j<buffer.size(); j++) dLogT[idx++] = buffer[j];
    T[i] = t;       // The i'th element in T it t
    Tidx[t] = i++;  // the index of t in T is i

    // increment buffer by one (in lexigoraphic order)
  } while (nextExpVector(buffer)); // until we cover all the group

  PhimX = Cyclotomic(m); // compute and store Phi_m(X)

  // initialize prods array
  long ndims = gens.size();
  prods.resize(ndims+1);
  prods[ndims] = 1;
  for (long j = ndims-1; j >= 0; j--) {
    prods[j] = OrderOf(j) * prods[j+1];
  } 
}

/***********************************************************************

  PAlgebraMod stuff....

************************************************************************/

PAlgebraModBase *buildPAlgebraMod(const PAlgebra& zMStar, long r)
{
  unsigned long p = zMStar.getP();
  assert(r > 0);

  if (p == 2 && r == 1) 
    return new PAlgebraModDerived<PA_GF2>(zMStar, r);
  else
    return new  PAlgebraModDerived<PA_zz_p>(zMStar, r);
}


template<class T> 
void PAlgebraLift(const ZZX& phimx, const T& lfactors, T& factors, T& crtc, long r);



// Missing NTL functionality

void EDF(vec_zz_pX& v, const zz_pX& f, long d)
{
   EDF(v, f, PowerXMod(zz_p::modulus(), f), d);
}

zz_pEX FrobeniusMap(const zz_pEXModulus& F)
{
  return PowerXMod(zz_pE::cardinality(), F);
}


template<class type> 
PAlgebraModDerived<type>::PAlgebraModDerived(const PAlgebra& _zMStar, long _r) 
  : zMStar(_zMStar), r(_r)

{
  long p = zMStar.getP();
  assert(r > 0);
  ZZ BigPPowR = power_ZZ(p, r);
  assert(BigPPowR.SinglePrecision());
  pPowR = to_long(BigPPowR);

  long m = zMStar.getM();
  long nSlots = zMStar.getNSlots();

  RBak bak; bak.save();
  SetModulus(p);

  // Compute the factors Ft of Phi_m(X) mod p, for all t \in T

  RX phimxmod;

  conv(phimxmod, zMStar.getPhimX()); // Phi_m(X) mod p

  vec_RX localFactors;

  EDF(localFactors, phimxmod, zMStar.getOrdP()); // equal-degree factorization
  
  RX* first = &localFactors[0];
  RX* last = first + localFactors.length();
  RX* smallest = min_element(first, last);
  swap(*first, *smallest);

  // We make the lexicographically smallest factor have index 0.
  // The remaining factors are ordered according to their representives.

  RXModulus F1(localFactors[0]); 
  for (long i=1; i<nSlots; i++) {
    unsigned long t =zMStar.ith_rep(i); // Ft is minimal polynomial of x^{1/t} mod F1
    unsigned long tInv = InvMod(t, m);  // tInv = t^{-1} mod m
    RX X2tInv = PowerXMod(tInv,F1);     // X2tInv = X^{1/t} mod F1
    IrredPolyMod(localFactors[i], X2tInv, F1);
  }
  /* Debugging sanity-check #1: we should have Ft= GCD(F1(X^t),Phi_m(X))
  for (i=1; i<nSlots; i++) {
    unsigned long t = T[i];
    RX X2t = PowerXMod(t,phimxmod);  // X2t = X^t mod Phi_m(X)
    RX Ft = GCD(CompMod(F1,X2t,phimxmod),phimxmod);
    if (Ft != localFactors[i]) {
      cout << "Ft != F1(X^t) mod Phi_m(X), t=" << t << endl;
      exit(0);
    }
  }*******************************************************************/

  if (r == 1) {
    build(PhimXMod, phimxmod);
    factors = localFactors;
    pPowRContext.save();

    // Compute the CRT coefficients for the Ft's
    crtCoeffs.SetLength(nSlots);
    for (long i=0; i<nSlots; i++) {
      RX te = phimxmod / factors[i]; // \prod_{j\ne i} Fj
      te %= factors[i];              // \prod_{j\ne i} Fj mod Fi
      InvMod(crtCoeffs[i], te, factors[i]); // \prod_{j\ne i} Fj^{-1} mod Fi
    }
  }
  else {
    PAlgebraLift(zMStar.getPhimX(), localFactors, factors, crtCoeffs, r);
    RX phimxmod1;
    conv(phimxmod1, zMStar.getPhimX());
    build(PhimXMod, phimxmod1);
    pPowRContext.save();
  }

  // set factorsOverZZ
  factorsOverZZ.resize(nSlots);
  for (long i = 0; i < nSlots; i++)
    conv(factorsOverZZ[i], factors[i]);
}

// Assumes current zz_p modulus is p^r
// computes S = F^{-1} mod G via Hensel lifting
void InvModpr(zz_pX& S, const zz_pX& F, const zz_pX& G, long p, long r)
{
  ZZX ff, gg, ss, tt;

  ff = to_ZZX(F); 
  gg = to_ZZX(G);

  zz_pBak bak;
  bak.save();
  zz_p::init(p);

  zz_pX f, g, s, t;
  f = to_zz_pX(ff);
  g = to_zz_pX(gg);
  s = InvMod(f, g);
  t = (1-s*f)/g;
  assert(s*f + t*g == 1);
  ss = to_ZZX(s);
  tt = to_ZZX(t);

  ZZ pk = to_ZZ(1);

  for (long k = 1; k < r; k++) {
    // lift from p^k to p^{k+1}
    pk = pk * p;

    assert(divide(ss*ff + tt*gg - 1, pk));

    zz_pX d = to_zz_pX( (1 - (ss*ff + tt*gg))/pk );
    zz_pX s1, t1;
    s1 = (s * d) % g;
    t1 = (d-s1*f)/g;
    ss = ss + pk*to_ZZX(s1);
    tt = tt + pk*to_ZZX(t1);
  }

  bak.restore();

  S = to_zz_pX(ss);

  assert((S*F) % G == 1);
}

template<class T> 
void PAlgebraLift(const ZZX& phimx, const T& lfactors, T& factors, T& crtc, long r)
{
   Error("uninstatiated version of PAlgebraLift");
}

// This specialized version of PAlgebraLift does the hensel
// lifting needed to finish off the initialization.
// It assumes the zz_p modulus is initialized to p
// when called, and leaves it set to p^r

template<> 
void PAlgebraLift(const ZZX& phimx, const vec_zz_pX& lfactors, vec_zz_pX& factors, vec_zz_pX& crtc, long r)
{
  long p = zz_p::modulus(); 
  long nSlots = lfactors.length();


  vec_ZZX vzz;             // need to go via ZZX

  // lift the factors of Phi_m(X) from mod-2 to mod-2^r
  if (lfactors.length() > 1)
    MultiLift(vzz, lfactors, phimx, r); // defined in NTL::ZZXFactoring
  else {
    vzz.SetLength(1);
    vzz[0] = phimx;
  }

  // Compute the zz_pContext object for mod p^r arithmetic
  zz_p::init(power_long(p, r));

  zz_pX phimxmod = to_zz_pX(phimx);
  factors.SetLength(nSlots);
  for (long i=0; i<nSlots; i++)             // Convert from ZZX to zz_pX
    conv(factors[i], vzz[i]);

  // Finally compute the CRT coefficients for the factors
  crtc.SetLength(nSlots);
  for (long i=0; i<nSlots; i++) {
    zz_pX& fct = factors[i];
    zz_pX te = phimxmod / fct; // \prod_{j\ne i} Fj
    te %= fct;                // \prod_{j\ne i} Fj mod Fi
    InvModpr(crtc[i], te, fct, p, r);// \prod_{j\ne i} Fj^{-1} mod Fi
  }

}

// Returns a vector crt[] such that crt[i] = p mod Ft (with t = T[i])
template<class type> 
void PAlgebraModDerived<type>::CRT_decompose(vector<RX>& crt, const RX& H) const
{
  unsigned long nSlots = zMStar.getNSlots();

  crt.resize(nSlots);
  for (unsigned long i=0; i<nSlots; i++)
    rem(crt[i], H, factors[i]); // crt[i] = H % factors[i]
}

template<class type>
void PAlgebraModDerived<type>::embedInAllSlots(RX& H, const RX& alpha, 
                                            const MappingData<type>& mappingData) const
{
  FHE_TIMER_START;
  long nSlots = zMStar.getNSlots();

  vector<RX> crt(nSlots); // alloate space for CRT components

  // The i'th CRT component is (H mod F_t) = alpha(maps[i]) mod F_t,
  // where with t=T[i].

  
  if (IsX(mappingData.G) || deg(alpha) <= 0) {
    // special case...no need for CompMod, which is
    // is not optimized for this case

    for (long i=0; i<nSlots; i++)   // crt[i] = alpha(maps[i]) mod Ft
      crt[i] = ConstTerm(alpha);
  }
  else {
    // general case...

    for (long i=0; i<nSlots; i++)   // crt[i] = alpha(maps[i]) mod Ft
      CompMod(crt[i], alpha, mappingData.maps[i], factors[i]);
  }

  CRT_reconstruct(H,crt); // interpolate to get H
  FHE_TIMER_START;
}

template<class type>
void PAlgebraModDerived<type>::embedInSlots(RX& H, const vector<RX>& alphas, 
                                         const MappingData<type>& mappingData) const
{
  FHE_TIMER_START;

  long nSlots = zMStar.getNSlots();
  assert(lsize(alphas) == nSlots);

  for (long i = 0; i < nSlots; i++) assert(deg(alphas[i]) < mappingData.degG); 
 
  vector<RX> crt(nSlots); // alloate space for CRT components

  // The i'th CRT component is (H mod F_t) = alphas[i](maps[i]) mod F_t,
  // where with t=T[i].

  if (IsX(mappingData.G)) {
    // special case...no need for CompMod, which is
    // is not optimized for this case

    for (long i=0; i<nSlots; i++)   // crt[i] = alpha(maps[i]) mod Ft
      crt[i] = ConstTerm(alphas[i]);
  }
  else {
    // general case...still try to avoid CompMod when possible,
    // which is the common case for encoding masks

    for (long i=0; i<nSlots; i++) {   // crt[i] = alpha(maps[i]) mod Ft
      if (deg(alphas[i]) <= 0) 
        crt[i] = alphas[i];
      else
        CompMod(crt[i], alphas[i], mappingData.maps[i], factors[i]);
    }
  }

  CRT_reconstruct(H,crt); // interpolate to get p

  FHE_TIMER_STOP;
}

template<class type>
void PAlgebraModDerived<type>::CRT_reconstruct(RX& H, vector<RX>& crt) const
{
  FHE_TIMER_START;
  long nSlots = zMStar.getNSlots();

  const vector<RX>& ctab = getCrtTable();

  clear(H);
  RX tmp;
  for (long i=0; i<nSlots; i++) {
    // optimize special cases 0 and 1, which are common
    if (!IsZero(crt[i])) {
      if (IsOne(crt[i]))
        H += ctab[i];
      else {
        MulMod(tmp, ctab[i], crt[i], PhimXMod);
        H += tmp;
      }
    }
  }
  FHE_TIMER_STOP;
}

template<class type>
void PAlgebraModDerived<type>::mapToFt(RX& w,
			     const RX& G,unsigned long t,const RX* rF1) const
{
  long i = zMStar.indexOfRep(t);
  if (i < 0) { clear(w); return; }


  if (rF1==NULL) {               // Compute the representation "from scratch"
    // special case
    if (G == factors[i]) {
      SetX(w);
      return;
    }

    //special case
    if (deg(G) == 1) {
      w = -ConstTerm(G);
      return;
    }

    // the general case: currently only works when r == 1
    assert(r == 1);  

    REBak bak; bak.save();
    RE::init(factors[i]);        // work with the extension field GF_p[X]/Ft(X)
    REX Ga;
    conv(Ga, G);                 // G as a polynomial over the extension field

    vec_RE roots;
    FindRoots(roots, Ga);        // Find roots of G in this field
    RE* first = &roots[0];
    RE* last = first + roots.length();
    RE* smallest = min_element(first, last);
                                // make a canonical choice
    w=rep(*smallest);         
    return;
  }
  // if rF1 is set, then use it instead, setting w = rF1(X^t) mod Ft(X)
  RXModulus Ft(factors[i]);
  //  long tInv = InvMod(t,m);
  RX X2t = PowerXMod(t,Ft);    // X2t = X^t mod Ft
  w = CompMod(*rF1,X2t,Ft);      // w = F1(X2t) mod Ft

  /* Debugging sanity-check: G(w)=0 in the extension field (Z/2Z)[X]/Ft(X)
  RE::init(factors[i]);
  REX Ga;
  conv(Ga, G); // G as a polynomial over the extension field
  RE ra;
  conv(ra, w);         // w is an element in the extension field
  eval(ra,Ga,ra);  // ra = Ga(ra)
  if (!IsZero(ra)) {// check that Ga(w)=0 in this extension field
    cout << "rF1(X^t) mod Ft(X) != root of G mod Ft, t=" << t << endl;
    exit(0);    
  }*******************************************************************/
}

template<class type> 
void PAlgebraModDerived<type>::mapToSlots(MappingData<type>& mappingData, const RX& G) const 
{
  assert(deg(G) > 0 && zMStar.getOrdP() % deg(G) == 0);
  assert(LeadCoeff(G) == 1);


  mappingData.G = G;
  mappingData.degG = deg(G);


  long nSlots = zMStar.getNSlots();
  long m = zMStar.getM();

  mappingData.maps.resize(nSlots);


  mapToF1(mappingData.maps[0],G); // mapping from base-G to base-F1
  for (long i=1; i<nSlots; i++)
    mapToFt(mappingData.maps[i], G, zMStar.ith_rep(i), &(mappingData.maps[0])); 

  if (deg(G)==1) return;

  REBak bak; bak.save(); 
  RE::init(G);
  mappingData.contextForG.save();

  mappingData.rmaps.resize(nSlots);

  if (G == factors[0]) {
    // an important special case

    for (long i = 0; i < nSlots; i++) {
        long t = zMStar.ith_rep(i);
        long tInv = InvMod(t, m);

        RX ct_rep;
        PowerXMod(ct_rep, tInv, G);
        
        RE ct;
        conv(ct, ct_rep);

        REX Qi;
        SetCoeff(Qi, 1, 1);
        SetCoeff(Qi, 0, -ct);

        mappingData.rmaps[i] = Qi;
    }
  }
  else
  {
    // the general case: currently only works when r == 1

    assert(r == 1);

    vec_REX FRts;
    for (long i=0; i<nSlots; i++) {
      // We need to lift Fi from R[Y] to (R[X]/G(X))[Y]
      REX  Qi;
      long t, tInv=0;

      if (i == 0) {
        conv(Qi,factors[i]);
        FRts=EDF(Qi, FrobeniusMap(Qi), deg(Qi)/deg(G)); 
        // factor Fi over GF(p)[X]/G(X)
      }
      else {
        t = zMStar.ith_rep(i);
        tInv = InvMod(t, m);
      }

      // need to choose the right factor, the one that gives us back X
      long j;
      for (j=0; j<FRts.length(); j++) { 
        // lift maps[i] to (R[X]/G(X))[Y] and reduce mod j'th factor of Fi

        REX FRtsj;
        if (i == 0) 
           FRtsj = FRts[j];
        else {
            REX X2tInv = PowerXMod(tInv, FRts[j]);
            IrredPolyMod(FRtsj, X2tInv, FRts[j]);
        }

        // FRtsj is the jth factor of factors[i] over the extension field.
        // For j > 0, we save some time by computing it from the jth factor 
        // of factors[0] via a minimal polynomial computation.
        
        REX GRti;
        conv(GRti, mappingData.maps[i]);
        GRti %= FRtsj;

        if (IsX(rep(ConstTerm(GRti)))) { // is GRti == X?
          Qi = FRtsj;                // If so, we found the right factor
          break;
        } // If this does not happen then move to the next factor of Fi
      }

      assert(j < FRts.length());
      mappingData.rmaps[i] = Qi;
    }
  }
}

template<class type> 
void PAlgebraModDerived<type>::decodePlaintext(
   vector<RX>& alphas, const RX& ptxt, const MappingData<type>& mappingData) const
{
  long nSlots = zMStar.getNSlots();

  // First decompose p into CRT components
  vector<RX> CRTcomps(nSlots); // allocate space for CRT component
  CRT_decompose(CRTcomps, ptxt);  // CRTcomps[i] = p mod facors[i]

  if (mappingData.degG==1) {
    alphas = CRTcomps;
    return;
  }

  alphas.resize(nSlots);

  REBak bak; bak.save(); mappingData.contextForG.restore();

  for (long i=0; i<nSlots; i++) {
    REX te; 
    conv(te, CRTcomps[i]);   // lift i'th CRT componnet to mod G(X)
    te %= mappingData.rmaps[i];  // reduce CRTcomps[i](Y) mod Qi(Y), over (Z_2[X]/G(X))

    // the free term (no Y component) should be our answer (as a poly(X))
    alphas[i] = rep(ConstTerm(te));
  }
}

template<class type> 
void PAlgebraModDerived<type>::
buildLinPolyCoeffs(vector<RX>& C, const vector<RX>& L,
                   const MappingData<type>& mappingData) const
{
  REBak bak; bak.save(); mappingData.contextForG.restore();

  long d = RE::degree();
  long p = zMStar.getP();

  assert(lsize(L) == d);

  vec_RE LL;
  LL.SetLength(d);

  for (long i = 0; i < d; i++)
    conv(LL[i], L[i]);

  vec_RE CC;
  ::buildLinPolyCoeffs(CC, LL, p, r);

  C.resize(d);
  for (long i = 0; i < d; i++)
    C[i] = rep(CC[i]);
}

// code for generating mask tables
// the tables are generated "on demand"

template<class type> 
void PAlgebraModDerived<type>::genMaskTable() const
{
  if (maskTable.size() > 0) return;

  RBak bak; bak.save(); restoreContext();

  // strip const
  vector< vector< RX > >& mtab = (vector< vector< RX > >&) maskTable;

  RX tmp1;
  
  mtab.resize(zMStar.numOfGens());
  for (long i = 0; i < (long)zMStar.numOfGens(); i++) {
    long ord = zMStar.OrderOf(i);
    mtab[i].resize(ord+1);
    mtab[i][ord] = 0;
    for (long j = ord-1; j >= 1; j--) {
      // initialize mask that is 1 whenever the ith coordinate is at least j
      // Note: mtab[i][0] = constant 1, mtab[i][ord] = constant 0
      mtab[i][j] = mtab[i][j+1];
      for (long k = 0; k < (long)zMStar.getNSlots(); k++) {
         if (zMStar.coordinate(i, k) == j) {
           div(tmp1, PhimXMod, factors[k]);
           mul(tmp1, tmp1, crtCoeffs[k]);
           add(mtab[i][j], mtab[i][j], tmp1);
         }
      }
    }
    mtab[i][0] = 1;
  }
}

// code for generating crt tables
// the tables are generated "on demand"

template<class type> 
void PAlgebraModDerived<type>::genCrtTable() const
{
  if (crtTable.size() > 0) return;

  RBak bak; bak.save(); restoreContext();

  // strip const
  vector< RX >& tab = (vector< RX >&) crtTable;

  
  long nslots = zMStar.getNSlots();
  tab.resize(nslots);
  for (long i = 0; i < nslots; i++) {
    RX allBut_i = PhimXMod / factors[i]; // = \prod_{j \ne i }Fj
    allBut_i *= crtCoeffs[i]; // = 1 mod Fi and = 0 mod Fj for j \ne i
    tab[i] = allBut_i;
  }
}


// Explicit instantiation

template class PAlgebraModDerived<PA_GF2>;
template class PAlgebraModDerived<PA_zz_p>;

