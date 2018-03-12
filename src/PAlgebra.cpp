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

#include <algorithm>   // defines count(...), min(...)

#include "PAlgebra.h"
#include "hypercube.h"
#include "timing.h"

#include <NTL/ZZXFactoring.h>
#include <NTL/GF2EXFactoring.h>
#include <NTL/lzz_pEXFactoring.h>

#include <NTL/BasicThreadPool.h>

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


bool PAlgebra::operator==(const PAlgebra& other) const
{
  if (m != other.m) return false;
  if (p != other.p) return false;

  return true;
}


unsigned long PAlgebra::exponentiate(const vector<unsigned long>& exps,
				bool onlySameOrd) const
{
  if (isDryRun()) return 1;
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
  cout << "m = " << m << ", p = " << p;
  if (isDryRun()) { cout << " (dry run)\n"; return; }
  cout << ", phi(m) = " << phiM << endl;
  cout << "  ord(p)=" << ordP << endl;

  unsigned long i;
  for (i=0; i<gens.size(); i++) if (gens[i]) {
      cout << "  generator " << gens[i] << " has order ("
           << (SameOrd(i)? "=":"!") << "= Z_m^*) of " 
	   << OrderOf(i) << endl;
  }
  if (cube.getSize()<100) {
    cout << "  T = [";
    for (i=0; i<T.size(); i++) cout << T[i] << " ";
    cout << "]\n";
  }
}


PAlgebra::PAlgebra(unsigned long mm, unsigned long pp,  
                   const vector<long>& _gens, const vector<long>& _ords )
{
  assert( ProbPrime(pp) );
  assert( (mm % pp) != 0 );
  assert( mm < NTL_SP_BOUND );
  assert( mm > 1 );

  cM  = 1.0; // default value for the ring constant
  m = mm;
  p = pp;

  long k = NextPowerOfTwo(m);
  if (mm == (1UL << k))
    pow2 = k;
  else
    pow2 = 0;

  // For dry-run, use a tiny m value for the PAlgebra tables
  if (isDryRun()) mm = (p==3)? 4 : 3;

  // Compute the generators for (Z/mZ)^* (defined in NumbTh.cpp)

  std::vector<long> tmpOrds;
  if (_gens.size()>0 && _gens.size() == _ords.size() && !isDryRun()) {
    // externally supplied generator,orders
    tmpOrds = _ords;
    this->gens = _gens;
    this->ordP = multOrd(pp, mm);
  }
  else
    // treat externally supplied generators (if any) as candidates
    this->ordP = findGenerators(this->gens, tmpOrds, mm, pp, _gens);

  // Record for each generator gi whether it has the same order in
  // ZM* as in Zm* /(p,g1,...,g_{i-1})
  resize(native, lsize(tmpOrds));
  for (long j=0; j<lsize(tmpOrds); j++) {
    native[j] = (tmpOrds[j]>0);
    tmpOrds[j] = abs(tmpOrds[j]);
  }
  cube.initSignature(tmpOrds); // set hypercume with these dimensions

  phiM = ordP * getNSlots();

  // Allocate space for the various arrays
  T.resize(getNSlots());
  Tidx.assign(mm,-1);    // allocate m slots, initialize them to -1
  zmsIdx.assign(mm,-1);  // allocate m slots, initialize them to -1
  zmsRep.resize(phiM);
  long i, idx;
  for (i=idx=0; i<(long)mm; i++) {
    if (GCD(i,mm)==1) {
      zmsIdx[i] = idx++;
      zmsRep[zmsIdx[i]] = i;
    }
  }

  // Now fill the Tidx translation table. We identify an element t \in T
  // with its representation t = \prod_{i=0}^n gi^{ei} mod m (where the
  // gi's are the generators in gens[]) , represent t by the vector of
  // exponents *in reverse order* (en,...,e1,e0), and order these vectors
  // in lexicographic order.

  // FIXME: is the comment above about reverse order true?
  // It doesn't seem like it to me, VJS.
  // The comment about reverse order is correct, SH.

  // buffer is initialized to all-zero, which represents 1=\prod_i gi^0
  vector<unsigned long> buffer(gens.size()); // temporaty holds exponents
  i = idx = 0;
  long ctr = 0;
  do {
    ctr++;
    unsigned long t = exponentiate(buffer);

    assert(GCD(t,mm) == 1); // sanity check for user-supplied gens
    assert(Tidx[t] == -1);

    T[i] = t;       // The i'th element in T it t
    Tidx[t] = i++;  // the index of t in T is i

    // increment buffer by one (in lexigoraphic order)
  } while (nextExpVector(buffer)); // until we cover all the group

  assert(ctr == getNSlots()); // sanity check for user-supplied gens

  PhimX = Cyclotomic(mm); // compute and store Phi_m(X)
  //  pp_factorize(mFactors,mm); // prime-power factorization from NumbTh.cpp
}


long PAlgebra::frobenuisPow(long j) const
{
  return PowerMod(mcMod(p, m), j, m);
  // Don't forget to reduce p mod m!!
}

long PAlgebra::genToPow(long i, long j) const
{
  assert(i >= -1 && i < LONG(gens.size()));

  long res;
  if (i == -1)
    res = frobenuisPow(j);
  else
    res = PowerMod(gens[i], j, m);

  return res;
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
  long m = zMStar.getM();

  // For dry-run, use a tiny m value for the PAlgebra tables
  if (isDryRun()) m = (p==3)? 4 : 3;

  assert(r > 0);

  ZZ BigPPowR = power_ZZ(p, r);
  assert(BigPPowR.SinglePrecision());
  pPowR = to_long(BigPPowR);

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
    unsigned long t =zMStar.ith_rep(i); // Ft is minimal poly of x^{1/t} mod F1
    unsigned long tInv = InvMod(t, m);  // tInv = t^{-1} mod m
    RX X2tInv = PowerXMod(tInv,F1);     // X2tInv = X^{1/t} mod F1
    NTL::IrredPolyMod(localFactors[i], X2tInv, F1);
          // IrredPolyMod(X,P,Q) returns in X the minimal polynomial of P mod Q
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

  genCrtTable();
  genMaskTable();
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

  if (isDryRun()) {
    crt.clear();
    return;
  }
  crt.resize(nSlots);
  for (unsigned long i=0; i<nSlots; i++)
    rem(crt[i], H, factors[i]); // crt[i] = H % factors[i]
}

template<class type>
void PAlgebraModDerived<type>::embedInAllSlots(RX& H, const RX& alpha, 
                                            const MappingData<type>& mappingData) const
{
  if (isDryRun()) {
    H = RX::zero();
    return;
  }
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

    // FIXME: should update this to use matrix_maps, but this routine
    // isn't actually used anywhere

    for (long i=0; i<nSlots; i++)   // crt[i] = alpha(maps[i]) mod Ft
      CompMod(crt[i], alpha, mappingData.maps[i], factors[i]);
  }

  CRT_reconstruct(H,crt); // interpolate to get H
  FHE_TIMER_STOP;
}

template<class type>
void PAlgebraModDerived<type>::embedInSlots(RX& H, const vector<RX>& alphas, 
                                         const MappingData<type>& mappingData) const
{
  if (isDryRun()) {
    H = RX::zero();
    return;
  }
  FHE_TIMER_START;

  long nSlots = zMStar.getNSlots();
  assert(lsize(alphas) == nSlots);

  long d = mappingData.degG;
  for (long i = 0; i < nSlots; i++) assert(deg(alphas[i]) < d); 
 
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

    FHE_NTIMER_START(CompMod);

#if 0
    for (long i: range(nSlots)) {
      if (deg(alphas[i]) <= 0) 
        crt[i] = alphas[i];
      else
        CompMod(crt[i], alphas[i], mappingData.maps[i], factors[i]);
    }
#else
    vec_R in, out;

    for (long i: range(nSlots)) {
      if (deg(alphas[i]) <= 0) 
        crt[i] = alphas[i];
      else {
        VectorCopy(in, alphas[i], d);
        mul(out, in, mappingData.matrix_maps[i]);
        conv(crt[i], out);
      }
    }
#endif
  }

  CRT_reconstruct(H,crt); // interpolate to get p

  FHE_TIMER_STOP;
}

template<class type>
void PAlgebraModDerived<type>::CRT_reconstruct(RX& H, vector<RX>& crt) const
{
  if (isDryRun()) {
    H = RX::zero();
    return;
  }
  FHE_TIMER_START;
  long nslots = zMStar.getNSlots();


  const vector<RX>& ctab = crtTable;

  clear(H);
  RX tmp1, tmp2;

  bool easy = true;
  for (long i = 0; i < nslots; i++) 
    if (!IsZero(crt[i]) && !IsOne(crt[i])) {
      easy = false;
      break;
    }
    
  if (easy) {
    for (long i=0; i<nslots; i++) 
      if (!IsZero(crt[i])) 
        H += ctab[i];
  }
  else {
    vector<RX> crt1;
    crt1.resize(nslots);
    for (long i = 0; i < nslots; i++)
       MulMod(crt1[i], crt[i], crtCoeffs[i], factors[i]);

    evalTree(H, crtTree, crt1, 0, nslots);
  }
  FHE_TIMER_STOP;
}

template<class type>
void PAlgebraModDerived<type>::mapToFt(RX& w,
			     const RX& G,unsigned long t,const RX* rF1) const
{
  if (isDryRun()) {
    w = RX::zero();
    return;
  }
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
  mappingData.degG = deg(mappingData.G);
  long d = deg(G);
  long ordp = zMStar.getOrdP();

  long nSlots = zMStar.getNSlots();
  long m = zMStar.getM();

  mappingData.maps.resize(nSlots);

  mapToF1(mappingData.maps[0],mappingData.G); // mapping from base-G to base-F1
  for (long i=1; i<nSlots; i++)
    mapToFt(mappingData.maps[i], mappingData.G, zMStar.ith_rep(i), &(mappingData.maps[0])); 


  // create matrices to streamline CompMod operations
  mappingData.matrix_maps.resize(nSlots);
  for (long i: range(nSlots)) {
    mat_R& mat = mappingData.matrix_maps[i];
    mat.SetDims(d, ordp);
    RX pow;
    pow = 1;
    for (long j: range(d)) {
      VectorCopy(mat[j], pow, ordp);
      if (j < d-1) MulMod(pow, pow, mappingData.maps[i], factors[i]);
    }
  }


  REBak bak; bak.save(); 
  RE::init(mappingData.G);
  mappingData.contextForG.save();

  if (deg(mappingData.G)==1) return;

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
  if (isDryRun()) {
    alphas.assign(nSlots, RX::zero());
    return;
  }

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
void PAlgebraModDerived<type>::decodeSomePlaintexts(std::vector<RX> &alphas,
                                                    const RX& ptxt,
                                                    const std::vector<long> &positions,
                                                    const MappingData<type> &mappingData) const
{
  long nSlots = zMStar.getNSlots();
  size_t nPositions = positions.size();
  if (isDryRun()) {
    alphas.assign(nPositions, RX::zero());
    return;
  }

  std::vector<RX> CRTcomps(nPositions);
  for (size_t i = 0; i < nPositions; i++)
      rem(CRTcomps[i], ptxt, factors[positions[i]]); // CRTcomp = H % factors[positions[i]]

  if (mappingData.degG == 1) {
    alphas = CRTcomps;
    return;
  }

  alphas.resize(nPositions);
  REBak bak; bak.save(); mappingData.contextForG.restore();
  for (size_t i = 0; i < nPositions; i++) {
    REX te;
    conv(te, CRTcomps[i]);
    te %= mappingData.rmaps[positions[i]];

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
void PAlgebraModDerived<type>::genMaskTable() 
{
  // This is only called by the constructor, which has already
  // set the zz_p context

  RX tmp1;
  
  maskTable.resize(zMStar.numOfGens());
  for (long i = 0; i < (long)zMStar.numOfGens(); i++) {
    long ord = zMStar.OrderOf(i);
    maskTable[i].resize(ord+1);
    maskTable[i][ord] = 0;
    for (long j = ord-1; j >= 1; j--) {
      // initialize mask that is 1 whenever the ith coordinate is at least j
      // Note: maskTable[i][0] = constant 1, maskTable[i][ord] = constant 0
      maskTable[i][j] = maskTable[i][j+1];
      for (long k = 0; k < (long)zMStar.getNSlots(); k++) {
         if (zMStar.coordinate(i, k) == j) {
           div(tmp1, PhimXMod, factors[k]);
           mul(tmp1, tmp1, crtCoeffs[k]);
           add(maskTable[i][j], maskTable[i][j], tmp1);
         }
      }
    }
    maskTable[i][0] = 1;
  }
}

// code for generating crt tables
// the tables are generated "on demand"

template<class type> 
void PAlgebraModDerived<type>::genCrtTable() 
{
  // This is only called by the constructor, which has already
  // set the zz_p context

  long nslots = zMStar.getNSlots();
  crtTable.resize(nslots);
  for (long i = 0; i < nslots; i++) {
    RX allBut_i = PhimXMod / factors[i]; // = \prod_{j \ne i }Fj
    allBut_i *= crtCoeffs[i]; // = 1 mod Fi and = 0 mod Fj for j \ne i
    crtTable[i] = allBut_i;
  }

  buildTree(crtTree, 0, nslots);
}

template<class type> 
void PAlgebraModDerived<type>::
  buildTree(shared_ptr< TNode<RX> >& res, long offset, long extent) const
{
  if (extent == 1)
    res = buildTNode<RX>(nullTNode<RX>(), nullTNode<RX>(), 
                            factors[offset]);
  else {
    long half = extent/2;
    shared_ptr< TNode<RX> > left, right;
    buildTree(left, offset, half);
    buildTree(right, offset+half, extent-half);
    RX data = left->data * right->data;
    res = buildTNode<RX>(left, right, data);
  }
}

template<class type> 
void PAlgebraModDerived<type>::evalTree(RX& res,
              shared_ptr< TNode<RX> > tree,
              const vector<RX>& crt1,
              long offset, long extent) const
{
  if (extent == 1) 
    res = crt1[offset];
  else {
    long half = extent/2;
    RX lres, rres;
    evalTree(lres, tree->left, crt1, offset, half);
    evalTree(rres, tree->right, crt1, offset+half, extent-half);
    RX tmp1, tmp2;
    mul(tmp1, lres, tree->right->data);
    mul(tmp2, rres, tree->left->data);
    add(tmp1, tmp1, tmp2);
    res = tmp1;
  }
}

// Explicit instantiation

template class PAlgebraModDerived<PA_GF2>;
template class PAlgebraModDerived<PA_zz_p>;
