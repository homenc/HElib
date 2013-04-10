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
/* PAlgebraMod.cpp - Implementation of the class template PAlgebraModTmpl
 *
 * This template implements the structure of the plaintext spaces
 * A_2 = R[X]/Phi_m(X) with R=(Z/2Z) or R=(Z/2^rZ). 
 *
 * Phi_m(X) is factored mod 2 as Phi_m(X)=\prod_{t\in T} F_t(X) mod 2, where
 * the F_t's are irreducible modulo 2. An arbitrarily factor is chosen as F_1,
 * then for each t \in T we associate with the index t the factor
 * F_t(X) = GCD(F_1(X^t), Phi_m(X)).
 *
 * Note that fixing a representation of the field K=(Z/2Z)[X]/F_1(X) and
 * letting z be a root of F_1 in K (which is a primitive m-th root of unity
 * in K), we get that F_t is the minimal polynomial of z^{1/t}.
 *
 * The case R=(Z/2^rZ) is implemented by lifting from R=(Z/2Z)
 */
#include <NTL/ZZ.h>
#include <NTL/GF2EX.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2EXFactoring.h>
#include <NTL/lzz_pEX.h>
#include <NTL/ZZXFactoring.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pEXFactoring.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEXFactoring.h>

NTL_CLIENT

#include <algorithm>   // defines count(...), min(...)
#include <iostream>
#include <cassert>

#include "NumbTh.h"
#include "PAlgebra.h"


void PAlgebraModTwo::init(const GF2X& factor)
{
  // Initialization depends on zmStar being initialized
  if (zmStar.M()<2 || zmStar.NSlots()<=0) {
    if (PhimXmod != NULL) delete PhimXmod;
    PhimXmod = NULL;
    factors.kill();
    crtCoeffs.kill();
    return;
  }

  // Compute the factors Ft of Phi_m(X) mod 2, for all t \in T

  if (PhimXmod != NULL) delete PhimXmod;
  PhimXmod = new GF2XModulus(to_GF2X(zmStar.PhimX())); // Phi_m(X) mod 2

  // If a factor is specified, use it, else factorize
  long nSlots = zmStar.NSlots();
  if (deg(factor)>0 && IsZero(*PhimXmod % factor)) {
    factors.SetLength(nSlots);
    factors[0] = factor;
  }
  else if (factors.length()<nSlots
	   || deg(factors[0])<=0 || !IsZero(*PhimXmod % factors[0]))
    EDF(factors, *PhimXmod, zmStar.OrdTwo()); // equal-degree factorization

  // It is left to order the factors according to their representatives

  GF2XModulus F1(factors[0]); // We arbitrarily choose factors[0] as F1
  for (long i=1; i<nSlots; i++) {
    unsigned t =zmStar.ith_rep(i); // Ft is minimal polynomial of x^{1/t} mod F1
    unsigned tInv = rep(inv(to_zz_p(t))); // tInv = t^{-1} mod m
    GF2X X2tInv = PowerXMod(tInv,F1);     // X2tInv = X^{1/t} mod F1
    IrredPolyMod(factors[i], X2tInv, F1);
  }
  /* Debugging sanity-check #1: we should have Ft= GCD(F1(X^t),Phi_m(X))
  GF2XModulus Pm2(*PhimXmod);
  for (i=1; i<nSlots; i++) {
    unsigned t = T[i];
    GF2X X2t = PowerXMod(t,*PhimXmod);  // X2t = X^t mod Phi_m(X)
    GF2X Ft = GCD(CompMod(F1,X2t,Pm2),Pm2);
    if (Ft != factors[i]) {
      cout << "Ft != F1(X^t) mod Phi_m(X), t=" << t << endl;
      exit(0);
    }
  }*******************************************************************/

  // Compute the CRT coefficients for the Ft's
  crtCoeffs.SetLength(nSlots);
  for (long i=0; i<nSlots; i++) {
    GF2X te = *PhimXmod / factors[i]; // \prod_{j\ne i} Fj
    te %= factors[i];              // \prod_{j\ne i} Fj mod Fi
    InvMod(crtCoeffs[i], te, factors[i]);// \prod_{j\ne i} Fj^{-1} mod Fi
  }
}


// Generate the representation of Z[X]/(Phi_m(X),2) for the odd integer m
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
  t = (s*f-1)/g;
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
    t1 = (s1*f - d)/g;
    ss = ss + pk*to_ZZX(s1);
    tt = tt + pk*to_ZZX(t1);
  }

  bak.restore();

  S = to_zz_pX(ss);

  assert((S*F) % G == 1);

}

void PAlgebraMod2r::init(unsigned r)
{
  zz_pBak bak;
  bak.save();

  assert (r>0 && r<NTL_SP_NBITS);   // sanity checks
  assert (&zmStar == &modTwo.zmStar);

  // Initialization depends on zmStar and modTwo being initialized
  if (zmStar.M()<2 || zmStar.NSlots()<=0 || modTwo.PhimXmod==NULL) {
    rr = 1; // Just so we can have something to return from getR()
    if (PhimXmod!=NULL) delete PhimXmod;
    PhimXmod = NULL;
    factors.kill();
    crtCoeffs.kill();
    return;
  }

  // Take the factors and their CRT coefficients mod 2 and lift them to mod 2^r
  zz_p::init(2);

  // convert the factors of Phi_m(X) from GF2X to zz_pX objects with p=2
  long nSlots = zmStar.NSlots();
  vec_zz_pX vzzp;          // no direct conversion from GF2X to zz_pX,
  vec_ZZX vzz;             // need to go via ZZX
  vzzp.SetLength(nSlots);
  vzz.SetLength(nSlots);
  for (long i=0; i<nSlots; i++) {
    vzz[i] = to_ZZX(modTwo.factors[i]);
    conv(vzzp[i], vzz[i]);
  }

  // lift the factors of Phi_m(X) from mod-2 to mod-2^r
  MultiLift(vzz, vzzp, zmStar.PhimX(), r); // defined in NTL::ZZXFactoring

  // Compute the zz_pContext object for mod-2^r arithmetic
  rr = r;
  unsigned two2r = 1UL << r; // compute 2^r
  zz_p::init(two2r);
  mod2rContext.save();

  if (PhimXmod!=NULL) delete PhimXmod;
  PhimXmod = new zz_pXModulus(to_zz_pX(zmStar.PhimX()));
  factors.SetLength(nSlots);
  for (long i=0; i<nSlots; i++)             // Convert from ZZX to zz_pX
    conv(factors[i],vzz[i]);

  // Finally compute the CRT coefficients for the factors
  crtCoeffs.SetLength(nSlots);
  for (long i=0; i<nSlots; i++) {
    zz_pX& fct = factors[i];
    zz_pX te = *PhimXmod / fct; // \prod_{j\ne i} Fj
    te %= fct;                // \prod_{j\ne i} Fj mod Fi
    InvModpr(crtCoeffs[i], te, fct, 2, r);// \prod_{j\ne i} Fj^{-1} mod Fi
  }

  // VJS: bak will restore the original zz_p modulus here, so
  // this function has no side effects.
}

// Returns a vector crt[] such that crt[i] = p mod Ft (with t = T[i])
template<class RX, class vec_RX, class RXM> void
PAlgebraModTmpl<RX,vec_RX,RXM>::CRT_decompose(vector<RX>&crt, const RX&p) const
{
  unsigned nSlots = zmStar.NSlots();
  if (nSlots==0) { crt.resize(0); return; }

  crt.resize(nSlots);
  for (unsigned i=0; i<nSlots; i++)
    rem(crt[i], p, factors[i]); // crt[i] = p % factors[i]
}
// explicit instantiations
template void PAlgebraModTmpl<GF2X,vec_GF2X,GF2XModulus>::CRT_decompose(
                                         vector<GF2X>&crt, const GF2X&p) const;
template void PAlgebraModTmpl<zz_pX,vec_zz_pX,zz_pXModulus>::CRT_decompose(
                                       vector<zz_pX>&crt, const zz_pX&p) const;


template<class RX, class vec_RX, class RXM>
bool PAlgebraModTmpl<RX,vec_RX,RXM>::operator==(const PAlgebraModTmpl<RX,vec_RX,RXM>& other) const
{
  if (zmStar.M() != other.zmStar.M()) return false;
  if (*PhimXmod != *other.PhimXmod) return false;
  if (factors != other.factors) return false;
  if (crtCoeffs  != other.crtCoeffs) return false;

  return true;
}
template bool PAlgebraModTmpl<GF2X,vec_GF2X,GF2XModulus>::operator==(const PAlgebraModTmpl<GF2X,vec_GF2X,GF2XModulus>& other) const;
template bool PAlgebraModTmpl<zz_pX,vec_zz_pX,zz_pXModulus>::operator==(const PAlgebraModTmpl<zz_pX,vec_zz_pX,zz_pXModulus>& other) const;

// Returns p \in R[X]/Phi_m(X) s.t. for every t \in T, the element
// pt = (p mod Ft) \in R[X]/Ft(X) represents the same element as
// alpha \in R[X]/G(X). Must have deg(alpha)<deg(G).
// The maps argument should contain the output of mapToAll(G).
template<class RX, class vec_RX, class RXM>
void PAlgebraModTmpl<RX,vec_RX,RXM>::embedInAllSlots(RX& p, const RX& alpha,
						 const vector<RX>& maps) const
{
  unsigned nSlots = zmStar.NSlots();
  if (nSlots==0 || maps.size()!=nSlots) { p=RX::zero(); return; }

  vector<RX> crt(nSlots); // alloate space for CRT components

  // The i'th CRT component is (p mod F_t) = alpha(maps[i]) mod F_t,
  // where with t=T[i].

  for (unsigned i=0; i<nSlots; i++)   // crt[i] = alpha(maps[i]) mod Ft
    CompMod(crt[i], alpha, maps[i], factors[i]);

  CRT_reconstruct(p,crt); // interpolate to get p
}
// explicit instantiations
template void PAlgebraModTmpl<GF2X,vec_GF2X,GF2XModulus>::embedInAllSlots(
                    GF2X& p, const GF2X& alpha, const vector<GF2X>&maps) const;
template void PAlgebraModTmpl<zz_pX,vec_zz_pX,zz_pXModulus>::embedInAllSlots(
                  zz_pX& p, const zz_pX& alpha, const vector<zz_pX>&maps)const;

// Returns p \in R[X]/Phi_m(X) s.t. for every i<nSlots and t=T[i], the
// element pt = (p mod Ft) \in R[X]/Ft(X) represents the same element
// as alphas[i] \in R[X]/G(X). Must have deg(alphas[i])<deg(G).
template<class RX, class vec_RX, class RXM>
void PAlgebraModTmpl<RX,vec_RX,RXM>::embedInSlots(RX& p, const vector<RX>& alphas,
					      const vector<RX>& maps) const
{
  unsigned nSlots = zmStar.NSlots();
  if (nSlots==0 || maps.size()!=nSlots || alphas.size()!= nSlots) {
    p=RX::zero(); return;
  }
  vector<RX> crt(nSlots); // alloate space for CRT components

  // The i'th CRT component is (p mod F_t) = alphas[i](maps[i]) mod F_t,
  // where with t=T[i].

  for (unsigned i=0; i<nSlots; i++)   // crt[i] = alpha(maps[i]) mod Ft
    CompMod(crt[i], alphas[i], maps[i], factors[i]);

  CRT_reconstruct(p,crt); // interpolate to get p
}
// explicit instantiations
template void PAlgebraModTmpl<GF2X,vec_GF2X,GF2XModulus>::embedInSlots(
                GF2X& p, const vector<GF2X>& alphas, const vector<GF2X>& maps) const;
template void PAlgebraModTmpl<zz_pX,vec_zz_pX,zz_pXModulus>::embedInSlots(
              zz_pX& p, const vector<zz_pX>& alphas, const vector<zz_pX>& maps)const;

// Returns p \in R[X]/Phi_m(X) s.t. for every i<nSlots and t=T[i],
// we have p == crt[i] (mod Ft)
template<class RX, class vec_RX, class RXM> 
void PAlgebraModTmpl<RX,vec_RX,RXM>::CRT_reconstruct(RX&p,vector<RX>&crt) const
{
  unsigned nSlots = zmStar.NSlots();
  if (nSlots==0) { p=RX::zero(); return; }

  // Recall that we have crtCoeffs[i] = \prod_{j \ne i} Fj^{-1} mod Fi
  p = RX::zero();
  for (unsigned i=0; i<nSlots; i++) {
    RX allBut_i = *PhimXmod/ factors[i]; // = \prod_{j \ne i} Fj
    allBut_i *= crtCoeffs[i];      // =1 mod Fi and =0 mod Fj for j \ne i
    MulMod(allBut_i, allBut_i, crt[i], *PhimXmod);
                              // =crt[i] mod Fi and =0 mod Fj for j \ne i
    p += allBut_i;
  }
}
// explicit instantiations
template void PAlgebraModTmpl<GF2X,vec_GF2X,GF2XModulus>::CRT_reconstruct(
                                               GF2X&p, vector<GF2X>&crt) const;
template void PAlgebraModTmpl<zz_pX,vec_zz_pX,zz_pXModulus>::CRT_reconstruct(
                                             zz_pX&p, vector<zz_pX>&crt) const;


/** The following functions are implemented as separate specializations of the
 * template becase they use additional types such as GF2E/GF2EX etc, and I
 * cannot figure out a way for C++ templates to do something like "if GF2X then
 * use GF2E and if zz_pX then use zz_pXE" without replicating all the code.
 ****************************************************************************/

// r \in (Z/2Z)[X]/Ft(X) represents the same as X \in (Z/2Z)[X]/G(X). The
// optional rF1 contains the output of mapToF1, to speed this operation.
template<> void 
PAlgebraModTmpl<GF2X,vec_GF2X,GF2XModulus>::mapToFt(GF2X& r,
			     const GF2X& G,unsigned t,const GF2X* rF1) const
{
  int i = zmStar.indexOfRep(t);
  if (i < 0) { r=GF2X::zero(); return; }

  if (rF1==NULL) {               // Compute the representation "from scratch"
    GF2E::init(factors[i]);      // work with the extension field GF_2[X]/Ft(X)
    GF2EX Ga=to_GF2EX((GF2X&)G); // G as a polynomial over the extension field
    r=rep(FindRoot(Ga));         // Find a root of G in this field
    return;
  }
  // if rF1 is set, then use it instead, setting r = rF1(X^t) mod Ft(X)
  GF2XModulus Ft(factors[i]);
  //  long tInv = InvMod(t,m);
  GF2X X2t = PowerXMod(t,Ft);    // X2t = X^t mod Ft
  r = CompMod(*rF1,X2t,Ft);      // r = F1(X2t) mod Ft

  /* Debugging sanity-check: G(r)=0 in the extension field (Z/2Z)[X]/Ft(X)
  GF2E::init(factors[i]);
  GF2EX Ga=to_GF2EX((GF2X&)G); // G as a polynomial over the extension field
  GF2E ra =to_GF2E(r);         // r is an element in the extension field
  eval(ra,Ga,ra);  // ra = Ga(ra)
  if (!IsZero(ra)) {// check that Ga(r)=0 in this extension field
    cout << "rF1(X^t) mod Ft(X) != root of G mod Ft, t=" << t << endl;
    exit(0);    
  }*******************************************************************/
}
template<> void 
PAlgebraModTmpl<zz_pX,vec_zz_pX,zz_pXModulus>::mapToFt(zz_pX& r,
			     const zz_pX& G,unsigned t,const zz_pX* rF1) const
{
  int i = zmStar.indexOfRep(t);
  if (i < 0) { r=zz_pX::zero(); return; }

  if (rF1==NULL) {              // Compute the representation "from scratch"
    zz_pE::init(factors[i]);    // work with the extension field GF_2[X]/Ft(X)
    zz_pEX Ga=to_zz_pEX((zz_pX&)G);// G is polynomial over the extension field
    r=rep(FindRoot(Ga));        // Find a root of G in this field
    return;
  }
  // if rF1 is set, then use it instead, setting r = rF1(X^t) mod Ft(X)
  zz_pXModulus Ft(factors[i]);
  //  long tInv = InvMod(t,m);
  zz_pX X2t = PowerXMod(t,Ft);    // X2t = X^t mod Ft
  r = CompMod(*rF1,X2t,Ft);      // r = F1(X2t) mod Ft

  /* Debugging sanity-check: G(r)=0 in the extension field (Z/2Z)[X]/Ft(X)
  zz_pE::init(factors[i]);
  zz_pEX Ga=to_zz_pEX((zz_pX&)G);// G as a polynomial over the extension field
  zz_pE ra =to_zz_pE(r);         // r is an element in the extension field
  eval(ra,Ga,ra);  // ra = Ga(ra)
  if (!IsZero(ra)) {// check that Ga(r)=0 in this extension field
    cout << "rF1(X^t) mod Ft(X) != root of G mod Ft, t=" << t << endl;
    exit(0);    
  }*******************************************************************/
}


// return an array such that alphas[i] \in R[X]/G(X) represent the
// same element as rt = (p mod Ft) \in R[X]/Ft(X) where t=T[i].
// If the optional map pointer is non-NULL, then (*map)[i] contains
// the output of mapToFt(...G,T[i]).
template<> void PAlgebraModTmpl<GF2X,vec_GF2X,GF2XModulus>::decodePlaintext(
   vector<GF2X>& alphas, const GF2X& ptxt, const GF2X &G, 
   const vector<GF2X>& maps) const
{
  unsigned nSlots = zmStar.NSlots();
  if (nSlots==0 || maps.size()!= nSlots ||deg(G)<1|| zmStar.OrdTwo()%deg(G)!=0){
    alphas.resize(0); return;
  }
  // First decompose p into CRT componsnt
  vector<GF2X> CRTcomps(nSlots); // allocate space for CRT component
  CRT_decompose(CRTcomps, ptxt);  // CRTcomps[i] = p mod facors[i]

  // Next convert each entry in CRTcomps[] back to base-G representation
  // (this is roughly the inverse of the embedInSlots operation)

  // maps contains all the base-G ==> base-Fi maps

  // SHAI: The code below is supposed to be Nigel's code, borrowed partly
  // from the constructor Algebra::Algebra(...) and partly from the
  // function AElement::to_type(2). I don't really unerstand that code,
  // so hopefully nothing was lost in translation

  alphas.resize(nSlots);
  GF2E::init(G); // the extension field R[X]/G(X)
  for (unsigned i=0; i<nSlots; i++) {
    // We need to lift Fi from R[Y] to (R[X]/G(X))[Y]
    GF2EX Qi; conv(Qi,(GF2X&)factors[i]);
    vec_GF2EX FRts=SFBerlekamp(Qi); // factor Fi over Z_2[X]/G(X)

    // need to choose the right factor, the one that gives us back X
    int j;
    for (j=0; j<FRts.length(); j++) { 
      // lift maps[i] to (R[X]/G(X))[Y] and reduce mod j'th factor of Fi
      GF2EX GRti = to_GF2EX((GF2X&)maps[i]);
      GRti %= FRts[j];

      if (IsX(rep(coeff(GRti,0)))) { // is GRti == X?
	Qi = FRts[j];                // If so, we found the right factor
	break;
      } // If this does not happen then move to the next factor of Fi
    }

    GF2EX te; conv(te,CRTcomps[i]);   // lift i'th CRT componnet to mod G(X)
    te %= Qi;       // reduce CRTcomps[i](Y) mod Qi(Y), over (Z_2[X]/G(X))

    // the free term (no Y component) should be our answer (as a poly(X))
    alphas[i] = rep(coeff(te,0));
  }
}

template<> void PAlgebraModTmpl<zz_pX,vec_zz_pX,zz_pXModulus>::decodePlaintext(
   vector<zz_pX>& alphas, const zz_pX& ptxt, const zz_pX &G,
   const vector<zz_pX>& maps) const
{
  unsigned nSlots = zmStar.NSlots();
  if (nSlots==0 || maps.size()!= nSlots ||deg(G)<1|| zmStar.OrdTwo()%deg(G)!=0){
    alphas.resize(0); return;
  }
  // First decompose p into CRT componsnt
  vector<zz_pX> CRTcomps(nSlots); // allocate space for CRT component
  CRT_decompose(CRTcomps, ptxt);  // CRTcomps[i] = p mod facors[i]

  if (IsX(G)) {
    // if modulus is not prime, this is the only case
    // that will actually work
    alphas = CRTcomps;
    return;
  }


  // Next convert each entry in CRTcomps[] back to base-G representation
  // (this is roughly the inverse of the embedInSlots operation)

  // maps contains all the base-G ==> base-Fi maps

  // SHAI: The code below is supposed to be Nigel's code, borrowed partly
  // from the constructor Algebra::Algebra(...) and partly from the
  // function AElement::to_type(2). I don't really unerstand that code,
  // so hopefully nothing was lost in translation

  alphas.resize(nSlots);
  zz_pE::init(G); // the extension field R[X]/G(X)
  for (unsigned i=0; i<nSlots; i++) {
    // We need to lift Fi from R[Y] to (R[X]/G(X))[Y]
    zz_pEX Qi; conv(Qi,(zz_pX&)factors[i]);
    vec_zz_pEX FRts=SFCanZass(Qi); // factor Fi over Z_2[X]/G(X)

    // need to choose the right factor, the one that gives us back X
    int j;
    for (j=0; j<FRts.length(); j++) { 
      // lift maps[i] to (R[X]/G(X))[Y] and reduce mod j'th factor of Fi
      zz_pEX GRti = to_zz_pEX((zz_pX&)maps[i]);
      GRti %= FRts[j];

      if (IsX(rep(coeff(GRti,0)))) { // is GRti == X?
	Qi = FRts[j];                // If so, we found the right factor
	break;
      } // If this does not happen then move to the next factor of Fi
    }

    zz_pEX te; conv(te,CRTcomps[i]);   // lift i'th CRT componnet to mod G(X)
    te %= Qi;       // reduce CRTcomps[i](Y) mod Qi(Y), over (Z_2[X]/G(X))

    // the free term (no Y component) should be our answer (as a poly(X))
    alphas[i] = rep(coeff(te,0));
  }
}
