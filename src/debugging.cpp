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
// debugging.cpp - debugging utilities
#include <NTL/ZZ.h>
NTL_CLIENT
#include "FHE.h"
#include "EncryptedArray.h"

#define FLAG_PRINT_ZZX  1
#define FLAG_PRINT_POLY 2
#define FLAG_PRINT_VEC  4

FHESecKey* dbgKey=NULL;
EncryptedArray* dbgEa=NULL;
ZZX dbg_ptxt;
Vec<ZZ> ptxt_pwr;

template<class T> ostream& printVec(ostream& s, const Vec<T>& v, long nCoeffs)
{
  long d = v.length();
  if (d<nCoeffs) return s << v; // just print the whole thing

  // otherwise print only 1st nCoeffs coefficiants
  s << '[';
  for (long i=0; i<nCoeffs-2; i++) s << v[i] << ' ';
  s << "... " << v[d-2] << ' ' << v[d-1] << ']';
  return s;
}

ostream& printZZX(ostream& s, const ZZX& poly, long nCoeffs=40)
{
  return printVec(s, poly.rep, nCoeffs);
  /*  long d = deg(poly);
  if (d<nCoeffs) return s << poly; // just print the whole thing

  // otherwise print only 1st nCoeffs coefficiants
  s << '[';
  for (long i=0; i<nCoeffs-2; i++) s << poly[i] << ' ';
  s << "... " << poly[d-1] << ' ' << poly[d] << ']';
  return s; */
}

/*void baseRep(Vec<long>& rep, long nDigits, ZZ num, long base)
{
  rep.SetLength(nDigits);
  for (long j=0; j<nDigits; j++) {
    rep[j] = rem(num, base);
    if (rep[j] > base/2)         rep[j] -= base;
    else if (rep[j] < -(base/2)) rep[j] += base;
    num = (num - rep[j]) / base;
  }
  }*/

void decryptAndPrint(ostream& s, const Ctxt& ctxt, const FHESecKey& sk,
		     const EncryptedArray& ea, long flags)
{
  const FHEcontext& context = ctxt.getContext();
  xdouble noiseEst = sqrt(ctxt.getNoiseVar());
  xdouble modulus = xexp(context.logOfProduct(ctxt.getPrimeSet()));
  vector<ZZX> ptxt;
  ZZX p, pp;
  sk.Decrypt(p, ctxt, pp);

  s << "plaintext space mod "<<ctxt.getPtxtSpace()
       << ", level="<<ctxt.findBaseLevel()
       << ", \n           |noise|=q*" << (coeffsL2Norm(pp)/modulus)
       << ", |noiseEst|=q*" << (noiseEst/modulus)
       <<endl;

  if (flags & FLAG_PRINT_ZZX) {
    s << "   before mod-p reduction=";
    printZZX(s,pp) <<endl;
  }
  if (flags & FLAG_PRINT_POLY) {
    s << "   after mod-p reduction=";
    printZZX(s,p) <<endl;
  }
  if (flags & FLAG_PRINT_VEC) {
    ea.decode(ptxt, p);
    if (ea.getAlMod().getTag() == PA_zz_p_tag
	&& ctxt.getPtxtSpace() != ea.getAlMod().getPPowR()) {
      long g = GCD(ctxt.getPtxtSpace(), ea.getAlMod().getPPowR());
      for (long i=0; i<ea.size(); i++)
	PolyRed(ptxt[i], g, true);
    }
    s << "   decoded to ";
    if (deg(p) < 40) // just pring the whole thing
      s << ptxt << endl;
    else if (ptxt.size()==1) // a single slot
      printZZX(s, ptxt[0]) <<endl;
    else { // print first and last slots
      printZZX(s, ptxt[0],20) << "--";
      printZZX(s, ptxt[ptxt.size()-1], 20) <<endl;      
    }
  }
}
