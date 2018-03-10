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
// debugging.cpp - debugging utilities
#include "debugging.h"
NTL_CLIENT

FHESecKey* dbgKey;
EncryptedArray* dbgEa;
NTL::ZZX dbg_ptxt;
NTL::Vec<NTL::ZZ> ptxt_pwr; // powerful basis

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

bool decryptAndCompare(const Ctxt& ctxt, const FHESecKey& sk,
		       const EncryptedArray& ea, const NewPlaintextArray& pa)
{
  NewPlaintextArray ppa(ea);
  ea.decrypt(ctxt, sk, ppa);

  return equals(ea, pa, ppa);
}
