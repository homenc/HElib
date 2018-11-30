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
#include <NTL/xdouble.h>
#include "debugging.h"
#include "norms.h"
#include "FHEContext.h"
#include "Ctxt.h"
#include "EncryptedArray.h"

NTL_CLIENT

FHESecKey* dbgKey = 0;
EncryptedArray* dbgEa = 0;
NTL::ZZX dbg_ptxt;
NTL::Vec<NTL::ZZ> ptxt_pwr; // powerful basis

double realToEstimatedNoise(const Ctxt& ctxt, const FHESecKey& sk)
{
  ZZX p, pp;

  xdouble noiseEst = ctxt.getNoiseBound();
  sk.Decrypt(p, ctxt, pp);
  xdouble actualNoise = coeffsL2Norm(pp);

  return conv<double>(actualNoise/noiseEst);
}

void checkNoise(const Ctxt& ctxt, const FHESecKey& sk, const std::string& msg, double thresh)
{
   double ratio;
   if ((ratio=realToEstimatedNoise(ctxt, sk)) > thresh) {
      std::cerr << "\n*** too much noise: " << msg << ": " << ratio << "\n";
   }
}

xdouble embeddingLargestCoeff(const Ctxt& ctxt, const FHESecKey& sk) 
{
  const FHEcontext& context = ctxt.getContext();
  ZZX p, pp;
  sk.Decrypt(p, ctxt, pp);
  return embeddingLargestCoeff(pp, context.zMStar);
  //return conv<xdouble>(largestCoeff(pp));
}





void decryptAndPrint(ostream& s, const Ctxt& ctxt, const FHESecKey& sk,
		     const EncryptedArray& ea, long flags)
{
  const FHEcontext& context = ctxt.getContext();
  xdouble noiseEst = ctxt.getNoiseBound();
  xdouble modulus = xexp(context.logOfProduct(ctxt.getPrimeSet()));
  vector<ZZX> ptxt;
  ZZX p, pp;
  sk.Decrypt(p, ctxt, pp);
  xdouble actualNoise = coeffsL2Norm(pp);

  s << "plaintext space mod "<<ctxt.getPtxtSpace()
    << ", bitCapacity="<<ctxt.bitCapacity()
    << ", \n           |noise|=q*" << (actualNoise/modulus)
    << ", |noiseEst|=q*" << (noiseEst/modulus);
//#if FFT_IMPL
// FIXME
#if 0
  xdouble embL2 = embeddingL2Norm(pp, ea.getPAlgebra());
  s << ", |noise|_canonical=q*"<< (embL2/modulus);
#endif
  s << endl;

  if (flags & FLAG_PRINT_ZZX) {
    s << "   before mod-p reduction=";
    printZZX(s,pp) <<endl;
  }
  if (flags & FLAG_PRINT_POLY) {
    s << "   after mod-p reduction=";
    printZZX(s,p) <<endl;
  }
  if (flags & FLAG_PRINT_VEC) { // decode to a vector of ZZX
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
  else if (flags & FLAG_PRINT_DVEC) { // decode to a vector of doubles
    const EncryptedArrayCx& eacx = ea.getCx();
    vector<double> v;
    eacx.decrypt(ctxt, sk, v);
    printVec(s<<"           ", v,20)<<endl;
  }
  else if (flags & FLAG_PRINT_XVEC) { // decode to a vector of complex
    const EncryptedArrayCx& eacx = ea.getCx();
    vector<cx_double> v;
    eacx.decrypt(ctxt, sk, v);
    printVec(s<<"           ", v,20)<<endl;
  }
}

bool decryptAndCompare(const Ctxt& ctxt, const FHESecKey& sk,
		       const EncryptedArray& ea, const PlaintextArray& pa)
{
  PlaintextArray ppa(ea);
  ea.decrypt(ctxt, sk, ppa);

  return equals(ea, pa, ppa);
}
