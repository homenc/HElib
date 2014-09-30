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
#include <NTL/ZZ.h>
NTL_CLIENT
#include "EncryptedArray.h"
//#include "EvalMap.h"
#include "AltEvalMap.h"
#include "powerful.h"

//#define DEBUG_PRINTOUT
#ifdef DEBUG_PRINTOUT /*********** Debugging utilities **************/
extern FHESecKey* dbgKey;
extern EncryptedArray* dbgEa;
extern ZZX dbg_ptxt;
extern Vec<ZZ> ptxt_pwr;

#define FLAG_PRINT_ZZX  1
#define FLAG_PRINT_POLY 2
#define FLAG_PRINT_VEC  4
extern void decryptAndPrint(ostream& s, const Ctxt& ctxt, const FHESecKey& sk,
			    const EncryptedArray& ea, long flags=0);
extern void baseRep(Vec<long>& rep, long nDigits, ZZ num, long base=2);
#endif                /********* End Debugging utilities **************/


/************* Some local functions *************/
static void x2iInSlots(ZZX& poly, long i,
		       vector<ZZX>& xVec, const EncryptedArray& ea);

static void initUnpackEncoding(vector< vector<ZZX> >& unpackSlotEncoding,
			       const EncryptedArray& ea);

// Make every entry of vec divisible by p^e by adding/subtracting multiples
// of p^r and q, while keeping the added multiples small. 
template<class VecInt>
long makeDivisible(VecInt& vec, long p2e, long p2r, long q, double alpha=0.42);

// Extract digits from fully packed slots
void extractDigitsPacked(Ctxt& ctxt, long botHigh, long r, long ePrime,
			 const vector< vector<ZZX> >& unpackSlotEncoding);
 
// bootstrap a ciphertext to reduce noise
void FHEPubKey::reCrypt(Ctxt &ctxt)
{
  FHE_TIMER_START;
  assert(bootstrapKeyID>=0); // check that we have bootstrapping data

  long p = getContext().zMStar.getP();
  long r = getContext().alMod.getR();
  long p2r = getContext().alMod.getPPowR();
  long p2e = bootstrapEkey.getPtxtSpace()/p2r;
  // the bootstrapping key is encrypted relative to plaintext space p^{e+r}.
  long e = round( log((double)p2e)/log((double)p) );

  // The bootstrapping data in the context is for plaintext space p^{e+r-e'}
  long ePrime = e + r - context.bootstrapPAM->getR();
  long p2ePrime = power_long(p,ePrime);
  long q = p2e+1;

  assert(e>=r);
#ifdef DEBUG_PRINTOUT
  cerr << "reCrypt: p="<<p<<", r="<<r<<", e="<<e<<" ePrime="<<ePrime
       << ", q="<<q<<endl;
#endif

  // can only bootstrap ciphertext with plaintext-space dividing p^r
  long ptxtSpace = ctxt.getPtxtSpace();
  assert(p2r % ptxtSpace == 0);

  FHE_NTIMER_START(preProcess);

  // Make sure that this ciphertxt is in canonical form
  if (!ctxt.inCanonicalForm()) ctxt.reLinearize();

  // Mod-switch down if needed
  const IndexSet& curPS = ctxt.getPrimeSet();
  if (curPS.card()>2) { // leave only bottom two primes
    long frst = curPS.first();
    long scnd = curPS.next(frst);
    IndexSet s(frst,scnd);
    ctxt.modDownToSet(s);
  }

  // key-switch to the bootstrapping key
  ctxt.reLinearize(bootstrapKeyID);

#ifdef DEBUG_PRINTOUT
  if (dbgKey && dbgEa) {
    ZZX poly;
    dbgKey->Decrypt(poly,ctxt);
    if (poly!=dbg_ptxt)
      cerr << "  Decryption error after key-switching to bootsrapping key\n";
    else {
      cerr << "  After key-switching to bootsrapping key: ";
      decryptAndPrint(cerr, ctxt, *dbgKey, *dbgEa, 0);
    }
  }
#endif

  // "raw mod-switch" to the bootstrapping mosulus q=p^e+1.
  vector<ZZX> zzParts; // the mod-switched parts, in ZZX format
  double noise = ctxt.rawModSwitch(zzParts, q);
  noise = sqrt(noise);

#ifdef DEBUG_PRINTOUT
  ZZX dbgPoly, skPoly, tmp1, tmp2;
  if (dbgKey && dbgEa) { // apply the decryption procedure to the zzParts
    dbgKey->sKeys[bootstrapKeyID].toPoly(skPoly);
    MulMod(dbgPoly, skPoly, zzParts[1], context.zMStar.getPhimX());
    dbgPoly += zzParts[0];

    Vec<ZZ> powerful;
    context.p2dConversion->ZZXtoPowerful(powerful,dbgPoly);

    Vec<long> expansion(INIT_SIZE, e+r); // base-p expansion of coefs
    for (long i=0; i<powerful.length(); i++) {
      long c_lo=0, c_hi=0;
      baseRep(expansion, e+r, powerful[i], p);
      for (long j=0; j<r; j++) {
	c_lo = (c_lo*p) + expansion[j];
	c_hi = (c_hi*p) + expansion[e+j];
      }
      if (p==2) c_hi += expansion[e-1];
      long delta = (c_lo - c_hi - ptxt_pwr[i]) % p2r;
      if (delta != 0) {
	cerr << " @error in pwrful(c1*s+c0)["<<i<<"]: "
	     << powerful[i] << "="<<expansion<<" (hi: "<<c_hi<<", lo: "<<c_lo
	     << ") does not yeild "<< ptxt_pwr[i]<<endl;
	break;
      }
    }
  }
#endif

  // Add multiples of p2r and q to make the zzParts divisible by p^{e'}
  long maxU=0;
  for (long i=0; i<(long)zzParts.size(); i++) {
    // make divisible by p^{e'}
    long newMax = makeDivisible(zzParts[i].rep, p2ePrime, p2r, q);
    zzParts[i].normalize();   // normalize after working directly on the rep
    if (maxU < newMax)  maxU = newMax;
  }

  // Check that the estimated noise is still low
  if (noise + maxU*p2r*(skHwts[bootstrapKeyID]+1) > q/2) 
    cerr << " * noise/q after makeDisivible = "
	 << ((noise + maxU*p2r*(skHwts[bootstrapKeyID]+1))/q) << endl;

#ifdef DEBUG_PRINTOUT
  if (dbgKey && dbgEa) { // apply the decryption procedure to the new zzParts
    MulMod(dbgPoly, skPoly, zzParts[1], context.zMStar.getPhimX());
    dbgPoly += zzParts[0];

    Vec<ZZ> powerful;
    context.p2dConversion->ZZXtoPowerful(powerful,dbgPoly);
    Vec<long> expansion(INIT_SIZE, e+r); // base-p expansion of coefs
    for (long i=0; i<powerful.length(); i++) {
      long c_lo=0, c_hi=0;
      baseRep(expansion, e+r, powerful[i], p);
      for (long j=0; j<r; j++) {
	c_lo = (c_lo*p) + expansion[j];
	c_hi = (c_hi*p) + expansion[e+j];
      }
      if (p==2) c_hi += expansion[e-1];
      long delta = (c_lo - c_hi - ptxt_pwr[i]) % p2r;
      if (delta != 0) {
	cerr << " @error in makeDivisible(pwrful(c1*s+c0))["<<i<<"]: "
	     << powerful[i] << "="<<expansion<<" (hi: "<<c_hi<<", lo: "<<c_lo
	     << ") does not yeild "<< ptxt_pwr[i]<<endl;
	break;
      }
    }
  }
#endif

  for (long i=0; i<(long)zzParts.size(); i++)
    zzParts[i] /= p2ePrime;   // divide by p^{e'}

  double p0size = to_double(coeffsL2Norm(zzParts[0]));
  double p1size = to_double(coeffsL2Norm(zzParts[1]));

  // Multiply the post-processed cipehrtext by the encrypted sKey
  ctxt = bootstrapEkey;
  e -= ePrime;
  p2e /= p2ePrime;     // reduce the plaintext space by p^{e'} factor
  ctxt.reducePtxtSpace(/*newPtxtSpace=*/p2e*p2r);

  cerr << "+ before recryption, level=" << ctxt.findBaseLevel()<<endl;
  cerr << "  constant size = " << p0size <<", "<<p1size << endl;
  ctxt.multByConstant(zzParts[1], p1size*p1size);
  ctxt.addConstant(zzParts[0], p0size*p0size);

  FHE_NTIMER_STOP(preProcess);
  cerr << "+ before linearTrans1, level=" << ctxt.findBaseLevel()<<endl;

#ifdef DEBUG_PRINTOUT
  CheckCtxt(ctxt, "After pre-processing");
  if (dbgKey && dbgEa) { // apply the decryption procedure to the new zzParts
    MulMod(dbgPoly, skPoly, zzParts[1], context.zMStar.getPhimX());
    dbgPoly += zzParts[0];

    Vec<ZZ> powerful;
    context.p2dConversion->ZZXtoPowerful(powerful,dbgPoly);

    Vec<long> expansion(INIT_SIZE, e+r); // base-p expansion of coefs
    for (long i=0; i<powerful.length(); i++) {
      long c_lo=0, c_hi=0;
      baseRep(expansion, e+r, powerful[i], p);
      for (long j=0; j<r; j++) {
	if (j>=ePrime) c_lo = (c_lo*p) + expansion[j-ePrime];
	c_hi = (c_hi*p) + expansion[e+j];
      }
      if (p==2) c_hi += expansion[e-1];
      long delta = (c_lo - c_hi - ptxt_pwr[i]) % p2r;
      if (delta != 0) {
	cerr << " @error pwrfl["<<i<<"]="
	     << powerful[i] << "="<<expansion<<" (hi: "<<c_hi<<", lo: "<<c_lo
	     << ") does not yeild "<< ptxt_pwr[i]<<endl;
	break;
      }
    }
  }
#endif

  // Move the powerful-basis coefficients to the plaintext slots
  FHE_NTIMER_START(LinearTransform1);
  ctxt.getContext().firstMap->apply(ctxt);
  FHE_NTIMER_STOP(LinearTransform1);
  cerr << "+ after linearTrans1, level=" << ctxt.findBaseLevel()<<endl;

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)

  // On 1st execution, initialize the linearized polynomials
  if (unpackSlotEncoding.size()==0)
    initUnpackEncoding(unpackSlotEncoding, *context.bootstrapEA);

  FHE_NTIMER_START(DigitExtraction);
  extractDigitsPacked(ctxt, e, r, ePrime, unpackSlotEncoding);
  FHE_NTIMER_STOP(DigitExtraction);

  // Move the slots back to powerful-basis coefficients
  cerr << "+ before linearTrans2, level=" << ctxt.findBaseLevel()<<endl;
  FHE_NTIMER_START(LinearTransform2);
  ctxt.getContext().secondMap->apply(ctxt);
  FHE_NTIMER_STOP(LinearTransform2);
}

/*********************************************************************/
/*********************************************************************/

// Return in poly a polynomial with X^i encoded in all the slots
static void x2iInSlots(ZZX& poly, long i,
		       vector<ZZX>& xVec, const EncryptedArray& ea)
{
  xVec.resize(ea.size());
  ZZX x2i = ZZX(i,1);
  for (long j=0; j<(long)xVec.size(); j++) xVec[j] = x2i;
  ea.encode(poly, xVec);
}

// Initialized linearized polynomials for unpacking the slots
// we should just use a 1-D array, but for now I'll use
// a 2-D array so as to not break the interfaces
//


static void initUnpackEncoding(vector< vector<ZZX> >& unpackSlotEncoding,
			       const EncryptedArray& ea)
{
  FHE_TIMER_START;

  zz_pBak bak; bak.save(); ea.getAlMod().restoreContext();

  long nslots = ea.size();
  long d = ea.getDegree();

  const Mat<zz_p>& CBi = ea.getDerived(PA_zz_p()).getNormalBasisMatrixInverse();

  vector<ZZX> LM;
  LM.resize(d);
  for (long i = 0; i < d; i++)
    LM[i] = rep(CBi[i][0]);

  vector<ZZX> C; 
  ea.buildLinPolyCoeffs(C, LM);
  

  unpackSlotEncoding.resize(1);
  unpackSlotEncoding[0].resize(d);

  for (long j = 0; j < d; j++) {
    vector<ZZX> v(nslots);
    for (long k = 0; k < nslots; k++) v[k] = C[j];
    ea.encode(unpackSlotEncoding[0][j], v);
  }
}
  

// Make every entry of vec divisible by p^e by adding/subtracting multiples
// of p^r and q, while keeping the added multiples small. Specifically, for
// e>=2, an integer z can be made divisible by p^e via
//   z' = z + u * p^r + v*p^r * q,
// with
//   |u|<=ceil(alpha p^{r(e-1)}/2) and |v|<=0.5+floor(beta p^{r(e-1)}/2),
// for any alpha+beta=1. We assume that r<e and that q-1 is divisible by p^e.
// Returns the largest absolute values of the u's and the new entries.
//
// This code is more general than we need, for bootstrapping we will always
// have e>r.
template<class VecInt>
long makeDivisible(VecInt& vec, long p2e, long p2r, long q, double alpha)
{
  assert(((p2e % p2r == 0) && (q % p2e == 1)) ||
	 ((p2r % p2e == 0) && (q % p2r == 1)));

  long maxU =0;
  ZZ maxZ;
  for (long i=0; i<vec.length(); i++) {
    ZZ z, z2; conv(z, vec[i]);
    long u=0, v=0;

    long zMod1=0, zMod2=0;
    if (p2r < p2e) {
      zMod1 = rem(z,p2r);
      if (zMod1 > p2r/2) zMod1 -= p2r; // map to the symmetric interval

      // make z divisible by p^r by adding a multiple of q
      z2 = z - to_ZZ(zMod1)*q;
      zMod2 = rem(z2,p2e); // z mod p^e, still divisible by p^r
      if (zMod2 > p2e/2) zMod2 -= p2e; // map to the symmetric interval
      zMod2 /= -p2r; // now z+ p^r*zMod2=0 (mod p^e) and |zMod2|<=p^{r(e-1)}/2

      u = ceil(alpha * zMod2);
      v = zMod2 - u; // = floor((1-alpha) * zMod2)
      z = z2 + u*p2r + to_ZZ(q)*v*p2r;
    }
    else { // r >= e, use only mulitples of q
      zMod1 = rem(z,p2e);
      if (zMod1 > p2e/2) zMod1 -= p2e; // map to the symmetric interval
      z -= to_ZZ(zMod1) * q;
    }
    if (abs(u) > maxU) maxU = abs(u);
    if (abs(z) > maxZ) maxZ = abs(z);

    if (rem(z,p2e) != 0) { // sanity check
      cerr << "**error: original z["<<i<<"]=" << vec[i]
	   << std::dec << ", p^r="<<p2r << ", p^e="<<p2e << endl;
      cerr << "z' = z - "<<zMod1<<"*q = "<< z2<<endl;
      cerr << "z''=z' +" <<u<<"*p^r +"<<v<<"*p^r*q = "<<z<<endl;
      exit(1);
    }
    conv(vec[i], z); // convert back to native format
  }
  return maxU;
}
// explicit instantiation for vec_ZZ and vec_long
template long makeDivisible<vec_ZZ>(vec_ZZ& v, long p2e,
				    long p2r, long q, double alpha);
// template long makeDivisible<vec_long>(vec_long& v, long p2e,
// 				      long p2r, long q, double alpha);


// Extract digits from fully packed slots
void extractDigitsPacked(Ctxt& ctxt, long botHigh, long r, long ePrime,
			 const vector< vector<ZZX> >& unpackSlotEncoding)
{
#ifdef DEBUG_PRINTOUT
  CheckCtxt(ctxt, "Before unpacking");
  double tm = -GetTime();
#endif


  // Step 1: unpack the slots of ctxt
  FHE_NTIMER_START(unpack);

  ctxt.cleanUp();

  // Apply the d automorphisms and store them in scratch area
  long d = ctxt.getContext().zMStar.getOrdP();


  vector<Ctxt> scratch; // used below 
  vector<Ctxt> unpacked(d, Ctxt(ZeroCtxtLike, ctxt));

  {
    // explicit scope to force all temporaries
    // to be released

    FHE_NTIMER_START(UnpackConvertConstants);

    vector< shared_ptr<DoubleCRT> > coeff_vector;
    coeff_vector.resize(d);
    for (long i = 0; i < d; i++)
      coeff_vector[i] = shared_ptr<DoubleCRT>(new 
        DoubleCRT(unpackSlotEncoding[0][i], ctxt.getContext(), ctxt.getPrimeSet()) );

    FHE_NTIMER_STOP(UnpackConvertConstants);


    FHE_NTIMER_START(UnpackMain);
    Ctxt tmp1(ZeroCtxtLike, ctxt);
    Ctxt tmp2(ZeroCtxtLike, ctxt);

    for (long j = 0; j < d; j++) { // process jth Frobenius 
      tmp1 = ctxt;
      tmp1.frobeniusAutomorph(j);
      tmp1.cleanUp();

      for (long i = 0; i < d; i++) {
        tmp2 = tmp1;
        tmp2.multByConstant(*coeff_vector[mcMod(i+j, d)]);
        unpacked[i] += tmp2;
      }
    }

    FHE_NTIMER_STOP(UnpackMain);
  }


  FHE_NTIMER_STOP(unpack);
  cerr << "+ after unpacking, level=" << unpacked[0].findBaseLevel()<<endl;

  // Step 2: extract the digits top-1,...,0 from the slots of unpacked[i]

  long p = ctxt.getContext().zMStar.getP();
  long p2r = power_long(p,r);
  long topHigh = botHigh + r-1;
  if (p==2 && r>2)
    topHigh--; // For p==2 we sometime get a bit for free

#ifdef DEBUG_PRINTOUT
  tm += GetTime();
  cerr << "  Unpacking in "<<tm<<" seconds.";
  CheckCtxt(unpacked[0], "After unpacking");
  tm = -GetTime();
  bool foundError = false;
  cerr << "  extracting "<<(topHigh+1)<<" digits:\n";
#endif

  for (long i=0; i<(long)unpacked.size(); i++) {
    FHE_NTIMER_START(extractDigits);
    if (topHigh<=0) { // extracting LSB = no-op
      scratch.assign(1,unpacked[i]);
    } else {          // extract digits topHigh...0, store them in scratch
      extractDigits(scratch, unpacked[i], topHigh+1);

#ifdef DEBUG_PRINTOUT      // check digit-extraction
      long nDigits = scratch.size();
      if (dbgKey && dbgEa && !foundError) {
	vector<ZZX> slots;
	dbgEa->decrypt(unpacked[i], *dbgKey, slots);

	vector< vector<long> > digits(nDigits);
	for (long j=0; j<nDigits; j++)
	  dbgEa->decrypt(scratch[j], *dbgKey, digits[j]);

	Vec<long> expansion;
	for (long k=0; k<(long)slots.size() && !foundError; k++) {
	  baseRep(expansion, nDigits, coeff(slots[k],0), p);
	  long p2j = 1;
	  for (long j=nDigits-1; j>=0; j--) {
	    p2j *= p;
	    if ( ((digits[j][k]-expansion[j]) % p2j) != 0) {
	      cerr << " @ Digit extraction error (p="<<p<<"):\n";
	      cerr << "   slots["<<k<<"]="<<slots[k]<<", rep="<<expansion<<endl;
	      cerr << "   but extract returns [";
	      for (long t=0; t<nDigits-1; t++) cerr << digits[t][k]<< " ";
	      cerr << digits[nDigits-1][k]<< "]\n";
	      cerr << "   "<<j<<"'th digits mismatch\n";
	      foundError = true;
	      break;
	    }
	  }
	}
      }
#endif
    }
    FHE_NTIMER_STOP(extractDigits);

    // set upacked[i] = -\sum_{j=botHigh}^{topHigh} scratch[j] * p^{j-botHigh}
    FHE_NTIMER_START(combiningDigits);
    if (topHigh >= (long)scratch.size()) {
      topHigh = scratch.size() -1;
      cerr << " @ suspect: not enough digits in extractDigitsPacked\n";
    }
    unpacked[i] = scratch[topHigh];

    for (long j=topHigh-1; j>=botHigh; --j) {
      unpacked[i].multByP();
      unpacked[i] += scratch[j];
    }
    if (p==2 && botHigh>0)   // For p==2, subtract also the previous bit
      unpacked[i] += scratch[botHigh-1];
    unpacked[i].negate();

    if (r>ePrime) {          // Add in digits from the bottom part, if any
      long topLow = r-1 - ePrime;
      Ctxt tmp = scratch[topLow];
      for (long j=topLow-1; j>=0; --j) {
	tmp.multByP();
	tmp += scratch[j];
      }
      if (ePrime>0)
	tmp.multByP(ePrime); // multiply by p^e'
      unpacked[i] += tmp;
    }

    // Our plaintext space is now mod p^r
    unpacked[i].reducePtxtSpace(p2r);
    FHE_NTIMER_STOP(combiningDigits);
  }
#ifdef DEBUG_PRINTOUT
  tm += GetTime();
  cerr << "  Extracting in "<<tm<<" seconds.";
  CheckCtxt(unpacked[0], "After extraction");
  tm = -GetTime();
#endif
  cerr << "+ before repacking, level=" << unpacked[0].findBaseLevel()<<endl;

  // Step 3: re-pack the slots
  FHE_NTIMER_START(repack);

  const EncryptedArray& ea2 = *ctxt.getContext().secondEA;
  ZZX xInSlots;
  vector<ZZX> xVec(ea2.size());
  ctxt = unpacked[0];
  for (long i=1; i<d; i++) {
    x2iInSlots(xInSlots, i, xVec, ea2);
    unpacked[i].multByConstant(xInSlots);
    ctxt += unpacked[i];
  }
  FHE_NTIMER_STOP(repack);
#ifdef DEBUG_PRINTOUT
  tm += GetTime();
  cerr << "  Repacking in "<<tm<<" seconds.";
  CheckCtxt(ctxt, "After re-packing");
#endif
}
