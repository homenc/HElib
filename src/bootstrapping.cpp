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
#include "EvalMap.h"
#include "powerful.h"

#define DEBUGGING

/*********** Debugging utilities **************/
FHESecKey* dbgKey=NULL;
EncryptedArray* dbgEa=NULL;
ZZX dbg_ptxt;
Vec<ZZ> ptxt_pwr;

#define FLAG_PRINT_ZZX  1
#define FLAG_PRINT_POLY 2
#define FLAG_PRINT_VEC  4

template<class T> ostream& printVec(ostream& s, const Vec<T>& v, long nCoeffs=40);
ostream& printZZX(ostream& s, const ZZX& poly, long nCoeffs=40);
void decryptAndPrint(const Ctxt& ctxt, const FHESecKey& sk,
		     const EncryptedArray& ea, long flags=0);
/********* End Debugging utilities **************/


// Initialized linearized polynomials for unpacking the slots
void initUnpackEncoding(vector< vector<ZZX> >& unpackSlotEncoding,
			const EncryptedArray& ea)
{
  long nslots = ea.size();
  long d = ea.getDegree();
  unpackSlotEncoding.resize(d);
  for (long i=0; i<d; i++) {
    // Specify transforms to extracts i'th coeff.: LM[i]=1 & LM[j]=0 for j!=i
    vector<ZZX> LM(d, ZZX::zero());
    conv(LM[i],1); //  FIXME: This may be the wrong transformation

    vector<ZZX> C; // "build" intermediate format for this transformation
    ea.buildLinPolyCoeffs(C, LM);

    // "encode" final format for this transformation
    unpackSlotEncoding[i].resize(d);
    for (long j = 0; j < d; j++) {
      vector<ZZX> v(nslots);
      for (long k = 0; k < nslots; k++) v[k] = C[j];
      ea.encode(unpackSlotEncoding[i][j], v);
    }
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
// have e>r/
template<class VecInt>
long makeDivisible(VecInt& vec, long p2e, long p2r, long q, double alpha=0.4)
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
  // if ((FHEcontext::bootstrapHwt+1)*maxZ/q > q/2) {
  //   double numer = to_double((FHEcontext::bootstrapHwt+1)*maxZ);
  //   cerr << "*make-divisible overflow: max-entry/q^2 = "
  // 	 << (numer/(q*(double)q)) << endl;
  // }
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
  long d = unpackSlotEncoding.size();
  CheckCtxt(ctxt, "Before unpacking");

  // Step 1: unpack the slots of ctxt

  // Apply the d automorphisms and store them in scratch area
  vector<Ctxt> scratch(d, ctxt);
  for (long i=1; i<d; i++) scratch[i].frobeniusAutomorph(i);

  // Now apply the d linear transformations to get the coefficients
  vector<Ctxt> unpacked(d, scratch[0]);
  for (long i=0; i<d; i++) { // compute the i'th linear transformation
    unpacked[i].multByConstant(unpackSlotEncoding[i][0]);
    for (long j = 1; j < d; j++) {
      Ctxt tmp(scratch[j]);
      tmp.multByConstant(unpackSlotEncoding[i][j]);
      unpacked[i] += tmp;
    }
  }
  CheckCtxt(unpacked[0], "After unpacking");

  // Step 2: extract the digits top-1,...,0 from the slots of unpacked[i]

  long p = ctxt.getContext().zMStar.getP();
  long p2r = power_long(p,r);
  long topHigh = botHigh + r-1;
  if (p==2 && r>2)
    topHigh--; // For p==2 we sometime get a bit for free

  cerr << "extracting "<<(topHigh+1)<<" digits\n";
  for (long i=0; i<(long)unpacked.size(); i++) {
    if (topHigh>0) // extract digits topHigh...0, store them in scratch
      extractDigits(scratch, unpacked[i], topHigh+1);
    else
      scratch.assign(1,unpacked[i]); // extracting LSB = no-op

    // set upacked[i] = -\sum_{j=botHigh}^{topHigh} scratch[j] * p^{j-botHigh}
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
  }
  CheckCtxt(unpacked[0], "After extraction");

  // Step 3: re-pack the slots
  EncryptedArray ea2(ctxt.getContext(), ctxt.getContext().alMod);
  long nSlots = ea2.size();
  vector<ZZX> xVec(nSlots, ZZX(1,1)); // X in all the slots
  ZZX xInSlots;
  ea2.encode(xInSlots, xVec);         // encode as ZZX
  ZZX x2iInSlots = xInSlots;          // X^i := X

  ctxt = unpacked[0];
  for (long i=1; i<d; i++) {
    unpacked[i].multByConstant(x2iInSlots);
    MulMod(x2iInSlots, x2iInSlots,xInSlots,ctxt.getContext().zMStar.getPhimX());
    PolyRed(x2iInSlots, conv<ZZ>(p2r));
    ctxt += unpacked[i];
  }
  CheckCtxt(ctxt, "After re-packing");
}
 
// Move the powerful-basis coefficients to the plaintext slots
void movePwrflCoefs2Slots(Ctxt& ctxt)
{
  ZZX ptxt1; // for debugging
  long p = ctxt.getPtxtSpace();
  if (dbgKey && dbgEa) dbgKey->Decrypt(ptxt1, ctxt);

  const FHEcontext& context = ctxt.getContext();
  context.firstMap->apply(ctxt);

  // For debugging: Check that we have powerful representation in the slots
  if (dbgKey && dbgEa) {
    zz_pBak bak; bak.save(); // backup NTL's current modulus
    zz_p::init(p);
    PowerfulConversion pConv(context.p2dConversion->getIndexTranslation());
    HyperCube<zz_p> powerful(pConv.getShortSig());
    zz_pX poly = conv<zz_pX>(ptxt1);
    pConv.polyToPowerful(powerful, poly);

    long nSlots = dbgEa->size();
    long d =  dbgEa->getDegree();
    vector<ZZX> v(nSlots);
    dbgEa->decrypt(ctxt,*dbgKey,v);

    long idx = 0;
    for (long i=0; i<nSlots; i++) for (long j=0; j<d; j++) {
	long coef1 = conv<long>(v[i][j]);
	long coef2 = conv<long>(powerful[idx++]);
	if (coef1 != coef2) {
	  cerr << " @ error: powerful coefficient is "<<coef2
	       << " but slot contains "<<coef1<<endl;
	  return;
	}
      }
    cerr << "LinTrans1 successful. ";
    decryptAndPrint(ctxt, *dbgKey, *dbgEa, 0);
  }
}

// Move the slots back to powerful-basis coefficients
void moveSlots2PwrflCoefs(Ctxt& ctxt)
{
  const FHEcontext& context = ctxt.getContext();
  context.secondMap->apply(ctxt);
  return;

  // For testing purposes: Decrypt, move slots to coefficients, then re-encryp
  if (dbgKey && dbgEa) {
    vector<ZZX> v;
    dbgEa->decrypt(ctxt, *dbgKey, v);

    ZZX ptxt; 
    long d =  dbgEa->getDegree();
    // copy the coefficients from ptxt to the slots
    for (long i=v.size()-1; i>=0; --i) for (long j=d-1; j>=0; --j) {
	const ZZ& coef = coeff(v[i], j);
	SetCoeff(ptxt, i*d + j, coef);
      }
    dbgKey->Encrypt(ctxt, ptxt, ctxt.getContext().alMod.getPPowR());
  }
}

// bootstrap a ciphertext to reduce noise
void FHEPubKey::reCrypt(Ctxt &ctxt)
{
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
  cerr << "reCrypt: p="<<p<<", r="<<r<<", e="<<e<<" ePrime="<<ePrime
       << ", q="<<q<<endl;

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

#ifdef DEBUGGING
  if (dbgKey && dbgEa) {
    ZZX poly;
    dbgKey->Decrypt(poly,ctxt);
    if (poly!=dbg_ptxt)
      cerr << " Decryption error after key-switching to bootsrapping key\n";
    else {
      cerr << " After key-switching to bootsrapping key: ";
      decryptAndPrint(ctxt, *dbgKey, *dbgEa, 0);
    }
  }
#endif

  // "raw mod-switch" to the bootstrapping mosulus q=p^e+1.
  vector<ZZX> zzParts; // the mod-switched parts, in ZZX format
  double noise = ctxt.rawModSwitch(zzParts, q);
  noise = sqrt(noise);

#ifdef DEBUGGING
  ZZX dbgPoly, skPoly, tmp1, tmp2;
  if (dbgKey && dbgEa) { // apply the decryption procedure to the zzParts
    dbgKey->sKeys[bootstrapKeyID].toPoly(skPoly);
    MulMod(dbgPoly, skPoly, zzParts[1], context.zMStar.getPhimX());
    dbgPoly += zzParts[0];

    Vec<ZZ> powerful;
    context.p2dConversion->ZZXtoPowerful(powerful,dbgPoly);
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

#ifdef DEBUGGING
  if (dbgKey && dbgEa) { // apply the decryption procedure to the new zzParts
    MulMod(dbgPoly, skPoly, zzParts[1], context.zMStar.getPhimX());
    dbgPoly += zzParts[0];

    Vec<ZZ> powerful;
    context.p2dConversion->ZZXtoPowerful(powerful,dbgPoly);
  }
#endif

  //  if (p==2) // add the constant 2^{e-1} * allOnes to the ciphertext
  //    zzParts[0] += (p2e/p) * context.allOnes;

  for (long i=0; i<(long)zzParts.size(); i++)
    zzParts[i] /= p2ePrime;   // divide by p^{e'}

  // Multiply the post-processed cipehrtext by the encrypted sKey
  ctxt = bootstrapEkey;
  ctxt.multByConstant(zzParts[1]);
  ctxt.addConstant(zzParts[0]);

  e -= ePrime;
  p2e /= p2ePrime; // reduce the plaintext space by p^{e'} factor
  ctxt.reducePtxtSpace(p2e*p2r);
  FHE_NTIMER_STOP(preProcess);

  // Move the powerful-basis coefficients to the plaintext slots
  FHE_NTIMER_START(LinearTransform1);
  movePwrflCoefs2Slots(ctxt);
  FHE_NTIMER_STOP(LinearTransform1);

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)

  // On 1st execution, initialize the linearized polynomials
  if (unpackSlotEncoding.size()==0)
    initUnpackEncoding(unpackSlotEncoding, *context.bootstrapEA);

  FHE_NTIMER_START(DigitExtraction);
  extractDigitsPacked(ctxt, e, r, ePrime, unpackSlotEncoding);
  FHE_NTIMER_STOP(DigitExtraction);

  // Move the slots back to powerful-basis coefficients
  FHE_NTIMER_START(LinearTransform2);
  moveSlots2PwrflCoefs(ctxt);
  FHE_NTIMER_STOP(LinearTransform2);
}


/********************************************************************
 ********************************************************************/
void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'p=2 e=4 q=257'\n\n";
  cerr << "  p is the plaintext base [default=2]\n";
  cerr << "  e is the exponent [default=4]\n";
  cerr << "  q is the modulus [default=p^e+1]\n";
  exit(0);
}

int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["L"] = "20";
  argmap["p"] = "2";
  argmap["m1"] = "5";
  argmap["m2"] = "7";
  argmap["m3"] = "0";
  argmap["r"] = "1";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long p = atoi(argmap["p"]);
  long m1 = atoi(argmap["m1"]);
  long m2 = atoi(argmap["m2"]);
  long m3 = atoi(argmap["m3"]);
  long r = atoi(argmap["r"]);
  long L =  atoi(argmap["L"]);

  if (m1 <= 1) {
    cerr << "m needs at least one factor\n";
    exit(0);
  }
  Vec<long> mvec(INIT_SIZE, 1);
  long m = mvec[0] = m1;
  if (m2 > 1) {
    m *= m2;
    append(mvec, m2);
  }
  if (m3 > 1) {
    m *= m3;
    append(mvec, m3);
  }

  cerr << "*** TestIt: p=" << p
       << ", r=" << r
       << ", L=" << L
       << ", m=" << m
       << " (=" << mvec << ")"
       << endl;

  setTimersOn();
  FHE_NTIMER_START(initialize);

  cerr << "Initializing context..." << std::flush;
  FHEcontext context(m,p,r);
  buildModChain(context, L, /*c=*/3);
  cerr << " done. Computing bootstrapping information..." << std::flush;
  context.makeBootstrappable(mvec);
  cerr << " done\n";
  long p2r = context.alMod.getPPowR();

  cerr << "Generating keys..." << std::flush;
  FHESecKey secretKey(context);
  FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(64);      // A Hamming-weight-64 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey);
  secretKey.genBootstrapData();
  cerr << " done\n";

  // the bootstrapping key is encrypted relative to plaintext space p^{e+r}.
  //  long ePr = FHEPubKey::ePlusR(p);

  //  PAlgebraMod almod2(context.zMStar, context.bootstrapPAM->getR());
  //  EncryptedArray ea2(context, almod2);
  //  EncryptedArray ea1(context, context.alMod);

  FHE_NTIMER_STOP(initialize);
  //  cerr << "****Initialization time:\n";
  //  printAllTimers();
  resetAllTimers();

  dbgKey = &secretKey; // debugging key and ea

  // The bootstrapping data in the context is for plaintext space p^{e+r-e'}
  dbgEa = context.bootstrapEA;

  zz_p::init(p2r);
  zz_pX poly_p = random_zz_pX(context.zMStar.getPhiM());
  PowerfulConversion pConv(context.p2dConversion->getIndexTranslation());
  HyperCube<zz_p> powerful(pConv.getShortSig());
  pConv.polyToPowerful(powerful, poly_p);
  
  conv(dbg_ptxt, poly_p);
  PolyRed(dbg_ptxt, p2r, true);
  context.p2dConversion->ZZXtoPowerful(ptxt_pwr, dbg_ptxt);

  ZZX poly2;
  Ctxt c1(publicKey);
  cerr << "encrypting (in powerful rep) "; // <<dbg_ptxt<<endl;
  printVec(cerr, powerful.getData())<<endl;

  secretKey.Encrypt(c1,dbg_ptxt,p2r);
  FHE_NTIMER_START(reCrypt);
  publicKey.reCrypt(c1);
  FHE_NTIMER_STOP(reCrypt);
  secretKey.Decrypt(poly2,c1);
  decryptAndPrint(c1, secretKey, *dbgEa, 0);

  if (dbg_ptxt != poly2) {
    conv(poly_p,poly2);
    HyperCube<zz_p> powerful2(pConv.getShortSig());
    cerr << "\ndecryption error, result is "; //<<poly2<<endl;
    pConv.polyToPowerful(powerful2, poly_p);
    printVec(cerr, powerful2.getData())<<endl;
    long numDiff = 0;
    for (long i=0; i<powerful.getSize(); i++) 
      if (powerful[i] != powerful2[i]) {
        numDiff++;
	cerr << i << ": " << powerful[i] << " != " << powerful2[i]<<", ";
	if (numDiff >5) {
	  cerr << endl;
	  break;
	}
      }
  }
  else cerr << "  decryption succeeds!!\n";

  //  cerr << "\n****Bootstrapping time:\n";
  //  printAllTimers();
  return 0;
}




/**********************************************/
/*********** Debugging utilities **************/
/**********************************************/

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

ostream& printZZX(ostream& s, const ZZX& poly, long nCoeffs)
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

void decryptAndPrint(const Ctxt& ctxt, const FHESecKey& sk,
		     const EncryptedArray& ea, long flags)
{
  const FHEcontext& context = ctxt.getContext();
  xdouble noiseEst = sqrt(ctxt.getNoiseVar());
  xdouble modulus = xexp(context.logOfProduct(ctxt.getPrimeSet()));
  vector<ZZX> ptxt;
  ZZX p, pp;
  sk.Decrypt(p, ctxt, pp);

  cerr << "plaintext space mod "<<ctxt.getPtxtSpace()
       << ", level="<<ctxt.findBaseLevel()
       << ", \n           |noise|=q*" << (coeffsL2Norm(pp)/modulus)
       << ", |noiseEst|=q*" << (noiseEst/modulus)
       <<endl;

  if (flags & FLAG_PRINT_ZZX) {
    cerr << "   before mod-p reduction=";
    printZZX(cerr,pp) <<endl;
  }
  if (flags & FLAG_PRINT_POLY) {
    cerr << "   after mod-p reduction=";
    printZZX(cerr,p) <<endl;
  }
  if (flags & FLAG_PRINT_VEC) {
    ea.decode(ptxt, p);
    if (ea.getAlMod().getTag() == PA_zz_p_tag
	&& ctxt.getPtxtSpace() != ea.getAlMod().getPPowR()) {
      long g = GCD(ctxt.getPtxtSpace(), ea.getAlMod().getPPowR());
      for (long i=0; i<ea.size(); i++)
	PolyRed(ptxt[i], g, true);
    }
    cerr << "   decoded to ";
    if (deg(p) < 40) // just pring the whole thing
      cerr << ptxt << endl;
    else if (ptxt.size()==1) // a single slot
      printZZX(cerr, ptxt[0]) <<endl;
    else { // print first and last slots
      printZZX(cerr, ptxt[0],20) << "--";
      printZZX(cerr, ptxt[ptxt.size()-1], 20) <<endl;      
    }
  }
}
