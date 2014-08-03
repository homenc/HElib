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


/*********** Debugging utilities **************/
FHESecKey* dbgKey=NULL;
EncryptedArray* dbgEa=NULL;

#define FLAG_PRINT_ZZX  1
#define FLAG_PRINT_POLY 2
#define FLAG_PRINT_VEC  4

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
    unpacked[i].negate();
    if (r>ePrime) { // Add in digits from the bottom part, if any
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
  // For testing purposes.
  // Decrypt a polynomial, put the coefficients in slots, then re-encryp
  if (dbgKey && dbgEa) {
    ZZX ptxt;
    dbgKey->Decrypt(ptxt, ctxt);

    long nSlots = dbgEa->size();
    long d =  dbgEa->getDegree();
    vector<ZZX> v(nSlots);
    // copy the coefficients from ptxt to the slots
    for (long i=nSlots-1; i>=0; --i) for (long j=d-1; j>=0; --j) {
	const ZZ& coef = coeff(ptxt, i*d + j);
	SetCoeff(v[i], j, coef);
      }
    dbgEa->skEncrypt(ctxt, *dbgKey, v);
  }
}

// Move the slots back to powerful-basis coefficients
void moveSlots2PwrflCoefs(Ctxt& ctxt)
{
  // For testing purposes.
  // Decrypt a polynomial, put the coefficients in slots, then re-encryp
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
  ZZX dbgPoly, skPoly, tmp1, tmp2;

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

  // "raw mod-switch" to the bootstrapping mosulus q=p^e+1.
  vector<ZZX> zzParts; // the mod-switched parts, in ZZX format
  double noise = ctxt.rawModSwitch(zzParts, q);
  noise = sqrt(noise);

  // if (noise>q/64.0)
  //   cerr << "  @ initial noise exceeds q/64, decryption error likely\n";

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

  // if (dbgKey) {
  //   MulMod(dbgPoly, skPoly, zzParts[1], context.zMStar.getPhimX());
  //   dbgPoly += zzParts[0];
  //   tmp1 = tmp2 = dbgPoly;
  //   if (p==2) tmp1 += ((p2e/p) * context.allOnes);
  //   PolyRed(dbgPoly, conv<ZZ>(q));
  //   PolyRed(dbgPoly, conv<ZZ>(p2r), true);
  //   ZZ maxCoef = ZZ::zero();
  //   ZZ maxNoise = ZZ::zero();
  //   for (long i=0; i<=deg(tmp2); i++) {
  //     ZZ coef2 = coeff(tmp2,i);
  //     ZZ coef1 = coeff(tmp1,i);
  //     ZZ cDiv2 = coef2 / p2e;
  //     ZZ cDiv1 = coef1 / p2e;
  //     ZZ cModQ2; rem(cModQ2,coef2,to_ZZ(q));
  //     if (cModQ2 > q/2) cModQ2 -= q;
  //     if (abs(coef2)>maxCoef)   maxCoef = abs(coef2);
  //     if (abs(cModQ2)>maxNoise) maxNoise = abs(cModQ2);
  //     if ((coef1 - cDiv1 - dbgPoly[i])% p2r != 0) {
  // 	cerr << "simple decryption formula failed\n";
  // 	cerr << " * original coef_"<<i<<"="<< coef2
  // 	     << "(~q^2*"<<(to_double(coef2)/(q*(double)q))<<")"
  // 	     << ", [orig-coef]_q="<<cModQ2
  // 	     << "(~q*"<<(to_double(cModQ2)/(double)q)<<")"
  // 	     << ", orig-coef/p^e="<<cDiv2<<endl;
  // 	cerr << " * coef="<< coef1
  // 	     << ", coef/p^e="<<cDiv1
  // 	     << " but [[orig-coef]_q]_"<<p2r<<"="<< coeff(dbgPoly,i)<<endl;
  // 	exit(0);
  //     }
  //     if (coef1 % p2r != 0)
  // 	cerr << "coef_"<<i<<"="<<coef1<<"!=0 (mod "<<p2r<<endl;
  //   } // end-for
  // }

  if (p==2) // add the constant 2^{e-1} * allOnes to the ciphertext
    zzParts[0] += (p2e/p) * context.allOnes;

  for (long i=0; i<(long)zzParts.size(); i++)
    zzParts[i] /= p2ePrime;   // divide by p^{e'}

  // Multiply the post-processed cipehrtext by the encrypted sKey
  ctxt = bootstrapEkey;
  ctxt.multByConstant(zzParts[1]);
  ctxt.addConstant(zzParts[0]);

  e -= ePrime;
  p2e /= p2ePrime; // reduce the plaintext space by p^{e'} factor
  ctxt.reducePtxtSpace(p2e*p2r);

  // if (dbgKey) {
  //   tmp1 /= p2ePrime;
  //   PolyRed(tmp1, conv<ZZ>(ctxt.getPtxtSpace()), true);
  //   dbgKey->Decrypt(tmp2, ctxt);
  //   if (tmp2 != tmp1) {
  //     cerr << "  @ decryption error afer mult-by-c, ciphretext ";
  //     decryptAndPrint(ctxt, *dbgKey, *dbgEa, FLAG_PRINT_POLY);
  //     cerr << "  but plaintext poly = ";
  //     printZZX(cerr, tmp1) << endl;
  //   }
  //   for (long i=0; i<=deg(tmp2); i++) {
  //     tmp2[i] /= p2e;
  //     tmp2[i] = -tmp2[i];
  //   }
  //   tmp2.normalize();
  //   if (ePrime<r)
  //     tmp2 += tmp1;
  //   PolyRed(tmp2, conv<ZZ>(p2r), true);
  //   if (dbgPoly != tmp2) {
  //     cerr << "  @ after division by p^{e'}, simple decryption fails:\n";
  //     cerr << "    "<<dbgPoly<<"!="<<tmp2<<endl;
  //   }
  // }

  // Move the powerful-basis coefficients to the plaintext slots
  movePwrflCoefs2Slots(ctxt);

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)

  // On 1st execution, initialize the linearized polynomials
  if (unpackSlotEncoding.size()==0)
    initUnpackEncoding(unpackSlotEncoding, *context.bootstrapEA);

  extractDigitsPacked(ctxt, e, r, ePrime, unpackSlotEncoding);
  // return in ctxt an encryption of the r-digit integer
  // <z_topLow,...,z_{topLow-r+1}> - <z_{botHigh+r-1},...,z_{botHigh}>
  // (modulo p^r), where digits with negative index are zero

  // Move the slots back to powerful-basis coefficients
  moveSlots2PwrflCoefs(ctxt);
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
  argmap["L"] = "7";
  argmap["p"] = "2";
  argmap["m"] = "17";
  argmap["r"] = "1";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long p = atoi(argmap["p"]);
  long m = atoi(argmap["m"]);
  long r = atoi(argmap["r"]);
  long L =  atoi(argmap["L"]);

  FHEcontext context(m,p,r,/*botstrappable=*/true);
  long p2r = context.alMod.getPPowR();
  buildModChain(context, L, /*c=*/3);

  FHESecKey secretKey(context);
  FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(64);      // A Hamming-weight-64 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey);
  secretKey.genBootstrapData();

  // the bootstrapping key is encrypted relative to plaintext space p^{e+r}.
  //  long ePr = FHEPubKey::ePlusR(p);

  //  PAlgebraMod almod2(context.zMStar, context.bootstrapPAM->getR());
  //  EncryptedArray ea2(context, almod2);
  //  EncryptedArray ea1(context, context.alMod);

  dbgKey = &secretKey; // debugging key and ea

  // The bootstrapping data in the context is for plaintext space p^{e+r-e'}
  dbgEa = context.bootstrapEA;

  ZZX poly = RandPoly(context.zMStar.getPhiM(), to_ZZ(p2r));
  PolyRed(poly, p2r, true);
  ZZX poly2;
  Ctxt c1(publicKey);
  cerr << "encrypting ";
  printZZX(cerr, poly)<<endl;

  secretKey.Encrypt(c1,poly,p2r);
  publicKey.reCrypt(c1);
  secretKey.Decrypt(poly2,c1);
  if (poly != poly2) {
    cerr << "\ndecryption error, result is ";
    printZZX(cerr, poly2)<<endl;
  }
  else cerr << "  decryption succeeds!!\n";
  return 0;
}




/**********************************************/
/*********** Debugging utilities **************/
/**********************************************/

ostream& printZZX(ostream& s, const ZZX& poly, long nCoeffs)
{
  long d = deg(poly);
  if (d<nCoeffs) return s << poly; // just print the whole thing

  // otherwise print only 1st nCoeffs coefficiants
  s << '[';
  for (long i=0; i<nCoeffs-2; i++) s << poly[i] << ' ';
  s << "... " << poly[d-1] << ' ' << poly[d] << ']';
  return s;
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

