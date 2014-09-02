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

//#define DEBUGGING

static long mValues[][13] = { 
//{ p, phi(m),  m,    d, m1,  m2, m3,   g1,    g2,    g3,ord1,ord2,ord3}
  {  2,   600,  1023, 10, 11,  93,  0,   838,   584,    0, 10,  6,   0}, // m=(3)*11*{31} m/phim(m)=1.7    C=24  D=2 E=1
  {  2,  1200,  1705, 20, 11, 155,  0,   156,   936,    0, 10,  6,   0}, // m=(5)*11*{31} m/phim(m)=1.42   C=34  D=2 E=2
  {  2, 12800, 17425, 40, 41, 425,  0,  5951,  8078,    0, 40, -8,   0}, // m=(5^2)*{17}*41 m/phim(m)=1.36 C=93  D=3 E=3
  {  2, 15004, 15709, 22, 23, 683,  0,  4099, 13663,    0, 22, 31,   0}, // m=23*(683) m/phim(m)=1.04      C=73  D=2 E=1
  {  2, 18000, 18631, 25, 31, 601,  0, 15627,  1334,    0, 30, 24,   0}, // m=31*(601) m/phim(m)=1.03      C=77  D=2 E=0
  {  2, 18816, 24295, 28, 43, 565,  0, 16386, 16427,    0, 42, 16,   0}, // m=(5)*43*{113} m/phim(m)=1.29  C=84  D=2 E=2
  {  2, 21168, 27305, 28, 43, 635,  0, 10796, 26059,    0, 42, 18,   0}, // m=(5)*43*{127} m/phim(m)=1.28  C=86  D=2 E=2
  {  2, 23040, 28679, 24, 17,  7, 241, 15184,  4098,28204, 16,  6, -10}, // m=7*17*(241) m/phim(m)=1.24    C=63  D=4 E=3
  {  2, 24000, 31775, 20, 41, 775,  0,  6976, 24806,    0, 40, 30,   0}, // m=(5^2)*{31}*41 m/phim(m)=1.32 C=88  D=2 E=2
  {  2, 26400, 27311, 55, 31, 881,  0, 21145,  1830,    0, 30, 16,   0}, // m=31*(881) m/phim(m)=1.03      C=99  D=2 E=0
  {  2, 31104, 35113, 36, 37, 949,  0, 16134,  8548,    0, 36, 24,   0}, // m=(13)*37*{73} m/phim(m)=1.12  C=94  D=2 E=2
  {  2, 34848, 45655, 44, 23, 1985, 0, 33746, 27831,    0, 22, 36,   0}, // m=(5)*23*{397} m/phim(m)=1.31  C=100 D=2 E=2
  {  2, 42336, 42799, 21, 127, 337, 0, 25276, 40133,    0,126, 16,   0}, // m=127*(337) m/phim(m)=1.01     C=161 D=2 E=0
  {  2, 45360, 46063, 45, 73, 631,  0, 35337, 20222,    0, 72, 14,   0}, // m=73*(631) m/phim(m)=1.01      C=129 D=2 E=0
  {  2, 46080, 53261, 24, 17, 13, 241, 43863, 28680,15913, 16, 12, -10}, // m=13*17*(241) m/phim(m)=1.15   C=69  D=4 E=3
  {  2, 49500, 49981, 30, 151, 331, 0,  6952, 28540,    0,150, 11,   0}, // m=151*(331) m/phim(m)=1        C=189 D=2 E=1
  {  2, 54000, 55831, 25, 31, 1801, 0, 19812, 50593,    0, 30, 72,   0}, // m=31*(1801) m/phim(m)=1.03     C=125 D=2 E=0
  {  2, 60016, 60787, 22, 89, 683,  0,  2050, 58741,    0, 88, 31,   0}, // m=89*(683) m/phim(m)=1.01      C=139 D=2 E=1

  { 17,   576,  1365, 12,  7,   3, 65,   976,   911,  463,  6,  2,   4}, // m=3*(5)*7*{13} m/phim(m)=2.36  C=22  D=3
  { 17, 18000, 21917, 30, 101, 217, 0,  5860,  5455,    0, 100, 6,  0}, // m=(7)*{31}*101 m/phim(m)=1.21  C=134 D=2 
  { 17, 30000, 34441, 30, 101, 341, 0,  2729, 31715,    0, 100, 10,  0}, // m=(11)*{31}*101 m/phim(m)=1.14 C=138 D=2
  { 17, 40000, 45551, 40, 101, 451, 0, 19394,  7677,    0, 100, 10,  0}, // m=(11)*{41}*101 m/phim(m)=1.13 C=148 D=2
  { 17, 46656, 52429, 36, 109, 481, 0, 46658,  5778,    0, 108, 12,  0}, // m=(13)*{37}*109 m/phim(m)=1.12 C=154 D=2
  { 17, 54208, 59363, 44, 23, 2581, 0, 25811,  5199,    0, 22, 56,   0}, // m=23*(29)*{89} m/phim(m)=1.09  C=120 D=2
  { 17, 70000, 78881, 10, 101, 781, 0, 67167, 58581,    0, 100, 70,  0}, // m=(11)*{71}*101 m/phim(m)=1.12 C=178 D=2

  {127,   576,  1365, 12,  7,   3, 65,   976,   911,  463,  6,  2,   4}, // m=3*(5)*7*{13} m/phim(m)=2.36   C=22  D=3
  {127,  1200,  1925, 20,  11, 175, 0,  1751,   199,    0, 10, 6,    0}, //  m=(5^2)*{7}*11 m/phim(m)=1.6   C=34 D=2
  {127,  2160,  2821, 30,  13, 217, 0,   652,   222,    0, 12, 6,    0}, // m=(7)*13*{31} m/phim(m)=1.3     C=46 D=2

  {127, 18816, 24295, 28, 43, 565,  0, 16386, 16427,    0, 42, 16,   0}, // m=(5)*43*{113} m/phim(m)=1.29   C=84  D=2
  {127, 26112, 30277, 24, 17, 1781, 0, 14249, 10694,    0, 16, 68,   0}, // m=(13)*17*{137} m/phim(m)=1.15  C=106 D=2
  {127, 39168, 62335, 24, 7,  5, 1781, 17811, 37402,49876,  6,  4,  68}, // m=5*7*(13)*{137} m/phim(m)=1.59 C=100 D=3
  {127, 49392, 61103, 28, 43, 1421, 0,  1422, 14234,    0, 42, 42,   0}, // m=(7^2)*{29}*43 m/phim(m)=1.23  C=110 D=2
  {127, 54400, 61787, 40, 41, 1507, 0, 30141, 46782,    0, 40, 34,   0}, // m=(11)*41*{137} m/phim(m)=1.13  C=112 D=2
  {127, 72000, 77531, 30, 61, 1271, 0,  7627, 34344,    0, 60, 40,   0}  // m=(31)*{41}*61 m/phim(m)=1.07   C=128 D=2
};
#define num_mValues (sizeof(mValues)/(13*sizeof(long)))


/*********** Debugging utilities **************/
FHESecKey* dbgKey=NULL;
EncryptedArray* dbgEa=NULL;
ZZX dbg_ptxt;
Vec<ZZ> ptxt_pwr;

#define FLAG_PRINT_ZZX  1
#define FLAG_PRINT_POLY 2
#define FLAG_PRINT_VEC  4

void baseRep(Vec<long>& rep, long nDigits, ZZ num, long base=2);
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
long makeDivisible(VecInt& vec, long p2e, long p2r, long q, double alpha=0.42)
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
  double tm = -GetTime();

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
  tm += GetTime();
  cerr << "  Unpacking in "<<tm<<" seconds.";
  CheckCtxt(unpacked[0], "After unpacking");

  // Step 2: extract the digits top-1,...,0 from the slots of unpacked[i]
  tm = -GetTime();
  long p = ctxt.getContext().zMStar.getP();
  long p2r = power_long(p,r);
  long topHigh = botHigh + r-1;
  if (p==2 && r>2)
    topHigh--; // For p==2 we sometime get a bit for free

  cerr << "  extracting "<<(topHigh+1)<<" digits:\n";
#ifdef DEBUGGING // check digit-extraction
  bool foundError = false;
#endif
  for (long i=0; i<(long)unpacked.size(); i++) {
    if (topHigh<=0) { // extracting LSB = no-op
      scratch.assign(1,unpacked[i]);
    } else {          // extract digits topHigh...0, store them in scratch
      FHE_NTIMER_START(extractDigits_method);
      extractDigits(scratch, unpacked[i], topHigh+1);
      FHE_NTIMER_STOP(extractDigits_method);
#ifdef DEBUGGING      // check digit-extraction
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
  tm += GetTime();
  cerr << "  Extracting in "<<tm<<" seconds.";
  CheckCtxt(unpacked[0], "After extraction");

  // Step 3: re-pack the slots
  tm = -GetTime();
  FHE_NTIMER_START(repackOverhead);
  EncryptedArray ea2(ctxt.getContext(), ctxt.getContext().alMod);
  long nSlots = ea2.size();
  vector<ZZX> xVec(nSlots, ZZX(1,1)); // X in all the slots
  ZZX xInSlots;
  ea2.encode(xInSlots, xVec);         // encode as ZZX
  FHE_NTIMER_STOP(repackOverhead);
  ZZX x2iInSlots = xInSlots;          // X^i := X

  ctxt = unpacked[0];
  for (long i=1; i<d; i++) {
    unpacked[i].multByConstant(x2iInSlots);
    FHE_NTIMER_START(repackOverhead2);
    MulMod(x2iInSlots, x2iInSlots,xInSlots,ctxt.getContext().zMStar.getPhimX());
    PolyRed(x2iInSlots, conv<ZZ>(p2r));
    FHE_NTIMER_STOP(repackOverhead2);
    ctxt += unpacked[i];
  }
  tm += GetTime();
  cerr << "  Repacking in "<<tm<<" seconds.";
  CheckCtxt(ctxt, "After re-packing");
}
 
// Move the powerful-basis coefficients to the plaintext slots
void movePwrflCoefs2Slots(Ctxt& ctxt)
{
#ifdef DEBUGGING
  ZZX ptxt1; // for debugging
  long p = ctxt.getPtxtSpace();
  if (dbgKey && dbgEa) dbgKey->Decrypt(ptxt1, ctxt);
#endif
  const FHEcontext& context = ctxt.getContext();
  context.firstMap->apply(ctxt);

#ifdef DEBUGGING
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
	long coef1 = conv<long>(coeff(v[i],j));
	long coef2 = conv<long>(powerful[idx++]);
	if (coef1 != coef2) {
	  cerr << " @ error: powerful coefficient "<<i<<" is "<<coef2
	       << " but slot contains "<<coef1<<endl;
	  return;
	}
      }
    cerr << "  LinTrans1 successful. ";
    decryptAndPrint(ctxt, *dbgKey, *dbgEa, 0);
  }
#endif
}

// Move the slots back to powerful-basis coefficients
void moveSlots2PwrflCoefs(Ctxt& ctxt)
{
#ifdef DEBUGGING
  vector<ZZX> v; // for debugging
  if (dbgKey && dbgEa) dbgEa->decrypt(ctxt,*dbgKey,v);
#endif

  const FHEcontext& context = ctxt.getContext();
  context.secondMap->apply(ctxt);

  // For debugging: Check that we have powerful representation in the slots
#ifdef DEBUGGING
  if (dbgKey && dbgEa) {
    ZZX ptxt1;
    dbgKey->Decrypt(ptxt1, ctxt);

    zz_pBak bak; bak.save(); // backup NTL's current modulus
    long p = ctxt.getPtxtSpace();
    zz_p::init(p);
    PowerfulConversion pConv(context.p2dConversion->getIndexTranslation());
    HyperCube<zz_p> powerful(pConv.getShortSig());
  
    zz_pX poly = conv<zz_pX>(ptxt1);
    pConv.polyToPowerful(powerful, poly);

    long idx = 0;
    long nSlots = dbgEa->size();
    long d =  dbgEa->getDegree();
    for (long i=0; i<nSlots; i++) for (long j=0; j<d; j++) {
	long coef1 = conv<long>(coeff(v[i],j) % p);
	long coef2 = conv<long>(powerful[idx++]);
	if (coef1 != coef2) {
	  cerr << " @ error: powerful coefficient "<<i<<" is "<<coef2
	       << " but slot contains "<<coef1<<endl;
	  return;
	}
      }
    cerr << "  LinTrans2 successful. ";
    decryptAndPrint(ctxt, *dbgKey, *dbgEa, 0);
  }
#endif
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
      cerr << "  Decryption error after key-switching to bootsrapping key\n";
    else {
      cerr << "  After key-switching to bootsrapping key: ";
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

#ifdef DEBUGGING
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

  // Multiply the post-processed cipehrtext by the encrypted sKey
  ctxt = bootstrapEkey;
  ctxt.multByConstant(zzParts[1]);
  ctxt.addConstant(zzParts[0]);
  CheckCtxt(ctxt, "After pre-processing");

  e -= ePrime;
  p2e /= p2ePrime;     // reduce the plaintext space by p^{e'} factor
  ctxt.reducePtxtSpace(/*newPtxtSpace=*/p2e*p2r);
  FHE_NTIMER_STOP(preProcess);

#ifdef DEBUGGING
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
  /*  cerr << "  e.g, 'p=2 e=4 q=257'\n\n";
      cerr << "  p is the plaintext base [default=2]\n";
      cerr << "  e is the exponent [default=4]\n";
  */
  cerr << "  r determines plaintext space mod 2^r [default=1]\n";
  cerr << "  L is # of primes in the chain [default=20]\n";
  cerr << "  N is a lower bound on phi(m) [default=0]\n";
  exit(0);
}

void TestIt(long idx, long p, long r, long L)
{
  Vec<long> mvec;
  vector<long> gens;
  vector<long> ords;

  long m = mValues[idx][1];
  assert(GCD(p, m) == 1);

  append(mvec, mValues[idx][3]);
  append(mvec, mValues[idx][4]);
  if (mValues[idx][5]>1) append(mvec, mValues[idx][5]);
  gens.push_back(mValues[idx][6]);
  gens.push_back(mValues[idx][7]);
  if (mValues[idx][8]>1) gens.push_back(mValues[idx][8]);
  ords.push_back(mValues[idx][9]);
  ords.push_back(mValues[idx][10]);
  if (abs(mValues[idx][11])>1) ords.push_back(mValues[idx][11]);

  cerr << "*** TestIt: p=" << p
       << ", r=" << r
       << ", L=" << L
       << ", m=" << m
       << " (=" << mvec << "), gens="<<gens<<", ords="<<ords
       << endl;

  setTimersOn();
  FHE_NTIMER_START(initialize);
  FHEcontext context(m, p, r, gens, ords);
  buildModChain(context, L, /*c=*/3);
  context.zMStar.printout();
  double t = -GetTime();
  cerr << "Computing key-independent bootstrapping tables..." << std::flush;
  context.makeBootstrappable(mvec);
  t += GetTime();
  cerr << " done in "<<t<<" seconds\n";
  long p2r = context.alMod.getPPowR();

  cerr << "Generating keys\n" << std::flush;
  FHESecKey secretKey(context);
  FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(64);      // A Hamming-weight-64 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey);
  cerr << "Computing key-dependent bootstrapping tables..." << std::flush;
  t = -GetTime();
  secretKey.genBootstrapData();
  t += GetTime();
  cerr << " done in "<<t<<" seconds\n";

  FHE_NTIMER_STOP(initialize);
  cerr << "****Initialization time:\n";
  printAllTimers();
  cerr << endl;
  resetAllTimers();

  dbgKey = &secretKey; // debugging key and ea
  dbgEa = context.bootstrapEA; // EA for plaintext space p^{e+r-e'}

  zz_p::init(p2r);
  zz_pX poly_p = random_zz_pX(context.zMStar.getPhiM());
  PowerfulConversion pConv(context.p2dConversion->getIndexTranslation());
  HyperCube<zz_p> powerful(pConv.getShortSig());
  pConv.polyToPowerful(powerful, poly_p);
  
  conv(dbg_ptxt, poly_p);
  PolyRed(dbg_ptxt, p2r, true);
  context.p2dConversion->ZZXtoPowerful(ptxt_pwr, dbg_ptxt);
  vecRed(ptxt_pwr, ptxt_pwr, p2r, true);

  ZZX poly2;
  Ctxt c1(publicKey);
  cerr << "encrypting (in powerful rep) "; // <<dbg_ptxt<<endl;
  printVec(cerr, powerful.getData())<<endl;

  secretKey.Encrypt(c1,dbg_ptxt,p2r);
  FHE_NTIMER_START(reCrypt);
  publicKey.reCrypt(c1);
  FHE_NTIMER_STOP(reCrypt);
  secretKey.Decrypt(poly2,c1);
  //  decryptAndPrint(c1, secretKey, *dbgEa, 0);

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
	if (numDiff >5) break;
      }
    cerr << endl<< endl;
  }
  else cerr << "  decryption succeeds!!\n\n";

  cerr << "****Bootstrapping time:\n";
  printAllTimers();
  cerr << endl<< endl;
}

int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  //  argmap["p"] = "2";
  argmap["r"] = "1";
  argmap["L"] = "20";
  argmap["N"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  //  long p = atoi(argmap["p"]);
  long p=2;
  long r = atoi(argmap["r"]);
  long L =  atoi(argmap["L"]);
  long N =  atoi(argmap["N"]);

  for (long i=0; i<(long)num_mValues; i++) if (mValues[i][0]>=N) {
      TestIt(i,p,r,L);
      break;
    }
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

void baseRep(Vec<long>& rep, long nDigits, ZZ num, long base)
{
  rep.SetLength(nDigits);
  for (long j=0; j<nDigits; j++) {
    rep[j] = rem(num, base);
    if (rep[j] > base/2)         rep[j] -= base;
    else if (rep[j] < -(base/2)) rep[j] += base;
    num = (num - rep[j]) / base;
  }
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
