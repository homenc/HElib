/** homAES.cpp - homomorphic AES using HElib
 */
namespace std {} using namespace std;
namespace NTL {} using namespace NTL;
#include <cstring>
#include "homAES.h"

#ifdef DEBUG_PRINTOUT
#define FLAG_PRINT_ZZX  1
#define FLAG_PRINT_POLY 2
#define FLAG_PRINT_VEC  4
extern void decryptAndPrint(ostream& s, const Ctxt& ctxt, const SecKey& sk,
			    const EncryptedArray& ea, long flags=0);
extern SecKey* dbgKey;
extern EncryptedArray* dbgEa;
#endif

// Decletration of local functions

static void buildAffine(vector<PolyType>& binMat, PolyType* binVec,
			const unsigned char cc[],
			const EncryptedArrayDerived<PA_GF2>& ea2);
inline void buildAffineEnc(vector<PolyType>& binMat, PolyType& binVec,
			   const EncryptedArrayDerived<PA_GF2>& ea2)
{
  /* The affine GF2-transformation on encryption is:
   *     [ 1 0 0 0 1 1 1 1 ] ( b0 )   ( 1 )
   *     [ 1 1 0 0 0 1 1 1 ] ( b1 )   ( 1 )
   *     [ 1 1 1 0 0 0 1 1 ] ( b2 )   ( 0 )
   *     [ 1 1 1 1 0 0 0 1 ] ( b3 )   ( 0 )
   *     [ 1 1 1 1 1 1 0 0 ]*( b4 ) + ( 0 )
   *     [ 0 1 1 1 1 1 0 0 ] ( b5 )   ( 1 )
   *     [ 0 0 1 1 1 1 1 0 ] ( b6 )   ( 1 )
   *     [ 0 0 0 1 1 1 1 1 ] ( b7 )   ( 0 )
   *
   * To encode this transformation, we use the binary representation of
   * the columns, which is 31,62,124,248,241,227,199,143, and 99
   */
  static unsigned char cc[] = { 31, 62, 124, 248, 241, 227, 199, 143, 99 };
  buildAffine(binMat, &binVec, cc, ea2);
}
inline void buildAffineDec(vector<PolyType>& binMat,
			   const EncryptedArrayDerived<PA_GF2>& ea2)
{
  /* The affine GF2-transformation on decryption is:
   *     [ 0 0 1 0 0 1 0 1 ] [( b0 )   ( 1 )]
   *     [ 1 0 0 1 0 0 1 0 ] [( b1 )   ( 1 )]
   *     [ 0 1 0 0 1 0 0 1 ] [( b2 )   ( 0 )]
   *     [ 1 0 1 0 0 1 0 0 ] [( b3 )   ( 0 )]
   *     [ 0 1 0 1 0 0 1 0 ]*[( b4 ) - ( 0 )]
   *     [ 0 0 1 0 1 0 0 1 ] [( b5 )   ( 1 )]
   *     [ 1 0 0 1 0 1 0 0 ] [( b6 )   ( 1 )]
   *     [ 0 1 0 0 1 0 1 0 ] [( b7 )   ( 0 )]
   *
   * To encode this transformation, we use the binary representation of
   * the columns, which is 74,148,41,82,164,73,146,37
   */
  static unsigned char cc[] = { 74, 148, 41, 82, 164, 73, 146, 37, 0};
  buildAffine(binMat, nullptr, cc, ea2);
}

// FIXME: The implementation below relies on the first generator in thw
//   PAlgebra corresponding to the MSB of the slot number. A more robust
//   implementation would use the Hypercube class for slot index arithmetic.


// Compute the constants for the shoftRow/mixCol transformations
static void buildLinEnc(vector<PolyType>& encLinTran,
			const EncryptedArrayDerived<PA_GF2>& ea2);
static void buildLinDec(vector<PolyType>& decLinTran,
			const EncryptedArrayDerived<PA_GF2>& ea2);

// Apply the shoftRow/mixCol transformations
static void encRowColTran(Ctxt& eData, const vector<PolyType>& encLinTran,
			  const EncryptedArrayDerived<PA_GF2>& ea2);
static void encRowShift(Ctxt& c, const vector<PolyType>& encLinTran,
			const EncryptedArrayDerived<PA_GF2>& ea2);
static void decRowColTran(Ctxt& eData, const vector<PolyType>& decLinTran,
			  const EncryptedArrayDerived<PA_GF2>& ea2);
static void decRowShift(Ctxt& c, const vector<PolyType>& decLinTran,
			const EncryptedArrayDerived<PA_GF2>& ea2);

static void invert(vector<Ctxt>& data); // Z -> Z^{-1} in GF(2^8)

// Pack the ciphertexts in c in as few "fully packed" cipehrtext as possible.
static void packCtxt(vector<Ctxt>& to, const vector<Ctxt>& from,
		     const GF2X& XinSlots);

// Unpack the fully-packed ciphertext in from into the vector to. If to.size()>0
// then do not unpack into more than to.size() ciphertexts. If 'from' does not
// have enough ciphertexts to fill all of 'to' then pad with zeros.
static void unackCtxt(vector<Ctxt>& to, const vector<Ctxt>& from,
		      const Mat<GF2X>& unpackConsts);

// Implementation of the class HomAES

static const uint8_t aesPolyBytes[] = { 0x1B, 0x1 }; // X^8+X^4+X^3+X+1
const GF2X HomAES::aesPoly = GF2XFromBytes(aesPolyBytes, 2);

HomAES::HomAES(const Context& context): ea2(context,aesPoly,context.alMod)
#ifndef USE_ZZX_POLY // initialize DoubleCRT using the context
, affVec(context)
#endif
{
  // Sanity-check: we need the first dimension to be divisible by 16.
  //OLD: assert( context.zMStar.OrderOf(0) % 16 == 0 );
  helib::assertEq(context.zMStar.OrderOf(0) % 16, 0l);

  // Compute the GF2-affine transformation constants
  buildAffineEnc(encAffMat, affVec, ea2);
  buildAffineDec(decAffMat, ea2);

  // Compute the rowShift/colMix constants
  buildLinEnc(encLinTran, ea2);
  buildLinDec(decLinTran, ea2);

  if (context.isBootstrappable())
    setPackingConstants();
}

// run the AES key-expansion and then encrypt the expanded key.
void HomAES::encryptAESkey(vector<Ctxt>& eKey, Vec<uint8_t>& aesKey,
			   const PubKey& hePK) const
{
  //OLD: assert(aesKey.length()==16 || aesKey.length()==24 || aesKey.length()==32);
  helib::assertTrue<helib::InvalidArgument>(aesKey.length()==16 || aesKey.length()==24 || aesKey.length()==32,
                                            "Key length must be 16, 24, or 32");


  // We rely on some implementation of this function
  extern long AESKeyExpansion(unsigned char roundKeySchedule[],
			      unsigned char key[], int keyBits);

  // Compute the key expansion
  uint8_t roundKeySchedule[240];
  long nRoundKeys =
    AESKeyExpansion(roundKeySchedule, aesKey.data(), aesKey.length()*8);

  long blocksPerCtxt = ea2.size() / 16;

  // Expand the key-schedule, copying each round key blocksPerCtxt times
  Vec<uint8_t> expanded(INIT_SIZE, nRoundKeys*blocksPerCtxt*16);
  for (long i=0; i<nRoundKeys; i++) {
    uint8_t* roundKey = &roundKeySchedule[16*i];
    for (long j=0; j<blocksPerCtxt; j++)
      memcpy(&expanded[16*(i*blocksPerCtxt +j)], roundKey, 16);
  }
  Vec<ZZX> encoded;
  encode4AES(encoded, expanded, ea2);      // encode as HE plaintext

  {Ctxt tmpCtxt(hePK);
  eKey.resize(encoded.length(), tmpCtxt);} // allocate space
  for (long i=0; i<(long)eKey.size(); i++) // encrypt the encoded key
    hePK.Encrypt(eKey[i], encoded[i]);
}


// Compute the packing/unpacking constants
void HomAES::setPackingConstants()
{
  // Get the context and the ea for "fully packed" polynomials
  const Context& context = ea2.getContext();
  const EncryptedArrayDerived<PA_GF2>& ea = context.ea->getDerived(PA_GF2());

  // Compute the packing constants, with X in all the slots
  {vector<GF2X> slots(ea.size(), GF2X(1,1)); // X in all the slots
  ZZX tmp; ea.encode(tmp, slots); // encode as ZZX
  conv(XinSlots, tmp);}           // convert to ZZ2X

  // Compute the unpacking constants

  long e = ea.getDegree() / 8; // the extension degree
  //OLD: assert(ea.getDegree()==e*8 && e<=(long) sizeof(long));
  helib::assertEq(ea.getDegree()==e*8, "ea must have degree divisible by 8");
  helib::assertTrue(e<=(long) sizeof(long), "extension degree must be at most 8 times sizeof(long)");

  GF2EBak bak; bak.save(); // save current modulus (if any)
  GF2XModulus F0(ea.getTab().getFactors()[0]);
  GF2E::init(F0);

  // Set a matrix for converting from z = sum_i zi*X^i in GF(2^d)
  // back to the vectors of zi's in GF(2^8)

  mat_GF2E Kinv(INIT_SIZE, e, e);
  for (long j=0; j<e; j++)
    conv(Kinv[0][j], GF2X(j,1)); // Kinv[0][j] = X^j
  for (long i=1; i<e; i++) for (long j=0; j<e; j++) {
    power(Kinv[i][j], Kinv[0][j], 1L<<(8*i)); // (X^j)^{2^{8*i}}
  }
  mat_GF2E K(INIT_SIZE, e, e);
  inv(K, Kinv); // invert Kinv to get K

  // Encode K in slots
  unpacking.SetDims(e,e);
  for (long i=0; i<e; i++) for (long j=0; j<e; j++) {
    vector<GF2X> slots(ea.size(), rep(K[i][j])); // K[i][j] in all the slots
    ZZX tmp; ea.encode(tmp, slots);
    conv(unpacking[i][j], tmp);
  }
}



// Perform in-place AES encryption on HE-encrypted plaintext bytes (ECB mode).
// The input and AES key are encrypted under HE. The output is doubly-encrypted,
// out=Enc_HE(Enc_AES(X)). The aesKey array contains an encryption of the
// expanded AES key, the number of AES rounds is aesKey.size() -1.
// It is assumed that all the input data cipehrtexts are at the same
// level, as they will be recrypted together.
void HomAES::homAESenc(vector<Ctxt>& eData, const vector<Ctxt>& aesKey) const
{
  if (1>(long)eData.size() || 1>(long)aesKey.size()) return; // no data/key
  //  long lvlBits = eData[0].getContext().bitsPerLevel;

  for (long j=0; j<(long)eData.size(); j++)
    eData[j] += aesKey[0];  // initial key addition

  for (long i=1; i<(long)aesKey.size(); i++) { // apply the AES rounds

    // ByteSub
    if (eData[0].findBaseLevel() < 4) batchRecrypt(eData);
    invert(eData);     // apply Z -> Z^{-1} to all elements of eData
#ifdef DEBUG_PRINTOUT
    CheckCtxt(eData[0], "+ After invert");
    //    cerr << " + After invert ";
    //    decryptAndPrint(cerr, eData[0], *dbgKey, *dbgEa);
#endif
    if (eData[0].findBaseLevel() < 2) batchRecrypt(eData);
    for (long j=0; j<(long)eData.size(); j++) { // GF2 affine transformation
      applyLinPolyLL(eData[j], encAffMat, ea2.getDegree());
      eData[j].addConstant(affVec);
    }
#ifdef DEBUG_PRINTOUT
    CheckCtxt(eData[0], "+ After affine");
    //    cerr << " + After affine ";
    //    decryptAndPrint(cerr, eData[0], *dbgKey, *dbgEa);
#endif

    // Apply RowShift/ColMix to each ciphertext
    if (eData[0].findBaseLevel() < 2) batchRecrypt(eData);
    if (i<(long)aesKey.size()-1) {
      for (long j=0; j<(long)eData.size(); j++)
	encRowColTran(eData[j], encLinTran, ea2);
    }
    else { // For the last round apply only RowShift, not ColMix
      for (long j=0; j<(long)eData.size(); j++)
	encRowShift(eData[j], encLinTran, ea2);
    }
#ifdef DEBUG_PRINTOUT
    CheckCtxt(eData[0], "+ After rowShift/colMix");
    //    cerr << " + After rowShift/colMix ";
    //    decryptAndPrint(cerr, eData[0], *dbgKey, *dbgEa);
#endif

    // Key addition
    for (long j=0; j<(long)eData.size(); j++) eData[j] += aesKey[i];
  }
}

// Perform AES encryption on plaintext bytes (ECB mode). The input are
// raw plaintext bytes, and the AES key encrypted under HE. The output
// is a doubly-encrypted ciphertext, out=Enc_HE(Enc_AES(X)). The aesKey
// array contains an encryption of the expanded AES key, the number of
// AES rounds is aesKey.size() -1.
// NOTE: This is a rather useless method, other than for benchmarking
void HomAES::homAESenc(vector<Ctxt>& eData, const vector<Ctxt>& aesKey,
		       const Vec<uint8_t> inBytes) const
{
  {Vec<ZZX> encodedBytes;
  encode4AES(encodedBytes, inBytes, ea2); // encode as HE plaintext

  // Allocate space for the output ciphertexts, initialized to zero
  eData.resize(encodedBytes.length(), Ctxt(ZeroCtxtLike,aesKey[0]));
  for (long i=0; i<(long)eData.size(); i++)   // encode ptxt as HE ctxt
    eData[i].DummyEncrypt(encodedBytes[i]);}

  homAESenc(eData, aesKey); // do the real work
}


// Perform in-place AES decryption on Doubly encrypted bytes (ECB mode).
// The input is a doubly encrypted cipehrtext, in=Enc_HE(Enc_AES(X)), and
// the HE-encrypted AES key. The output is a "plaintext" (but encrypted
// under the HE scheme), out=Enc_HE(X). The aesKey array containsencryption
// of the expanded AES key, the number of AES rounds is aesKey.length()-1.
// It is assumed that all the input data cipehrtexts are at the same
// level, as they will be recrypted together.
void HomAES::homAESdec(vector<Ctxt>& eData, const vector<Ctxt>& aesKey) const
{
  if (1>(long)eData.size() || 1>(long)aesKey.size()) return; // no data/key
  //  long lvlBits = eData[0].getContext().bitsPerLevel;

  for (long i=aesKey.size()-1; i>0; i--) { // apply the AES rounds
    // Key addition
    for (long j=0; j<(long)eData.size(); j++) eData[j] -= aesKey[i];

    // Apply RowShift/ColMix to each ciphertext
    if (eData[0].findBaseLevel() < 2) batchRecrypt(eData);
    //    if (eData[0].log_of_ratio() > (-lvlBits)) batchRecrypt(eData);
    if (i<(long)aesKey.size()-1)
      for (long j=0; j<(long)eData.size(); j++)
	decRowColTran(eData[j], decLinTran, ea2);

    else // For the first round apply only RowShift, not ColMix
      for (long j=0; j<(long)eData.size(); j++)
	decRowShift(eData[j], decLinTran, ea2);
#ifdef DEBUG_PRINTOUT
    CheckCtxt(eData[0], "+ After rowShift/colMix");
    //    cerr << " + After rowShift/colMix ";
    //    decryptAndPrint(cerr, eData[0], *dbgKey, *dbgEa);
#endif

    // ByteSub
    if (eData[0].findBaseLevel() < 2) batchRecrypt(eData);
    for (long j=0; j<(long)eData.size(); j++) { // GF2 affine transformation
      eData[j].addConstant(affVec);
      applyLinPolyLL(eData[j], decAffMat, ea2.getDegree());
    }
#ifdef DEBUG_PRINTOUT
    CheckCtxt(eData[0], "+ After affine");
    //    cerr << " + After affine ";
    //    decryptAndPrint(cerr, eData[0], *dbgKey, *dbgEa);
#endif
    if (eData[0].findBaseLevel() < 4) batchRecrypt(eData);
     invert(eData); // apply Z -> Z^{-1} to all elements of eData
#ifdef DEBUG_PRINTOUT
    CheckCtxt(eData[0], "+ After invert");
    //    cerr << " + After invert ";
    //    decryptAndPrint(cerr, eData[0], *dbgKey, *dbgEa);
#endif
  }

  for (long j=0; j<(long)eData.size(); j++)
    eData[j] -= aesKey[0];  // final key addition
}

// Perform AES decryption on AES ciphertext bytes (ECB mode). The input
// are "raw AES-encrypted bytes", in=Enc_AES(X), and the HE-encrypted AES key.
// The output is encrypted under HE but not AES, out=Enc_HE(X). The aesKey
// array containsencryption of the expanded AES key, the number of AES rounds
// is aesKey.size()-1.
void HomAES::homAESdec(vector<Ctxt>& eData, const vector<Ctxt>& aesKey,
		       const Vec<uint8_t> inBytes) const
{
  {Vec<ZZX> encodedBytes;
  encode4AES(encodedBytes, inBytes, ea2); // encode as HE plaintext

  // Allocate space for the output ciphertexts, initialized to zero
  eData.resize(encodedBytes.length(), Ctxt(ZeroCtxtLike,aesKey[0]));
  for (long i=0; i<(long)eData.size(); i++)   // encode ptxt as HE ctxt
    eData[i].DummyEncrypt(encodedBytes[i]);}

  homAESdec(eData, aesKey); // do the real work
}

// Implementation of local functions


void HomAES::batchRecrypt(vector<Ctxt>& data) const
{
  FHE_TIMER_START;
  PubKey& pk = (PubKey&) data[0].getPubKey();
  if (!pk.isBootstrappable()) return;

  if (data.size()>1 && unpacking.NumRows()==0) // lazy initialization
    return; // setPackingConstants();

  vector<Ctxt>* pData = &data;
  vector<Ctxt> fullyPacked; // empty at first
  if (data.size()>1) {      // pack to save on recryption operations
    packCtxt(fullyPacked, data, XinSlots);
    pData = &fullyPacked;
  }

#ifdef DEBUG_PRINTOUT
  long dSize = data.size();
  Vec<ZZX> ptxt(INIT_SIZE, dSize);
  for (long i=0; i<dSize; i++)
    dbgKey->Decrypt(ptxt[i], data[i]);

  cerr << "   > Before recrypt ";
  decryptAndPrint(cerr, (*pData)[0], *dbgKey, *dbgEa);
#endif

  // recrypt each ciphertext in the vector
  FHE_NTIMER_START(recryption);
  for (long i=0; i<(long)pData->size(); i++)
    pk.reCrypt((*pData)[i]);
  FHE_NTIMER_STOP(recryption);

  // unpack back to the original vector, if needed
  if (fullyPacked.size()>0) {
    unackCtxt(data, fullyPacked, unpacking);
  }

#ifdef DEBUG_PRINTOUT
  if (dSize != (long)data.size()) {
    cerr << " size mismatch after pack/unpack\n";
    exit(0);
  }
  for (long i=0; i<(long)data.size(); i++) {
    ZZX ptxt2;
    dbgKey->Decrypt(ptxt2, data[i]);
    if (ptxt[i] != ptxt2) {
      cerr << " ciphertext "<<i<<" mismatch after pack/unpack\n";
      exit(0);
    }
  }
  cerr << "   > After batchRecrypt ";
  decryptAndPrint(cerr, data[0], *dbgKey, *dbgEa);
#endif
}



// Buils a GF2-affine transformation from the constants in cc
static void buildAffine(vector<PolyType>& binMat, PolyType* binVec,
			const unsigned char cc[],
			const EncryptedArrayDerived<PA_GF2>& ea2)
{
  vector<GF2X> scratch(8); // Specify the different columns
  for (long j = 0; j < 8; j++) // convert from byte to degree-7 polynomial
    GF2XFromBytes(scratch[j], &cc[j], 1);

  // "building" the linearized-polynomial coefficients
  vector<GF2X> C;
  ea2.buildLinPolyCoeffs(C, scratch);

  // "encoding" the coefficients
#ifdef USE_ZZX_POLY
  vector<ZZX>& zzxMat=binMat;
#else
  vector<ZZX> zzxMat;
#endif
  zzxMat.resize(8);
  scratch.resize(ea2.size());
  for (long j = 0; j < 8; j++) {
    for (long i = 0; i < ea2.size(); i++) // set all slots to C[j]
      scratch[i] = C[j];
    ea2.encode(zzxMat[j], scratch);       // encode these slots
  }
#ifndef USE_ZZX_POLY
  binMat.resize(8,DoubleCRT(ea2.getContext()));
  for (long j=0; j<8; j++) binMat[j] = zzxMat[j]; // convert to DoubleCRT
#endif

  if (binVec != nullptr) {
    GF2X cc8;
    GF2XFromBytes(cc8, &cc[8], 1);
    for (long i = 0; i < ea2.size(); i++)  // set all slots to cc8
      scratch[i] = cc8;
#ifdef USE_ZZX_POLY
    ea2.encode(*binVec, scratch);       // encode these slots
#else
    ZZX tmpZZX;
    ea2.encode(tmpZZX, scratch);        // encode these slots
    *binVec = tmpZZX;                   // conveer to DoubleCRT
#endif
  }
}


// Compute the constants for the shoftRow/mixCol transformations
static void buildLinEnc(vector<PolyType>& encLinTran,
			const EncryptedArrayDerived<PA_GF2>& ea2)
{
  // The constants only have nonzero entires in their slots corresponding
  // to bytes 0,4,8,12 of each blocks. Since the blocks are interleaved, it
  // means that the non-zero slots are nBlocks*{0,4,8,12}+{0,1,...,nBlocks-1}.
  // The three constants have in these slots 1, X, and X+1, corresponding
  // to hex values 1, 2, 3 respectively.
  uint8_t hexVals[3] = { 1, 2, 3 };

  Vec<uint8_t> bytes(INIT_SIZE, ea2.size());
  memset(bytes.data(), 0, bytes.length());
  long blocksPerCtxt = ea2.size() / 16;
  Vec<ZZX> tmp;
#ifdef USE_ZZX_POLY
  encLinTran.resize(6);
#else
  encLinTran.resize(6,DoubleCRT(ea2.getContext()));
#endif
  for (long i=0; i<3; i++) { // constants for the RowShift/ColMix trans
    for (long j=0; j<blocksPerCtxt; j++) {
      uint8_t* bptr = &bytes[16*j];
      bptr[0] = bptr[4] = bptr[8] = bptr[12] = hexVals[i];
    }
    encode4AES(tmp, bytes, ea2);
    encLinTran[i] = tmp[0];
  }
  for (long i=0; i<3; i++) { // masks for only RowShift (last round)
    memset(bytes.data(), 0, bytes.length());
    for (long j=0; j<blocksPerCtxt; j++) {
      uint8_t* bptr = &bytes[16*j +i+1];
      bptr[0] = bptr[4] = bptr[8] = bptr[12] = 1;
    }
    encode4AES(tmp, bytes, ea2);
    encLinTran[i+3] = tmp[0];
  }
}

// Apply the shoftRow/mixCol transformations
static void encRowColTran(Ctxt& c, const vector<PolyType>& encLinTran,
			  const EncryptedArrayDerived<PA_GF2>& ea2)
{
  // The basic rotation amount along the 1st dimension
  long rotAmount = ea2.getContext().zMStar.OrderOf(0) / 16;

  c.cleanUp();
  Ctxt c1(c), c6(c), c11(c);
  ea2.rotate1D(c1, 0, rotAmount);   // rotation + re-linearization
  ea2.rotate1D(c6, 0, 6*rotAmount); // rotation + re-linearization
  ea2.rotate1D(c11,0,11*rotAmount); // rotation + re-linearization
  c1.cleanUp();  c6.cleanUp();  c1.cleanUp();

  const PolyType& p1 = encLinTran[0]; // 1
  const PolyType& p2 = encLinTran[1]; // X
  const PolyType& p3 = encLinTran[2]; // X+1

  /* return ( c*p2 + c1*p1 + c6*p1 + c11*p3
   *         + (c*p1 + c1*p1 + c6*p3 + c11*p2)>>1
   *         + (c*p1 + c1*p3 + c6*p2 + c11*p1)>>2
   *         + (c*p3 + c1*p2 + c6*p1 + c11*p1)>>3 )
   */

  // Compute the top row of the AES matrix
  Ctxt sum(c), tmpA(c1), tmpB(c11);
  tmpA += c6;
  sum.multByConstant(p2);
  tmpA.multByConstant(p1);
  tmpB.multByConstant(p3);
  sum += tmpA; sum += tmpB;

  // Compute the 2nd row of the AES matrix
  tmpA = c; tmpA += c1; tmpA.multByConstant(p1);
  tmpB = c6;            tmpB.multByConstant(p3);
  Ctxt tmpC(c11);       tmpC.multByConstant(p2);
  tmpA += tmpB; tmpA += tmpC;
  ea2.rotate1D(tmpA, 0, rotAmount); // rotation + re-linearization
  sum += tmpA;

  // Compute the 3rd row of the AES matrix
  tmpA = c;  tmpA += c11;  tmpA.multByConstant(p1);
  tmpB = c1;               tmpB.multByConstant(p3);
  tmpC = c6;               tmpC.multByConstant(p2);
  tmpA += tmpB; tmpA += tmpC;
  ea2.rotate1D(tmpA, 0, 2*rotAmount); // rotation + re-linearization
  sum += tmpA;

  // Compute the bottom row of the AES matrix
  c.multByConstant(p3);
  c1.multByConstant(p2);
  c6 += c11; c6.multByConstant(p1);
  c += c1; c += c6;
  ea2.rotate1D(c, 0, 3*rotAmount); // rotation + re-linearization
  c += sum;

  c.cleanUp();
}

static void encRowShift(Ctxt& c, const vector<PolyType>& encLinTran,
			const EncryptedArrayDerived<PA_GF2>& ea2)
{
  // The basic rotation amount along the 1st dimension
  long rotAmount = ea2.getContext().zMStar.OrderOf(0) / 16;

  c.cleanUp();
  Ctxt c4(c), c8(c), c12(c);
  ea2.rotate1D(c4, 0, 4*rotAmount);
  ea2.rotate1D(c8, 0, 8*rotAmount);
  ea2.rotate1D(c12,0, 12*rotAmount);
  c4.cleanUp();  c8.cleanUp();  c12.cleanUp();

  const PolyType& p0 = encLinTran[0]; // 1
  const PolyType& p4 = encLinTran[3]; // 1>>4
  const PolyType& p8 = encLinTran[4]; // 1>>8
  const PolyType& p12 = encLinTran[5]; // 1>>12

  // return c*p0 + c4*p12 + c8*p8 + c12*p4

  c.multByConstant(p0);
  c4.multByConstant(p12); c += c4;
  c8.multByConstant(p8);  c += c8;
  c12.multByConstant(p4); c += c12;
}

static void buildLinDec(vector<PolyType>& decLinTran,
			const EncryptedArrayDerived<PA_GF2>& ea2)
{
  // The constants only have nonzero entires in their slots corresponding
  // to bytes 0,4,8,12 of each blocks. Since the blocks are interleaved, it
  // means that the non-zero slots are nBlocks*{0,4,8,12}+{0,1,...,nBlocks-1}.
  // The three constants have in these slots the polynomials corresponding
  // to hex values ox9, 0xB, 0xD, and 0xE respectively.

  uint8_t hexVals[4] = { 0x9, 0xB, 0xD, 0xE };

  Vec<uint8_t> bytes(INIT_SIZE, ea2.size());
  memset(bytes.data(), 0, bytes.length());
  long blocksPerCtxt = ea2.size() / 16;

  Vec<ZZX> tmp;
#ifdef USE_ZZX_POLY
  decLinTran.resize(8);
#else
  decLinTran.resize(8,DoubleCRT(ea2.getContext()));
#endif
  for (long i=0; i<4; i++) { // constants for the RowShift/ColMix trans
    for (long j=0; j<blocksPerCtxt; j++) {
      uint8_t* bptr = &bytes[16*j];
      bptr[0] = bptr[4] = bptr[8] = bptr[12] = hexVals[i];
    }
    encode4AES(tmp, bytes, ea2);
    decLinTran[i] = tmp[0];
  }
  for (long i=0; i<4; i++) { // masks for only RowShift (first round)
    memset(bytes.data(), 0, bytes.length());
    for (long j=0; j<blocksPerCtxt; j++) {
      uint8_t* bptr = &bytes[16*j +i];
      bptr[0] = bptr[4] = bptr[8] = bptr[12] = 1;
    }
    encode4AES(tmp, bytes, ea2);
    decLinTran[i+4] = tmp[0];
  }
}


static void decRowColTran(Ctxt& c, const vector<PolyType>& decLinTran,
			  const EncryptedArrayDerived<PA_GF2>& ea2)
{
  // The basic rotation amount along the 1st dimension
  long rotAmount = ea2.getContext().zMStar.OrderOf(0) / 16;

  c.cleanUp();
  Ctxt cf(c), ce(c), cd(c);
  ea2.rotate1D(cf, 0, 15*rotAmount); // rotation + re-linearization
  ea2.rotate1D(ce, 0, 14*rotAmount); // rotation + re-linearization
  ea2.rotate1D(cd, 0, 13*rotAmount); // rotation + re-linearization
  cf.cleanUp();  ce.cleanUp();  cd.cleanUp();

  const PolyType& p9 = decLinTran[0]; // poly represented by 0x9
  const PolyType& pB = decLinTran[1]; // poly represented by 0xE
  const PolyType& pD = decLinTran[2]; // poly represented by 0xD
  const PolyType& pE = decLinTran[3]; // poly represented by 0xE

  /* return ( c*pE + cf*pB + ce*pD + cd*p9
   *         + (c*p9 + cf*pE + ce*pB + cd*pD)>>1
   *         + (c*pD + cf*p9 + ce*pE + cd*pB)>>2
   *         + (c*pB + cf*pD + ce*p9 + cd*pE)>>3 )
   */

  // Compute the top row of the AES matrix
  Ctxt sum(c), tmpF(cf), tmpE(ce), tmpD(cd);
  sum.multByConstant(pE);
  tmpF.multByConstant(pB);
  tmpE.multByConstant(pD);
  tmpD.multByConstant(p9);
  sum += tmpF; sum += tmpE; sum += tmpD;

  // Compute the 2nd row of the AES matrix
  Ctxt tmp(c); tmpF = cf; tmpE = ce; tmpD = cd;
  tmp.multByConstant(p9);
  tmpF.multByConstant(pE);
  tmpE.multByConstant(pB);
  tmpD.multByConstant(pD);
  tmp += tmpF; tmp += tmpE; tmp += tmpD;
  ea2.rotate1D(tmp, 0, 5*rotAmount); // rotation + re-linearization
  sum += tmp;

  // Compute the 3rd row of the AES matrix
  tmp = c; tmpF = cf; tmpE = ce; tmpD = cd;
  tmp.multByConstant(pD);
  tmpF.multByConstant(p9);
  tmpE.multByConstant(pE);
  tmpD.multByConstant(pB);
  tmp += tmpF; tmp += tmpE; tmp += tmpD;
  ea2.rotate1D(tmp, 0, 10*rotAmount); // rotation + re-linearization
  sum += tmp;

  // Compute the bottom row of the AES matrix
  c.multByConstant(pB);
  cf.multByConstant(pD);
  ce.multByConstant(p9);
  cd.multByConstant(pE);
  c += cf; c += ce; c += cd;
  ea2.rotate1D(c, 0, 15*rotAmount); // rotation + re-linearization
  c += sum;

  c.cleanUp();
}

static void decRowShift(Ctxt& c, const vector<PolyType>& decLinTran,
			const EncryptedArrayDerived<PA_GF2>& ea2)
{
  // The basic rotation amount along the 1st dimension
  long rotAmount = ea2.getContext().zMStar.OrderOf(0) / 16;

  c.cleanUp();
  Ctxt c4(c), c8(c), c12(c);
  ea2.rotate1D(c4, 0, 4*rotAmount);
  ea2.rotate1D(c8, 0, 8*rotAmount);
  ea2.rotate1D(c12,0, 12*rotAmount);
  c4.cleanUp();  c8.cleanUp();  c12.cleanUp();

  const PolyType& p0 = decLinTran[4]; // 1
  const PolyType& p4 = decLinTran[5]; // 1>>4
  const PolyType& p8 = decLinTran[6]; // 1>>8
  const PolyType& p12 = decLinTran[7]; // 1>>12

  // return c*p0 + c4*p4 + c8*p8 + c12*p12

  c.multByConstant(p0);
  c4.multByConstant(p4); c += c4;
  c8.multByConstant(p8);  c += c8;
  c12.multByConstant(p12); c += c12;
}

// Encode AES plaintext/ciphertext bytes as native HE plaintext
void encode4AES(Vec<ZZX>& encData, const Vec<uint8_t>& data,
		const EncryptedArrayDerived<PA_GF2>& ea2)
{
  long nBlocks = divc(data.length(),16); // ceil( data.length()/16 )
  long blocksPerCtxt = ea2.size() / 16;  // = nSlots/16
  long nCtxt = divc(nBlocks, blocksPerCtxt);

  // We encode blocksPerCtxt = n/16 blocks in the slots of one ctxt.
  encData.SetLength(nCtxt);

  for (long i=0; i<nCtxt; i++) {         // i is the cipehrtext number
    // Copy the bytes into Hypercube<GF2X>'es to be used for encoding
    vector<GF2X> slots(ea2.size(), GF2X::zero());
    for (long j=0; j<blocksPerCtxt; j++) { // j is the block number in this ctxt
      long blockShift = (i*blocksPerCtxt +j)*16;  // point to block
      for (long k=0; k<16; k++) {         // k is the byte number in this block
	long byteIdx= blockShift+ k;      // column orded within block
	if (byteIdx < data.length()) {
	  long slotIdx = j + k*blocksPerCtxt;
	  GF2XFromBytes(slots[slotIdx], &data[byteIdx], 1);// copy byte as poly
	}
      }
    }
    ea2.encode(encData[i], slots);
  }
}

// Decode native HE plaintext as AES plaintext/ciphertext bytes
void decode4AES(Vec<uint8_t>& data, const Vec<ZZX>& encData,
		const EncryptedArrayDerived<PA_GF2>& ea2)
{
  // Check the size of the data array
  long nBytes = encData.length() * ea2.size(); // total number of input bytes
  if (data.length()<=0 || data.length()>nBytes)
    data.SetLength(nBytes);
  long nBlocks = divc(data.length(),16);       // ceil( data.length()/16 )
  long blocksPerCtxt = ea2.size() / 16;        // = nSlots/16
  long nCtxt = divc(nBlocks, blocksPerCtxt);   // <= encData.length()

  // We encode blocksPerCtxt = n/16 blocks in the slots of one ctxt.

  vector<GF2X> slots;
  for (long i=0; i<nCtxt; i++) {         // i is the cipehrtext number
    ea2.decode(slots, encData[i]);
    for (long j=0; j<blocksPerCtxt; j++) { // j is the block number in this ctxt
      long blockShift = (i*blocksPerCtxt +j)*16;  // point to block
      for (long k=0; k<16; k++) {         // k is the byte number in this block
	long byteIdx= blockShift +k;      // column orded within block
	if (byteIdx < data.length()) {
	  long slotIdx = j + k*blocksPerCtxt;
	  BytesFromGF2X(&data[byteIdx], slots[slotIdx], 1);// copy poly as byte
	}
      }
    }
  }
}

// the transformation X -> X^{-1} in GF(2^8)
static void invert(vector<Ctxt>& data)
{
  for (long i=0; i<(long)data.size(); i++){ // compute X -> X^{254} on i'th ctxt
    Ctxt tmp1(data[i]);           // tmp1   = data[i] = X
    tmp1.frobeniusAutomorph(1);   // tmp1   = X^2   after Z -> Z^2
    data[i].multiplyBy(tmp1);     // data[i]= X^3
    Ctxt tmp2(data[i]);           // tmp2   = X^3
    tmp2.frobeniusAutomorph(2);   // tmp2   = X^12  after Z -> Z^4
    tmp1.multiplyBy(tmp2);        // tmp1   = X^14
    data[i].multiplyBy(tmp2);     // data[i]= X^15
    data[i].frobeniusAutomorph(4);// data[i]= X^240 after Z -> Z^16
    data[i].multiplyBy(tmp1);     // data[i]= X^254
  }
}

// Pack the ciphertexts in c in as few "fully packed" cipehrtext as possible.
static void packCtxt(vector<Ctxt>& to, const vector<Ctxt>& from,
		     const GF2X& XinSlots)
{
  FHE_TIMER_START;
  if (from.size() <= 1) { // nothing to do here
    to = from; return;
  }

  // Get the context and the ea for "fully packed" polynomials
  const Context& context = from[0].getContext();
  const EncryptedArrayDerived<PA_GF2>& ea = context.ea->getDerived(PA_GF2());
  const GF2XModulus& PhimX = ea.getTab().getPhimXMod();
  long e = ea.getDegree() / 8; // the extension degree
  long nPacked = divc(from.size(), e); // How many fully-packed ciphertexts

  // Initialize the vector 'to' with empty cipehrtexts
  to.assign(nPacked, Ctxt(ZeroCtxtLike, from[0]));

  // Each ctxt in 'to' is the result of packing <= e ctxts from 'from'
  for (long i=0; i<(long) to.size(); i++) {
    to[i] = from[e*i];
    if (e*i == (long)from.size()-1) break; // check if we are done packing

    Ctxt tmp = from[e*i +1];
    tmp.multByConstant(conv<ZZX>(XinSlots));
    to[i] += tmp;
    if (e*i == (long)from.size()-2) break; // check if we are done packing

    GF2X X2i = XinSlots; // X^i, initially i=1
    for (long j= e*i +2; j<e*(i+1) && j<(long)from.size(); j++) {
      MulMod(X2i, X2i, XinSlots, PhimX);  // X^i
      tmp = from[j];
      tmp.multByConstant(conv<ZZX>(X2i));
      to[i] += tmp;
    }
  }
}

// Unpack the fully-packed ciphertext in from into the vector to. If to.size()>0
// then do not unpack into more than to.size() ciphertexts. If 'from' does not
// have enough ciphertexts to fill all of 'to' then pad with zeros.
static void unackCtxt(vector<Ctxt>& to, const vector<Ctxt>& from,
		      const Mat<GF2X>& unpackConsts)
{
  FHE_TIMER_START;
  // Get the context and the ea for "fully packed" polynomials
  const Context& context = from[0].getContext();
  const EncryptedArrayDerived<PA_GF2>& ea = context.ea->getDerived(PA_GF2());
  long e = ea.getDegree() / 8; // the extension degree
  long nUnpacked = from.size()*e; // How many lightly-packed ciphertexts
  if (to.size()==0) to.resize(nUnpacked, Ctxt(ZeroCtxtLike, from[0]));
  else {
    if (nUnpacked > (long) to.size()) nUnpacked = to.size();
    for (long i=0; i<(long)to.size(); i++) to[i].clear();
  }
  // At this point 'to' contains empty (zero) ciphertexts

  long nPacked = divc(nUnpacked, 8);
  for (long idx=0; idx<nPacked; idx++) {
    vector<Ctxt> conjugates(e, from[idx]); // Compute the conjugates, Z^{2^{8j}}
    for (long j=1; j<e; j++)
      conjugates[j].frobeniusAutomorph(8*j);

    for (long i=0; i<e && (idx*e +i)<(long)to.size(); i++) {
      // Recall that to[idx*e +i] was initialize to zero
      for (long j=0; j<e; j++) {
	Ctxt tmp = conjugates[j];
	tmp.multByConstant(conv<ZZX>(unpackConsts[i][j]));
	to[idx*e +i] += tmp;
      }
    }
  }
}
