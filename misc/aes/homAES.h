/** homAES.h - homomorphic AES using HElib
 */
#include <stdint.h>
#include <NTL/ZZX.h>
#include <NTL/GF2X.h>
#include "EncryptedArray.h"
#include "hypercube.h"

#ifdef USE_ZZX_POLY
#define PolyType ZZX
#else
#if (ALT_CRT)
#define PolyType AltCRT
#else
#define PolyType DoubleCRT
#endif
#endif

class HomAES {
  const EncryptedArrayDerived<PA_GF2> ea2;

  vector<PolyType> encAffMat, decAffMat; // The GF2 affine map constants
  PolyType affVec;

  vector<PolyType> encLinTran, decLinTran; // The rowShift/colMix constants

  GF2X XinSlots; // "Fully packed" poly with X in all the slots, for packing
  Mat<GF2X> unpacking; // constants for unpacking after recryption

  void batchRecrypt(vector<Ctxt>& data) const; // recryption during AES computation

public:
  static const GF2X aesPoly;  // The AES polynomial: X^8+X^4+X^3+X+1

  //! Constructor. If context is bootstrappable then also
  //! the packing/unpacking constants are computed.
  explicit HomAES(const Context& context);

  //! Method for copmuting packing/unpacking constants after initialization
  void setPackingConstants();

  //! run the AES key-expansion and then encrypt the expanded key
  void encryptAESkey(vector<Ctxt>& eKey, Vec<uint8_t>& aesKey,
		     const PubKey& hePK) const;

  //! Perform AES encryption/decryption on "raw bytes" (ECB mode)
  //! The input bytes are either plaintext or AES-encrypted ciphertext
  void homAESenc(vector<Ctxt>& eData, const vector<Ctxt>& eKey,
		 const Vec<uint8_t> inBytes) const;
  void homAESdec(vector<Ctxt>& eData, const vector<Ctxt>& eKey,
		 const Vec<uint8_t> inBytes) const;

  //! In-place AES encryption/decryption on HE encrypted bytes (ECB mode)
  void homAESenc(vector<Ctxt>& eData, const vector<Ctxt>& eKey) const;
  void homAESdec(vector<Ctxt>& eData, const vector<Ctxt>& eKey) const;

  // utility functions
  const EncryptedArrayDerived<PA_GF2>& getEA() const { return ea2; }
};


// Encode/decode AES plaintext/ciphertext bytes as native HE plaintext
void encode4AES(Vec<ZZX>& encData, const Vec<uint8_t>& data,
		const EncryptedArrayDerived<PA_GF2>& ea2);
void decode4AES(Vec<uint8_t>& data, const Vec<ZZX>& encData,
		const EncryptedArrayDerived<PA_GF2>& ea2);
