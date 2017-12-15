// intraSlot.h - Packing/unpacking of mod-p integers in GF(p^d) slots.
#ifndef _INTRASLOT_H_
#define _INTRASLOT_H_

#include <NTL/BasicThreadPool.h>

#include "FHEContext.h"
#include "EncryptedArray.h"
#include "Ctxt.h"
#include "CtPtrs.h"

// Prepare the constants for unpacking
void buildUnpackSlotEncoding(std::vector<zzX>& unpackSlotEncoding,
                             const EncryptedArray& ea);

// Low-level unpack of one ciphertext using pre-copmuted constants
void unpack(const CtPtrs& unpacked, const Ctxt& packed, 
            const EncryptedArray& ea,
            const std::vector<zzX>& unpackSlotEncoding);

// unpack many ciphertexts, returns the number of unpacked ciphertexts
long unpack(const CtPtrs& unpacked, const CtPtrs& packed,
            const EncryptedArray& ea, 
            const std::vector<zzX>& unpackSlotEncoding);

// Low-level (re)pack in slots of one ciphertext
void repack(Ctxt& packed, const CtPtrs& unpacked, const EncryptedArray& ea);

// pack many ciphertexts, returns the number of packed ciphertexts
long repack(const CtPtrs& packed, const CtPtrs& unpacked, const EncryptedArray& ea);

// Returns in 'value' a vector of slot values. The bits of value[i] are
// the content of the i'th slot pa[i], as represented on the normal basis
void unpackSlots(std::vector<unsigned long>& value,
                 NewPlaintextArray& pa, const EncryptedArray& ea);


// Packs the low-order nbits of the integer 'data' into each slot,
// where the bits get mapped to coefficients on the normal basis
void packConstant(zzX& result, unsigned long data, long nbits,
                  const EncryptedArray& ea);

// Packs the low-order nbits of the entries of 'data' into each slot,
// where the bits get mapped to coefficients on the normal basis
void packConstants(zzX& result, const std::vector<unsigned long>& data,
                   long nbits, const EncryptedArray& ea);
#endif // #ifndef _INTRASLOT_H_
