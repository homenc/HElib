#ifndef _BINARY_ARITH_H_
#define _BINARY_ARITH_H_
/* binaryArith.h
 * Implementing integer addition and multiplication in binary representation.
 */
#include "EncryptedArray.h"
#include "CtPtrs.h" //  defines CtPtrs, CtPtrMat

void addTwoNumbers(CtPtrs& sum, const CtPtrs& a, const CtPtrs& b,
                   long limitSize=0,
                   std::vector<zzX>* unpackSlotEncoding=nullptr);

// Adding fifteen input bits, getting a 4-bit counter. Some of the
// input pointers may be null, but output pointers must point to
// allocated Ctxt objects. If limitSize<4, only that many LSBs are
// computed.
// Returns number of output bits that are not identically zero.
long fifteenOrLess4Four(const CtPtrs& out, const CtPtrs& in, long limitSize=4);

// Calculates the sum of many numbers using the 3-for-2 method
void addManyNumbers(CtPtrs& sum, CtPtrMat& numbers, long sizeLimit=0,
                    std::vector<zzX>* unpackSlotEncoding=nullptr);

// Multiply two integers (i.e. an array of bits) a, b.
// Computes the pairwise products x_{i,j} = a_i * b_j
// then sums the prodcuts using the 3-for-2 method.
void multTwoNumbers(CtPtrs& product, const CtPtrs& a, const CtPtrs& b, bool bNegative);

#ifdef DEBUG_PRINTOUT
void decryptBinaryNums(vector<long>& pNums, const CtPtrs& eNums,
                  const FHESecKey& sKey, const EncryptedArray& ea,
                  bool negative=false, bool allSlots=true);
#endif

#endif // ifndef _BINARY_ARITH_H_
