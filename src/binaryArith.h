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
#ifndef _BINARY_ARITH_H_
#define _BINARY_ARITH_H_
/**
 * @file binaryArith.h
 * @brief Implementing integer addition, multiplication in binary representation
 **/
#include "EncryptedArray.h"
#include "CtPtrs.h" //  defines CtPtrs, CtPtrMat

//! Add two integers (i.e. two array of bits) a, b.
void
addTwoNumbers(CtPtrs& sum, const CtPtrs& a, const CtPtrs& b,
              long sizeLimit=0, std::vector<zzX>* unpackSlotEncoding=nullptr);

//! Adding fifteen input bits, getting a 4-bit counter. Some of the
//! input pointers may be null, but output pointers must point to
//! allocated Ctxt objects. If sizeLimit<4, only that many LSBs are
//! computed.
//! Returns number of output bits that are not identically zero.
long fifteenOrLess4Four(const CtPtrs& out, const CtPtrs& in, long sizeLimit=4);

//! Calculate the sum of many numbers using the 3-for-2 method
void addManyNumbers(CtPtrs& sum, CtPtrMat& numbers, long sizeLimit=0,
                    std::vector<zzX>* unpackSlotEncoding=nullptr);

//! Multiply two integers (i.e. two array of bits) a, b.
void multTwoNumbers(CtPtrs& product, const CtPtrs& a, const CtPtrs& b,
                    bool bNegative=false, long sizeLimit=0,
                    std::vector<zzX>* unpackSlotEncoding=nullptr);

//! Decrypt the binary numbers that are encrypted in eNums.
void decryptBinaryNums(vector<long>& pNums, const CtPtrs& eNums,
                  const FHESecKey& sKey, const EncryptedArray& ea,
                  bool negative=false, bool allSlots=true);
// The bits are encrypted in a bit-sliced manner. Namely, encNums[0]
// contains the LSB of all the numbers, encNums[1] the next bits from
// all, etc. If negative==true then the number is interpreted as a
// signed integer in 2's-complement representation.
// If allSlots==false then we only return the subcube with index i=0
// in the last dimension within each ciphertext. Namely, the bit for
// the j'th counter is found in slot of index j*sizeOf(lastDim).

//! Use packed bootstrapping, so we can bootstrap all in just one go.
void packedRecrypt(const CtPtrs& a, const CtPtrs& b,
                   std::vector<zzX>* unpackSlotEncoding);

#endif // ifndef _BINARY_ARITH_H_
