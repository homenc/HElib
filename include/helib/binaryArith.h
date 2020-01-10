/* Copyright (C) 2012-2019 IBM Corp.
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
#ifndef HELIB_BINARYARITH_H
#define HELIB_BINARYARITH_H
/**
 * @file binaryArith.h
 * @brief Implementing integer addition, multiplication in binary representation
 **/
#include <helib/EncryptedArray.h>
#include <helib/CtPtrs.h> //  defines CtPtrs, CtPtrMat

namespace helib {

/**
 * @brief Adds two numbers in binary representation where each ciphertext of the input vector contains a bit.
 * @param sum result of the addition operation.
 * @param lhs left hand side of the addition.
 * @param rhs right hand side of the addition.
 * @param sizeLimit number of bits to compute on, taken from the least significant end.
 * @param unpackSlotEncoding vector of constants for unpacking, as used in bootstrapping.
 **/
void
addTwoNumbers(CtPtrs& sum, const CtPtrs& lhs, const CtPtrs& rhs,
              long sizeLimit=0, std::vector<zzX>* unpackSlotEncoding=nullptr);

/**
 * @brief Add together up to fifteen {0,1} integers, producing a 4-bit counter.
 * @param out 4-bit counter to be outputed.
 * @param in bits to be counted.
 * @param sizeLimit number of bits to compute on, taken from the least significant end.
 * @return number of output bits that are not identically zero (i.e. != null).
 *
 * Adding fifteen input bits, getting a 4-bit counter. Some of the
 * input pointers may be null, but output pointers must point to
 * allocated Ctxt objects. If sizeLimit<4, only that many bits are
 * computed (taken from the least significant end).
 **/
long fifteenOrLess4Four(const CtPtrs& out, const CtPtrs& in, long sizeLimit=4);

/**
 * @brief Sum an arbitrary amount of numbers in binary representation.
 * @param sum result of the summation.
 * @param numbers values of which to sum.
 * @param sizeLimit number of bits to compute on, taken from the least significant end.
 * @param unpackSlotEncoding vector of constants for unpacking, as used in bootstrapping.
 * 
 * Calculates the sum of many numbers using the 3-for-2 method.
 **/
void addManyNumbers(CtPtrs& sum, CtPtrMat& numbers, long sizeLimit=0,
                    std::vector<zzX>* unpackSlotEncoding=nullptr);

/**
 * @brief Multiply two numbers in binary representation where each ciphertext of the input vector contains a bit.
 * @param product result of the multiplication operation.
 * @param lhs left hand side of the multiplication.
 * @param rhs right hand side of the multiplication.
 * @param rhsTwosComplement flag to state the multiplier is potentially negative.
 * @param sizeLimit number of bits to compute on, taken from the least significant end.
 * @param unpackSlotEncoding vector of constants for unpacking, as used in bootstrapping.
 **/
void multTwoNumbers(CtPtrs& product, const CtPtrs& lhs, const CtPtrs& rhs,
                    bool rhsTwosComplement=false, long sizeLimit=0,
                    std::vector<zzX>* unpackSlotEncoding=nullptr);

/**
 * @brief Decrypt the binary numbers that are encrypted in eNums.
 * @param pNums vector to decrypt the binary numbers into.
 * @param eNums encrypted binary numbers of which to be decrypted.
 * @param sKey secret key used for decryption.
 * @param ea encrypted array that holds neccessary information for decryption.
 * @param twosComplement when set to true, the number to decrypt is a signed integer in 2's complement.
 * @param allSlots when set to false, return only the sub-cube with index=0 in the last
 * dimension within each ciphertext.
 *
 * The bits are encrypted in a bit-sliced manner. Namely, encNums[0]
 * contains the LSB of all the numbers, encNums[1] the next bits from
 * all, etc. If twosComplement==true then the number is interpreted as a
 * signed integer in 2's-complement representation.
 * If allSlots==false then we only return the subcube with index i=0
 * in the last dimension within each ciphertext. Namely, the bit for
 * the j'th counter is found in slot of index j*sizeOf(lastDim).
 **/
void decryptBinaryNums(std::vector<long>& pNums, const CtPtrs& eNums,
                  const SecKey& sKey, const EncryptedArray& ea,
                  bool twosComplement=false, bool allSlots=true);

/**
 * @brief Function for packed recryption to recrypt multiple numbers.
 * @param a first input of which to recrypt.
 * @param b second input of which to recrypt.
 * @param unpackSlotEncoding vector of constants for unpacking, as used in bootstrapping.
 **/
void packedRecrypt(const CtPtrs& a, const CtPtrs& b,
                   std::vector<zzX>* unpackSlotEncoding);

}
#endif // ifndef HELIB_BINARYARITH_H
