/* Copyright (C) 2012-2020 IBM Corp.
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
 * @brief Returns a number as a vector of bits with LSB on the left.
 * @param num Number to be converted.
 * @param bitSize Number of bits of the input and output.
 * @return Bit vector representation of num.
 * @note `bitSize` must be non-negative.
 **/
std::vector<long> longToBitVector(long num, long bitSize);

/**
 * @brief Implementation of `output = cond * trueValue + (1 - cond) *
 * falseValue`.
 * @param output Equal to `trueValue` in slots where `cond` is one and
 * `falseValue` in slots where `cond` is zero.
 * @param cond The condition, namely a `Ctxt` containing elements of {0,1} in
 * each slot.
 * @param trueValue Value of `output` wherever `cond` is one.
 * @param falseValue Value of `output` wherever `cond` is zero.
 * @note `trueValue`, `falseValue` and `output` must have the same size.
 **/
void binaryCond(CtPtrs& output,
                const Ctxt& cond,
                const CtPtrs& trueValue,
                const CtPtrs& falseValue);

/**
 * @brief Zeroes the slots of `binaryNums` where the corresponding slot of
 * `mask` is 0.
 * @param binaryNums Input bits on which to mask (this is done in place).
 * @param mask Encrypted mask indicating desired slots.
 **/
void binaryMask(CtPtrs& binaryNums, const Ctxt& mask);

/**
 * @brief Concatenates two binary numbers into a single `CtPtrs` object.
 * E.g. If `a=10111`, `b=00101` then `output = 1011100101`.
 * @param output Equal to the concatenation of `a` and `b`.
 * @param a First number to copy into `output`.
 * @param b Second number to concatenate to `a`.
 * @note The size of `output` must be of size `a.size() + b.size()`.
 **/
void concatBinaryNums(CtPtrs& output, const CtPtrs& a, const CtPtrs& b);

/**
 * @brief Splits a single binary number into two binary numbers `leftSplit` and
 * `rightSplit`.
 * @param leftSplit Left hand side of the split.
 * @param rightSplit Right hand side of the split.
 * @param input Binary number to be split.
 * @note The size of `leftSplit` and `rightSplit` must sum to the size of
 * `input`.
 * @note The location of the split is defined by the sizes of `leftSplit` and
 * `rightSplit`.
 **/
void splitBinaryNums(CtPtrs& leftSplit,
                     CtPtrs& rightSplit,
                     const CtPtrs& input);

/**
 * @brief Left shift `input` by `shamt`.
 * @param output Shifted result.
 * @param input The number to be shifted.
 * @param shamt The number to bits to shift by.
 * @note This is a left shift only, i.e. the bits are moved to the
 * most-significant end.
 * @note `shamt` must be positive.
 * @note The size of `output` and `input` must be the same.
 **/
void leftBitwiseShift(CtPtrs& output, const CtPtrs& input, const long shamt);

/**
 * @brief Rotate `input` by `rotamt`.
 * @param output Rotated result.
 * @param input The number to be bitwise-rotated.
 * @param rotamt The amount by which to rotate `input`.  May be negative for
 * opposite-direction rotations.
 * @note For positive `rotamt` arguments, this rotates towards the
 * most-significant end (i.e. the same direction as leftBitwiseShift).
 * @note The size of `output` and `input` must be the same.
 **/
void bitwiseRotate(CtPtrs& output, const CtPtrs& input, long rotamt);

/**
 * @brief Compute a bitwise XOR between `lhs` and `rhs`.
 * @param output Result of bitwise `lhs` XOR `rhs`.
 * @param lhs Left operand to the XOR operation.
 * @param rhs Right operand to the XOR operation.
 * @note `output`, `lhs` and `rhs` must all have the same size.
 **/
void bitwiseXOR(CtPtrs& output, const CtPtrs& lhs, const CtPtrs& rhs);

/**
 * @brief Compute a bitwise OR between `lhs` and `rhs`.
 * @param output Result of bitwise `lhs` OR `rhs`.
 * @param lhs Left operand to the OR operation.
 * @param rhs Right operand to the OR operation.
 * @note `output`, `lhs` and `rhs` must all have the same size.
 **/
void bitwiseOr(CtPtrs& output, const CtPtrs& lhs, const CtPtrs& rhs);

/**
 * @brief Compute a bitwise AND between `lhs` and `rhs`.
 * @param output Result of bitwise `lhs` AND `rhs`.
 * @param lhs Left operand to the AND operation.
 * @param rhs Right operand to the AND operation.
 * @note `output`, `lhs` and `rhs` must all have the same size.
 **/
void bitwiseAnd(CtPtrs& output, const CtPtrs& lhs, const CtPtrs& rhs);

/**
 * @brief Compute a bitwise AND between `input` and a `std::vector<long>`.
 * @param output Equal to the output of the AND operation.
 * @param input Number to AND.
 * @param mask Number to AND with `input`. This should be a vector of elements
 * of {0,1}.
 * @note The size of `output` and `input` must be the same.
 **/
void bitwiseAnd(CtPtrs& output,
                const CtPtrs& input,
                const std::vector<long> mask);

/**
 * @brief Compute a bitwise NOT of `input`.
 * @param output Result of bit-flipping `input`.
 * @param input Binary number to be bit-flipped.
 * @note The size of `output` and `input` must be the same.
 **/
void bitwiseNot(CtPtrs& output, const CtPtrs& input);

/**
 * @brief Adds two numbers in binary representation where each ciphertext of the
 * input vector contains a bit.
 * @param sum result of the addition operation.
 * @param lhs left hand side of the addition.
 * @param rhs right hand side of the addition.
 * @param sizeLimit number of bits to compute on, taken from the least
 * significant end.
 * @param unpackSlotEncoding vector of constants for unpacking, as used in
 * bootstrapping.
 **/
void addTwoNumbers(CtPtrs& sum,
                   const CtPtrs& lhs,
                   const CtPtrs& rhs,
                   long sizeLimit = 0,
                   std::vector<zzX>* unpackSlotEncoding = nullptr);

/**
 * @brief Negates a number in binary 2's complement representation.
 * @param negation Reference to the negated number that will be populated.
 * @param input Number to be negated.
 * @note `input` will be treated as a number in 2's complement.
 * @note `input` must not alias negation.
 **/
void negateBinary(CtPtrs& negation, const CtPtrs& input);

/**
 * @brief Subtracts `rhs` from `lhs` where `lhs`, `rhs` are in 2's complement.
 * @param difference Reference to the difference post subtraction.
 * @param lhs Left hand side of subtraction.
 * @param rhs Right hand side of subtraction.
 * @param unpackSlotEncoding vector of constants for unpacking, as used in
 * bootstrapping.
 * @note `lhs` and `rhs` must have the same size.
 **/
void subtractBinary(CtPtrs& difference,
                    const CtPtrs& lhs,
                    const CtPtrs& rhs,
                    std::vector<zzX>* unpackSlotEncoding = nullptr);
/**
 * @brief Add together up to fifteen {0,1} integers, producing a 4-bit counter.
 * @param out 4-bit counter to be outputted.
 * @param in bits to be counted.
 * @param sizeLimit number of bits to compute on, taken from the least
 * significant end.
 * @return number of output bits that are not identically zero (i.e. != null).
 *
 * Adding fifteen input bits, getting a 4-bit counter. Some of the
 * input pointers may be null, but output pointers must point to
 * allocated Ctxt objects. If sizeLimit<4, only that many bits are
 * computed (taken from the least significant end).
 * @note This function is currently not thread safe.
 **/
long fifteenOrLess4Four(const CtPtrs& out,
                        const CtPtrs& in,
                        long sizeLimit = 4);

/**
 * @brief Sum an arbitrary amount of numbers in binary representation.
 * @param sum result of the summation.
 * @param numbers values of which to sum.
 * @param sizeLimit number of bits to compute on, taken from the least
 * significant end.
 * @param unpackSlotEncoding vector of constants for unpacking, as used in
 * bootstrapping.
 *
 * Calculates the sum of many numbers using the 3-for-2 method.
 **/
void addManyNumbers(CtPtrs& sum,
                    CtPtrMat& numbers,
                    long sizeLimit = 0,
                    std::vector<zzX>* unpackSlotEncoding = nullptr);

/**
 * @brief Multiply two numbers in binary representation where each ciphertext of
 * the input vector contains a bit.
 * @param product result of the multiplication operation.
 * @param lhs left hand side of the multiplication.
 * @param rhs right hand side of the multiplication.
 * @param rhsTwosComplement flag to state the multiplier is potentially
 * negative.
 * @param sizeLimit number of bits to compute on, taken from the least
 * significant end.
 * @param unpackSlotEncoding vector of constants for unpacking, as used in
 * bootstrapping.
 **/
void multTwoNumbers(CtPtrs& product,
                    const CtPtrs& lhs,
                    const CtPtrs& rhs,
                    bool rhsTwosComplement = false,
                    long sizeLimit = 0,
                    std::vector<zzX>* unpackSlotEncoding = nullptr);

/**
 * @brief Decrypt the binary numbers that are encrypted in eNums.
 * @param pNums vector to decrypt the binary numbers into.
 * @param eNums encrypted binary numbers of which to be decrypted.
 * @param sKey secret key used for decryption.
 * @param ea encrypted array that holds necessary information for decryption.
 * @param twosComplement when set to true, the number to decrypt is a signed
 * integer in 2's complement.
 * @param allSlots when set to false, return only the sub-cube with index=0 in
 * the last dimension within each ciphertext.
 *
 * The bits are encrypted in a bit-sliced manner. Namely, encNums[0]
 * contains the LSB of all the numbers, encNums[1] the next bits from
 * all, etc. If twosComplement==true then the number is interpreted as a
 * signed integer in 2's-complement representation.
 * If allSlots==false then we only return the subcube with index i=0
 * in the last dimension within each ciphertext. Namely, the bit for
 * the j'th counter is found in slot of index j*sizeOf(lastDim).
 **/
void decryptBinaryNums(std::vector<long>& pNums,
                       const CtPtrs& eNums,
                       const SecKey& sKey,
                       const EncryptedArray& ea,
                       bool twosComplement = false,
                       bool allSlots = true);

/**
 * @brief Function for packed recryption to recrypt multiple numbers.
 * @param a first input of which to recrypt.
 * @param b second input of which to recrypt.
 * @param unpackSlotEncoding vector of constants for unpacking, as used in
 * bootstrapping.
 **/
void packedRecrypt(const CtPtrs& a,
                   const CtPtrs& b,
                   std::vector<zzX>* unpackSlotEncoding);

} // namespace helib
#endif // ifndef HELIB_BINARYARITH_H
