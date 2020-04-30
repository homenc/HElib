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
#ifndef HELIB_TABLELOOKUP_H
#define HELIB_TABLELOOKUP_H
/**
 * @file tableLookup.h
 * @brief Code for homomorphic table lookup and fixed-point functions
 **/
#include <functional>
#include <helib/EncryptedArray.h>
#include <helib/CtPtrs.h>

namespace helib {

// inline void convert(Ctxt& a, const Ctxt& b) { a=b; }

//! For an n-size array, compute the 2^n products
//!     products[j] = \prod_{i s.t. j_i=1} array[i]
//!                   \times \prod_{i s.t. j_i=0}(a-array[i])
void computeAllProducts(CtPtrs& products,
                        const CtPtrs& array,
                        std::vector<zzX>* unpackSlotEncoding = nullptr);

//! The input is a plaintext table T[] and an array of encrypted bits
//! I[], holding the binary representation of an index i into T.
//! The output is the encrypted value T[i].
void tableLookup(Ctxt& out,
                 const std::vector<zzX>& table,
                 const CtPtrs& idx,
                 std::vector<zzX>* unpackSlotEncoding = nullptr);

//! The input is an encrypted table T[] and an array of encrypted bits
//! I[], holding the binary representation of an index i into T.
//! This function increments by one the entry T[i].
void tableWriteIn(const CtPtrs& table,
                  const CtPtrs& idx,
                  std::vector<zzX>* unpackSlotEncoding = nullptr);

/**
 * @function buildLookupTable
 * @brief Built a table-lookup for a function in fixed-point representation
 *
 * FIXED-POINT CONVENTIONS:
 * Fixed-point numbers are specified by a triple (nbits,scale,signed).
 * Such a number is represented as an integer x with nbits bits. If
 * signed == 1, then x is treated as a signed integer in 2's compliment;
 * otherwise it is as an unsigned integer. The value represented by x
 * is x*2^{scale}.
 *
 * The buildLookupTable function builds a lookup table T, which can be
 * used in conjunction with the tableLookup function above. The size of
 * T will be 2^{nbits_in}. For every signed integer x with bit-size
 * 'nbits_in', we will have
 *     T[x] = f(x * 2^{scale_in}) * 2^{-scale_out}),
 * rounded to the nearest integer and truncated to 'nbits_out' bits.
 * The bits are packed inside the slots, so it is assumed that each
 * slot has enough room to fit these many bits. (Otherwise we only
 * keep as many low-order bits as fit in a slot.)
 *
 * SATURATED ARITHMETIC:
 * Applications of f that return a result that is too large to represent
 * in the output format will be converted to the maximum representable
 * value. Similarly, Applications of f that return a result that is too
 * small will be converted to the minimal representable value. (This
 * applies also to applications of f that return infinites, NaNs will
 * just be mapped to zero.) For this to work correctly, you should be
 * working with standard IEEE arithmetic...which will be the case on
 * almost all platforms.
 *
 * EXAMPLE:
 *
 * buildLookupTable(T, [](double x){ return 1/x;},
 *                  nbits_in, scale_in, nbits_out, scale_out, sign_out, ea)
 *
 * will build a lookup table for inversion.
 **/
void buildLookupTable(std::vector<zzX>& T, // encoded result is returned in T
                      std::function<double(double)> f,
                      long nbits_in, // number of precision bits
                      long scale_in, // scaling factor
                      long sign_in,  // 1: 2's complement signed, 0: unsigned
                      long nbits_out,
                      long scale_out,
                      long sign_out,
                      const EncryptedArray& ea);

} // namespace helib

#endif // ifndef HELIB_TABLELOOKUP_H
