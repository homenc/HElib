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
#ifndef HELIB_BINARYCOMPARE_H
#define HELIB_BINARYCOMPARE_H
/**
 * @file binaryCompare.h
 * @brief Implementing integer comparison in binary representation.
 **/
#include <helib/EncryptedArray.h>
#include <helib/CtPtrs.h> //  defines CtPtrs, CtPtrMat

namespace helib {

/**
 * @brief Compares two integers in binary `a`, `b`. Returns `max(a, b)`, `min(a,
 *b)` and indicator bits `mu`=(`a`>`b`) and `ni`=(`a`<`b`)
 * @param max Maximum of `a` and `b`.
 * @param min Minimum of `a` and `b`.
 * @param mu Indicator bits `mu`=(`a`>`b`).
 * @param ni Indicator bits `ni`=(`a`<`b`).
 * @param a First number to compare.
 * @param b Second number to compare.
 * @param twosComplement When set to `true`, the inputs are signed integers in
 *2's complement. If set to `false` (default), unsigned comparison is performed.
 * @param unpackSlotEncoding Vector of constants for unpacking, as used in
 *bootstrapping.
 * @note If `a`=`b` then `mu`=`ni`=`0`
 **/
void compareTwoNumbers(CtPtrs& max,
                       CtPtrs& min,
                       Ctxt& mu,
                       Ctxt& ni,
                       const CtPtrs& a,
                       const CtPtrs& b,
                       bool twosComplement = false,
                       std::vector<zzX>* unpackSlotEncoding = nullptr);

/**
 * @brief Compares two integers in binary `a`, `b`. Returns only indicator bits
 * `mu`=(`a`>`b`) and `ni`=(`a`<`b`).
 * @param mu Indicator bits `mu`=(`a`>`b`).
 * @param ni Indicator bits `ni`=(`a`<`b`).
 * @param a First number to compare.
 * @param b Second number to compare.
 * @param twosComplement When set to `true`, the inputs are signed integers in
 *2's complement. If set to `false` (default), unsigned comparison is performed.
 * @param unpackSlotEncoding Vector of constants for unpacking, as used in
 *bootstrapping.
 * @note If `a`=`b` then `mu`=`ni`=`0`
 **/
void compareTwoNumbers(Ctxt& mu,
                       Ctxt& ni,
                       const CtPtrs& a,
                       const CtPtrs& b,
                       bool twosComplement = false,
                       std::vector<zzX>* unpackSlotEncoding = nullptr);

} // namespace helib
#endif // ifndef HELIB_BINARYCOMPARE_H
