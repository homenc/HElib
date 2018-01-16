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
#ifndef _BINARY_COMPARE_H_
#define _BINARY_COMPARE_H_
/**
 * @file binaryCompare.h
 * @brief Implementing integer comparison in binary representation.
 **/
#include "EncryptedArray.h"
#include "CtPtrs.h" //  defines CtPtrs, CtPtrMat

//! Compares two integers in binary a,b.
//! Returns max(a,b), min(a,b) and indicator bits mu=(a>b) and ni=(a<b)
void compareTwoNumbers(CtPtrs& max, CtPtrs& min, Ctxt& mu, Ctxt& ni,
                       const CtPtrs& a, const CtPtrs& b,
                       std::vector<zzX>* unpackSlotEncoding=nullptr);

#endif // ifdef _BINARY_COMPARE_H_
