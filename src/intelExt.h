/* Copyright (C) 2021 Intel Corporation
 * SPDX-License-Identifier: Apache-2.0
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

#ifndef HELIB_INTELEXT_H
#define HELIB_INTELEXT_H

#include <stdint.h>

namespace intel {

// Simple wrappers around HEXL 
void FFTFwd(long* output, const long* input, long n, long q, long root);
void FFTRev1(long* output, const long* input, long n, long q, long root);

void EltwiseAddMod(long* result, const long* operand1,
                    const long* operand2, long n, long modulus);
void EltwiseAddMod(long* result, const long* operand,
                    long scalar, long n, long modulus);

void EltwiseSubMod(long* result, const long* operand1,
                    const long* operand2, long n, long modulus);
void EltwiseSubMod(long* result, const long* operand,
                    long scalar, long n, long modulus);

void EltwiseMultMod(long* result, const long* operand1,
                    const long* operand2, long n, long modulus);
void EltwiseMultMod(long* result, const long* operand,
                    long scalar, long n, long modulus);

} // namespace intel

#endif // HELIB_INTELEXT_H
