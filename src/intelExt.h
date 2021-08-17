// Copyright (C) 2021 Intel Corporation
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// You initialize a PrimeGenerator as follows:
//    PrimeGenerator gen(len, m);
// Each call to gen.next() generates a prime p with
// (1-1/2^B)*2^len <= p < 2^len and p = 2^k*t*m + 1,
// where t is odd and k is as large as possible
// and B is a small constant (typically, B in {2,3,4}).
// If no such prime is found, then an error is raised.


#ifndef HELIB_INTELEXT_H
#define HELIB_INTELEXT_H

#include <stdint.h>

namespace intel {

// Simple wrappers around HEXL
//void FFTFwd(long* output, const long* input, long n, long q, long root);
//void FFTRev1(long* output, const long* input, long n, long q, long root);
void FFTFwd(long* output, const long* input, long n, long q);
void FFTRev1(long* output, const long* input, long n, long q);

void EltwiseAddMod(long* result,
                   const long* operand1,
                   const long* operand2,
                   long n,
                   long modulus);
void EltwiseAddMod(long* result,
                   const long* operand,
                   long scalar,
                   long n,
                   long modulus);

void EltwiseSubMod(long* result,
                   const long* operand1,
                   const long* operand2,
                   long n,
                   long modulus);
void EltwiseSubMod(long* result,
                   const long* operand,
                   long scalar,
                   long n,
                   long modulus);

void EltwiseMultMod(long* result,
                    const long* operand1,
                    const long* operand2,
                    long n,
                    long modulus);
void EltwiseMultMod(long* result,
                    const long* operand,
                    long scalar,
                    long n,
                    long modulus);

} // namespace intel

#endif // HELIB_INTELEXT_H
