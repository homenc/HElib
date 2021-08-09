// Copyright (C) 2021 Intel Corporation
// SPDX-License-Identifier: Apache-2.0
    

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
