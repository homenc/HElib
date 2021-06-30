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

#ifdef USE_INTEL_HEXL

#include "intelExt.h"

#include <hexl/hexl.hpp>

#include <iostream>

namespace intel {

// TODO: Create a lookup table to avoid re-creating previously created NTTs?
// Use a hash map?
intel::hexl::NTT initNTT(uint64_t degree, uint64_t q, uint64_t root)
{
  return intel::hexl::NTT(degree, q, root);
}

void FFTFwd(long* output,
               const long* input,
               long n,
               long q, 
               long root)//, const helib::FFTPrimeInfo& info)
{
  initNTT(/*degree=*/n, /*modulus=*/q, root)
      .ComputeForward(reinterpret_cast<uint64_t*>(output),
                      reinterpret_cast<const uint64_t*>(input),
                      /*input_mod_factor=*/1,
                      /*output_mod_factor=*/1);
  return;
}

void FFTRev1(long* output,
                const long* input,
                long n,
                long q, 
                long root)//, const helib::FFTPrimeInfo& info)
{
  initNTT(/*degree=*/n, /*modulus=*/q, root)
      .ComputeInverse(reinterpret_cast<uint64_t*>(output),
                      reinterpret_cast<const uint64_t*>(input),
                      /*input_mod_factor=*/1,
                      /*output_mod_factor=*/1);
  return;
}

void BitReverseCopy(long* B, const long* A, long bit_width) 
{   
  // API as close to current HElib equiv API.
  for(long i = 0; i < (1L << bit_width); ++i){
    B[intel::hexl::ReverseBits(i, bit_width)] = A[i];
  }
}

} // namespace intel

#endif // USE_INTEL_HEXL