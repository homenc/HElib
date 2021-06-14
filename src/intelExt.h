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

#include <hexl/hexl.hpp>

#include <iostream>

namespace intel {

  // TODO: Create a lookup table to avoid re-creating previously created NTTs?
  // Use a hash map?
  inline intel::hexl::NTT initNTT(uint64_t degree, uint64_t q, uint64_t root)
  {
    return intel::hexl::NTT(degree, q, root);
  }

  inline void AltFFTFwd(long* output, const long* input, long n, long q, long root)//, const helib::FFTPrimeInfo& info)
  {
    std::cout << "HEXL Fwd\n";
    initNTT(/*degree=*/n, /*modulus=*/q, /*root=*/root)
      .ComputeForward(reinterpret_cast<uint64_t*>(output), reinterpret_cast<const uint64_t*>(input), /*input_mod_factor=*/4, /*output_mod_factor=*/4);
    return;
  }

  inline void AltFFTRev1(long* output, const long* input, long n, long q, long root)//, const helib::FFTPrimeInfo& info)
  {
    std::cout << "HEXL inverse\n";
    initNTT(/*degree=*/n, /*modulus=*/q, /*root=*/root)
      .ComputeInverse(reinterpret_cast<uint64_t*>(output), reinterpret_cast<const uint64_t*>(input), /*input_mod_factor=*/2, /*output_mod_factor=*/2);
    return;
  }

} // namespace intel

#endif // ifndef HELIB_INTELEXT_H
