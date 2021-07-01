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

/*
  Some simple wrappers to enable HElib to use HEXL.
  For HEXL's NTTs we make use of roots that HElib have already computated for
  the qs.
*/

#ifdef USE_INTEL_HEXL

#include "intelExt.h"

#include <hexl/hexl.hpp>

#include <iostream>
#include <unordered_map>
#include <functional>
#include <mutex>

namespace intel {

using intel::hexl::NTT;

std::mutex table_mutex;

struct Key {
  uint64_t degree;
  uint64_t q;
  uint64_t root;

  bool operator==(const Key& other) const
  {
    return this->degree == other.degree &&
           this->q == other.q &&
           this->root == other.root;
  }
};

struct Hash {      
  size_t operator()(const Key& key) const 
  {
    return std::hash<uint64_t>{}(key.degree) ^
           (std::hash<uint64_t>{}(key.q) ^
           (std::hash<uint64_t>{}(key.root) << 1) << 1 );
  }
};

// Lookup table to avoid re-creating previously created NTTs
// For life of program.
static std::unordered_map<Key, NTT, Hash> table;

static NTT& initNTT(uint64_t degree, uint64_t q, uint64_t root)
{       
  Key key = {degree, q, root};
  auto it = table.find(key);
  if(it != table.end()) {
    return it->second;
  }
  else {
    // Lock the table for writing
    std::scoped_lock table_lock(table_mutex);
    auto ret = table.emplace(key, NTT(degree, q, root));
    return (ret.first)->second; // The NTT object just created.
  }
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
                      /*input_mod_factor=*/4,
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
                      /*input_mod_factor=*/2,
                      /*output_mod_factor=*/1);
  return;
}

void EltwiseMultMod(long* result, const long* operand1,
                    const long* operand2, long n, long modulus,
                    long input_mod_factor) 
{
  intel::hexl::EltwiseMultMod(reinterpret_cast<uint64_t*>(result), 
                        reinterpret_cast<const uint64_t*>(operand1),
                        reinterpret_cast<const uint64_t*>(operand2), 
                        n, 
                        modulus, 
                        input_mod_factor);
}


} // namespace intel

#endif // USE_INTEL_HEXL
