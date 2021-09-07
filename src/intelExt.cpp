// Copyright (C) 2021 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

/*
 * Some simple wrappers to enable HElib to use HEXL.
 */

#ifdef USE_INTEL_HEXL

#include "intelExt.h"

#include <hexl/hexl.hpp>

#include <unordered_map>
#include <functional>
#include <shared_mutex>
#include <mutex>

namespace intel {

using intel::hexl::NTT;

static std::shared_mutex table_mutex;

struct Key
{
  uint64_t degree;
  uint64_t q;

  bool operator==(const Key& other) const
  {
    return this->degree == other.degree && this->q == other.q;
  }
};

struct Hash
{
  size_t operator()(const Key& key) const
  {
    // Using the popular hash_combine method of boost.
    return key.degree ^
           (key.q + 0x9e3779b9 + (key.degree << 6) + (key.degree >> 2));
  }
};

// Lookup table to avoid re-creating previously created NTTs
// For life of program.
static std::unordered_map<Key, NTT, Hash> table;

static NTT& initNTT(uint64_t degree, uint64_t q)
{
  Key key = {degree, q};

  { // First Read to see if there is an NTT
    std::shared_lock read_lock(table_mutex);
    auto it = table.find(key);
    if (it != table.end()) { // Found
      return it->second;
    }
  }

  { // Didn't find an NTT, lock the table for writing
    std::scoped_lock write_lock(table_mutex);
    // Check again to see another thread snuck it in.
    auto it = table.find(key);
    if (it != table.end()) { // Found
      return it->second;
    } else {
      auto ret = table.emplace(key, NTT(degree, q));
      return (ret.first)->second; // The NTT object just created.
    }
  }
}

void FFTFwd(long* output, const long* input, long n, long q)
{
  initNTT(/*degree=*/n, /*modulus=*/q)
      .ComputeForward(reinterpret_cast<uint64_t*>(output),
                      reinterpret_cast<const uint64_t*>(input),
                      /*input_mod_factor=*/1,
                      /*output_mod_factor=*/1);
  return;
}

void FFTRev1(long* output, const long* input, long n, long q)
{
  initNTT(/*degree=*/n, /*modulus=*/q)
      .ComputeInverse(reinterpret_cast<uint64_t*>(output),
                      reinterpret_cast<const uint64_t*>(input),
                      /*input_mod_factor=*/1,
                      /*output_mod_factor=*/1);
  return;
}

void EltwiseAddMod(long* result,
                   const long* operand1,
                   const long* operand2,
                   long n,
                   long modulus)
{
  intel::hexl::EltwiseAddMod(reinterpret_cast<uint64_t*>(result),
                             reinterpret_cast<const uint64_t*>(operand1),
                             reinterpret_cast<const uint64_t*>(operand2),
                             n,
                             modulus);
}

void EltwiseAddMod(long* result,
                   const long* operand,
                   long scalar,
                   long n,
                   long modulus)
{
  intel::hexl::EltwiseAddMod(reinterpret_cast<uint64_t*>(result),
                             reinterpret_cast<const uint64_t*>(operand),
                             scalar,
                             n,
                             modulus);
}

void EltwiseSubMod(long* result,
                   const long* operand1,
                   const long* operand2,
                   long n,
                   long modulus)
{
  intel::hexl::EltwiseSubMod(reinterpret_cast<uint64_t*>(result),
                             reinterpret_cast<const uint64_t*>(operand1),
                             reinterpret_cast<const uint64_t*>(operand2),
                             n,
                             modulus);
}

void EltwiseSubMod(long* result,
                   const long* operand,
                   long scalar,
                   long n,
                   long modulus)
{
  intel::hexl::EltwiseSubMod(reinterpret_cast<uint64_t*>(result),
                             reinterpret_cast<const uint64_t*>(operand),
                             scalar,
                             n,
                             modulus);
}

void EltwiseMultMod(long* result,
                    const long* operand1,
                    const long* operand2,
                    long n,
                    long modulus)
{
  intel::hexl::EltwiseMultMod(reinterpret_cast<uint64_t*>(result),
                              reinterpret_cast<const uint64_t*>(operand1),
                              reinterpret_cast<const uint64_t*>(operand2),
                              n,
                              modulus,
                              /*input_mod_factor=*/1);
}

void EltwiseMultMod(long* result,
                    const long* operand,
                    long scalar,
                    long n,
                    long modulus)
{

  intel::hexl::EltwiseFMAMod(reinterpret_cast<uint64_t*>(result),
                             reinterpret_cast<const uint64_t*>(operand),
                             scalar,
                             nullptr,
                             n,
                             modulus,
                             /*input_mod_factor=*/1);
}

} // namespace intel

#endif // USE_INTEL_HEXL
