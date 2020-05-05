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
#ifndef HELIB_DEBUGGING_H
#define HELIB_DEBUGGING_H
//! @file debugging.h
//! @brief debugging utilities
#include <iostream>
#include <string>
#include <NTL/ZZX.h>
#include <helib/NumbTh.h>

#define FLAG_PRINT_ZZX 1
#define FLAG_PRINT_POLY 2
#define FLAG_PRINT_VEC 4   /* decode to ZZX */
#define FLAG_PRINT_DVEC 8  /* decode to double float */
#define FLAG_PRINT_XVEC 16 /* decode to complex numbers */

namespace helib {

// forward declarations
class Ctxt;
class SecKey;
class EncryptedArray;
class PlaintextArray;
class DoubleCRT;

extern SecKey* dbgKey;
extern std::shared_ptr<const EncryptedArray> dbgEa;
extern NTL::ZZX dbg_ptxt;

/**
 * @brief Setup function for setting up the global debug variables.
 * @note Works only if `HELIB_DEBUG` is defined. It does not do anything
 * otherwise
 */
inline void setupDebugGlobals(
    SecKey* debug_key,
    const std::shared_ptr<const EncryptedArray>& debug_ea,
    NTL::ZZX debug_ptxt = NTL::ZZX{})
{
#ifdef HELIB_DEBUG
  dbgKey = debug_key;
  dbgEa = debug_ea;
  dbg_ptxt = debug_ptxt;
#else
  (void)debug_key;
  (void)debug_ea;
  (void)debug_ptxt;
#endif
}

/**
 * @brief Cleanup function for clearing the global debug variables.
 */
inline void cleanupDebugGlobals()
{
  dbgKey = nullptr;
  dbgEa = nullptr;
  dbg_ptxt = NTL::ZZX{};
}

void decryptAndPrint(std::ostream& s,
                     const Ctxt& ctxt,
                     const SecKey& sk,
                     const EncryptedArray& ea,
                     long flags = 0);

NTL::xdouble embeddingLargestCoeff(const Ctxt& ctxt, const SecKey& sk);

double realToEstimatedNoise(const Ctxt& ctxt, const SecKey& sk);

void checkNoise(const Ctxt& ctxt,
                const SecKey& sk,
                const std::string& msg,
                double thresh = 10.0);

bool decryptAndCompare(const Ctxt& ctxt,
                       const SecKey& sk,
                       const EncryptedArray& ea,
                       const PlaintextArray& pa);

void rawDecrypt(NTL::ZZX& plaintxt,
                const std::vector<NTL::ZZX>& zzParts,
                const DoubleCRT& sKey,
                long q = 0);

// Debug printing routines for vectors, ZZX'es, print only a few entries

template <typename VEC>
std::ostream& printVec(std::ostream& s, const VEC& v, long nCoeffs = 40)
{
  long d = lsize(v);
  if (d < nCoeffs)
    return s << v; // just print the whole thing

  // otherwise print only 1st nCoeffs coefficients
  s << '[';
  for (long i = 0; i < nCoeffs - 2; i++)
    s << v[i] << ' ';
  s << "... " << v[d - 2] << ' ' << v[d - 1] << ']';
  return s;
}

inline std::ostream& printZZX(std::ostream& s,
                              const NTL::ZZX& poly,
                              long nCoeffs = 40)
{
  return printVec(s, poly.rep, nCoeffs);
}

} // namespace helib

#endif // ifndef HELIB_DEBUGGING_H
