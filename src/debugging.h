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
#ifndef HELIB_DEBUGGING_H
#define HELIB_DEBUGGING_H
//! @file debugging.h
//! @brief debugging utilities
#include <iostream>
#include <string>
#include <NTL/ZZX.h>
#include "NumbTh.h"

#define FLAG_PRINT_ZZX  1
#define FLAG_PRINT_POLY 2
#define FLAG_PRINT_VEC  4 /* decode to ZZX */
#define FLAG_PRINT_DVEC 8 /* decode to double float */
#define FLAG_PRINT_XVEC 16 /* decode to complex numbers */

// forward declerations
class Ctxt;
class FHESecKey;
class EncryptedArray;
class PlaintextArray;
class DoubleCRT;

extern FHESecKey* dbgKey;
extern EncryptedArray* dbgEa;
extern NTL::ZZX dbg_ptxt;
extern NTL::Vec<NTL::ZZ> ptxt_pwr; // powerful basis

// Cleanup method for using the above debug variables several times in one process
inline void cleanupGlobals()
{
  dbgKey = nullptr;
  dbgEa = nullptr;
  dbg_ptxt = NTL::ZZX{};
  ptxt_pwr = NTL::Vec<NTL::ZZ>{};
}

void decryptAndPrint(std::ostream& s, const Ctxt& ctxt, const FHESecKey& sk,
		     const EncryptedArray& ea, long flags=0);

NTL::xdouble embeddingLargestCoeff(const Ctxt& ctxt, const FHESecKey& sk);

double realToEstimatedNoise(const Ctxt& ctxt, const FHESecKey& sk);

void checkNoise(const Ctxt& ctxt, const FHESecKey& sk, const std::string& msg, double thresh=10.0);

bool decryptAndCompare(const Ctxt& ctxt, const FHESecKey& sk,
		       const EncryptedArray& ea, const PlaintextArray& pa);

void rawDecrypt(NTL::ZZX& plaintxt, const std::vector<NTL::ZZX>& zzParts,
                const DoubleCRT& sKey, long q=0);

// Debug printing routines for vectors, ZZX'es, print only a few entries

template<class VEC> std::ostream&
printVec(std::ostream& s, const VEC& v, long nCoeffs=40)
{
  long d = lsize(v);
  if (d<nCoeffs) return s << v; // just print the whole thing

  // otherwise print only 1st nCoeffs coefficiants
  s << '[';
  for (long i=0; i<nCoeffs-2; i++) s << v[i] << ' ';
  s << "... " << v[d-2] << ' ' << v[d-1] << ']';
  return s;
}

inline std::ostream&
printZZX(std::ostream& s, const NTL::ZZX& poly, long nCoeffs=40)
{
  return printVec(s, poly.rep, nCoeffs);
}


#endif // ifndef HELIB_DEBUGGING_H
