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
#ifndef _DEBUGGING_H_
#define _DEBUGGING_H_
//! @file debugging.h
//! @brief debugging utilities
#include <iostream>
#include <NTL/ZZ.h>
#include "FHE.h"
#include "EncryptedArray.h"

#define FLAG_PRINT_ZZX  1
#define FLAG_PRINT_POLY 2
#define FLAG_PRINT_VEC  4

extern FHESecKey* dbgKey;
extern EncryptedArray* dbgEa;
extern NTL::ZZX dbg_ptxt;
extern NTL::Vec<NTL::ZZ> ptxt_pwr; // powerful basis

void decryptAndPrint(std::ostream& s, const Ctxt& ctxt, const FHESecKey& sk,
		     const EncryptedArray& ea, long flags=0);

bool decryptAndCompare(const Ctxt& ctxt, const FHESecKey& sk,
		       const EncryptedArray& ea, const NewPlaintextArray& pa);
#endif // _DEBUGGING_H_
