/* Copyright (C) 2021 Intel Corporation
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *  http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef HELIB_MACRO_H
#define HELIB_MACRO_H

// HEXL works best with primes of fewer than 50 bits.
#ifndef HELIB_SP_NBITS
#ifdef USE_INTEL_HEXL
#define HELIB_SP_NBITS (49)
#else
#define HELIB_SP_NBITS NTL_SP_NBITS
#endif // USE_INTEL_HEXL
#endif // HELIB_SP_NBITS

#endif // HELIB_MACRO_H
