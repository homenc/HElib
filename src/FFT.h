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
#include "NumbTh.h"
#include <NTL/FFT.h>


struct AltFFTPrimeInfo : FFTPrimeInfo  {
   bool initialized;

   AltFFTPrimeInfo() : initialized(false) { }

};

void InitAltFFTPrimeInfo(AltFFTPrimeInfo& alt_info, const FFTPrimeInfo& info, long k);

void AltFFT(long* A, const long* a, long k, const AltFFTPrimeInfo& info, long dir);


static inline
void AltFFTFwd(long* A, const long *a, long k, const AltFFTPrimeInfo& info)
{
   AltFFT(A, a, k, info, 0);
}

static inline
void AltFFTRev(long* A, const long *a, long k, const AltFFTPrimeInfo& info)
{
   AltFFT(A, a, k, info, 1);
}

static inline
void AltFFTMulTwoInv(long* A, const long *a, long k, const AltFFTPrimeInfo& info)
{
   VectorMulModPrecon(1L << k, A, a, info.TwoInvTable[k], info.q,
                      info.TwoInvPreconTable[k]);
}

static inline
void AltFFTRev1(long* A, const long *a, long k, const AltFFTPrimeInfo& info)
{
   AltFFTRev(A, a, k, info);
   AltFFTMulTwoInv(A, A, k, info);
}


