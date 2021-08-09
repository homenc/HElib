// Copyright (C) 2021 Intel Corporation
// SPDX-License-Identifier: Apache-2.0

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
