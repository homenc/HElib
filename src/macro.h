#ifndef HELIB_MACRO_H
#define HELIB_MACRO_H

// HEXL works best with primes of fewer than 50 bits.
#ifdef USE_INTEL_HEXL 
#define HELIB_SP_NBITS (49)
#else 
#define HELIB_SP_NBITS NTL_SP_NBITS
#endif

#endif // HELIB_MACRO_H
