#include "NumbTh.h"
#include <NTL/FFT.h>


struct AltFFTPrimeInfo : FFTPrimeInfo  {
   bool initialized;

   AltFFTPrimeInfo() : initialized(false) { }

};

void InitAltFFTPrimeInfo(AltFFTPrimeInfo& alt_info, const FFTPrimeInfo& info, long k);

void AltFFT(long* A, const long* a, long k, const AltFFTPrimeInfo& info, long dir);
