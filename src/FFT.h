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


