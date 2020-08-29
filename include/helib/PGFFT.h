
/**************************************************************************


PGFFT: Pretty Good FFT (v1.8)

Copyright (C) 2019, Victor Shoup

See below for more details.

**************************************************************************/

#if __cplusplus <  201103L
#error "C++11 required to compile PGFFT"
#endif



#ifndef PGFFT_H
#define PGFFT_H

#include <vector>
#include <complex>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"

namespace helib {

class PGFFT {
public:

   // initialize data strctures for n-point FFT
   // REQUIREMENT: n > 0
   explicit PGFFT(long n);

   void apply(const std::complex<double>* src, std::complex<double>* dst) const;
   // Apply n-point FFT to src[0..n-1], storing result in dst[0..n-1].
   // That is,
   //   dst[i] = \sum_{j=0}^{n-1} src[j] W^{ij}
   // and where W is the nth root of unity polar(1, -2*pi/n).
   // src and dst may be equal, but should not otherwise overlap

   void apply(std::complex<double>* v) const { apply(v, v); }
   // same as apply(v, v)

   // Copy/move constructors/assignment ops deleted, as future implementations
   // may not support them.
   PGFFT(const PGFFT&) = delete;
   PGFFT(PGFFT&&) = delete;
   PGFFT& operator=(const PGFFT&) = delete;
   PGFFT& operator=(PGFFT&&) = delete;

   static bool
   simd_enabled();

   // Define aligned vectors.
   // This stuff should probably be private, but I'm not sure
   // of the rules for accessing these in the implementation file.

   static void *
   aligned_allocate(std::size_t n, std::size_t nelts);

   static void
   aligned_deallocate(void *p);

   template <class T>
   class aligned_allocator
   {
   public:
       using value_type    = T;

       aligned_allocator() noexcept {}
       template <class U> aligned_allocator(aligned_allocator<U> const&) noexcept {}

       value_type*
       allocate(std::size_t n)
       {
           void *p = aligned_allocate(n, sizeof(T));
           if (p) {
	      return static_cast<value_type*>(p);
           }

           throw std::bad_alloc();
       }

       void
       deallocate(value_type* p, std::size_t) noexcept
       {
	   aligned_deallocate(p);
       }

       template <class U>
       bool
       operator==(aligned_allocator<U> const&) noexcept
       {
	   return true;
       }

       template <class U>
       bool
       operator!=(aligned_allocator<U> const& y) noexcept
       {
	   return false;
       }
   };

   template<class T>
   using aligned_vector= std::vector<T,aligned_allocator<T>>;



private:

   long n;
   long k;

   long strategy;

   // holds all of the twiddle factors
   std::vector<aligned_vector<std::complex<double>>> tab;

   // additional data structures needed for Bluestein
   aligned_vector<std::complex<double>> powers;
   aligned_vector<std::complex<double>> Rb;

   // additonal data structures needed for 2^k-point FFT
   std::vector<long> rev, rev1;





};

}

#endif

/**************************************************************************

PGFFT: Pretty Good FFT (v1.8)

Copyright (C) 2019, Victor Shoup

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

****************************************************************************

The logic of this code is derived from code originally developed by David Harvey,
even though the code itself has been essentially rewritten from scratch.
Here is David Harvey's original copyright notice.

fft62: a library for number-theoretic transforms

Copyright (C) 2013, David Harvey

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

****************************************************************************/

#pragma GCC diagnostic pop
