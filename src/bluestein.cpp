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
/* bluestein.cpp -
 * An implementation of non-power-of-two FFT using Bluestein's trick
 *
 */

#include <helib/bluestein.h>
#include <helib/timing.h>
#include <helib/CModulus.h>
#include <helib/apiAttributes.h>

#define NEW_BLUE (1)

namespace helib {

/************************************************************

Victor says: I really need to document the optimizations I made
to the Bluestein logic in relation to the truncated FFT.

In the mean time, here are two emails I wrote that explain the
situation (from June 7, 2018):


What we are really computing is f*g mod x^m-1, where f and g have degree less
than m.

The way it's done now is that we replace f by a larger polynomial F, and then
compute F*g mod x^N-1, where N is the next power of 2 after 2*m-1. This is done
implicitly in NTL's FFT routine.

With the truncated FFT, a better way is as follows. Just compute the polynomial
h = f*g, which we can do with the same size FFT, but truncated to 2*m-1 terms.
Then compute h mod x^m-1 separately, which takes negligible time.


..........................................

There are some complications with the idea that I had, because of the way we
currently implement Bluestein.  For my idea to work, I need convolutions modulo
m, but that is not what we currently have.

It relates to the fact that we are working with roots of order 2*m, rather than
m.  In fact, if m is odd (which it almost always is, unless it's a power of 2),
there is no need to do this.

First, you can look here to refresh your memory on Bluestein:

https://www.dsprelated.com/freebooks/mdft/Bluestein_s_FFT_Algorithm.html

Now, the problem is that in these formulas, we work with W^{1/2}, which is the
root of order 2*m.  But if W has order m and m is itself odd, then 2 has an
inverse mod m, and so we can just work with W.  In all of the computations, we
are just computing W^{(1/2)*i} for various values of i, and so in the exponent
we're just doing arithmetic mod m.

To get the speedup using the truncated FFT, I really need these two be circular
convolutions modulo m, not modulo 2*m.

I never really understood why we needed to work with a root of order 2*m, and
now I see that we don't, at least when m is odd.


*************************************************************/

void BluesteinInit(long n,
                   const NTL::zz_p& root,
                   NTL::zz_pX& powers,
                   NTL::Vec<NTL::mulmod_precon_t>& powers_aux,
                   NTL::fftRep& Rb)
{
  long p = NTL::zz_p::modulus();

  NTL::zz_p one;
  one = 1;
  powers.SetMaxLength(n);

  long e;
  if (n % 2 == 0)
    e = 2 * n;
  else
    e = n;

  SetCoeff(powers, 0, one);
  for (long i = 1; i < n; i++) {
    long iSqr = NTL::MulMod(i, i, e);       // i^2 mod 2n
    SetCoeff(powers, i, power(root, iSqr)); // powers[i] = root^{i^2}
  }

  // powers_aux tracks powers
  powers_aux.SetLength(n);
  for (long i = 0; i < n; i++)
    powers_aux[i] = NTL::PrepMulModPrecon(rep(powers[i]), p);

  long k = NTL::NextPowerOfTwo(2 * n - 1);
  long k2 = 1L << k; // k2 = 2^k

  Rb.SetSize(k);

  NTL::zz_pX b(NTL::INIT_SIZE, k2);

  if (NEW_BLUE && n == e) {
    NTL::zz_p rInv = inv(root);
    for (long i = 0; i < n; i++) {
      long iSqr = NTL::MulMod(i, i, e); // i^2 mod 2n
      NTL::zz_p bi = power(rInv, iSqr);
      SetCoeff(b, i, bi);
    }
  } else {
    NTL::zz_p rInv = inv(root);
    SetCoeff(b, n - 1, one); // b[n-1] = 1
    for (long i = 1; i < n; i++) {
      long iSqr = NTL::MulMod(i, i, e); // i^2 mod 2n
      NTL::zz_p bi = power(rInv, iSqr);
      // b[n-1+i] = b[n-1-i] = root^{-i^2}
      SetCoeff(b, n - 1 + i, bi);
      SetCoeff(b, n - 1 - i, bi);
    }
  }

  TofftRep(Rb, b, k);
}

void BluesteinFFT(NTL::zz_pX& x,
                  long n,
                  UNUSED const NTL::zz_p& root,
                  const NTL::zz_pX& powers,
                  const NTL::Vec<NTL::mulmod_precon_t>& powers_aux,
                  const NTL::fftRep& Rb)
{
  HELIB_TIMER_START;

  if (IsZero(x))
    return;
  if (n <= 0) {
    clear(x);
    return;
  }

  long p = NTL::zz_p::modulus();

  long dx = deg(x);
  for (long i = 0; i <= dx; i++) {
    x[i].LoopHole() =
        NTL::MulModPrecon(rep(x[i]), rep(powers[i]), p, powers_aux[i]);
  }
  x.normalize();

  long k = NTL::NextPowerOfTwo(2 * n - 1);
  NTL::fftRep& Ra = Cmodulus::getScratch_fftRep(k);

  // Careful! we are multiplying polys of degrees 2*(n-1)
  // and (n-1) modulo x^k-1.  This gives us some
  // truncation in certain cases.

  if (NEW_BLUE && n % 2 != 0) {
    TofftRep_trunc(Ra, x, k, 2 * n - 1);

    mul(Ra, Ra, Rb); // multiply in FFT representation

    FromfftRep(x, Ra, 0, 2 * (n - 1)); // then convert back
    dx = deg(x);
    if (dx >= n) {
      // reduce mod x^n-1
      for (long i = n; i <= dx; i++) {
        x[i - n].LoopHole() = NTL::AddMod(rep(x[i - n]), rep(x[i]), p);
      }
      x.SetLength(n);
      x.normalize();
      dx = deg(x);
    }

    for (long i = 0; i <= dx; i++) {
      x[i].LoopHole() =
          NTL::MulModPrecon(rep(x[i]), rep(powers[i]), p, powers_aux[i]);
    }
    x.normalize();
  } else {
    TofftRep_trunc(Ra, x, k, 3 * (n - 1) + 1);

    mul(Ra, Ra, Rb); // multiply in FFT representation

    FromfftRep(x, Ra, n - 1, 2 * (n - 1)); // then convert back
    dx = deg(x);
    for (long i = 0; i <= dx; i++) {
      x[i].LoopHole() =
          NTL::MulModPrecon(rep(x[i]), rep(powers[i]), p, powers_aux[i]);
    }
    x.normalize();
  }
}

} // namespace helib
