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
/* EaCx.cpp - Encoding/decoding and data-movement for encrypted complex data
 */
#include <algorithm>
#include "zzX.h"
#include "EncryptedArray.h"

#include "timing.h"
#include "cloned_ptr.h"
#include "norms.h"

NTL_CLIENT

void EncryptedArrayCx::decrypt(const Ctxt& ctxt,
                               const FHESecKey& sKey, vector<cx_double>& ptxt) const
{
  assert(&getContext() == &ctxt.getContext());
  NTL::ZZX pp;
  sKey.Decrypt(pp, ctxt);

  // convert to zzX, if the pp is too big, scale it down
  long nBits = NTL::MaxBits(pp);
  zzX zpp(INIT_SIZE, deg(pp)+1);
  for (long i=0; i<lsize(zpp); i++) {
    if (nBits<=NTL_SP_NBITS)
      conv(zpp[i], pp[i]);
    else
      conv(zpp[i], pp[i]>>(nBits-NTL_SP_NBITS));
  }
  canonicalEmbedding(ptxt, zpp, getPAlgebra()); // decode without scaling

  // divide by some factor
  double factor;
  if (nBits <= NTL_SP_NBITS)
    factor = to_double(ctxt.getRatFactor());
  else
    factor = to_double(ctxt.getRatFactor()/power2_xdouble(nBits-NTL_SP_NBITS));

  for (cx_double& cx : ptxt)
    cx /= factor;
}

// rotate ciphertext in dimension 0 by amt
void EncryptedArrayCx::rotate1D(Ctxt& ctxt, long i, long amt, bool dc) const
{
  assert(&getContext() == &ctxt.getContext());
  assert(nativeDimension(i));

  const PAlgebra& palg = getPAlgebra();
  long ord = sizeOfDimension(i);
  amt %= ord;// DIRT: assumes division w/ remainder follows C++11 and C99 rules
  if (amt == 0) return;

  ctxt.smartAutomorph(palg.genToPow(i, amt));
}

// Shift k positions along the i'th dimension with zero fill.
// Negative shift amount denotes shift in the opposite direction.
void EncryptedArrayCx::shift1D(Ctxt& ctxt, long i, long k) const
{
  throw std::logic_error("EncryptedArrayCx::shift1D not implemented");
}

// We only support linear arrays for approximate numbers,
// so rotate,shift are the same as rotate1D, shift1D
void EncryptedArrayCx::rotate(Ctxt& ctxt, long amt) const
{
  rotate1D(ctxt, 0, amt, true);
}
void EncryptedArrayCx::shift(Ctxt& ctxt, long amt) const
{
  shift1D(ctxt, 0, amt);
}

void EncryptedArrayCx::encode(zzX& ptxt, const vector<cx_double>& array,
                              long precision) const
{
  // This factor ensures that encode/decode introduce less than 1/precision
  // error. If precision=0 then the error bound defaults to 2^{-almod.getR()}.  
  double factor = alMod.encodeScalingFactor(precision);
  embedInSlots(ptxt, array, getPAlgebra(), factor);
}

void EncryptedArrayCx::decode(vector<cx_double>& array, const zzX& ptxt) const
{
  canonicalEmbedding(array, ptxt, getPAlgebra());
  double factor = alMod.encodeScalingFactor();
  for (auto& x: array) x /= factor;
}

// return an array of random complex numbers in the unit circle
void EncryptedArrayCx::random(vector<cx_double>& array) const
{
  const double twoPi = 8 * std::atan(1);

  resize(array, size()); // allocate space
  for (auto& x : array) {
    long bits = NTL::RandomLen_long(32); // 32 random bits
    double r = std::sqrt(bits & 0xffff)/256.0; // sqrt(uniform[0,1])
    double theta = twoPi * ((bits>>16)& 0xffff) / 65536.0; // uniform(0,2pi)
    x = std::polar(r,theta);
  }
}
