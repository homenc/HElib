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

static constexpr cx_double the_imaginary_i = cx_double(0.0, 1.0);

void EncryptedArrayCx::decrypt(const Ctxt& ctxt,
                               const FHESecKey& sKey, vector<cx_double>& ptxt) const
{
  assert(&getContext() == &ctxt.getContext());
  NTL::ZZX pp;
  sKey.Decrypt(pp, ctxt);

  // convert to zzX, if the pp is too big, scale it down
  long nBits = NTL::MaxBits(pp) - NTL_SP_NBITS;
  zzX zpp(INIT_SIZE, deg(pp)+1);
  double factor;
  if (nBits<=0) { // convert to zzX, double
    for (long i=0; i<lsize(zpp); i++)
      conv(zpp[i], pp[i]);
    factor = to_double(ctxt.getRatFactor());
  } else { // scale and then convert to zzX, double
    for (long i=0; i<lsize(zpp); i++)
      conv(zpp[i], pp[i]>>nBits);
    factor = to_double(ctxt.getRatFactor()/power2_xdouble(nBits));
  }
  canonicalEmbedding(ptxt, zpp, getPAlgebra()); // decode without scaling
  for (cx_double& cx : ptxt)  // divide by the factor
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
         // if precision==0 use the default PAlgebraCx::encodeScalingFactor()
  embedInSlots(ptxt, array, getPAlgebra(), factor);
}

void EncryptedArrayCx::encode(zzX& ptxt, double num, long precision) const
{
  // This factor ensures that encode/decode introduce less than
  // 1/precision error. If precision=0 then the scaling factor defaults
  // to PAlgebraCx::encodeScalingFactor(), corresponding to precision
  // error bound of 2^{-almod.getR()}
  num *= alMod.encodeScalingFactor(precision);

  resize(ptxt, 1, long(round(num))); // Constant polynomial
}

void EncryptedArrayCx::encodei(zzX& ptxt, long precision) const
{
  vector<cx_double> v(size(), the_imaginary_i); // i in all the slots
  this->encode(ptxt, v, precision);
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

void EncryptedArrayCx::extractRealPart(Ctxt& c) const
{
  Ctxt tmp = c;
  tmp.complexConj(); // the complex conjugate of c
  c += tmp;          // c + conj(c) = 2*real(c)
  c.multByConstantCKKS(0.5); // divide by two
}

// Note: If called with dcrt==nullptr, it will perform FFT's when
// encoding i as a DoubleCRT object. If called with dcrt!=nullptr,
// it assumes that dcrt points to an object that encodes i. If the
// primeSet of the given DoubleCRT is missing some of the moduli in
// c.getPrimeSet(), many extra FFTs/iFFTs will be called.
void EncryptedArrayCx::extractImPart(Ctxt& c, DoubleCRT* iDcrtPtr) const
{
  DoubleCRT tmpDcrt(getContext(), IndexSet::emptySet());
  {Ctxt tmp = c;
  c.complexConj(); // the complex conjugate of c
  c -= tmp;        // conj(c) - c = -2*i*imaginary(c)
  }
  if (iDcrtPtr==nullptr) { // Need to encode i in a DoubleCRt object
    tmpDcrt.addPrimes(c.getPrimeSet());
    tmpDcrt.FFT(getiEncoded(), c.getPrimeSet());
    // FFT is a low-level DoubleCRT procedure to initialize an
    // existing object with a given PrimeSet and a given polynomial
    iDcrtPtr = &tmpDcrt;
  }
  c.multByConstantCKKS(*iDcrtPtr); // multiply by i
  c.multByConstantCKKS(0.5);       // divide by two
}

void EncryptedArrayCx::buildLinPolyCoeffs(vector<zzX>& C,
              const cx_double& oneImage, const cx_double& iImage) const
{
  resize(C,2); // allocate space

  // Compute the constants x,y such that L(z) = x*z + y*conjugate(z)
  cx_double x = (oneImage - the_imaginary_i*iImage)*0.5;
  cx_double y = (oneImage + the_imaginary_i*iImage)*0.5;

  // Encode x,y in zzX objects
  long n = size();
  vector<cx_double> v(n, x); // x in all the slots
  encode(C[0], v);
  v.assign(n, y);            // y in all the slots
  encode(C[1], v);
}

void EncryptedArrayCx::buildLinPolyCoeffs(vector<zzX>& C,
     const vector<cx_double>&oneImages, const vector<cx_double>&iImages) const
{
  resize(C,2); // allocate space

  // Compute the constants x,y such that L(z) = x*z + y*conjugate(z)
  vector<cx_double> x(size());
  vector<cx_double> y(size());
  for (long j=0; j<size(); j++) {
    x[j] = (oneImages[j] - the_imaginary_i*iImages[j])*0.5;
    y[j] = (oneImages[j] + the_imaginary_i*iImages[j])*0.5;
  }
  // Encode x,y in zzX objects
  encode(C[0], x);
  encode(C[1], y);
}
