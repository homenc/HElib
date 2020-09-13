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
/* EaCx.cpp - Encoding/decoding and data-movement for encrypted complex data
 */
#include <algorithm>
#include <type_traits>

#include <helib/zzX.h>
#include <helib/EncryptedArray.h>

#include <helib/timing.h>
#include <helib/clonedPtr.h>
#include <helib/norms.h>
#include <helib/debugging.h>
#include <helib/apiAttributes.h>

namespace helib {

static constexpr cx_double the_imaginary_i = cx_double(0.0, 1.0);

void EncryptedArrayCx::decrypt(const Ctxt& ctxt,
                               const SecKey& sKey,
                               std::vector<cx_double>& ptxt) const
{
  assertEq(&getContext(),
           &ctxt.getContext(),
           "Cannot decrypt with non-matching context");
  NTL::ZZX pp;
  sKey.Decrypt(pp, ctxt);

#if 0

  // convert to zzX, if the pp is too big, scale it down
  long nBits = NTL::MaxBits(pp) - NTL_SP_NBITS;
  zzX zpp(NTL::INIT_SIZE, deg(pp)+1);
  double factor;
  if (nBits<=0) { // convert to zzX, double
    for (long i=0; i<lsize(zpp); i++)
      conv(zpp[i], pp[i]);
    factor = NTL::to_double(ctxt.getRatFactor());
  } else { // scale and then convert to zzX, double
    for (long i=0; i<lsize(zpp); i++)
      conv(zpp[i], pp[i]>>nBits);
    factor = NTL::to_double(ctxt.getRatFactor()/NTL::power2_xdouble(nBits));
  }
  CKKS_canonicalEmbedding(ptxt, zpp, getPAlgebra()); // decode without scaling
  for (cx_double& cx : ptxt)  // divide by the factor
    cx /= factor;

#else

  // NOTE: I changed the code so that we convert to a
  // vector<double> instead of a zzX. This is more
  // efficient and more precise.  It should not affect overflow,
  // as far as I can tell. The old code is above.
  //   --Victor

  const long MAX_BITS = 400;
  long nBits = NTL::MaxBits(pp) - MAX_BITS;
  double factor;
  if (nBits <= 0) { // convert to zzX, double
    CKKS_canonicalEmbedding(ptxt, pp, getPAlgebra());
    factor = NTL::to_double(ctxt.getRatFactor());
  } else {
    long dpp = deg(pp);
    std::vector<double> pp_scaled(dpp + 1);
    NTL::ZZ tmp;
    for (long i : range(dpp + 1)) {
      RightShift(tmp, pp.rep[i], nBits);
      pp_scaled[i] = NTL::to_double(tmp);
    }
    CKKS_canonicalEmbedding(ptxt, pp_scaled, getPAlgebra());
    factor = NTL::to_double(ctxt.getRatFactor() / NTL::power2_xdouble(nBits));
  }
  for (cx_double& cx : ptxt) // divide by the factor
    cx /= factor;

#endif
}

// rotate ciphertext in dimension 0 by amt
void EncryptedArrayCx::rotate1D(Ctxt& ctxt,
                                long i,
                                long amt,
                                UNUSED bool dc) const
{
  assertEq(&getContext(),
           &ctxt.getContext(),
           "Cannot decrypt with non-matching context");
  assertTrue(nativeDimension(i),
             "Rotation in " + std::to_string(i) + " is not a native operation");

  const PAlgebra& palg = getPAlgebra();
  long ord = sizeOfDimension(i);
  amt %= ord; // DIRT: assumes division w/ remainder follows C++11 and C99 rules
  if (amt == 0)
    return;

  ctxt.smartAutomorph(palg.genToPow(i, amt));
}

// TODO: Shift k positions along the i'th dimension with zero fill.
// Negative shift amount denotes shift in the opposite direction.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void EncryptedArrayCx::shift1D(Ctxt& ctxt, long i, long k) const
{
  throw LogicError("EncryptedArrayCx::shift1D not implemented");
}
#pragma GCC diagnostic pop

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

double EncryptedArrayCx::encode(zzX& ptxt,
                                const std::vector<cx_double>& array,
                                double useThisSize,
                                long precision) const
{
  // VJS-FIXME: does it really make sense to use the *largest*
  // size in determining the factor?
  // It might make sense to use the *smallest* size.

  if (useThisSize < 0)
    for (auto& x : array) {
      if (useThisSize < std::abs(x))
        useThisSize = std::abs(x);
    }
  if (useThisSize <= 0)
    useThisSize = 1.0;

  // This factor ensures that encode/decode introduce less than 1/precision
  // error. If precision=0 then the error bound defaults to 2^{-almod.getR()}
  double factor = encodeScalingFactor(precision) / useThisSize;
  CKKS_embedInSlots(ptxt, array, getPAlgebra(), factor);
  return factor;
}

double EncryptedArrayCx::encode(zzX& ptxt,
                                double num,
                                double useThisSize,
                                long precision) const
{
  // This factor ensures that encode/decode introduce less than
  // 1/precision error. If precision=0 then the scaling factor
  // defaults to encodeScalingFactor(), corresponding to precision
  // error bound of 2^{-almod.getR()}
  if (useThisSize <= 0)
    useThisSize = roundedSize(num);
  double factor = encodeScalingFactor(precision) / useThisSize;

  resize(ptxt, 1, long(round(num * factor))); // Constant polynomial
  return factor;
}

double EncryptedArrayCx::encodei(zzX& ptxt, long precision) const
{
  std::vector<cx_double> v(size(), the_imaginary_i); // i in all the slots
  return this->encode(ptxt, v, /*size=*/1.0, precision);
}

const zzX& EncryptedArrayCx::getiEncoded() const
{
  if (lsize(iEncoded) <= 0)              // encoded-i not yet initialized
    encodei(const_cast<zzX&>(iEncoded)); // temporarily suspend cont-ness
  return iEncoded;
}

void EncryptedArrayCx::decode(std::vector<cx_double>& array,
                              const zzX& ptxt,
                              double scaling) const
{
  assertTrue<InvalidArgument>(scaling > 0,
                              "Scaling must be positive to decode");
  CKKS_canonicalEmbedding(array, ptxt, getPAlgebra());
  for (auto& x : array)
    x /= scaling;
}

// return an array of random complex numbers in a circle of radius rad
void EncryptedArrayCx::random(std::vector<cx_double>& array, double rad) const
{
  if (rad == 0)
    rad = 1.0; // radius

  resize(array, size()); // allocate space
  for (auto& x : array) {
    // VJS-FIXME: on 32-bit machines, this probably does not work correctly
    long bits = NTL::RandomLen_long(32);         // 32 random bits
    double r = std::sqrt(bits & 0xffff) / 256.0; // sqrt(uniform[0,1])
    // VJS-FIXME: Should use RandomBits_long?? should be unsigned long?
    // VJS-FIXME: why use 16 bits and take a square root?
    double theta =
        2.0L * PI * ((bits >> 16) & 0xffff) / 65536.0; // uniform(0,2pi)
    x = std::polar(rad * r, theta);
  }
}

void EncryptedArrayCx::extractRealPart(Ctxt& c) const
{
  Ctxt tmp = c;
  tmp.complexConj();         // the complex conjugate of c
  c += tmp;                  // c + conj(c) = 2*real(c)
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
  {
    Ctxt tmp = c;
    c.complexConj(); // the complex conjugate of c
    c -= tmp;        // conj(c) - c = -2*i*imaginary(c)
  }
  if (iDcrtPtr == nullptr) { // Need to encode i in a DoubleCRt object
    tmpDcrt.addPrimes(c.getPrimeSet());
    const zzX& iEncoded = getiEncoded();
    tmpDcrt.FFT(iEncoded, c.getPrimeSet());
    // FFT is a low-level DoubleCRT procedure to initialize an
    // existing object with a given PrimeSet and a given polynomial
    iDcrtPtr = &tmpDcrt;
  }
  c.multByConstantCKKS(*iDcrtPtr); // multiply by i
  c.multByConstantCKKS(0.5);       // divide by two
}

void EncryptedArrayCx::buildLinPolyCoeffs(std::vector<zzX>& C,
                                          const cx_double& oneImage,
                                          const cx_double& iImage,
                                          long precision) const
{
  resize(C, 2); // allocate space

  // Compute the constants x,y such that L(z) = x*z + y*conjugate(z)
  cx_double x = (oneImage - the_imaginary_i * iImage) * 0.5;
  cx_double y = (oneImage + the_imaginary_i * iImage) * 0.5;
  double sizex = std::abs(x);
  double sizey = std::abs(y);
  double msize = roundedSize(std::max(sizex, sizey));

  // Encode x,y in zzX objects
  long n = size();
  std::vector<cx_double> v(n, x); // x in all the slots
  encode(C[0], v, msize, precision);
  v.assign(n, y); // y in all the slots
  encode(C[1], v, msize, precision);
}

void EncryptedArrayCx::buildLinPolyCoeffs(
    std::vector<zzX>& C,
    const std::vector<cx_double>& oneImages,
    const std::vector<cx_double>& iImages,
    long precision) const
{
  resize(C, 2); // allocate space

  // Compute the constants x,y such that L(z) = x*z + y*conjugate(z)
  std::vector<cx_double> x(size());
  std::vector<cx_double> y(size());
  double msize = 0.0;
  for (long j = 0; j < size(); j++) {
    x[j] = (oneImages[j] - the_imaginary_i * iImages[j]) * 0.5;
    y[j] = (oneImages[j] + the_imaginary_i * iImages[j]) * 0.5;
    if (msize < std::abs(x[j]))
      msize = std::abs(x[j]);
    if (msize < std::abs(y[j]))
      msize = std::abs(y[j]);
  }
  // Encode x,y in zzX objects
  msize = roundedSize(msize);
  encode(C[0], x, msize, precision);
  encode(C[1], y, msize, precision);
}

} // namespace helib
