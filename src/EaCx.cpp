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
#include <helib/log.h>
#include <helib/fhe_stats.h>

namespace helib {

static constexpr cx_double the_imaginary_i = cx_double(0.0, 1.0);

void EncryptedArrayCx::rawDecrypt(const Ctxt& ctxt,
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

void EncryptedArrayCx::rawDecrypt(const Ctxt& ctxt,
                                  const SecKey& sKey,
                                  std::vector<double>& ptxt) const
{
  std::vector<cx_double> v;
  rawDecrypt(ctxt, sKey, v);
  project(ptxt, v);
}

void EncryptedArrayCx::decrypt(const Ctxt& ctxt,
                               const SecKey& sKey,
                               std::vector<cx_double>& ptxt,
                               OptLong prec) const
{
  // This mitigates against the attack in
  // "On the Security of Homomorphic Encryption on Approximate Numbers",
  // by Baiyu Li and Daniele Micciancio.

  // As a preprocessing step, we add noise so that
  // the scaled error increases by at most eps (with some
  // futher adjustments made in addNoiseForCKKSDecryption
  // to maintain a certain level of security as the cost
  // of accuracy.

  double eps;

  if (prec.isDefined()) {
    // use eps = 2^{-prec}

    long p = prec;
    // some sanity checking/correcting
    const int p_bound = 100;
    if (p > p_bound) p = p_bound;
    if (p < -p_bound) p = -p_bound;
    eps = std::pow(2.0, -p);
  }
  else {
    // use eps = err/2^{adjust}
    const int adjust = 4;
    
    eps = ctxt.errorBound() * std::pow(2.0, -adjust);
  }

  // now add noise to a copy of ctxt
  Ctxt ctxt1 = ctxt;
  ctxt1.addNoiseForCKKSDecryption(sKey, eps);


  // finally, perform the decryption
  rawDecrypt(ctxt1, sKey, ptxt);
}

void EncryptedArrayCx::decrypt(const Ctxt& ctxt,
                               const SecKey& sKey,
                               std::vector<double>& ptxt,
                               OptLong prec) const
{
  std::vector<cx_double> v;
  decrypt(ctxt, sKey, v, prec);
  project(ptxt, v);
}

// rotate ciphertext in dimension 0 by amt
void EncryptedArrayCx::rotate1D(Ctxt& ctxt,
                                long i,
                                long amt,
                                UNUSED bool dc) const
{
  helib::assertEq(&context, &ctxt.getContext(), "Context mismatch");
  helib::assertInRange(i,
                       0l,
                       dimension(),
                       "i must be between 0 and dimension()");
  assertTrue(nativeDimension(i),
             "Rotation in " + std::to_string(i) + " is not a native operation");

  const PAlgebra& palg = getPAlgebra();
  long ord = sizeOfDimension(i);
  amt %= ord; // DIRT: assumes division w/ remainder follows C++11 and C99 rules
  if (amt == 0)
    return;
  if (amt < 0)
    amt += ord; // Make sure amt is in the range [1,ord-1]

  ctxt.smartAutomorph(palg.genToPow(i, amt));
}

// TODO: Shift k positions along the i'th dimension with zero fill.
// Negative shift amount denotes shift in the opposite direction.
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
void EncryptedArrayCx::shift1D(Ctxt& ctxt, long i, long k) const
{
  helib::assertEq(&context, &ctxt.getContext(), "Context mismatch");
  helib::assertInRange(i,
                       0l,
                       dimension(),
                       "i must be between 0 and dimension()");

  // NOTE: this works even in a non-native dimension, but
  // this is a bit academic

  const PAlgebra& al = getPAlgebra();

  long ord = al.OrderOf(i);

  if (k <= -ord || k >= ord) {
    ctxt.clear();
    return;
  }

  // Make sure amt is in the range [1,ord-1]
  long amt = k % ord;
  if (amt == 0)
    return;
  if (amt < 0)
    amt += ord;

  long val;
  if (k < 0)
    val = al.genToPow(i, amt - ord);
  else
    val = al.genToPow(i, amt);

  long n = size();
  std::vector<bool> maskArray(n);

  for (long j : range(n)) {
    long c = coordinate(i, j);
    if (c + k >= ord || c + k < 0)
      maskArray[j] = false;
    else
      maskArray[j] = true;
  }

  EncodedPtxt mask;
  encode(mask, maskArray);
  ctxt.multByConstant(mask);
  ctxt.smartAutomorph(val);
}
#pragma GCC diagnostic pop

// We only support linear arrays for approximate numbers,
// so rotate,shift are the same as rotate1D, shift1D
void EncryptedArrayCx::rotate(Ctxt& ctxt, long amt) const
{
  assertTrue(dimension() == 1,
             "CKKS rotation not supported in multi-dimensional hypercube");
  rotate1D(ctxt, 0, amt, true);
}
void EncryptedArrayCx::shift(Ctxt& ctxt, long amt) const
{
  assertTrue(dimension() == 1,
             "CKKS rotation not supported in multi-dimensional hypercube");
  shift1D(ctxt, 0, amt);
}

//====== New Encoding Functions ====

void EncryptedArrayCx::encode(EncodedPtxt& eptxt,
                              const std::vector<cx_double>& array,
                              double mag,
                              OptLong prec) const
{
  double actual_mag = Norm(array);
  if (mag < 0)
    mag = actual_mag;
  else {
    if (actual_mag > mag)
      Warning(
          "EncryptedArrayCx::encode: actual magnitude exceeds mag parameter");
  }

  double err = defaultErr(); 
  // For now, we use defaultErr().  We may want to eventually
  // allow APIs that use a different err value (such as the *actual*
  // err value).  However, if we encrypt this encoding, we
  // should not use a data-dependent err value. Moreover, I did
  // not want to have yet another esteric parameter for the user
  // to worry about.  We can revisit this later.

  double  scale = defaultScale(err, prec); // default scale

  zzX poly;
  CKKS_embedInSlots(poly, array, getPAlgebra(), scale);
  eptxt.resetCKKS(poly, mag, scale, err, getContext());

  // Check that error is actually bounded.
  // If this is too costly, we can consider only
  // running it in "debug mode".
  std::vector<cx_double> array1;
  decode(array1, poly, scale);
  double dist = Distance(array1, array);
  double scaled_err = err / scale;
  double ratio = dist / scaled_err;
  if (ratio > 1) {
    Warning("CKKS encode: error exceeds bound");
  }
  HELIB_STATS_UPDATE("CKKS_encode_ratio", ratio);
}

void EncryptedArrayCx::encode(EncodedPtxt& eptxt,
                              const PlaintextArray& array,
                              double mag,
                              OptLong prec) const
{
  encode(eptxt, array.getData<PA_cx>(), mag, prec);
}

void EncryptedArrayCx::decryptComplex(const Ctxt& ctxt,
                                      const SecKey& sKey,
                                      PlaintextArray& ptxt,
                                      OptLong prec) const
{
  decrypt(ctxt, sKey, ptxt.getData<PA_cx>(), prec);
}

void EncryptedArrayCx::rawDecryptComplex(const Ctxt& ctxt,
                                         const SecKey& sKey,
                                         PlaintextArray& ptxt) const
{
  rawDecrypt(ctxt, sKey, ptxt.getData<PA_cx>());
}

void EncryptedArrayCx::decryptReal(const Ctxt& ctxt,
                                   const SecKey& sKey,
                                   PlaintextArray& ptxt,
                                   OptLong prec) const
{
  std::vector<double> v;
  decrypt(ctxt, sKey, v, prec);
  convert(ptxt.getData<PA_cx>(), v);
}

void EncryptedArrayCx::rawDecryptReal(const Ctxt& ctxt,
                                      const SecKey& sKey,
                                      PlaintextArray& ptxt) const
{
  std::vector<double> v;
  rawDecrypt(ctxt, sKey, v);
  convert(ptxt.getData<PA_cx>(), v);
}

//======================

double EncryptedArrayCx::encode(zzX& ptxt,
                                const std::vector<cx_double>& array,
                                double useThisSize,
                                long precision) const
{
  // VJS-FIXME: this routine has a number of issues and should
  // be deprecated in favor of the new EncodedPtxt-based routines

  // VJS-NOTE: does it really make sense to use the *largest*
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
  // VJS-FIXME: this is NOT thread-safe
  // It also seems like it is not used anywhere, so I suggest we
  // get rid of it...

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
  // VJS-FIXME: this routine has a number of issues and should
  // be deprecated in favor of either the RandomComplex() routine
  // in NumbTh.h or PtxtArray::random().

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
