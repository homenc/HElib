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
/**
 * @file zzX.cpp - manipulating polynomials with single-precision coefficient
 *               It is assumed that the result is also single-precision
 **/
#include <mutex>
#include <map>

#include <helib/PAlgebra.h>
#include <helib/timing.h>
#include <helib/zzX.h>
#include <helib/range.h>

namespace helib {

void MulMod(zzX& res, const zzX& a, const zzX& b, const PAlgebra& palg)
{
  HELIB_TIMER_START;
  NTL::zz_pPush push; // backup the NTL current modulus
  const NTL::zz_pXModulus& phimX = getPhimXMod(palg);

  NTL::zz_pX aa, bb;
  convert(aa, a); // convert to zz_pX
  convert(bb, b); // convert to zz_pX
  NTL::MulMod(aa, aa, bb, phimX);
  convert(res, aa);
}

void add(zzX& res, const zzX& a, const zzX& b)
{
  HELIB_TIMER_START;
  // The shorter one is aa, the longer is bb
  const zzX& aa = (lsize(a) < lsize(b)) ? a : b;
  const zzX& bb = (lsize(a) < lsize(b)) ? b : a;
  res.SetLength(lsize(bb));
  for (long i = 0; i < lsize(aa); i++)
    res[i] = aa[i] + bb[i];
  for (long i = lsize(aa); i < lsize(bb); i++)
    res[i] = bb[i];
}
void mul(zzX& res, const zzX& a, long b)
{
  res.SetLength(lsize(a));
  for (long i = 0; i < lsize(a); i++)
    res[i] *= a[i] * b;
}

void div(zzX& res, const zzX& a, long b)
{
  res.SetLength(lsize(a));
  for (long i = 0; i < lsize(a); i++)
    res[i] = a[i] / b;
}

void normalize(zzX& f)
{
  long deg = f.length() - 1;
  if (deg < 0)
    return; // empty vector = zero polynomial
  while (deg >= 0 && f[deg] == 0)
    deg--;
  f.SetLength(deg + 1);
}

// Helper functions, returns a zz_pXModulus object, modulo Phi_m(X)
// and a single 60-bit prime. Can be used to get faster operation
// modulo Phi_m(X), where we know apriori that the numbers do not wrap.
// This function changes the NTL current zz_p modulus.
const NTL::zz_pXModulus& getPhimXMod(const PAlgebra& palg)
{
  static std::map<long, NTL::zz_pXModulus*> moduli; // pointer per value of m
  static std::mutex pt_mtx; // control access to modifying the map

  NTL::zz_p::FFTInit(0); // set "the best FFT prime" as NTL's current modulus

  long m = palg.getM();
  auto it = moduli.find(m); // check if we already have zz_pXModulus for m

  if (it == moduli.end()) { // init a new zz_pXModulus for this value of m
    std::unique_lock<std::mutex> lck(pt_mtx); // try to get a lock

    // Got the lock, insert a new entry for this value of m into the map
    NTL::zz_pX phimX = NTL::conv<NTL::zz_pX>(palg.getPhimX());
    NTL::zz_pXModulus* ptr =
        new NTL::zz_pXModulus(phimX); // will "never" be deleted

    // insert returns a pair (iterator, bool)
    auto ret = moduli.insert(std::pair<long, NTL::zz_pXModulus*>(m, ptr));
    if (ret.second == false) // Another thread inserted it, delete your copy
      delete ptr;
    // FIXME: Could leak memory if insert throws an exception
    //        without inserting the element (but who cares)

    it = ret.first; // point to the entry in the map
  }
  return *(it->second);
}

// DIRT: We use modular arithmetic mod p \approx 2^{60} as a
//       substitute for computing on rational numbers
void reduceModPhimX(zzX& poly, const PAlgebra& palg)
{
  NTL::zz_pPush push; // backup the NTL current modulus
  const NTL::zz_pXModulus& phimX = getPhimXMod(palg);

  NTL::zz_pX pp;
  convert(pp, poly); // convert to zz_pX
  rem(pp, pp, phimX);
  convert(poly, pp);
}

zzX balanced_zzX(const NTL::zz_pX& f)
{
  long p = NTL::zz_p::modulus();
  long len = deg(f) + 1;
  zzX out;
  out.SetLength(len);
  for (long i : range(len)) {
    long coef = NTL::conv<long>(f[i]);
    if (coef > p / 2 || (p % 2 == 0 && coef == p / 2 && NTL::RandomBnd(2)))
      coef -= p;

    out[i] = coef;
  }

  return out;
}

zzX balanced_zzX(const NTL::GF2X& f)
{
  long len = deg(f) + 1;
  zzX out;
  out.SetLength(len);
  for (long i : range(len)) {
    if (f[i] == 0)
      out[i] = 0;
    else if (NTL::RandomBnd(2))
      out[i] = -1;
    else
      out[i] = 1;
  }

  return out;
}

} // namespace helib
