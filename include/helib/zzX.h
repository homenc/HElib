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
#ifndef HELIB_ZZX_H
#define HELIB_ZZX_H
/**
 * @file zzX.h - manipulating polynomials with single-precision coefficient
 *               It is assumed that the result is also single-precision
 **/
#include <NTL/vector.h>
#include <NTL/lzz_pX.h>
#include <NTL/GF2X.h>

namespace helib {

class PAlgebra;
typedef NTL::Vec<long> zzX;

inline bool IsZero(const zzX& a) { return a.length() == 0; }
inline void clear(zzX& a) { a.SetLength(0); }

inline void convert(NTL::zz_pX& x, const zzX& a)
{
  NTL::conv(x.rep, a);
  x.normalize();
}

void add(zzX& res, const zzX& a, const zzX& b);
inline zzX operator+(const zzX& a, const zzX& b)
{
  zzX tmp;
  add(tmp, a, b);
  return tmp;
}
inline zzX& operator+=(zzX& a, const zzX& b)
{
  add(a, a, b);
  return a;
}

void div(zzX& res, const zzX& a, long b);
inline zzX operator/(const zzX& a, long b)
{
  zzX tmp;
  div(tmp, a, b);
  return tmp;
}
inline zzX& operator/=(zzX& a, long b)
{
  div(a, a, b);
  return a;
}

void mul(zzX& res, const zzX& a, long b);
inline zzX operator*(const zzX& a, long b)
{
  zzX tmp;
  mul(tmp, a, b);
  return tmp;
}
inline zzX& operator*=(zzX& a, long b)
{
  mul(a, a, b);
  return a;
}

void normalize(zzX& f);

const NTL::zz_pXModulus& getPhimXMod(const PAlgebra& palg);
void reduceModPhimX(zzX& poly, const PAlgebra& palg);

void MulMod(zzX& res, const zzX& a, const zzX& b, const PAlgebra& palg);
inline zzX MulMod(const zzX& a, const zzX& b, const PAlgebra& palg)
{
  zzX tmp;
  MulMod(tmp, a, b, palg);
  return tmp;
}

// these produce properly balanced residues, with randomization
// if necessary
zzX balanced_zzX(const NTL::zz_pX& f);
zzX balanced_zzX(const NTL::GF2X& f);

} // namespace helib

#endif // ifndef HELIB_ZZX_H
