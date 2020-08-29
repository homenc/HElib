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
 * @file eqtesting.cpp
 * @brief Useful functions for equality testing...
 */
#include <NTL/lzz_pXFactoring.h>
#include <helib/timing.h>
#include <helib/EncryptedArray.h>
#include <helib/Ptxt.h>

#include <cstdio>

namespace helib {

// Map all non-zero slots to 1, leaving zero slots as zero.
// Assumes that r=1, and that all the slot contain elements from GF(p^d).
//
// We compute x^{p^d-1} = x^{(1+p+...+p^{d-1})*(p-1)} by setting y=x^{p-1}
// and then outputting y * y^p * ... * y^{p^{d-1}}, with exponentiation to
// powers of p done via Frobenius.

// FIXME: the computation of the "norm" y * y^p * ... * y^{p^{d-1}}
// can be done using O(log d) automorphisms, rather than O(d).

void mapTo01(const EncryptedArray& ea, Ctxt& ctxt)
{
  long p = ctxt.getPtxtSpace();
  if (p != ea.getPAlgebra().getP()) // ptxt space is p^r for r>1
    throw LogicError("mapTo01 not implemented for r>1");

  if (p > 2)
    ctxt.power(p - 1); // set y = x^{p-1}

  long d = ea.getDegree();
  if (d > 1) { // compute the product of the d automorphisms
    std::vector<Ctxt> v(d, ctxt);
    for (long i = 1; i < d; i++)
      v[i].frobeniusAutomorph(i);
    totalProduct(ctxt, v);
  }
}

template <typename Scheme>
void mapTo01(const EncryptedArray&, Ptxt<Scheme>& ptxt)
{
  ptxt.mapTo01();
}

template void mapTo01(const EncryptedArray&, Ptxt<BGV>& ptxt);
template void mapTo01(const EncryptedArray&, Ptxt<CKKS>& ptxt);

// computes ctxt^{2^d-1} using a method that takes
// O(log d) automorphisms and multiplications
void fastPower(Ctxt& ctxt, long d)
{
  assertEq(ctxt.getPtxtSpace(), 2l, "ptxtSpace must be 2");
  if (d <= 1)
    return;

  Ctxt orig = ctxt;

  long k = NTL::NumBits(d);
  long e = 1;

  for (long i = k - 2; i >= 0; i--) {
    Ctxt tmp1 = ctxt;
    tmp1.smartAutomorph(1L << e);
    ctxt.multiplyBy(tmp1);
    e = 2 * e;

    if (NTL::bit(d, i)) {
      ctxt.smartAutomorph(2);
      ctxt.multiplyBy(orig);
      e += 1;
    }
  }
}

// ===> This function only works for p=2, r=1 <===
// Test if prefixes of bits in slots are all zero: Set slot j of res[i] to 0
// if bits 0..i of j'th slot in ctxt are all zero, else it is set to 1
// It is assumed that res and the res[i]'s are initialized by the caller.
// Complexity: O(d + n log d) smart automorphisms
//             O(n d)
void incrementalZeroTest(Ctxt* res[],
                         const EncryptedArray& ea,
                         const Ctxt& ctxt,
                         long n)
{
  HELIB_TIMER_START;
  long nslots = ea.size();
  long d = ea.getDegree();

  // compute linearized polynomial coefficients

  std::vector<std::vector<NTL::ZZX>> Coeff;
  Coeff.resize(n);

  for (long i = 0; i < n; i++) {
    // coefficients for mask on bits 0..i
    // L[j] = X^j for j = 0..i, L[j] = 0 for j = i+1..d-1

    std::vector<NTL::ZZX> L;
    L.resize(d);

    for (long j = 0; j <= i; j++)
      SetCoeff(L[j], j);

    std::vector<NTL::ZZX> C;

    ea.buildLinPolyCoeffs(C, L);

    Coeff[i].resize(d);
    for (long j = 0; j < d; j++) {
      // Coeff[i][j] = to the encoding that has C[j] in all slots
      // FIXME: maybe encrypted array should have this functionality
      //        built in
      std::vector<NTL::ZZX> T;
      T.resize(nslots);
      for (long s = 0; s < nslots; s++)
        T[s] = C[j];
      ea.encode(Coeff[i][j], T);
    }
  }

  std::vector<Ctxt> Conj(d, ctxt);
  // initialize Cong[j] to ctxt^{2^j}
  for (long j = 0; j < d; j++) {
    Conj[j].smartAutomorph(1L << j);
  }

  for (long i = 0; i < n; i++) {
    res[i]->clear();
    for (long j = 0; j < d; j++) {
      Ctxt tmp = Conj[j];
      tmp.multByConstant(Coeff[i][j]);
      *res[i] += tmp;
    }

    // *res[i] now has 0..i in each slot
    // next, we raise to the power 2^d-1

    fastPower(*res[i], d);
  }
  HELIB_TIMER_STOP;
}

} // namespace helib
