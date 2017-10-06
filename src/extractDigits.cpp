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
/* EncryptedArray.cpp - Data-movement operations on arrays of slots
 */
#include <NTL/ZZ.h>
NTL_CLIENT
#include "EncryptedArray.h"
#include "polyEval.h"

static void buildDigitPolynomial(ZZX& result, long p, long e);

// extractDigits assumes that the slots of *this contains integers mod p^r
// i.e., that only the free terms are nonzero. (If that assumptions does
// not hold then the result will not be a valid ciphertext anymore.)
// 
// It returns in the slots of digits[j] the j'th-lowest gigits from the
// integers in the slots of the input. Namely, the i'th slot of digits[j]
// contains the j'th digit in the p-base expansion of the integer in the
// i'th slot of the *this. The plaintext space of digits[j] is mod p^{r-j},
// and all the digits are at the same level.
void extractDigits(vector<Ctxt>& digits, const Ctxt& c, long r)
{
  FHEcontext& context = (FHEcontext&) c.getContext();
  long rr = c.effectiveR();
  if (r<=0 || r>rr) r = rr; // how many digits to extract

  long p = context.zMStar.getP();
  ZZX x2p;
  if (p>3) { 
    buildDigitPolynomial(x2p, p, r);
  }

  Ctxt tmp(c.getPubKey(), c.getPtxtSpace());
  digits.resize(r, tmp);      // allocate space
  for (long i=0; i<r; i++) {
    tmp = c;
    for (long j=0; j<i; j++) {
      if (p==2) digits[j].square();
      else if (p==3) digits[j].cube();
      else polyEval(digits[j], x2p, digits[j]); // "in spirit" digits[j] = digits[j]^p
      tmp -= digits[j];
      tmp.divideByP();
    }
    digits[i] = tmp; // needed in the next round
  }
}



// Compute a degree-p polynomial poly(x) s.t. for any t<e and integr z of the
// form z = z0 + p^t*z1 (with 0<=z0<p), we have poly(z) = z0 (mod p^{t+1}).
//
// We get poly(x) by interpolating a degree-(p-1) polynomial poly'(x)
// s.t. poly'(z0)=z0 - z0^p (mod p^e) for all 0<=z0<p, and then setting
// poly(x) = x^p + poly'(x).
static void buildDigitPolynomial(ZZX& result, long p, long e)
{
  if (p<2 || e<=1) return; // nothing to do
  FHE_TIMER_START;
  long p2e = power_long(p,e); // the integer p^e

  // Compute x - x^p (mod p^e), for x=0,1,...,p-1
  vec_long x(INIT_SIZE, p);
  vec_long y(INIT_SIZE, p);
  long bottom = -(p/2);
  for (long j=0; j<p; j++) {
    long z = bottom+j;
    x[j] = z;
    y[j] = z-PowerMod((z < 0 ? z + p2e : z), p, p2e);  // x - x^p (mod p^e)

    while (y[j] > p2e/2)         y[j] -= p2e;
    while (y[j] < -(p2e/2))      y[j] += p2e;
  }
  interpolateMod(result, x, y, p, e);
  assert(deg(result)<p); // interpolating p points, should get deg<=p-1
  SetCoeff(result, p);   // return result = x^p + poly'(x)
  //  cerr << "# digitExt mod "<<p<<"^"<<e<<"="<<result<<endl;
  FHE_TIMER_STOP;
}
