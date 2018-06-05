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
/**
 * @file norms.cpp - computing various norms of ring elements
 **/
#include "norms.h"
NTL_CLIENT

ZZ largestCoeff(const ZZX& f)
{
  ZZ mx = ZZ::zero();
  for (long i=0; i<=deg(f); i++) {
    if (mx < abs(coeff(f,i)))
      mx = abs(coeff(f,i));
  }
  return mx;
}

ZZ sumOfCoeffs(const ZZX& f) // = f(1)
{
  ZZ sum = ZZ::zero();
  for (long i=0; i<=deg(f); i++) sum += coeff(f,i);
  return sum;
}

xdouble coeffsL2Norm(const ZZX& f) // l_2 norm
{
  xdouble s = to_xdouble(0.0);
  for (long i=0; i<=deg(f); i++) {
    xdouble coef = to_xdouble(coeff(f,i));
    s += coef * coef;
  }
  return sqrt(s);
}
