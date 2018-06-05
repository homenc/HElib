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
#ifndef _NORMS_H_
#define _NORMS_H_
/**
 * @file norms.h - computing various norms of ring elements
 **/
#include <NTL/ZZX.h>
#include <NTL/xdouble.h>

NTL::ZZ sumOfCoeffs(const NTL::ZZX& f);  // = f(1)
NTL::ZZ largestCoeff(const NTL::ZZX& f); // l_infty norm
NTL::xdouble coeffsL2Norm(const NTL::ZZX& f); // l_2 norm

#endif // _NORMS_H_
