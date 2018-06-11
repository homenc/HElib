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
#include <complex>
#include <NTL/ZZX.h>
#include <NTL/xdouble.h>

long sumOfCoeffs(const zzX& f);         // = f(1)
NTL::ZZ sumOfCoeffs(const NTL::ZZX& f); // = f(1)
NTL::ZZ sumOfCoeffs(const DoubleCRT& f); // somewhat lame implementation

//! The L-infinity norm of an element (in coefficient representation)
long largestCoeff(const zzX& f);         // l_infty norm
NTL::ZZ largestCoeff(const NTL::ZZX& f); // l_infty norm
NTL::ZZ largestCoeff(const DoubleCRT& f);

//! The L2-norm of an element (in coefficient representation)
double coeffsL2NormSquared(const zzX& f); // l2 norm^2
NTL::xdouble coeffsL2NormSquared(const NTL::ZZX& f); // l2 norm^2
NTL::xdouble coeffsL2NormSquared(const DoubleCRT& f); // l2 norm^2

inline double coeffsL2Norm(const zzX& f) // l2 norm
{ return sqrt(coeffsL2NormSquared(f)); }
inline NTL::xdouble coeffsL2Norm(const NTL::ZZX& f) // l2 norm
{ return sqrt(coeffsL2NormSquared(f)); }
inline NTL::xdouble coeffsL2Norm(const DoubleCRT& f) // l2 norm
{ return sqrt(coeffsL2NormSquared(f)); }


//! The L2-norm of the canonical embedding of an element
double embeddingL2NormSquared(const zzX& f, long m);
NTL::xdouble embeddingL2NormSquared(const ZZX& f, long m);
NTL::xdouble embeddingL2NormSquared(const DoubleCRT& f);

typedef complex<double> dcomplex;
typedef complex<NTL::xdouble> xcomplex;

inline double embeddingL2Norm(const zzX& f, long m)
{ return sqrt(embeddingL2NormSquared(f,m)); }
inline NTL::xdouble embeddingL2Norm(const NTL::ZZX& f, long m)
{ return sqrt(embeddingL2NormSquared(f,m)); }
inline NTL::xdouble embeddingL2Norm(const DoubleCRT& f)
{ return sqrt(embeddingL2NormSquared(f)); }

//! Computing the canonical embedding, returning only the first half
//! of the entries, the others are v[phi(m)-i] = conj(v[i])
void canonicalEmbedding(std::vector<dcomplex>& v, const zzX& f, long m);
void canonicalEmbedding(std::vector<xcomplex>& v, const NTL::ZZX& f, long m);
void canonicalEmbedding(std::vector<xcomplex>& v, const DoubleCRT& f);
#endif // _NORMS_H_
