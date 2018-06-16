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
#ifndef _SAMPLE_H_
#define _SAMPLE_H_
/**
 * @file sample.h - implementing various sampling routines 
 **/
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include "zzX.h"

//! Sample a degree-(n-1) poly, with -1/0/+1 coefficients.
//! Each coefficients is +-1 with probability prob/2 each,
//! and 0 with probability 1-prob. By default, pr[nonzero]=1/2.
void sampleSmall(zzX &poly, long n, double prob=0.5);
void sampleSmall(NTL::ZZX &poly, long n, double prob=0.5);

//! Sample a degree-(n-1) poly as bove, with only Hwt nonzero coefficients
void sampleHWt(zzX &poly, long n, long Hwt=100);
void sampleHWt(NTL::ZZX &poly, long n, long Hwt=100);

//! Sample polynomials with Gaussian coefficients.
void sampleGaussian(zzX &poly, long n, double stdev);
void sampleGaussian(NTL::ZZX &poly, long n, double stdev);

//! Sample a degree-(n-1) ZZX, with coefficients uniform in [-B,B]
void sampleUniform(zzX& poly, long n, long B=100);
void sampleUniform(NTL::ZZX& poly, long n, const NTL::ZZ& B=NTL::ZZ(100L));

//! Choose a vector of continuous Gaussians
void sampleGaussian(std::vector<double> &dvec, long n, double stdev);

class PAlgebra;
// Same as above, but sample mod X^m-1 and then reduce mod Phi_m(X)
void sampleHWt(zzX &poly, const PAlgebra& palg, long Hwt=100);
void sampleSmall(zzX &poly, const PAlgebra& palg);
void sampleGaussian(zzX &poly, const PAlgebra& palg, double stdev);
void sampleUniform(zzX &poly, const PAlgebra& palg, long B=100);
void sampleUniform(NTL::ZZX&poly, const PAlgebra& palg,
                   const NTL::ZZ& B=NTL::ZZ(100L));

//! Helper functions, return a bound B such that for random noise
//! terms we have Pr[|canonicalEmbed(noise)|_{\infty} > B] < epsilon.
//! (The default is epsilon = 2^{-40}.)
///@{
double boundFreshNoise(long m, long phim, long p2r, double epsilon=9e-13);
double boundRoundingNoise(long m, long phim, long p2r, double epsilon=9e-13);
///@}

void reduceModPhimX(zzX& poly, const PAlgebra& palg);
const NTL::zz_pXModulus& getPhimXMod(const PAlgebra& palg);

#endif // _SAMPLE_H_
