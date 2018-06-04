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

// Sample polynomials with entries {-1,0,1}. These functions are similar to
// the SampleSmall class from v1, but without a class around it.

// In sampleSmall, 
// sampleHWt, min(Hwt,n) random coefficients are chosen at random in {-1,+1}
// and the others are set to zero. If n=0 then n=poly.deg()+1 is used. 

//! @brief Sample polynomials with entries {-1,0,1}. Each coefficient is 0 with probability 1/2 and +-1 with probability 1/4.
void sampleSmall(NTL::ZZX &poly, long n=0);

//! @brief Sample polynomials with entries {-1,0,1} with a given HAming weight.
//!
//! Choose min(Hwt,n) coefficients at random in {-1,+1} and the others are set
//! to zero. If n=0 then n=poly.deg()+1 is used. 
void sampleHWt(NTL::ZZX &poly, long Hwt, long n=0);

//! Sample polynomials with Gaussian coefficients.
void sampleGaussian(NTL::ZZX &poly, long n=0, double stdev=1.0);

//! Sample polynomials with coefficients sampled uniformy
//! over [-B..B]
void sampleUniform(NTL::ZZX& poly, const NTL::ZZ& B, long n=0);

#endif // _SAMPLE_H_
