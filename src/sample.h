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
typedef NTL::Vec<long> zzX; // same as in NumbTh.h

//! Sample a degree-(n-1) poly, with -1/0/+1 coefficients. Each
//! coefficient is 0 with probability 1/2 and +-1 with probability 1/4.
void sampleSmall(zzX &poly, long n);
void sampleSmall(NTL::ZZX &poly, long n=0);

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

//! Implementing the Ducas-Durmus error procedure
void sampleErrorDD(zzX& err, const PAlgebra& palg, double stdev);
void sampleErrorDD(NTL::ZZX& err, const PAlgebra& palg, double stdev);

//! Helper function, returns a bound B such that for terms
//! of the form f = SampleSmall*SampleUniform(p), we have
//! Pr[|canonicalEmbed(f)|_{\infty} > B/3] < epsilon.
//! (The default is epsilon = 2^{-40}.)
double boundCanonEmb(long m, long phim, long p, double epsilon=9e-13);

#endif // _SAMPLE_H_
