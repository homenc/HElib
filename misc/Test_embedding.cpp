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
#include <cassert>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <cmath>
#include <atomic>
#include <NTL/BasicThreadPool.h>
#include <NTL/ZZ_pX.h>
#include <NTL/xdouble.h>
#include "NumbTh.h"
#include "PAlgebra.h"
#include "DoubleCRT.h"
#include "Context.h"
#include "sample.h"
#include "timing.h"
#include "norms.h"

NTL_CLIENT

void printHistogram(const vector<double>& data,
                    const cx_double& mean, double step)
{
  vector<long> hist(8,0); // histogram: hist[i]= # of x'es with |x|<i*stdev
  double max = 0.0;
  double sum=0.0, sumSqr=0.0;
  for (double x: data) {
    sum += x;
    sumSqr += x*x;
    long j = std::floor(x/step);
    if (j >= lsize(hist)) hist.resize(j+1, 0);
    hist[j]++;
    if (x > max) max = x;
  }
  sum /= data.size();    // E[x]
  sumSqr /= data.size(); // E[x^2]
  double stdev = sqrt(sumSqr - (sum*sum));

  for (long i=hist.size()-1; i>0; i--) // cumulative
    hist[i-1] += hist[i];

  vector<double> dhist(hist.size(), 1.0);
  vector<double> ratio(hist.size(), 0.0);
  for (long i=1; i<lsize(dhist); i++) {
    dhist[i] = hist[i]/double(hist[0]);
    if (dhist[i]>0.0) ratio[i] = dhist[i-1]/dhist[i];
  }
  cout << data.size() << " points, mean="<<mean<<endl;
  cout << "size mean="<<sum<<", stdev="<<stdev
       << ", max="<<max<<".  Step="<<step<<endl;
  cout << "  histogram ="<< hist<<endl
       << " probability="<<dhist<<endl
       << "      ratio ="<<ratio<<endl<<endl;
}

void freshCtxtNoise(zzX& f, const Context& context,
                    double sigma, bool modPhimX)
{
  zzX s, r, e1, e2, e3;
  const PAlgebra& palg = context.zMStar;
  sampleSmallBounded(s, context);
  if (modPhimX) {
    sampleSmall(r, palg.getPhiM()-1);
    sampleGaussian(e1, palg.getPhiM()-1, sigma);
    sampleGaussianBounded(e2, context, sigma);
    sampleGaussian(e3, palg.getPhiM()-1, sigma);
  } else {
    sampleSmall(r, context);
    sampleGaussian(e1, context, sigma);
    sampleGaussianBounded(e2, context, sigma);
    sampleGaussian(e3, context, sigma);
  }
  f = MulMod(s, e1, palg) + MulMod(r, e2, palg) + e3;
}

void roundingNoise(zzX& f, const Context& context,
                   long p2r, bool modPhimX)
{
  zzX s, e1, e2;
  const PAlgebra& palg = context.zMStar;
  sampleSmallBounded(s, context);
  if (modPhimX) {
    sampleUniform(e1, palg.getPhiM()-1, p2r);
    sampleUniform(e2, palg.getPhiM()-1, p2r);
  } else {
    sampleUniform(e1, context, p2r);
    sampleUniform(e2, context, p2r);
  }
  f = MulMod(s, e1, palg) + e2;
}

int main(int argc, char **argv)
{
  FHE_NTIMER_START(init);
  // get parameters from the command line
  ArgMapping amap;

  // long noPrint = 1;
  // amap.arg("noPrint", noPrint, "suppress printouts");

  long m = 15;
  amap.arg("m", m, "the cyclotomic index");
  double sigma = 0.0;
  amap.arg("sigma", sigma, "nomral standard deviation", "heristic");
  long p = 2;
  amap.arg("p", p, "plaintext base");
  long r = 1;
  amap.arg("r", r, "lifting");
  long N = 1000;
  amap.arg("N", N, "# of samples to use");
  long seed=0;
  amap.arg("seed", seed, "PRG seed");
  long nt=4;
  amap.arg("nt", nt, "number of threads");
  amap.parse(argc, argv);

  if (seed != 0) NTL::SetSeed(ZZ(seed));
  if (nt > 1)    NTL::SetNumThreads(nt);
  long p2r = power_long(p, r);
  if (sigma<=0.0) {
    if (m&1) // odd m
      sigma = 3.2*sqrt(m);
    else
      sigma=3.2;
  }

  Context context(m, p, r);
  const PAlgebra& palg = context.zMStar;
  buildModChain(context, /*L=*/5, /*c=*/3);
  long phim = palg.getPhiM();
  FHE_NTIMER_STOP(init);

  cout << "m="<<m<<", phi(m)="<<palg.getPhiM()
       << ", sigma="<<std::setprecision(3)<<sigma
       << ", p^r="<<p2r << endl;

  zzX f;
  vector<double> data;
  cx_double sum, mean;
  double step;

  step = (sigma+0.1)*0.54*(1+((m&1)? sqrt(phim*m): phim));
  // fresh ciphertext, sampling mod X^m-1
  data.resize(0);
  sum = 0;
  for (long i=0; i<N; i++) {
    std::vector<cx_double> cemb;
    freshCtxtNoise(f, context, sigma, false);
    canonicalEmbedding(cemb, f, palg); // Canonical embedding of f
    for (auto& entry : cemb) {
      double sz2 = conv<double>(std::norm<double>(entry));
      data.push_back(std::sqrt(sz2));
      sum += entry/double(lsize(cemb));
    }
  }
  mean = sum / double(N);
  cout << "fresh ctxt noise, sample mod X^m-1: ";
  printHistogram(data, mean, step);
  data.resize(0);
  sum = 0;
  // rounding error, mod Phi_m(X)

  step = (2*p2r+1)*(phim-2)/8.0;
  for (long i=0; i<N; i++) {
    std::vector<cx_double> cemb;
    roundingNoise(f, context, p2r, true);
    canonicalEmbedding(cemb, f, palg); // Canonical embedding of f
    for (auto& entry : cemb) {
      double sz2 = conv<double>(std::norm<double>(entry));
      data.push_back(std::sqrt(sz2));
      sum += entry/double(lsize(cemb));
    }
  }
  mean = sum / double(N);
  cout << "\nrounding noise, sample mod Phi_m(X): ";
  printHistogram(data, mean, step);
  data.resize(0);
  sum = 0;

  cout << endl;

  //  printAllTimers();
  return 0;
}

/********************************************************************/
#if 0 // OLD CODE
  vector<xdouble> l2(8, xdouble(0.0));
  vector<xdouble> ratio(8, xdouble(1.0));

  cout << "*** m="<<m<<", sampling "<<N<<" different f's ***\n";
  for (long i=0; i<N; i++) {
    ZZX f;
    sampleSmall(f, m);
    rem(f, f, phimX);
    if (IsZero(f)) continue;
    xdouble ll = embeddingL2NormSquared(f,m)/phim;
    l2[0] += ll;
    for (long j=1; j<lsize(l2); j++) {
      NTL::SqrMod(f, f, phimX); // fd -> fd^2
      ll *= ll;

      xdouble tt = embeddingL2NormSquared(f, m)/phim;
      l2[j] += tt;

      tt /= ll;
      if (tt>ratio[j]) ratio[j] = tt;
    }
  }
  xdouble base = l2[0] / N;
  long e = 1;
  xdouble factorial(1.0);
  for (long j=0; j<lsize(l2); j++) {
    l2[j] /= N;
    for (long i=e/2 +1; i<=e; i++) factorial *= i;
    cout << "E[|f^"<<e<<"|^2]="<<l2[j]
         << ",\t= "<<(l2[j]/base)<<"*E[|f^2|]^{"<<e<<"} (vs. "
         << e<<"! ="<<factorial<<")\n";
    cout << "\t\t\t max ratio = "<<ratio[j]<<endl;
    base *= base;
    e *= 2;
  }
#endif
#if 0
  step = (1+sqrt(phim*log(phim)))*0.85;
  for (long i=0; i<N; i++) {
    sampleSmall(f, palg);
    data.push_back(embeddingLargestCoeff(f, palg));
  }
  cout << "sampleSmall noise: ";
  printHistogram(data, mean, step);
  data.resize(0);
  step = sigma*(2+((m&1)? sqrt(m*log(phim)) : sqrt(phim*log(phim))))*1.2;
  for (long i=0; i<N; i++) {
    sampleGaussian(f, palg, sigma);
    data.push_back(embeddingLargestCoeff(f, palg));
  }
  cout << "sampleGaussian noise: ";
  printHistogram(data, mean, step);
  data.resize(0);
  exit(0);
#endif
#if 0
  // First test: fresh ciphertext, sampling mod Phi_m(X)
  for (long i=0; i<N; i++) {
    std::vector<cx_double> cemb;
    freshCtxtNoise(f, palg, sigma, true);
    canonicalEmbedding(cemb, f, palg); // Canonical embedding of f
    for (auto& entry : cemb) {
      double sz2 = conv<double>(std::norm<double>(entry));
      data.push_back(std::sqrt(sz2));
      sum += entry/double(lsize(cemb));
    }
  }
  mean = sum / double(N);
  cout << "fresh ctxt noise, sample mod Phi_m(X): ";
  printHistogram(data, mean, step);
#endif
#if 0
  for (long i=0; i<N; i++) {
    std::vector<cx_double> cemb;
    roundingNoise(f, palg, p2r, false);
    canonicalEmbedding(cemb, f, palg); // Canonical embedding of f
    for (auto& entry : cemb) {
      double sz2 = conv<double>(std::norm<double>(entry));
      data.push_back(std::sqrt(sz2));
      sum += entry/double(lsize(cemb));
    }
  }
  mean = sum / double(N);
  cout << "rounding noise, sample mod X^m-1: ";
  printHistogram(data, mean, step);
#endif
