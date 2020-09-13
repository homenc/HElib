/* Copyright (C) 2012-2018 IBM Corp.
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


#include <NTL/BasicThreadPool.h>
#include "NumbTh.h"
#include "powerful.h"
#include "fhe_stats.h"
#include "ArgMap.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <cstdio>

NTL_CLIENT

#ifdef NTL_THREADS
typedef atomic<long> ATOMIC_LONG;
#else
typedef long ATOMIC_LONG;
#endif

static string v_values_name = "";



static void
dump_v_values()
{
  const vector<double> *v_values = fetch_saved_values("v_values");
  if (v_values && v_values_name != "") {
    // write v_values to a file

    cerr << "writing v_values to " << v_values_name << "\n";

    ofstream F;
    F.open(v_values_name.c_str());
    for (long i: range(v_values->size()))
      F << (*v_values)[i] << "\n";
  }
}

static void
anderson_darling(const vector<double>& X, double& AD, double& p_val)
{
  long N = X.size();

  if (N < 2) {
    AD = 0;
    p_val = 1;
    return;
  }

  vector<double> Y(X);
  sort(Y.begin(), Y.end());

  // compute the sample mean
  double SM = 0;
  for (long i: range(N)) SM += Y[i];
  SM /= N;


  // compute the sample variance
  double SV = 0;
  for (long i: range(N)) SV += (Y[i]-SM)*(Y[i]-SM);
  SV /= (N-1);

  // replace Y[i] by CDF of Y[i]
  for (long i: range(N))
    Y[i] = 0.5*(1 + erf((Y[i]-SM)/sqrt(2*SV)));

  double S = 0;
  for (long i: range(N)) {
    S += (2*i+1)*(log(Y[i]) + log1p(-Y[N-1-i]));
  }
  AD = -N - S/N;

  AD *= (1 + 0.75/N + 2.25/N/N);
  // This adjustment and the p-values below come from:
  // R.B. D'Augostino and M.A. Stephens, Eds., 1986,
  // Goodness-of-Fit Techniques, Marcel Dekker.

  p_val;
  if (AD >= 0.6)
    p_val = exp(1.2937 - 5.709*(AD)+ 0.0186*fsquare(AD));
  else if (AD > 0.34)
    p_val = exp(0.9177 - 4.279*(AD) - 1.38*fsquare(AD));
  else if (AD > 0.2)
    p_val = 1 - exp(-8.318 + 42.796*(AD)- 59.938*fsquare(AD));
  else
    p_val = 1 - exp(-13.436 + 101.14*(AD)- 223.73*fsquare(AD));
}

static void
print_anderson_darling()
{
  const vector<double> *v_values = fetch_saved_values("v_values");
  if (v_values) {
    double AD, p_val;
    anderson_darling(*v_values, AD, p_val);
    printf("AD=%6.4f, p_val=%8.6f", AD, p_val);
  }
}

static void
print_sigma_info()
{
  const vector<double> *v_values = fetch_saved_values("v_values");
  if (v_values) {
    double max_sigma = 0;
    for (double val: *v_values) {
      if (fabs(val) > max_sigma) max_sigma = fabs(val);
    }
    double max_sigma_prob = 1 - exp(double(v_values->size())*log(erf(max_sigma/sqrt(2))));
    // max_sigma_prob is the probability of seeing max_sigma this big,
    // assuming independent gaussians
    printf("max_sigma=%6.4f, max_sigma_prob=%8.6f\n", max_sigma, max_sigma_prob);
  }
}




int main(int argc, char *argv[])
{
  ArgMap amap;

  long t=120;
  long nthreads=1;
  long seed=0;
  Vec<long> mvec;
  long m = 0;
  long iter;
  long nfacs = 0;

  amap.arg("t", t, "Hamming weight of recryption secret key");
  amap.arg("nthreads", nthreads, "number of threads");
  amap.arg("seed", seed, "random number seed");
  amap.arg("mvec", mvec);
  amap.arg("m", m);
  amap.arg("nfacs", nfacs);
  amap.arg("iter", iter);
  amap.arg("v_values", v_values_name);

  amap.parse(argc, argv);

  if (seed)
    SetSeed(ZZ(seed));

  SetNumThreads(nthreads);

  fhe_stats = true;

  const long N = (1L << 20) + 1;

  if (mvec.length() == 0) {
    if (m == 0) {
      if (nfacs == 0) nfacs = RandomBnd(5)+1;
      vector<long> mmvec;
      do {
	m = RandomBnd(15000) + 25000;
	if (m % 2 == 0) m++;
        pp_factorize(mmvec, m);
      } while (mmvec.size() != nfacs);
      convert(mvec, mmvec);
    }
    else {
      vector<long> mmvec;
      pp_factorize(mmvec, m);
      convert(mvec, mmvec);
    }
  }
  else {
    m = 1;
    for (long fac: mvec) m *= fac;
  }

  long phim = phi_N(m);

  ZZX phimX_ZZX = Cyclotomic(m);

  long k = mvec.length();
  double mrat = double(phim)/double(m);
  double sigma = double(N) * sqrt( mrat * double(t) * double(1L << k)  / 3.0 ) * 0.5;

  ATOMIC_LONG counter(iter);

  zz_p::FFTInit(0);
  long p = zz_p::modulus();

  PowerfulTranslationIndexes ind(mvec);
  PowerfulConversion pConv(ind);

  zz_pX phimX = conv<zz_pX>(phimX_ZZX);
  zz_pXModulus PhimX(phimX);

  zz_pContext context;
  context.save();


  NTL_EXEC_INDEX(nthreads, index)

  context.restore();

  HyperCube<zz_p> pwrfl(pConv.getShortSig());
  Vec<zz_p>& pwrfl_data = pwrfl.getData();

  zz_pX poly1, poly2;

  while (--counter >= 0) {
    // initialize poly1 to a polynomial mod X^m-1 with t non-zero
    // +/- 1 coeffs
    poly1.rep.SetLength(m);
    for (long i: range(m)) poly1.rep[i] = 0;
    for (long i: range(t)) {
      long idx;
      do {
	idx = RandomBnd(m);
      } while (poly1.rep[idx] != 0);
      if (RandomBnd(2))
        poly1.rep[idx] = 1;
      else
        poly1.rep[idx] = -1;
    }
    poly1.normalize();
    rem(poly1, poly1, PhimX);

    // initialize pwrfl1 to random mod N (balanced)
    for (zz_p& data: pwrfl_data) {
      long val = RandomBnd(N);
      data = balRem(val, N);
    }

    pConv.powerfulToPoly(poly2, pwrfl);

    MulMod(poly2, poly1, poly2, PhimX);

    pConv.polyToPowerful(pwrfl, poly2);

    double max_pwrfl = 0;
    for (long i: range(phim)) {
       long max_pwrfl0 = rep(pwrfl_data[i]);
       max_pwrfl0 = balRem(max_pwrfl0, p);
       double max_pwrfl1 = max_pwrfl0;
       double std_devs = fabs(max_pwrfl1)/sigma;
       if (std_devs > max_pwrfl) max_pwrfl = std_devs;
    }

    FHE_STATS_SAVE("max_pwrfl_saved", max_pwrfl);

    long ran_pwrfl0 = rep(pwrfl_data[RandomBnd(phim)]);
    ran_pwrfl0 = balRem(ran_pwrfl0, p);
    double ran_pwrfl = ran_pwrfl0;
    double std_devs = fabs(ran_pwrfl)/sigma;

    // update various indicator variables
    FHE_STATS_UPDATE("sigma_1_0", double(std_devs <= 1.0)); // 0.683
    FHE_STATS_UPDATE("sigma_2_0", double(std_devs <= 2.0)); // 0.954
    FHE_STATS_UPDATE("sigma_3_0", double(std_devs <= 3.0)); // 0.997, 1 in 370

    // compute sample variance, and scale by the variance we expect
    FHE_STATS_UPDATE("var_calc", fsquare(ran_pwrfl)/fsquare(sigma));

    // save the scaled value for application of other tests
    FHE_STATS_SAVE("v_values", ran_pwrfl/sigma);
  }

  NTL_EXEC_INDEX_END

  const vector<double>& max_pwrfl_saved = *fetch_saved_values("max_pwrfl_saved");
  double max_pwrfl = 0;
  for (double one_pwrfl: max_pwrfl_saved) {
    if (one_pwrfl > max_pwrfl) max_pwrfl = one_pwrfl;
  }

  double max_pwrfl_prob = 1 - exp(double(iter)*double(phim)*log(erf(max_pwrfl/sqrt(2))));
  // max_pwrfl_prob is the probability of seeing max_pwrfl this big,
  // assuming independent gaussians

  double max_pwrfl_prob0, max_pwrfl_prob1;
  // max_pwrfl_prob1 is the probability of seeing max_pwrfl this big,
  // using union bound

  max_pwrfl_prob0 = double(phim)*erfc(max_pwrfl/sqrt(2));
  if (max_pwrfl_prob0 > 1)
    max_pwrfl_prob1 = 1;
  else
    max_pwrfl_prob1 = 1 - exp(double(iter)*log(1-max_pwrfl_prob0));


  cerr << "\n";

  print_stats(cout);
  print_anderson_darling();
  cout << ", mvec=" << mvec << "\n";
  print_sigma_info();
  printf("max_pwrfl=%6.4f, max_pwrfl_prob=%8.6f, max_pwrfl_prob1=%8.6f\n", max_pwrfl, max_pwrfl_prob, max_pwrfl_prob1);
  dump_v_values();

  return 0;
}
