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
// g++ -o cgauss -I/usr/local/include cgauss.cpp -L/usr/local/lib -lntl -lm

// So this confirms that E(s^k)/H^k = k!

#include <NTL/ZZ.h>

#include <complex>


NTL_CLIENT



//evaluate f at e^{2 pi i/m}, returning a complex number

complex<double> evalPoly(double *f, long i, long m)
{
  complex<double> t(0.0, 2*M_PI*i/((double) m));
  complex<double> x = exp(t);

  complex<double> res = 0.0;
  for (long j = m-1; j >= 0; j--)
    res = res*x + f[j];

  return res;
}



int main()
{
   long m = 1001;
   long h = 100;

   long maxpow = 7;
   long iter = 100000;

   double *s = new double[m];
   complex<double> *pow = new complex<double> [m];
   double *sum = new double[maxpow+1];

   for (long k = 1; k <= maxpow; k++) sum[k] = 0;

   complex<double> t(0.0, 2*M_PI*1/((double) m));
   complex<double> x = exp(t);
   complex<double> xi = 1.0;
   for (long i = 0; i < m; i++) {
      pow[i] = xi;
      xi *= x;
   }

   for (long u = 0; u < iter; u++) {
      complex<double> v = 0.0;
      for (long i = 0; i < m; i++) {
         if (RandomBnd(m) < h) { // true w/ prob. h/m
            if (RandomBnd(2) == 0)
               v += pow[i];
            else
               v -= pow[i];
         }
      }

      double nv = norm(v);
      double nvk = nv;

      for (long k = 1; k <= maxpow; k++) {
         sum[k] += nvk;
         nvk *= nv;
      }
   }

   for (long k = 1; k <= maxpow; k++) {
      double ave = sum[k]/iter;
      cout << ave/exp(k*log(h)) << "\n";
   }
}
