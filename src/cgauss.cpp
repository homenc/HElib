/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
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
