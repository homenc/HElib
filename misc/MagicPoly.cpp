
#include <NTL/lzz_pX.h>
#include <helib/NumbTh.h>
#include <helib/ArgMap.h>

NTL_CLIENT

void compute_a_vals(Vec<long>& a, long p, long e)
// computes a[m] = a(m)/m! for m = p..(e-1)(p-1)+1,
// as defined by Chen and Han.
// a.length() is set to (e-1)(p-1)+2

{
   long p_to_e = power_long(p, e);
   long p_to_2e = power_long(p, 2*e);

   long len = (e-1)*(p-1)+2;

   zz_pPush push(p_to_2e);

   zz_pX x_plus_1_to_p = power(zz_pX(INIT_MONO, 1) + 1, p);
   zz_pX denom = InvTrunc(x_plus_1_to_p - zz_pX(INIT_MONO, p), len);
   zz_pX poly = MulTrunc(x_plus_1_to_p, denom, len);
   poly *= p;

   a.SetLength(len);

   long m_fac = 1;
   for (long m = 2; m < p; m++) {
      m_fac = MulMod(m_fac, m, p_to_2e);
   }

   for (long m = p; m < len; m++) {
      m_fac = MulMod(m_fac, m, p_to_2e);
      long c = rep(coeff(poly, m));
      long d = GCD(m_fac, p_to_2e);
      if (d == 0 || d > p_to_e || c % d != 0) Error("cannot divide");
      long m_fac_deflated = (m_fac / d) % p_to_e;
      long c_deflated = (c / d) % p_to_e;
      a[m] = MulMod(c_deflated, InvMod(m_fac_deflated, p_to_e), p_to_e);
   }

}

int main(int argc, char **argv)
{
   ArgMap amap;

   long p = 2;
   amap.arg("p", p);

   long e = 2;
   amap.arg("e", e);

   amap.parse(argc, argv);

   Vec<long> a;

   compute_a_vals(a, p, e);

   long p_to_e = power_long(p, e);
   long len = (e-1)*(p-1)+2;

   zz_pPush push(p_to_e);

   zz_pX poly(0);
   zz_pX term(1);
   zz_pX X(INIT_MONO, 1);

   poly = 0;
   term = 1;

   for (long m = 0; m < p; m++) {
      term *= (X-m);
   }

   for (long m = p; m < len; m++) {
      poly += term * a[m];
      term *= (X-m);
   }

   cout << poly << "\n";

   for (long i = 0; i < 1000; i++) {
      zz_p x;
      random(x);
      zz_p fx = eval(poly, x);
      if (fx + (rep(x) % p) != x) Error("bad eval");
   }
}
