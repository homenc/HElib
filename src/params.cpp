
namespace std {}
namespace NTL {}

using namespace std;
using namespace NTL;

#include "NumbTh.h"
#include "PAlgebra.h"
#include <iomanip>

const long MaxOrd = 100;
const long MaxGens = 2;

int main()
{
   long total = 0;
   long ord_skip = 0;
   long gen_skip = 0;
   long power_skip = 0;
   long good_m = 0;

   
   for (long m = 1001; m <= 100000; m += 2) {

#if 0
      if (m % 1000 == 1 && m > 1001) {
         cout << "*** " << (m-1) << "\n";
         cout.flush();
      }
#endif


      total++;

      long d = multOrd(2, m);

      if (d > MaxOrd) {
         ord_skip++;
         continue;
      }

#if 0
      PAlgebra pal(m);

      if (long(pal.numOfGens()) > MaxGens) {
         gen_skip++;
         continue;
      }

      if (long(pal.numOfGens()) > 1 && !pal.SameOrd(0)) {
         gen_skip++;
         continue;
      }
#endif

      long phim = phi_N(m);

      Vec< Pair<long, long> > fac;
      factorize(fac, m);

      long k = fac.length();

      Vec<long> fac1;
      fac1.SetLength(k);

      for (long i = 0; i < k; i++)
         fac1[i] = power_long(fac[i].a, fac[i].b);


      long good_guy = -1;
      for (long i = 0; i < k; i++) {
         long m1 = fac1[i];
         long phim1 = phi_N(m1);
         if (phim1 <= 5*sqrt(phim) && multOrd(2, m1) == d) {
            // test if Z_{m1}*/<2> has a single good generator
            if (phim1 == d) { good_guy = i; break; }
            PAlgebra pal1(m1);
            if (pal1.numOfGens() == 1 && pal1.SameOrd(0))
               { good_guy = i; break; }
         }
      }

      long alt_good_guy = -1;

      for (long i = 0; i < k; i++) {
         long m1 = m/fac1[i];
         if (m1 == 1) continue;
         long phim1 = phi_N(m1);
         if (phim1 <= 5*sqrt(phim) && multOrd(2, m1) == d) {
            // test if Z_{m1}*/<2> has a single good generator
            if (phim1 == d) { alt_good_guy = i; break; }
            PAlgebra pal1(m1);
            if (pal1.numOfGens() == 1 && pal1.SameOrd(0))
               { alt_good_guy = i; break; }
         }
      }


      if (good_guy == -1 && alt_good_guy == -1) {
         power_skip++;
         continue;
      }

      good_m++;

      cout << setw(6) << phim << "  ";
      cout << setw(4) << d << "  ";
      cout << setw(6) << m << "  ";

      cout << "m=";
      for (long i = 0; i < k; i++) {
        if (i > 0) cout << "*";
        if (i == good_guy) cout << "(";
        if (i == alt_good_guy) cout << "[";
        cout << fac[i].a;
        if (fac[i].b > 1) cout << "^" << fac[i].b;
        if (i == alt_good_guy) cout << "]";
        if (i == good_guy) cout << ")";
      }
      cout << " ";


      cout << "m/phim(m)=" << double(long(100*double(m)/double(phim)))/100.0 << " ";
#if 0
      cout << "ngens=" << pal.numOfGens() << " ";

      for (long i = 0; i < long(pal.numOfGens()); i++)
         cout << pal.OrderOf(i) << (pal.SameOrd(i) ? "" : "!") << " ";
#endif



      cout << "\n";
      cout.flush();
   }

#if 0

   cout << "*************\n";
   cout << "total=" << total << " "
        << "ord_skip=" << ord_skip << " "
        << "gen_skip=" << gen_skip << " "
        << "power_skip=" << power_skip << " "
        << "good_m=" << good_m << " ";

   cout << "\n";
#endif

}
