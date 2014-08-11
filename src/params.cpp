
namespace std {}
namespace NTL {}

using namespace std;
using namespace NTL;

#include "NumbTh.h"
#include "PAlgebra.h"

const long MaxOrd = 100;
const long MaxGens = 2;

int main()
{
   long total = 0;
   long ord_skip = 0;
   long gen_skip = 0;
   long fac_skip = 0;
   long good_m = 0;

   
   for (long m = 1001; m <= 100000; m += 2) {

      if (m % 1000 == 1 && m > 1001) {
         cout << "*** " << (m-1) << "\n";
         cout.flush();
      }

      total++;

      if (multOrd(2, m) > MaxOrd) {
         ord_skip++;
         continue;
      }

      PAlgebra pal(m);

      if (long(pal.numOfGens()) > MaxGens) {
         gen_skip++;
         continue;
      }

      if (long(pal.numOfGens()) > 1 && !pal.SameOrd(0)) {
         gen_skip++;
         continue;
      }


      Vec< Pair<long, long> > fac;
      factorize(fac, m);

      if (fac[fac.length()-1].a > 5*sqrt(m)) {
        fac_skip++;
        continue;
      }

      good_m++;

      cout << "m=" << m << "=";
      for (long i = 0; i < fac.length(); i++) {
        if (i > 0) cout << "*";
        cout << fac[i].a;
        if (fac[i].b > 1) cout << "^" << fac[i].b;
      }
      cout << " ";

      cout << "phi(m)=" << pal.getPhiM() << " "
           << "d=" << pal.getOrdP() << " "
           << "exp=" << double(long(100*double(m)/double(pal.getPhiM())))/100.0 << " "
           << "ngens=" << pal.numOfGens() << " ";

      for (long i = 0; i < long(pal.numOfGens()); i++)
         cout << pal.OrderOf(i) << (pal.SameOrd(i) ? "" : "!") << " ";



      cout << "\n";
   }

   cout << "*************\n";
   cout << "total=" << total << " "
        << "ord_skip=" << ord_skip << " "
        << "gen_skip=" << gen_skip << " "
        << "fac_skip=" << fac_skip << " "
        << "good_m=" << good_m << " ";

   cout << "\n";

}
