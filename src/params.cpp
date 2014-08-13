
namespace std {}
namespace NTL {}

using namespace std;
using namespace NTL;

#include "NumbTh.h"
#include "PAlgebra.h"
#include <iomanip>
#include <algorithm>

const long MaxOrd = 100;

long computePhi(const Pair<long,long>& x)
{
   return power_long(x.a, x.b-1) * (x.a-1);
}

bool comparePhi(const Pair<long,long>& x, const Pair<long,long>& y)
{
   return computePhi(x) < computePhi(y);
} 


int main()
{
   for (long m = 1001; m <= 100000; m += 2) {

      long d = multOrd(2, m);

      if (d > MaxOrd) continue;

      long phim = phi_N(m);

      Vec< Pair<long, long> > fac;
      factorize(fac, m);

      long k = fac.length();

      if (k == 1) continue;

      sort(fac.begin(), fac.end(), comparePhi);

      Vec<long> fac1;
      fac1.SetLength(k);

      for (long i = 0; i < k; i++)
         fac1[i] = power_long(fac[i].a, fac[i].b);

      Vec<long> phivec;
      phivec.SetLength(k);
      for (long i = 0; i < k; i++)
         phivec[i] = computePhi(fac[i]);

      if (phivec[k-1] >= 5*sqrt(phim)) continue;

      long gen_index = -1;
      bool good_gen = false;
      long first_gen = 0;


      for (long i = 0; i < k; i++) {
         long m1 = fac1[i];
         if (multOrd(2, m1) != d) continue;

         PAlgebra pal1(m1);
         bool good = (pal1.numOfGens() == 0 || 
                      (pal1.numOfGens() == 1 && pal1.SameOrd(0)));

         if (gen_index == -1 || good) {
            gen_index = i;
            good_gen = good;
            first_gen = pal1.ZmStarGen(0);
         }

         if (gen_index != -1 && good_gen) break;
      }

      bool cofactor = false;

      if (gen_index == -1 || !good_gen) {

         for (long i = 0; i < k; i++) {
            long m1 = m/fac1[i];
            long phim1 = phim/phivec[i];
            if (phim1 > 5*sqrt(phim) || multOrd(2, m1) != d) continue;

            PAlgebra pal1(m1);
            bool good = (pal1.numOfGens() == 0 || 
                         (pal1.numOfGens() == 1 && pal1.SameOrd(0)));

            if (gen_index == -1 || good) {
               gen_index = i;
               good_gen = good;
               first_gen = pal1.ZmStarGen(0);
               cofactor = true;
            }

            if (gen_index != -1 && good_gen) break;
         }
      }


      if (gen_index == -1) continue;


      Vec<long> fac2;

      if (cofactor) {
         fac2.SetLength(2);
         fac2[0] = m/fac1[gen_index];
         fac2[1] = fac1[gen_index];
      }
      else {
         fac2.SetLength(k);
         fac2 = fac1;

         for (long i = gen_index-1; i >= 0; i--)
           swap(fac2[i], fac2[i+1]);
      }

      long k2 = fac2.length();

      Vec<long> genvec, ordvec;
      genvec.SetLength(k2);
      ordvec.SetLength(k2);

      if (first_gen == 0) first_gen = 1;
      genvec[0] = first_gen;
      ordvec[0] = phi_N(fac2[0])/d;
      if (!good_gen) ordvec[0] = -ordvec[0];


      for (long i = 1; i < k2; i++) {
         vector<long> gens, ords;
         findGenerators(gens, ords, fac2[i], 1);
         assert(gens.size() == 1 && ords.size() == 1 &&
                ords[0] > 0);
         genvec[i] = gens[0];
         ordvec[i] = ords[0];
      }

      Vec<long> global_gen, crtvec;
      global_gen.SetLength(k2);
      crtvec.SetLength(k2);

      long all_1 = 0;

      for (long i = 0; i < k2; i++) {
        crtvec[i] = (m/fac2[i]) * InvMod((m/fac2[i]) % fac2[i], fac2[i]);
        all_1 += crtvec[i];
      }

      for (long i = 0; i < k2; i++) {
         global_gen[i] = (all_1 - crtvec[i] + crtvec[i]*genvec[i]) % m;
      }

      cout << setw(6) << phim << "  ";
      cout << setw(4) << d << "  ";
      cout << setw(6) << m << "  ";

      cout << "(" << NumTwos(conv<ZZ>(d)) << ") ";

      cout << "m=";
      for (long i = 0; i < k; i++) {
        if (i > 0) cout << "*";
        if (i == gen_index) { cout << (cofactor ? "[" : "("); }
        cout << fac[i].a;
        if (fac[i].b > 1) cout << "^" << fac[i].b;
        if (i == gen_index) { cout << (cofactor ? "]" : ")"); }
        if (i == gen_index && !good_gen) cout << "!";
      }
      cout << " ";


      cout << "m/phim(m)=" 
           << double(long(100*double(m)/double(phim)))/100.0 << " ";



      cout << "fac=" << fac2 << " ";
      cout << "gen=" << global_gen << " "; 
      cout << "ord=" << ordvec << " "; 


      cout << "\n";
      cout.flush();
   }


}




// params_x | sort -snk3 | sort -snk2 | sort -snk1 > params.txt
