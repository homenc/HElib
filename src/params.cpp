
namespace std {}
namespace NTL {}

using namespace std;
using namespace NTL;

#include "NumbTh.h"
#include "PAlgebra.h"
#include <iomanip>

Vec<long> rev(const Vec<long>& v) 
{
   long n = v.length();
   Vec<long> w;
   w.SetLength(n);
   for (long i = 0; i < n; i++) w[i] = v[n-1-i];
   return w;
}

Vec<long> trunc(const Vec<long>& v)
{
   long n = v.length();
   Vec<long> w = v;
   if (n > 0 && w[n-1] == 1) w.SetLength(n-1);
   return w;
}

const long MaxOrd = 100;

long computePhi(const Pair<long,long>& x)
{
   return power_long(x.a, x.b-1) * (x.a-1);
}

bool comparePhi(const Pair<long,long>& x, const Pair<long,long>& y)
{
   return computePhi(x) < computePhi(y);
} 


int main(int argc, char *argv[])
{
   argmap_t argmap;
   argmap["gens"] = "0";
   if (!parseArgs(argc, argv, argmap)) Error("bad args");

   long gens_flag = atoi(argmap["gens"]);


   for (long m = 1001; m <= 30000; m += 2) {

      long d = multOrd(2, m);

      if (d > MaxOrd) continue;

      long phim = phi_N(m);

      Vec< Pair<long, long> > fac;
      factorize(fac, m);

      long k = fac.length();

      if (k == 1) continue;

      Vec<long> fac1;
      fac1.SetLength(k);

      for (long i = 0; i < k; i++)
         fac1[i] = power_long(fac[i].a, fac[i].b);

      Vec<long> phivec;
      phivec.SetLength(k);
      for (long i = 0; i < k; i++)
         phivec[i] = computePhi(fac[i]);

      long phisum = 0;
      for (long i = 0; i < k; i++)
         phisum += phivec[i];

      long gen_index = -1;
      bool good_gen = false;
      long first_gen = 0;

      long best_weighted_cost = NTL_MAX_LONG;
      long best_cost = NTL_MAX_LONG;
      long best_depth = NTL_MAX_LONG;

      for (long i = 0; i < k; i++) {
         long m1 = fac1[i];
         long phim1 = phivec[i];
         if (multOrd(2, m1) != d) continue;

         PAlgebra pal1(m1);
         if (pal1.numOfGens() > 1) continue;

         bool good = (pal1.numOfGens() == 0 || 
                      (pal1.numOfGens() == 1 && pal1.SameOrd(0)));

         long cost = phisum - phivec[i] + d-1;
         long depth = k-1;
         cost += (2-long(good))*(phim1/d-1);
         depth += (2-long(good));

         if (cost*depth < best_weighted_cost) {
            gen_index = i;
            good_gen = good;
            first_gen = pal1.ZmStarGen(0);
            best_weighted_cost = cost*depth;
            best_cost = cost;
            best_depth = depth;
         }
      }

      long gen_index2 = -1;

      // search for two-generator solution

      for (long i = 0; i < k; i++) {
         for (long j = i+1; j < k; j++) {
            long m1 = fac1[i]*fac1[j];
            long phim1 = phivec[i]*phivec[j];
            if (multOrd(2, m1) != d) continue;

            PAlgebra pal1(m1);
            if (pal1.numOfGens() > 1) continue;

            bool good = (pal1.numOfGens() == 0 || 
                         (pal1.numOfGens() == 1 && pal1.SameOrd(0)));

            long cost = phisum - (phivec[i] + phivec[j]) + d-1;
            long depth = k-2;
            cost += (2-long(good))*(phim1/d-1);
            depth += (2-long(good));

            if (cost*depth < best_weighted_cost) {
               gen_index = i;
               gen_index2 = j;
               good_gen = good;
               first_gen = pal1.ZmStarGen(0);
               best_weighted_cost = cost*depth;
               best_cost = cost;
               best_depth = depth;
            }
         }
      }


      if (gen_index == -1) continue;

      if (best_cost > sqrt(phim)) continue;


      Vec<long> fac2;

      fac2 = fac1;

      for (long i = gen_index-1; i >= 0; i--)
        swap(fac2[i], fac2[i+1]);

      if (gen_index2 != -1) {
         for (long i = gen_index2; i < k-1; i++)
            swap(fac2[i], fac2[i+1]);
         fac2[0] *= fac2[k-1];
         fac2.SetLength(k-1);
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


      cout << "m=";
      for (long i = 0; i < k; i++) {
        if (i > 0) cout << "*";
        if (i == gen_index) cout << "(";
        if (i == gen_index2) cout << "{";
        cout << fac[i].a;
        if (fac[i].b > 1) cout << "^" << fac[i].b;
        if (i == gen_index) cout << ")";
        if (i == gen_index2) cout << "}";
      }
      cout << " ";

      if (!good_gen) 
         cout << ":-( ";

      cout << "m/phim(m)=" 
           << double(long(100*double(m)/double(phim)))/100.0 << " ";

      cout << "C=" << best_cost << " ";
      cout << "D=" << best_depth << " ";
      cout << "E=" << NumTwos(conv<ZZ>(d)) << " ";

      if (gens_flag) {
         cout << "mvec=" << "\"" << rev(fac2) << "\"" << " ";
         cout << "gens=" << "\"" << trunc(rev(global_gen)) << "\"" << " "; 
         cout << "ords=" << "\"" << trunc(rev(ordvec)) << "\"" << " "; 
      }


      cout << "\n";
      cout.flush();
   }


}




// params_x | sort -snk3 | sort -snk2 | sort -snk1 > params.txt
