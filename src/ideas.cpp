
#include "NumbTh.h"
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2EXFactoring.h>

// 2^9 - 1 = 7*73, 2 has order 9 mod 73


// sorting stuff...stolen from PAlgebra.cpp

template<class RX> bool poly_comp(const RX& a, const RX& b) 
{
  long na = deg(a) + 1;
  long nb = deg(b) + 1;

  long i = 0;
  while (i < na && i < nb && coeff(a, i) == coeff(b, i)) i++;

  if (i < na && i < nb)
    return coeff(a, i) < coeff(b, i);
  else 
    return na < nb;
}

namespace NTL {

// for some weird reason, these need to be in either the std or NTL 
// namespace; otherwise, the compiler won't find them...

bool operator<(GF2 a, GF2 b) { return rep(a) < rep(b); }

bool operator<(const GF2X& a, const GF2X& b) { return poly_comp(a, b); }


}




int main()
{

   long m1 = 73;
   long m2 = 7;
   long m = m1*m2;

   GF2X Phi_m1 = conv<GF2X>(Cyclotomic(m1));

   Vec<GF2X> fac1 = SFCanZass(Phi_m1);

   cout << fac1 << "\n";
   long l1 = fac1.length();

   // Phi_m1 factors as 8 factors, each of degree 9

   GF2X Phi_m2 = conv<GF2X>(Cyclotomic(m2));
   Vec<GF2X> roots(INIT_SIZE, l1);
   
   for (long i = 0; i < l1; i++) {
      GF2E::init(fac1[i]);
      roots[i] = rep(FindRoot(conv<GF2EX>(Phi_m2)));
   }

   cout << roots << "\n";

   // roots[i] is a root of Phi_m2 modulo fac1[i] 


   // now we construct a random element of the m-th cyclotomic
   // ring, represented on the powerful basis.

   Mat<GF2> M(INIT_SIZE, m2-1, m1-1); // matrix with m2-1 rows, m1-1 cols

   for (long i = 0; i < m2-1; i++)
      for (long j = 0; j < m1-1; j++) {
         random(M[i][j]);
         random(M[i][j]);
         random(M[i][j]);
         random(M[i][j]);
      }

   // next, we implement the isomorphism from the powerful basis
   // to the traditional basis:

   GF2X Phi_m = conv<GF2X>(Cyclotomic(m));

   GF2X g;
   for (long i = 0; i < m2-1; i++)
      for (long j = 0; j < m1-1; j++) {
         g += GF2X(i*m1 + j*m2, M[i][j]);
      }

   g = g % Phi_m;

   Vec<GF2X> fac = SFCanZass(Phi_m);
   Vec<GF2X> slots(INIT_SIZE, fac.length());
   for (long i = 0; i < fac.length(); i++)
      slots[i] = g % fac[i];

   // there are 48 slots, each of degree 9 

   cout << slots.length() << "\n\n";

   
   // for now, we sort them lexicographically, because I'm not
   // sure of the correspondence

   sort(slots.begin(), slots.end());
   cout << slots << "\n\n";

   // OK, now we compute the same slots, but using the powerful 
   // basis directly

   Vec<GF2X> slots1(INIT_SIZE, l1*(m2-1));

   for (long j = 0; j < l1; j++) {
      GF2X f = fac1[j];
      GF2X r = roots[j];

      for (long k = 0; k < m2-1; k++) {

         // r = roots[j]^{k+1}

         // Horner eval
         GF2X acc;
         for (long i = m2-2; i >= 0; i--) {
            GF2X c = conv<GF2X>(M[i]) % f;
            acc = (acc * r + c) % f;
         }

         slots1[j*(m2-1) + k] = acc;

         r = (r * roots[j]) % f;
      }
   }

   sort(slots1.begin(), slots1.end());
   cout << slots1 << "\n";
   
   if (slots == slots1) 
      cout << "yes!!\n";
   else
      cout << "NOOO!!!\n";

   for (long i = 0; i < slots.length(); i++) 
      for (long j = 0; j < slots1.length(); j++)
         if (slots[i] == slots1[j])
            cout << j << " " << slots[i] << "\n";
}
