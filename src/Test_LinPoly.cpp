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
#include <NTL/mat_lzz_pE.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/ZZXFactoring.h>
#include "NumbTh.h"
#include <cassert>

NTL_CLIENT


void usage(char *prog)
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'p=5 r=1 m=101'\n\n";
  cerr << "  p is the plaintext base [default=5]" << endl;
  cerr << "  r is the lifting [default=1]" << endl;
  cerr << "  m defines the cyclotomic polynomial Phi_m(X) [default=101]"<< endl;
  exit(0);
}

int main(int argc, char *argv[])
{
   argmap_t argmap;
   argmap["p"] = "5";
   argmap["m"] = "101";
   argmap["r"] = "1";

   if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

   long p = atoi(argmap["p"]);
   long m = atoi(argmap["m"]);
   long r = atoi(argmap["r"]);

   cout << "p=" << p << ", m=" << m << ", r=" << r << "\n";

   ZZX phimx = Cyclotomic(m);

   zz_p::init(p);
   zz_pX phimx_modp = to_zz_pX(phimx);

   

   vec_zz_pX factors = SFCanZass(phimx_modp);


   vec_ZZX FFactors;
   MultiLift(FFactors, factors, phimx, r);

   zz_p::init(power_long(p, r));

   vec_zz_pX Factors;
   Factors.SetLength(FFactors.length());
   for (long i = 0; i < Factors.length(); i++)
      conv(Factors[i], FFactors[i]);

   zz_pX G = Factors[0];
   long d = deg(G);
   cout << "d=" << d << "\n";

   zz_pE::init(G);

   // L selects the even coefficients
   vec_zz_pE L;
   L.SetLength(d);
   for (long j = 0; j < d; j++) {
      if (j % 2 == 0)
         L[j] = to_zz_pE(zz_pX(j, 1));
   }

   vec_zz_pE C;
   buildLinPolyCoeffs(C, L, p, r);

   zz_pE alpha, beta;
   random(alpha);

   applyLinPoly(beta, C, alpha, p);

   cout << alpha << "\n";
   cout << beta << "\n";

}
