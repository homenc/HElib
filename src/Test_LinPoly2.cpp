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
#include <NTL/mat_GF2E.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/ZZXFactoring.h>
#include "NumbTh.h"
#include <cassert>

NTL_CLIENT


void usage(char *prog)
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  m defines the cyclotomic polynomial Phi_m(X) [default=105]"<< endl;
  exit(0);
}

int main(int argc, char *argv[])
{
   argmap_t argmap;
   argmap["m"] = "105";

   if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

   long m = atoi(argmap["m"]);

   cout << "m=" << m << "\n";

   ZZX phimx = Cyclotomic(m);

   GF2X phimx_modp = conv<GF2X>(phimx);

   

   vec_GF2X factors = SFCanZass(phimx_modp);



   GF2X G = factors[0];
   long d = deg(G);
   cout << "d=" << d << "\n";

   GF2E::init(G);

   // L selects the even coefficients
   vec_GF2E L;
   L.SetLength(d);
   L[0] = 1;
/*
   for (long j = 0; j < d; j++) {
      if (j % 2 == 0)
         L[j] = to_GF2E(GF2X(j, 1));
   }
*/

   vec_GF2E C;
   buildLinPolyCoeffs(C, L, 2, 1);
  
   cout << C << "\n";

   GF2E alpha, beta;
   random(alpha);

   applyLinPoly(beta, C, alpha, 2);

   cout << alpha << "\n";
   cout << beta << "\n";

}
