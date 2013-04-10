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
#include "NumbTh.h"
#include "CModulus.h"

int main(int argc, char *argv[]) 
{
  if (argc<3) {
    cout << "\nUsage: " << argv[0] << " m |q|\n\n";
    cout << "  m is an integer that determines Z_m^*\n\n";
    cout << "  |q|<31 is the modulus bit-length, a prime q of this length\n"
	 << "      is selected such that Z_q^* has 2m-th roots of unity\n\n"; 
    exit(0);
  }
  unsigned m = atoi(argv[1]);
  unsigned l = atoi(argv[2]);
  if (l>30) Error("q bit-length cannot be more than 30");

  PAlgebra al(m);

  unsigned Qmax = (1UL<<(l+1))/(2*m);
  long q;
  for (int i=0; i<1000; i++) {
    q = RandomBnd(Qmax)*2*m +1;
    if (ProbPrime(q)) break;
  }
  Cmodulus mod(al,q,0);

  cout << "m=" << mod.getM()
       << ", q=" << mod.getQ()
       << ", root=" << mod.getRoot() << endl;

  vec_long v1,v2;
  ZZX poly1, poly2;
  { mod.restoreModulus();  // set poly1 = a random polynomial mod q
    zz_pX rpoly; random(rpoly, al.phiM());
    conv(poly1, rpoly);
  }

  mod.FFT(v1,poly1);  // v1 = FFT(poly1)
  mod.iFFT(poly2,v1); // poly2 = FFT^{-1}(v1)
  cout << "iFFT(FFT(poly)) is "
       << ((poly1==poly2)? "successful" : "UNsuccessful") << endl;

  mod.iFFT(poly1,v1); // poly2 = FFT^{-1}(v1)
  mod.FFT(v2,poly1);   // v1 = FFT(poly1)
  cout << "FFT(iFFT(v)) is "
       << ((v1==v2)? "successful" : "UNsuccessful") << endl;
}

