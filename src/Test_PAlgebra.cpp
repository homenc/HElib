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
#include <cassert>
#include <string>
#include <sstream>
#include <NTL/ZZ.h>
NTL_CLIENT
#include "NumbTh.h"
#include "FHEContext.h"

void usage() 
{
  cout << "Usage: Test_PAlgebra_x m=<int> [ p=<int> ] [ r=<int> ]" << endl;
  cout << "  m is an integer that determines Z_m^*" << endl;
  cout << "  p is an integer that determines the plaintext base [default=2]" << endl;
  cout << "  r is an integer that determines the lifting [default=1]" << endl;
  exit(0);
}

int main(int argc, char *argv[]) 
{

  argmap_t argmap;
  argmap["m"] = "0"; 
  argmap["p"] = "2";
  argmap["r"] = "1";

  if (!parseArgs(argc, argv, argmap)) usage();

  unsigned long m = atoi(argmap["m"]);

  if (!m) usage();

  unsigned long p = atoi(argmap["p"]);

  unsigned long r = atoi(argmap["r"]);

  cout << "m = " << m << ", p = " << p <<  ", r = " << r << endl;

  vector<long> f;
  factorize(f,m);
  cout << "factoring "<<m<<" gives [";
  for (unsigned long i=0; i<f.size(); i++)
    cout << f[i] << " ";
  cout << "]\n";

  PAlgebra al(m, p);
  al.printout();
  cout << "\n";

  PAlgebraMod almod(al, r);

  FHEcontext context(m, p, r);
  buildModChain(context, 5, 2);

  stringstream s1;
  writeContextBase(s1, context);
  s1 << context;

  string s2 = s1.str();

  cout << s2 << endl;

  stringstream s3(s2);

  unsigned long m1, p1, r1;
  vector<long> gens, ords;
  readContextBase(s3, m1, p1, r1, gens, ords);

  FHEcontext c1(m1, p1, r1, gens, ords);
  s3 >> c1;

  if (context == c1)
    cout << "equal\n";
  else
    cout << "not equal\n";

  return 0;
}
