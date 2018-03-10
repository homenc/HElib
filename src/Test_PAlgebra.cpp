/* Copyright (C) 2012-2017 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
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
  ArgMapping amap;

  long m=17;
  amap.arg("m", m, "cyclotomic index");
  amap.note("e.g., m=2047");

  long p=2;
  amap.arg("p", p, "plaintext base");

  long r=1;
  amap.arg("r", r,  "lifting");

  Vec<long> gens0;
  amap.arg("gens", gens0, "use specified vector of generators", NULL);
  amap.note("e.g., gens='[562 1871 751]'");

  Vec<long> ords0;
  amap.arg("ords", ords0, "use specified vector of orders", NULL);
  amap.note("e.g., ords='[4 2 -4]', negative means 'bad'");

  amap.parse(argc, argv);

  cout << "m = " << m << ", p = " << p <<  ", r = " << r << endl;

  vector<long> f;
  factorize(f,m);
  cout << "factoring "<<m<<" gives [";
  for (unsigned long i=0; i<f.size(); i++)
    cout << f[i] << " ";
  cout << "]\n";

  vector<long> gens1, ords1;
  convert(gens1, gens0);
  convert(ords1, ords0);

  PAlgebra al(m, p, gens1, ords1);
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
