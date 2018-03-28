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
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"

void printBits(const vector<ZZX>& v, long n) 
{
  long len = v.size();
  if (n>50 || len>32) return;
  for (long i = 0; i < len; i++) {
    for (long j = n-1; j >= 0; j--)
      cout << coeff(v[i], j);
    cout << " ";
  }
  cout << "\n";
}


void  TestIt(long c, long k, long w, long L, long m, long n)
{
  FHEcontext context(m, 2, 1); // p = 2, r = 1
  long d = context.zMStar.getOrdP(); 

  buildModChain(context, L, c);

  context.zMStar.printout();
  cerr << endl;
#ifdef DEBUG
  cerr << context << endl;
#endif

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(w); // A Hamming-weight-w secret key


  ZZX G;

  G = makeIrredPoly(2, d); 
  // G = context.alMod.getFactorsOverZZ()[0];

  cerr << "generating key-switching matrices... ";
  addFrbMatrices(secretKey);
  addSome1DMatrices(secretKey);
  cerr << "done\n";

  cerr << "computing masks and tables for rotation...";
  EncryptedArray ea(context, G);
  cerr << "done\n";

  long nslots = ea.size();

  if (n <= 0 || n > d) n = d;

  vector<ZZX> v;
  v.resize(nslots);
  for (long i = 0; i < nslots; i++) {
    GF2X f;
    random(f, n);
    conv(v[i], f);
  }

  printBits(v, n);

  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, v);
  // ctxt encrypts a vector where each slots is a random
  // polynomial of degree < n

  Ctxt* res[n];
  for (long j = 0; j < n; j++) res[j] = new Ctxt(publicKey); // allocate

  resetAllTimers();

  incrementalZeroTest(res, ea, ctxt, n);

  for (long j = 0; j < n; j++) {
    vector<ZZX> v1;
    ea.decrypt(*res[j], secretKey, v1); 
    printBits(v1, n);
  }

  for (long j = 0; j < n; j++) delete res[j]; // cleanup
}

void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'R=4 L=9 k=80'\n\n";
  cerr << "  c is number of columns in the key-switching matrices [default=2]\n";
  cerr << "  k is the security parameter [default=80]\n";
  cerr << "  L is the # of primes in the modulus chain [default=20]\n";
  cerr << "  m is a specific modulus\n";
  cerr << "  n is the number of masks (defaults to all)\n";
  exit(0);
}

int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["c"] = "2";
  argmap["k"] = "80";
  argmap["L"] = "20";
  argmap["m"] = "0";
  argmap["n"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long c = atoi(argmap["c"]);
  long k = atoi(argmap["k"]);
  long L = atoi(argmap["L"]);
  long chosen_m = atoi(argmap["m"]);
  long n = atoi(argmap["n"]);

  long w = 64; // Hamming weight of secret key
  //  long L = z*R; // number of levels

  long m = FindM(k, L, c, 2, 1, 0, chosen_m, true);

  setTimersOn();
  TestIt(c, k, w, L, m, n);

  cerr << endl;
  printAllTimers();
  cerr << endl;

}
// call to get our running test case:
// Test_eqtesting_x m=20485 
// Test_eqtesting_x m=105 for quick testing
