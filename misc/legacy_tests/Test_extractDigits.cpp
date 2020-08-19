/* Copyright (C) 2012-2020 IBM Corp.
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
/* Test_extractDigits.cpp - extracting digits.
 *   For a plaintext space modulo a prime-power $p^e$, extracting
 *   the base-$p$ representation of an encrypted values.
 */
#include <NTL/ZZ.h>
NTL_CLIENT
#include <helib/EncryptedArray.h>
#include <helib/polyEval.h>
#include <helib/debugging.h>
#include <helib/ArgMap.h>

using namespace helib;

static bool noPrint = false;
static bool debug = 0;


int main(int argc, char *argv[])
{
  ArgMap amap;
  long p = 5;
  long r=0;
  long m = 0;
  bool dry = false;

  amap.arg("p", p, "plaintext base");
  amap.arg("r", r, "lifting");
  amap.arg("m", m, "the cyclotomic ring", "heuristic");
  amap.arg("noPrint", noPrint, "suppress printouts");
  amap.arg("debug", debug, "extra debugging");
  amap.arg("dry", dry, "dry=1 for a dry-run");

  // get parameters from the command line
  amap.parse(argc, argv);

  if (p<2) exit(0);
  double lBound = 30.0;
  long bound = floor(lBound/log2((double)p));
  if (r<2 || r>bound) r = bound;
  long p2r = power_long(p,r); // p^r

  long ll = NextPowerOfTwo(p);
  long L = (r*ll*3 +2)*30; // how many levels do we need
  if (m==0)
    m = p+1; // FindM(/*secparam=*/80, L, /*c=*/4, p, /*d=*/1, 0, m);
  setDryRun(dry);

  if (!noPrint) {
    if (dry) cout << "dry run: ";
    cout << "m="<<m<<", p="<<p<<", r="<<r<<", L="<<L<<endl;
  }
  Context context(m, p, r);
  buildModChain(context, L, /*c=*/4);

  SecKey secretKey(context);
  const PubKey& publicKey = secretKey;
  secretKey.GenSecKey(); // A +-1/0 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need

  if (debug) {
      dbgKey = &secretKey; // debugging key and ea
      dbgEa = context.ea;
  }
#ifdef HELIB_DEBUG
  dbgKey = &secretKey; // debugging key and ea
  dbgEa = context.ea;
#endif //HELIB_DEBUG

  EncryptedArray ea(context);
  vector<long> v;
  vector<long> pDigits;
  ea.random(v); // random values in the slots

  Ctxt c(publicKey);
  ea.encrypt(c, publicKey, v);
  ea.decrypt(c, secretKey, pDigits);
  if (ea.size()<=20 && !noPrint)
    cout << "plaintext="<<pDigits<<endl;


  if (!noPrint)
    cout << "extracting "<<r<<" digits..." << std::flush;
  vector<Ctxt> digits;
  extractDigits(digits, c);
  if (!noPrint)
    cout << " done\n" << std::flush;

  vector<long> tmp = v;
  long pp = p2r;
  for (long i=0; i<(long)digits.size(); i++) {
    if (!digits[i].isCorrect()) {
      cout << " BAD, potential decryption error for "<<i<<"th digit ";
      CheckCtxt(digits[i], "");
      exit(0);
    }
    ea.decrypt(digits[i], secretKey, pDigits);
    if (ea.size()<=20 && !noPrint)
      cout << i << "th digit="<<pDigits<<endl;

    // extract the next digit from the plaintext, compare to pDigits
    for (long j=0; j<(long)v.size(); j++) {
      long digit = tmp[j] % p;
      if (digit > p/2) digit -= p;
      else if (digit < -p/2) digit += p;

      // assert ((pDigits[j]-digit) % p == 0);
      if ((pDigits[j]-digit) % pp != 0) {
        cout << "BAD\n";
        if (!noPrint)
          cout << " error: v["<<j<<"]="<<v[j]
               << " but "<<i<<"th digit comes "<< pDigits[j]
               << " rather than "<<digit<<endl<<endl;
	exit(0);
      }
      tmp[j] -= digit;
      tmp[j] /= p;
    }
    pp /= p;
  }
  cout << "GOOD\n";
}
