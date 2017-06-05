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
/* Test_matmul.cpp - Testing the functionality of multiplying an encrypted
 * vector by a plaintext matrix, either over the extension- or the
 * base-field/ring.
 */
#include <cassert>
#include <NTL/lzz_pXFactoring.h>
#include "multiAutomorph.h"
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"

// defined in debugging.cpp
void decryptAndPrint(ostream& s, const Ctxt& ctxt, const FHESecKey& sk,
		     const EncryptedArray& ea, long flags);

static void checkAuto(const Ctxt& cOrig, const Ctxt& cAuto, long amt,
                      EncryptedArray& ea, const FHESecKey& sk)
{
  NewPlaintextArray v1(ea), v2(ea);
  Ctxt cTmp = cOrig;
  cTmp.smartAutomorph(amt);

  ea.decrypt(cTmp, sk, v1);
  ea.decrypt(cAuto, sk, v2);

  if (!equals(ea, v1, v2)) { // check that we've got the right answer
    cout << " k = "<<amt<<"failed, Grrr@*\n";
    exit(0);
  }
}

class AutoTester: public AutomorphHandler {
public:
  const Ctxt& cOrig;
  FHESecKey& sk;
  EncryptedArray& ea;

    AutoTester(const Ctxt& c, FHESecKey& k, EncryptedArray& e, bool v=false):
    cOrig(c), sk(k), ea(e) {}

  // cPtr points to the original ciphertext, ctxtx is after automorphism
  bool handle(std::unique_ptr<Ctxt>& ctxt, long amt) override {
    // check that ctxt is indeed the original ctxt after automorphism
    checkAuto(cOrig, *ctxt, amt, ea, sk);

    return true;
  }
};

void  TestIt1(FHESecKey& secretKey, EncryptedArray& ea, bool verbose=false)
{
  const FHEcontext& context = ea.getContext();
  const FHEPubKey& publicKey = secretKey;

  // choose a random plaintext vector
  NewPlaintextArray v(ea);
  random(ea, v);

  // encrypt the random vector
  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, v);
  ctxt.square();
  ctxt.cube();

  AutoTester test(ctxt, secretKey, ea);
  for (long i=0; i<=ea.dimension(); i++) {
    const AutGraph& tree = publicKey.getTree4dim(i);
    multiAutomorph(ctxt, tree, test);
  }
  cout << "  All tests using handler passed successfully\n";
}

void  TestIt2(FHESecKey& secretKey, EncryptedArray& ea, bool verbose=false)
{
  const FHEcontext& context = ea.getContext();
  const FHEPubKey& publicKey = secretKey;

  // choose a random plaintext vector
  NewPlaintextArray v(ea);
  random(ea, v);

  // encrypt the random vector
  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, v);
  ctxt.square();
  ctxt.cube();

  Ctxt tmp(ZeroCtxtLike, ctxt);
  for (long i=0; i<=ea.dimension(); i++) {
    const AutGraph& tree = publicKey.getTree4dim(i);
    std::unique_ptr<AutoIterator> it(AutoIterator::build(ctxt, tree));
    while (long val = it->next(tmp))
      checkAuto(ctxt, tmp, val, ea, secretKey);
  }
  cout << "  All tests using iterator passed successfully\n";
}



/* Testing the new automorphism

 * Usage: Test_newAutomorph_x [optional params]
 *
 *  m defines the cyclotomic polynomial Phi_m(X) [default=2047]
 *    another useful setting to test is m=4369
 *  p is the plaintext base [default=2]
 *  L is the # of primes in the modulus chain [default=4]
 *  verbose print extra info [default=0]
 */
int main(int argc, char *argv[]) 
{
  ArgMapping amap;

  long p=2;
  amap.arg("p", p, "plaintext base");
  long m=2047;
  amap.arg("m", m, "defines the cyclotomic polynomial Phi_m(X)");
  amap.note("another useful setting to test is m=4369, p=2");
  long L=15;
  amap.arg("L", L, "# of levels in the modulus chain");
  bool verbose=false;
  amap.arg("verbose", verbose, "print extra information");
  amap.parse(argc, argv);

  cout << "*** "<<argv[0]
       << ": m=" << m
       << ", p=" << p
       << ", L=" << L
       << endl;

  FHEcontext context(m, p, 1);
  buildModChain(context, L, /*c=*/3);
    
  FHESecKey secretKey(context);
  secretKey.GenSecKey(/*w=*/64); // A Hamming-weight-w secret key

  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need
  EncryptedArray ea(context, context.alMod);

  if (verbose) {
    context.zMStar.printout();
    cout << endl;
    for (long i=0; i<=ea.dimension(); i++) {
      cout << "Tree("<<i<<") =\n";
      const AutGraph& tree = secretKey.getTree4dim(i);
      for (auto x: tree) {
        cout << "  "<< x.first<<": ";
        cout << x.second << endl;
      }
    }
  }

  resetAllTimers();  
  TestIt1(secretKey, ea, verbose);
  if (verbose) {
    printAllTimers();
    cout << endl;
  }

  resetAllTimers();
  TestIt2(secretKey, ea, verbose);
  if (verbose) {
    printAllTimers();
    cout << endl;
  }
}
