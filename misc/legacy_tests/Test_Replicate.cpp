/* Copyright (C) 2012-2019 IBM Corp.
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

/* Test_Replicate.cpp - Testing the functionality of replicating one
 * slot from a vector acress the whole vector (or replicating each slot
 * to a full cipehrtext)
 */

#include <cassert>
#include <NTL/lzz_pXFactoring.h>
NTL_CLIENT

#include <helib/helib.h>
#include <helib/replicate.h>
#include <helib/ArgMap.h>

using namespace helib;

static bool noPrint = false;

static bool check_replicate(const Ctxt& c1, const Ctxt& c0, long i,
			    const SecKey& sKey, const EncryptedArray& ea)
{
  PlaintextArray pa0(ea), pa1(ea);
  ea.decrypt(c0, sKey, pa0);
  ea.decrypt(c1, sKey, pa1);
  replicate(ea, pa0, i);

  return equals(ea, pa1, pa0); // returns true if replication succeeded
}


class StopReplicate { };

// A class that handles the replicated ciphertexts one at a time
class ReplicateTester : public ReplicateHandler {
public:
  const SecKey& sKey;
  const EncryptedArray& ea;
  const PlaintextArray& pa;
  long B;

  double t_last, t_total;
  long pos;
  bool error;

  ReplicateTester(const SecKey& _sKey, const EncryptedArray& _ea,
                  const PlaintextArray& _pa, long _B)
  : sKey(_sKey), ea(_ea), pa(_pa), B(_B)
  {
    t_last = GetTime();
    t_total = 0.0;
    pos = 0;
    error = false;
  }

  // This method is called for every replicated ciphertext: in the i'th time
  // that it is called, the cipehrtext will have in all the slots the content
  // of the i'th input slot. In this test program we only decrypt and check
  // the result, in a real program it will do something with the cipehrtext.
  virtual void handle(const Ctxt& ctxt) {

    double t_new = GetTime();
    double t_elapsed = t_new - t_last;
    t_total += t_elapsed;

    // Decrypt and check
    PlaintextArray pa1 = pa;
    replicate(ea, pa1, pos);
    PlaintextArray pa2(ea);

    if (pos==0 && !noPrint) CheckCtxt(ctxt, "replicateAll");

    ea.decrypt(ctxt, sKey, pa2);
    if (!equals(ea, pa1, pa2)) error = true; // record the error, if any
    t_last = GetTime();

    pos++;
    if (B > 0 && pos >= B) throw StopReplicate();
  }
};



void  TestIt(long m, long p, long r, long d, long L, long bnd, long B)
{
  if (!noPrint)
    std::cout << "*** TestIt" << (isDryRun()? "(dry run):" : ":")
	 << " m=" << m
	 << ", p=" << p
	 << ", r=" << r
	 << ", d=" << d
	 << ", L=" << L
	 << ", bnd=" << bnd
	 << ", B=" << B
	 << endl;

  setTimersOn();
  Context context(m, p, r);
  buildModChain(context, L, /*c=*/2);

  ZZX G;
  if (d == 0)
    G = context.alMod.getFactorsOverZZ()[0];
  else
    G = makeIrredPoly(p, d);

  if (!noPrint) {
    context.zMStar.printout();
    cout << endl;
    cout << "G = " << G << "\n";
  }

  SecKey secretKey(context);
  const PubKey& publicKey = secretKey;
  secretKey.GenSecKey(); // A +-1/0 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need

  EncryptedArray ea(context, G);
  PlaintextArray xp0(ea), xp1(ea);
  random(ea, xp0);
  random(ea, xp1);

  Ctxt xc0(publicKey);
  ea.encrypt(xc0, publicKey, xp0);

  ZZX poly_xp1;
  ea.encode(poly_xp1, xp1);

  if (!noPrint)  cout << "** Testing replicate():\n";
  bool error = false;
  Ctxt xc1 = xc0;
  if (!noPrint) CheckCtxt(xc1, "before replicate");
  replicate(ea, xc1, ea.size()/2);
  if (!check_replicate(xc1, xc0, ea.size()/2, secretKey, ea)) error = true;
  if (!noPrint) CheckCtxt(xc1, "after replicate");

  // Get some timing results
  for (long i=0; i<20 && i<ea.size(); i++) {
    xc1 = xc0;
    HELIB_NTIMER_START(replicate);
    replicate(ea, xc1, i);
    if (!check_replicate(xc1, xc0, i, secretKey, ea)) error = true;
    HELIB_NTIMER_STOP(replicate);
  }
  cout << (error? "BAD" : "GOOD") << endl;

  if (!noPrint) {
    printAllTimers();
    cout << "\n** Testing replicateAll()... " << std::flush;
  }
#ifdef HELIB_DEBUG
  replicateVerboseFlag = true;
#else
  replicateVerboseFlag = false;
#endif

  error = false;
  ReplicateTester *handler = new ReplicateTester(secretKey, ea, xp0, B);
  try {
    HELIB_NTIMER_START(replicateAll);
    replicateAll(ea, xc0, handler, bnd);
  }
  catch (StopReplicate) {
  }
  std::cout << (handler->error? "BAD" : "GOOD") << endl;
  if (!noPrint)
    cout << "  total time=" << handler->t_total << " ("
         << ((B>0)? B : ea.size()) << " vectors)\n";
  delete handler;
}

int main(int argc, char *argv[])
{
  ArgMap amap;

  bool dry=false;
  amap.arg("dry", dry, "dry=1 for a dry-run");

  long m=2047;
  amap.arg("m", m, "cyclotomic ring");

  long p=2;
  amap.arg("p", p, "plaintext base");

  long r=1;
  amap.arg("r", r,  "lifting");

  long d=1;
  amap.arg("d", d, "degree of the field extension");
  amap.note("d == 0 => factors[0] defines extension");

  long L=250;
  amap.arg("L", L, "# of bits in the modulus chain");

  long bnd = 64;
  amap.arg("bnd", bnd, "recursion bound for replication");

  long B = 0;
  amap.arg("B", B, "bound for # of replications", "all");

  amap.arg("noPrint", noPrint, "suppress printouts");

  amap.parse(argc, argv);
  setDryRun(dry);

  TestIt(m, p, r, d, L, bnd, B);
  cout << endl;
}
