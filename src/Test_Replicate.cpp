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

/* Test_Replicate.cpp - Testing the functionality of replicating one
 * slot from a vector acress the whole vector (or replicating each slot
 * to a full cipehrtext)
 */

#include <cassert>
#include <NTL/lzz_pXFactoring.h>
NTL_CLIENT

#include "FHE.h"
#include "replicate.h"
#include "timing.h"

static bool check_replicate(const Ctxt& c1, const Ctxt& c0, long i,
			    const FHESecKey& sKey, const EncryptedArray& ea)
{
  PlaintextArray pa0(ea), pa1(ea);
  ea.decrypt(c0, sKey, pa0);
  ea.decrypt(c1, sKey, pa1);
  pa0.replicate(i);
  
  return pa1.equals(pa0); // returns true if replication succeeded
}


class StopReplicate { };

// A class that handles the replicated ciphertexts one at a time
class ReplicateTester : public ReplicateHandler {
public:
  const FHESecKey& sKey;
  const EncryptedArray& ea;
  const PlaintextArray& pa;
  long B;

  double t_last, t_total;
  long pos;
  bool error;

  ReplicateTester(const FHESecKey& _sKey, const EncryptedArray& _ea, 
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
    pa1.replicate(pos);
    PlaintextArray pa2(ea);

    ea.decrypt(ctxt, sKey, pa2);
    if (!pa1.equals(pa2)) error = true; // record the error, if any
    t_last = GetTime();

    pos++;
    if (B > 0 && pos >= B) throw StopReplicate();
  }
};



void  TestIt(long m, long p, long r, long d, long L, long bnd, long B)
{
  cout << "*** TestIt" << (isDryRun()? "(dry run):" : ":")
       << " m=" << m
       << ", p=" << p
       << ", r=" << r
       << ", d=" << d
       << ", L=" << L
       << ", bnd=" << bnd
       << ", B=" << B
       << endl;

  setTimersOn();
  FHEcontext context(m, p, r);
  buildModChain(context, L, /*c=*/2);

  context.zMStar.printout();
  cout << endl;

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(/*w=*/64); // A Hamming-weight-w secret key

  ZZX G;
  if (d == 0)
    G = context.alMod.getFactorsOverZZ()[0];
  else
    G = makeIrredPoly(p, d); 

  cout << "G = " << G << "\n";
  cout << "generating key-switching matrices... ";
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  cout << "done\n";

  cout << "computing masks and tables for rotation...";
  EncryptedArray ea(context, G);
  cout << "done\n";

  PlaintextArray xp0(ea), xp1(ea);
  xp0.random();
  xp1.random();

  Ctxt xc0(publicKey);
  ea.encrypt(xc0, publicKey, xp0);

  ZZX poly_xp1;
  ea.encode(poly_xp1, xp1);

  cout << "** Testing replicate():\n";
  bool error = false;
  Ctxt xc1 = xc0;
  CheckCtxt(xc1, "before replicate");
  replicate(ea, xc1, ea.size()/2);
  if (!check_replicate(xc1, xc0, ea.size()/2, secretKey, ea)) error = true;
  CheckCtxt(xc1, "after replicate");

  // Get some timing results
  for (long i=0; i<20 && i<ea.size(); i++) {
    xc1 = xc0;
    FHE_NTIMER_START(replicate);
    replicate(ea, xc1, i);
    if (!check_replicate(xc1, xc0, i, secretKey, ea)) error = true;
    FHE_NTIMER_STOP(replicate);
  }
  cout << "  Replicate test " << (error? "failed :(\n" : "succeeded :)")
       << endl<< endl;
  printAllTimers();

  cout << "\n** Testing replicateAll()... " << std::flush;
#ifdef DEBUG_PRINTOUT
  replicateVerboseFlag = true;
#else
  replicateVerboseFlag = false;
#endif

  error = false;
  ReplicateTester *handler = new ReplicateTester(secretKey, ea, xp0, B);
  try {
    FHE_NTIMER_START(replicateAll);
    replicateAll(ea, xc0, handler, bnd);
  }
  catch (StopReplicate) {
  }
  cout << (handler->error? "failed :(\n" : "succeeded :)")
       << ", total time=" << handler->t_total << " ("
       << ((B>0)? B : ea.size())
       << " vectors)\n";
  delete handler;
}


int main(int argc, char *argv[]) 
{
  ArgMapping amap;

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

  long L=3;
  amap.arg("L", L, "# of levels in the modulus chain",  "heuristic");

  long bnd = 64;
  amap.arg("bnd", bnd, "recursion bound for replication");

  long B = 0;
  amap.arg("B", B, "bound for # of replications", "all");

  amap.parse(argc, argv);
  setDryRun(dry);

  TestIt(m, p, r, d, L, bnd, B);
  cout << endl;
}
