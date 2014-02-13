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

#include "FHE.h"
#include "replicate.h"
#include "timing.h"

#include <NTL/lzz_pXFactoring.h>

#include <cassert>

NTL_CLIENT

static bool decryptAndTest=true;

// k=10 p=5 r=2 -- 3 dimensions: 12, 3, !2
// p=11 -- 2 dimensions: 16, !4
// Test_General_x p=23 s=100 z=20  ## a large, realistic example, 48 x 4
// Test_General_x p=2 r=4 k=10 z=10 m=4369 ## 16 x !16
// Test_General_x p=2 r=4 k=10 z=10 m=1247 ## 42
// Test_General_x p=2 r=4 k=10 z=10 m=8191 ## 630
// Test_General_x p=2 r=4 k=10 z=10 m=3133 ## 60 x !2


class StopReplicate { };

class ReplicateTester : public ReplicateHandler {
public:
  const FHESecKey& sKey;
  const EncryptedArray& ea;
  const PlaintextArray& pa;
  long M;

  double t_last, t_total;
  long pos;

  ReplicateTester(const FHESecKey& _sKey, const EncryptedArray& _ea, 
                  const PlaintextArray& _pa, long _M)
  : sKey(_sKey), ea(_ea), pa(_pa), M(_M)
  {
    t_last = GetTime();
    t_total = 0.0;
    pos = 0;
  }

  virtual void handle(const Ctxt& ctxt) {

    if (decryptAndTest) {
      double t_new = GetTime();
      double t_elapsed = t_new - t_last;

      t_total += t_elapsed;

      // cerr << "*** " << pos << " t=" << t_elapsed 
      // 	   << ", t_total=" << t_total 
      // 	   << ", level=" << ctxt.findBaseLevel() 
      // 	   << ", log(noise/modulus)~" << ctxt.log_of_ratio() 
      // 	   << "\n";
    
      PlaintextArray pa1 = pa;
      pa1.replicate(pos);
      PlaintextArray pa2(ea);

      ea.decrypt(ctxt, sKey, pa2);
      if (!pa1.equals(pa2)) {
	cerr << "error:\n";
	pa2.print(cerr); cerr << "\n";
      }

      t_last = GetTime();
    }
    pos++;

    if (M > 0 && pos >= M) throw StopReplicate();
  }
};



void  TestIt(long p, long r, long d, long c, long k, long w, 
               long L, long m, long bnd, long M, long v)
{
  cerr << "*** TestIt: "
       << "  p=" << p
       << ", r=" << r
       << ", d=" << d
       << ", c=" << c
       << ", k=" << k
       << ", w=" << w
       << ", L=" << L
       << ", m=" << m
       << ", bnd=" << bnd
       << ", M=" << bnd
       << ", v=" << v
       << ", t=" << decryptAndTest
       << endl;

  FHEcontext context(m, p, r);
  buildModChain(context, L, c);

  context.zMStar.printout();
  cerr << endl;

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(w); // A Hamming-weight-w secret key


  ZZX G;

  if (d == 0)
    G = context.alMod.getFactorsOverZZ()[0];
  else
    G = makeIrredPoly(p, d); 



  cerr << "G = " << G << "\n";
  cerr << "generating key-switching matrices... ";
  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  cerr << "done\n";

  printAllTimers();
  resetAllTimers();

  cerr << "computing masks and tables for rotation...";
  EncryptedArray ea(context, G);
  cerr << "done\n";



  PlaintextArray xp0(ea), xp1(ea);


  xp0.random();
  xp1.random();


  Ctxt xc0(publicKey);
  ea.encrypt(xc0, publicKey, xp0);

  ZZX poly_xp1;
  ea.encode(poly_xp1, xp1);

#if 0
  double t;
  cerr << "multiplication test:\n";
  t = GetTime();
  for (long i = 0; i < ea.size(); i++) {
    Ctxt xc2 = xc0;
    xc2.multByConstant(poly_xp1);
  }
  t = GetTime()-t;
  cerr << "time = " << t << "\n";
#endif  

  cerr << "** Testing replicate():\n";
  Ctxt xc1 = xc0;
  CheckCtxt(xc1, "before replicate");
  replicate(ea, xc1, ea.size()/2);
  CheckCtxt(xc1, "after replicate");

  // Get some timing results
  for (long i=0; i<20 && i<ea.size(); i++) {
    xc1 = xc0;
    startFHEtimer("replicate");
    replicate(ea, xc1, i);
    stopFHEtimer("replicate");    
  }

  cerr << "** Testing replicateAll():\n";
  replicateVerboseFlag = v;
  ReplicateHandler *handler = new ReplicateTester(secretKey, ea, xp0, M);
  try {
    startFHEtimer("replicateAll");
    replicateAll(ea, xc0, handler, bnd);
  }
  catch (StopReplicate) {
  }
  stopFHEtimer("replicateAll");

  delete handler;

}


void usage(char *prog) 
{
  cerr << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cerr << "  p is the plaintext base [default=2]" << endl;
  cerr << "  r is the lifting [default=1]" << endl;
  cerr << "  d is the degree of the field extension [default==1]\n";
  cerr << "    (d == 0 => factors[0] defined the extension)\n";
  cerr << "  c is number of columns in the key-switching matrices [default=2]\n";
  cerr << "  k is the security parameter [default=80]\n";
  cerr << "  z is the # of primes [default=4]\n";
  cerr << "  s is the minimum number of slots [default=4]\n";
  cerr << "  m is a specific modulus\n";
  cerr << "  bnd is a recursion bound for replication\n";
  cerr << "  M is a bound for # of replications [default=0 => all]\n";
  cerr << "  v for verbose [default=0]\n";
  cerr << "  t for decrypt-and-compare [default=1], otherwise just timing\n";
  exit(0);
}


int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["p"] = "2";
  argmap["r"] = "1";
  argmap["d"] = "1";
  argmap["c"] = "2";
  argmap["k"] = "80";
  argmap["z"] = "4";
  argmap["s"] = "0";
  argmap["m"] = "0";
  argmap["bnd"] = "64";
  argmap["M"] = "0";
  argmap["v"] = "0";
  argmap["t"] = "1";

  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long d = atoi(argmap["d"]);
  long c = atoi(argmap["c"]);
  long k = atoi(argmap["k"]);
  long z = atoi(argmap["z"]);
  long s = atoi(argmap["s"]);
  long chosen_m = atoi(argmap["m"]);
  long bnd = atoi(argmap["bnd"]);
  long M = atoi(argmap["M"]);
  long v = atoi(argmap["v"]);

  long w = 64; // Hamming weight of secret key
  long L = z; // number of levels
  long m = FindM(k, L, c, p, d, s, chosen_m, true);

  decryptAndTest = atoi(argmap["t"]);

  //  setTimersOn();

  TestIt(p, r, d, c, k, w, L, m, bnd, M, v);

  cerr << endl;
  printAllTimers();
  cerr << endl;

}

