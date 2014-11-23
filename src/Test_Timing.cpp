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

/* Test_Timing.cpp - A program that tests the timing of various operations, outputs the results in a comma-separate-value (csv) format.
 */
#include <cassert>
#include <cstdio>
#include <memory>
#include <NTL/ZZ.h>
NTL_CLIENT

#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include "replicate.h"
#include "permutations.h"

// We measure low-level timing at all levels
class LowLvlTimingData {
public:
  long lvl;
  double addConst;
  double add;
  double multConst;
  double mult;
  double multBy2;
  double autoNative;
  double autoTypical;
  double innerProd;

  explicit LowLvlTimingData(long _lvl=-1) {lvl=_lvl;}
};

// We only measure high-level timing at one level (lvl=8 or lower)
class HighLvlTimingData {
public:
  long lvl;
  double rotate;
  double shift;
  double permute;
  double matmul;
  double replicate;
  double replAll;

  explicit HighLvlTimingData(long _lvl=-1) {lvl=_lvl;}
};

class OtherTimingData {
public:
  double init2;
  double init4;
  double keyGen;
  double encode2;
  double encode2d;
  double encode4;
  double encode4d;
  double encrypt;
  double decode2;
  double decode4;
  double decrypt;
};

class TimingData {
public:
  long m;
  long phim;
  long nSlots;
  long p;
  vector<LowLvlTimingData> lowLvl;
  HighLvlTimingData highLvl;
  OtherTimingData other;

  explicit TimingData(long _m=-1) {m=_m;}
};

static FHEtimer* _init_timer_2 = 0;
static FHEtimer* _init_timer_4 = 0;

void timeInit(long m, long p, long r, long d, long L, long nTests)
{
  for (long i=0; i<nTests; i++) {
    cerr << "." << std::flush;

    // Complicated mumbo-jumbo to get one of two timers, depending on r
    auto_timer _init_timer((r>1)? _init_timer_4 : _init_timer_2);

    FHEcontext context(m, p, r);

    ZZX G;
    if (d==1) SetX(G); // set G(X)=X
    else G = context.alMod.getFactorsOverZZ()[0];

    buildModChain(context, L, /*c=*/3);
    EncryptedArray ea(context, G);
    _init_timer.stop();

    FHE_NTIMER_START(keyGen);
    FHESecKey secretKey(context);
    const FHEPubKey& publicKey = secretKey;
    secretKey.GenSecKey(64); // A Hamming-weight-64 secret key
    addSome1DMatrices(secretKey); // compute key-switching matrices
    FHE_NTIMER_STOP(keyGen);

    ZZX poly;
    PlaintextArray pp(ea);
    pp.random();

    Ctxt cc(publicKey);

    if (r==1 && d==1) {
      FHE_NTIMER_START(encode2);
      ea.encode(poly, pp);
      FHE_NTIMER_STOP(encode2);
    }
    else if (r==1 && d!=1) {
      FHE_NTIMER_START(encode2d);
      ea.encode(poly, pp);
      FHE_NTIMER_STOP(encode2d);
    }
    else if (r!=1 && d==1) {
      FHE_NTIMER_START(encode4);
      ea.encode(poly, pp);
      FHE_NTIMER_STOP(encode4);
    }
    else {
      FHE_NTIMER_START(encode4d);
      ea.encode(poly, pp);
      FHE_NTIMER_STOP(encode4d);
    }

    FHE_NTIMER_START(encrypt);
    publicKey.Encrypt(cc, poly);
    FHE_NTIMER_STOP(encrypt);

    FHE_NTIMER_START(decrypt);
    secretKey.Decrypt(poly, cc);
    FHE_NTIMER_STOP(decrypt);

    if (r>1) {
      FHE_NTIMER_START(decode4);
      ea.decode(pp, poly);
      FHE_NTIMER_STOP(decode4);
    } else {
      FHE_NTIMER_START(decode2);
      ea.decode(pp, poly);
      FHE_NTIMER_STOP(decode2);
    }
  }
}


// Returns either a random automorphism amount or an amount
// for which we have a key-switching matrix s^k -> s.
long rotationAmount(const EncryptedArray& ea, const FHEPubKey& publicKey,
	       bool onlyWithMatrix)
{
  const PAlgebra& pa = ea.getContext().zMStar;
  long nSlots = pa.getNSlots();
  long r = RandomBnd(nSlots);
  long k = pa.ith_rep(r);
  if (onlyWithMatrix) { // return the 1st step in the path to k
    const KeySwitch& matrix = publicKey.getNextKSWmatrix(k,0);
    k = matrix.fromKey.getPowerOfX();
  }
  return k;
}

void timeOps(const EncryptedArray& ea, const FHEPubKey& publicKey, Ctxt& ret,
	     const vector<Ctxt>& c, ZZX& p, long nTests, LowLvlTimingData& td)
{
  assert(c.size()>=3);
  vector<Ctxt> cc = c;
  // perform operations at a lower level
  long level = td.lvl;
  if (level>0) for (long i=0; i<(long)cc.size(); i++)
    cc[i].modDownToLevel(level);
  else {
    level = cc[0].findBaseLevel();
    td.lvl = level;
  }

  // inner-product of vectors
  cerr << "." << std::flush;
  FHE_NTIMER_START(innerProduct);
  innerProduct(ret,cc,cc);
  ret.modDownToLevel(ret.findBaseLevel());
  FHE_NTIMER_STOP(innerProduct);

  // Multiplication with 2,3 arguments
  cerr << "." << std::flush;
  for (long i=0; i<nTests; i++) {
    Ctxt c0 = cc[0];
    FHE_NTIMER_START(multiplyBy);
    c0.multiplyBy(cc[1]);
    c0.modDownToLevel(c0.findBaseLevel());
    FHE_NTIMER_STOP(multiplyBy);
    ret += c0; // Just so the compiler doesn't optimize it away
  }

  if (level > 2) {
    cerr << "." << std::flush;
    for (long i=0; i<nTests; i++) {
      Ctxt c0 = cc[0];
      FHE_NTIMER_START(multiplyBy2);
      c0.multiplyBy2(cc[1],cc[2]);
      c0.modDownToLevel(c0.findBaseLevel()); // mod-down if needed
      FHE_NTIMER_STOP(multiplyBy2);
      ret += c0; // Just so the compiler doesn't optimize it away
    }
  }

  // Multiply by constant
  cerr << "." << std::flush;
  for (long i=0; i<4*nTests; i++) {
    Ctxt c0 = cc[0];
    FHE_NTIMER_START(multByConstant);
    c0.multByConstant(p);
    FHE_NTIMER_STOP(multByConstant);
    ret -= c0; // Just so the compiler doesn't optimize it away
  }

  // Add constant
  cerr << "." << std::flush;
  for (long i=0; i<10*nTests; i++) {
    Ctxt c0 = cc[0];
    FHE_NTIMER_START(addConstant);
    c0.addConstant(p);
    FHE_NTIMER_STOP(addConstant);
    ret += c0; // Just so the compiler doesn't optimize it away
  }

  // Addition
  cerr << "." << std::flush;
  for (long i=0; i<10*nTests; i++) {
    FHE_NTIMER_START(add);
    ret += cc[0];
    FHE_NTIMER_STOP(add);
  }

  // Rotation by an amount k for which we have a key-switching matrix
  cerr << "." << std::flush;
  for (long i=0; i<nTests; i++) {
    Ctxt c0 = cc[0];
    long k = rotationAmount(ea,publicKey,/*withMatrix=*/true);
    FHE_NTIMER_START(nativeAutomorph);
    c0.smartAutomorph(k);
    c0.modDownToLevel(c0.findBaseLevel());
    FHE_NTIMER_STOP(nativeAutomorph);    
    ret += c0; // Just so the compiler doesn't optimize it away
  }
  // Rotation by a random amount k
  cerr << "." << std::flush;
  for (long i=0; i<nTests; i++) {
    Ctxt c0 = cc[0];
    long k = rotationAmount(ea,publicKey,/*withMatrix=*/false);
    FHE_NTIMER_START(automorph);
    c0.smartAutomorph(k);
    c0.modDownToLevel(c0.findBaseLevel()); // mod-down if needed
    FHE_NTIMER_STOP(automorph);
    ret += c0; // Just so the compiler doesn't optimize it away
  }

  // record the results
  FHEtimer tmp(NULL,NULL);

  getTimerByName(tmp, "addConstant");
  td.addConst = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "add");
  td.add = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "multByConstant");
  td.multConst = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "multiplyBy");
  td.mult = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "multiplyBy2");
  if (tmp.getNumCalls() > 0)
    td.multBy2 = tmp.getTime() / tmp.getNumCalls();
  else td.multBy2 = 0;

  getTimerByName(tmp, "nativeAutomorph");
  td.autoNative = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "automorph");
  td.autoTypical = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "innerProduct");
  td.innerProd = tmp.getTime() / tmp.getNumCalls();
  resetAllTimers();
}


template<class type> 
class RandomMatrix : public  PlaintextMatrixInterface<type> {
public:
  PA_INJECT(type) 

private:
  const EncryptedArray& ea;

  vector< vector< RX > > data;

public:
  //  ~RandomMatrix() { cerr << "destructor: random matrix\n"; }

  RandomMatrix(const EncryptedArray& _ea) : ea(_ea) { 
    long n = ea.size();
    long d = ea.getDegree();

    RBak bak; bak.save(); ea.getContext().alMod.restoreContext();

    data.resize(n);
    for (long i = 0; i < n; i++) {
      data[i].resize(n);
      for (long j = 0; j < n; j++)
        random(data[i][j], d);
    }
  }

  virtual const EncryptedArray& getEA() const {
    return ea;
  }

  virtual bool get(RX& out, long i, long j) const {
    assert(i >= 0 && i < ea.size());
    assert(j >= 0 && j < ea.size());
    out = data[i][j];
    return false;
  }
};

PlaintextMatrixBaseInterface *buildRandomMatrix(const EncryptedArray& ea)
{
  switch (ea.getContext().alMod.getTag()) {
    case PA_GF2_tag: {
      return new RandomMatrix<PA_GF2>(ea);
    }

    case PA_zz_p_tag: {
      return new RandomMatrix<PA_zz_p>(ea);
    }

    default: return 0;
  }
}

class ReplicateDummy : public ReplicateHandler {
public:
  ReplicateDummy() {}
  virtual void handle(const Ctxt& ctxt) {}
};

void timeHighLvl(const EncryptedArray& ea, const FHEPubKey& publicKey,
		 Ctxt& ret, const vector<Ctxt>& c, GeneratorTrees& trees,
		 long nTests, HighLvlTimingData& td)
{
  Ctxt tmp = c[0];
  tmp.modDownToLevel(td.lvl);
  cerr << "." << std::flush;
  {
  std::unique_ptr<PlaintextMatrixBaseInterface> ptr(buildRandomMatrix(ea));
  FHE_NTIMER_START(MatMul);
  ea.mat_mul(tmp, *ptr);      // multiply the ciphertext vector
  FHE_NTIMER_STOP(MatMul);
  } // free the pointer
  ret = tmp;

  for (long i=0; i<nTests; i++) {
    cerr << "." << std::flush;
    long nSlots = ea.size();
    long r = RandomBnd(nSlots);
    tmp = c[i % c.size()];
    tmp.modDownToLevel(td.lvl);
    // time rotation
    FHE_NTIMER_START(rotate);
    ea.rotate(tmp, r);
    FHE_NTIMER_STOP(rotate);
    // time shift, amount between -nSlots/2 and +nSlots/2
    if (r>nSlots/2) r -= nSlots;
    FHE_NTIMER_START(shift);
    ea.shift(tmp, r);
    FHE_NTIMER_STOP(shift);
    ret += tmp; // just so the compiler will not optimize it out
  }

  cerr << "." << std::flush;
  for (long i=0; i<nTests && i<ea.size(); i++) {
    tmp = c[i % c.size()];
    tmp.modDownToLevel(td.lvl);
    FHE_NTIMER_START(replicate);
    replicate(ea, tmp, i);
    FHE_NTIMER_STOP(replicate);    
    ret += tmp; // just so the compiler will not optimize it out
  }

  cerr << "." << std::flush;
  ReplicateDummy handler;
  tmp = c[1];
  tmp.modDownToLevel(td.lvl);
  FHE_NTIMER_START(replicateAll);
  replicateAll(ea, tmp, &handler);
  FHE_NTIMER_STOP(replicateAll);
  ret += tmp;

  cerr << "." << std::flush;
  Permut pi;
  randomPerm(pi, trees.getSize());
  tmp = c[2];
  tmp.modDownToLevel(td.lvl);

  PermNetwork net;
  FHE_NTIMER_START(permutation);
  net.buildNetwork(pi, trees);  // Build a permutation network for pi
  net.applyToCtxt(tmp, ea);         // Apply permutation netwrok
  FHE_NTIMER_STOP(permutation);
  ret +=tmp; // just so the compiler will not optimize it out

  // record the timing data
  FHEtimer ttmp(NULL,NULL);

  getTimerByName(ttmp, "rotate");
  td.rotate = ttmp.getTime() / ttmp.getNumCalls();

  getTimerByName(ttmp, "shift");
  td.shift = ttmp.getTime() / ttmp.getNumCalls();

  getTimerByName(ttmp, "permutation");
  td.permute = ttmp.getTime() / ttmp.getNumCalls();

  getTimerByName(ttmp, "MatMul");
  td.matmul = ttmp.getTime() / ttmp.getNumCalls();

  getTimerByName(ttmp, "replicate");
  td.replicate = ttmp.getTime() / ttmp.getNumCalls();

  getTimerByName(ttmp, "replicateAll");
  td.replAll = ttmp.getTime() / ttmp.getNumCalls();
  resetAllTimers();
}


void  TimeIt(long m, long p, TimingData& data, bool high=false)
{
  setTimersOn();
  resetAllTimers();
  long phim = phi_N(m);
  long L = floor((7.2*phim)/(FHE_pSize* /*cc*/1.33* (110+/*k*/80)));
  if (L<2) L=2; // Make sure we have at least a few primes

  buildTimer(_init_timer_2, "init2", FHE_AT);
  buildTimer(_init_timer_4, "init4", FHE_AT);

  // Initialize a context with r=2,d=1
  auto_timer _init_timer(_init_timer_4);
  FHEcontext context(m, 2, 2);
  buildModChain(context, L, /*c=*/3);

  ZZX G; SetX(G); // G(X) = X
  EncryptedArray ea(context, G);
  _init_timer.stop();

  FHE_NTIMER_START(keyGen);
  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(64); // A Hamming-weight-64 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices
  FHE_NTIMER_STOP(keyGen);

  // record timing results for init/geygen, etc.
  data.m = m;
  data.phim = phim;
  data.nSlots = ea.size();
  data.p = p;

  timeInit(m, p, /*r=*/1, /*d=*/1, L, /*nTests=*/4);
  timeInit(m, p, /*r=*/1, /*d=*/0, L, /*nTests=*/4);
  timeInit(m, p, /*r=*/2, /*d=*/1, L, /*nTests=*/4);
  timeInit(m, p, /*r=*/2, /*d=*/0, L, /*nTests=*/4);

  //  printAllTimers();
  FHEtimer tmp(NULL,NULL);

  getTimerByName(tmp, "init2");
  data.other.init2 = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "init4");
  data.other.init4 = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "keyGen");
  data.other.keyGen = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "encode2");
  data.other.encode2 = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "encode2d");
  data.other.encode2d = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "encode4");
  data.other.encode4 = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "encode4d");
  data.other.encode4d = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "encrypt");
  data.other.encrypt = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "decrypt");
  data.other.decrypt = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "decode2");
  data.other.decode2 = tmp.getTime() / tmp.getNumCalls();

  getTimerByName(tmp, "decode4");
  data.other.decode4 = tmp.getTime() / tmp.getNumCalls();
  resetAllTimers();

  // time low-level operations
  cerr << "#" << std::flush;

  ZZX poly;
  PlaintextArray pp(ea);
  pp.random();
  ea.encode(poly, pp);

  Ctxt cc(publicKey);
  long nTests = 10;
  vector<Ctxt> vc(nTests,cc);
  for (long i=0; i<nTests; i++) {
    pp.random();
    ea.encrypt(vc[i], publicKey, pp);
  }

  LowLvlTimingData td;
  data.lowLvl.clear();
  for (long i=2; i<2*L-1; i*=2) {
    td.lvl=i;
    timeOps(ea, publicKey, cc,vc, poly, nTests, td);
    data.lowLvl.push_back(td);
  }
  td.lvl = 2*L -1;
  timeOps(ea, publicKey, cc,vc, poly, nTests, td);
  data.lowLvl.push_back(td);

  if (high && L > 4) { /// cannot test high-level routines at level <8
    cerr << "#" << std::flush;
    // Setup generator-descriptors for the PAlgebra generators
    Vec<GenDescriptor> vec(INIT_SIZE, ea.dimension());
    for (long i=0; i<ea.dimension(); i++)
      vec[i] = GenDescriptor(/*order=*/ea.sizeOfDimension(i),
			     /*good=*/ ea.nativeDimension(i), /*genIdx=*/i);

    // Some default for the width-bound, if not provided
    //  long widthBound = log2((double)ea.size()) -1;

    // Get the generator-tree structures and the corresponding hypercube
    GeneratorTrees trees;
    trees.buildOptimalTrees(vec, /*widthBound=*/7);
    //  cout << " cost =" << cost << endl;

    // build network for a random permutation, for the sole purpose
    // of adding key-switching matrices
    Permut pi;
    randomPerm(pi, trees.getSize());
    PermNetwork net;
    net.buildNetwork(pi, trees);
    addMatrices4Network(secretKey, net);

    data.highLvl.lvl = 8;
    timeHighLvl(ea, publicKey, cc, vc, trees, nTests, data.highLvl);
  }
  cerr << "!" << std::flush;
}

void printTimeData(TimingData& td)
{
  if (td.m <=0) return;

  OtherTimingData& otd = td.other;
  HighLvlTimingData& hl = td.highLvl;
  for (long i=0; i<(long)td.lowLvl.size(); i++) {
    LowLvlTimingData& ll = td.lowLvl[i];
    long lvl = ll.lvl;
    cout <<td.p<<","<<td.m<<","<<td.phim<<","<<td.nSlots<<","
	 <<otd.init2<<","<<otd.init4<<","<<otd.keyGen<<","
	 <<otd.encode2<<","<<otd.encode2d<<","
	 <<otd.encode4<<","<<otd.encode4d<<","
	 <<otd.encrypt<<","<<otd.decrypt<<","
	 <<otd.decode2<<","<<otd.decode4<<","
	 <<lvl<<","<<ll.addConst<<","<<ll.add<<","<<ll.multConst<<","
	 <<ll.mult<<",";
    if (ll.multBy2>0.0) cout << ll.multBy2;
    cout <<","<<ll.autoNative<<","
	 <<ll.autoTypical<<","<<ll.innerProd;
    if (lvl==hl.lvl)
      cout <<","<<hl.rotate<<","<<hl.shift<<","<<hl.permute<<","
	   <<hl.matmul<<","<<hl.replicate<<","<<hl.replAll;
    cout << endl;
  }
}

void usage(char *prog) 
{
  cerr << "A program that tests the timing of various operations,\n  outputs the results in a comma-separate-value (csv) format.\n"
  cerr << "Usage: "<<prog<<" [ optional parameters ]... 2> logfile > results-file\n";
  cerr << "results on stdout in comma-separated-value format, ";
  cerr << "progress printed on stderr\n";
  cerr << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cerr << "  e.g, 'm=11441 p=2 high=1'\n\n";
  cerr << "  m determines the cyclotomic ring, defaults to all the set\n";
  cerr << "    m in { 4051, 4369, 4859, 10261,11023,11441,\n";
  cerr << "          18631,20485,21845, 49981,53261       }\n";
  cerr << "  p is the plaintext base [default=2]" << endl;
  cerr << "  high=1 will time also high-level procedures [default==0]\n";
  exit(0);
}

int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["p"] = "2";
  argmap["m"] = "0";
  argmap["high"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long p = atoi(argmap["p"]);
  long m = atoi(argmap["m"]);
  long high = atoi(argmap["high"]);

#define numTests 11
  long ms[numTests] = { 4051, 4369, 4859, 10261,11023,11441,
		 18631,20485,21845, 49981,53261};

  cout << "p,m,phim,nSlots,init(p),init(p^2),keyGen,encode(F_p),encode(F_p^d),encode(Z_p^r),encode(R_(p^r)^d),encrypt,decrypt,decode(p),decode(p^2),level,addConst,add,multConst,mult,mult2,autoNative,autoTypical,inProd10";
  if (high) cout << ",rotate,shift,permute,matmul,replicate,replAll\n";
  else cout << endl;

  TimingData td;
  if (m>0) {
    cerr << "\nTesting m="<<m; 
    if (high) cout << " (including high-level)";
    TimeIt(m, p, td, high);
    printTimeData(td);
  }
  else for (long i=0; i<numTests; i++) {
      cerr << "\nTesting m="<<ms[i];
      if (high) cout << " (including high-level)";
      TimeIt(ms[i], p, td, /*timeHighLvl=*/high);
      printTimeData(td);
    }
}
