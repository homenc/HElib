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
#include <cstdio>

#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include "replicate.h"
#include "permutations.h"

#define pSize (NTL_SP_NBITS/2) /* The size of levels in the chain */

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


void timeInit(long m, long p, long r, long d, long L, long nTests)
{
  for (long i=0; i<nTests; i++) {
    cerr << "." << std::flush;
    startFHEtimer((r>1)? "init4" : "init2");
    FHEcontext context(m, p, r);

    ZZX G;
    if (d==1) SetX(G); // set G(X)=X
    else G = context.alMod.getFactorsOverZZ()[0];

    buildModChain(context, L, /*c=*/3);
    EncryptedArray ea(context, G);
    stopFHEtimer((r>1)? "init4" : "init2");

    startFHEtimer("keyGen");
    FHESecKey secretKey(context);
    const FHEPubKey& publicKey = secretKey;
    secretKey.GenSecKey(64); // A Hamming-weight-64 secret key
    addSome1DMatrices(secretKey); // compute key-switching matrices
    stopFHEtimer("keyGen");

    ZZX poly;
    PlaintextArray pp(ea);
    pp.random();

    Ctxt cc(publicKey);
    if (r==1 && d==1)      startFHEtimer("encode2");
    else if (r==1 && d!=1) startFHEtimer("encode2d");
    else if (r!=1 && d==1) startFHEtimer("encode4");
    else                   startFHEtimer("encode4d");

    ea.encode(poly, pp);

    if (r==1 && d==1)      stopFHEtimer("encode2");
    else if (r==1 && d!=1) stopFHEtimer("encode2d");
    else if (r!=1 && d==1) stopFHEtimer("encode4");
    else                   stopFHEtimer("encode4d");

    startFHEtimer("encrypt");
    publicKey.Encrypt(cc, poly);
    stopFHEtimer("encrypt");

    startFHEtimer("decrypt");
    secretKey.Decrypt(poly, cc);
    stopFHEtimer("decrypt");
    startFHEtimer((r>1)? "decode4" : "decode2");
    ea.decode(pp, poly);
    stopFHEtimer((r>1)? "decode4" : "decode2");
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
  startFHEtimer("innerProduct");
  innerProduct(ret,cc,cc);
  ret.modDownToLevel(ret.findBaseLevel());
  stopFHEtimer("innerProduct");

  // Multiplication with 2,3 arguments
  cerr << "." << std::flush;
  for (long i=0; i<nTests; i++) {
    Ctxt c0 = cc[0];
    startFHEtimer("multiplyBy");
    c0.multiplyBy(cc[1]);
    c0.modDownToLevel(c0.findBaseLevel());
    stopFHEtimer("multiplyBy");
    ret += c0; // Just so the compiler doesn't optimize it away
  }

  if (level > 2) {
    cerr << "." << std::flush;
    for (long i=0; i<nTests; i++) {
      Ctxt c0 = cc[0];
      startFHEtimer("multiplyBy2");
      c0.multiplyBy2(cc[1],cc[2]);
      c0.modDownToLevel(c0.findBaseLevel()); // mod-down if needed
      stopFHEtimer("multiplyBy2");
      ret += c0; // Just so the compiler doesn't optimize it away
    }
  }

  // Multiply by constant
  cerr << "." << std::flush;
  for (long i=0; i<4*nTests; i++) {
    Ctxt c0 = cc[0];
    startFHEtimer("multByConstant");
    c0.multByConstant(p);
    stopFHEtimer("multByConstant");
    ret -= c0; // Just so the compiler doesn't optimize it away
  }

  // Add constant
  cerr << "." << std::flush;
  for (long i=0; i<10*nTests; i++) {
    Ctxt c0 = cc[0];
    startFHEtimer("addConstant");
    c0.addConstant(p);
    stopFHEtimer("addConstant");
    ret += c0; // Just so the compiler doesn't optimize it away
  }

  // Addition
  cerr << "." << std::flush;
  for (long i=0; i<10*nTests; i++) {
    startFHEtimer("add");
    ret += cc[0];
    stopFHEtimer("add");
  }

  // Rotation by an amount k for which we have a key-switching matrix
  cerr << "." << std::flush;
  for (long i=0; i<nTests; i++) {
    Ctxt c0 = cc[0];
    long k = rotationAmount(ea,publicKey,/*withMatrix=*/true);
    startFHEtimer("nativeAutomorph");
    c0.smartAutomorph(k);
    c0.modDownToLevel(c0.findBaseLevel());
    stopFHEtimer("nativeAutomorph");    
    ret += c0; // Just so the compiler doesn't optimize it away
  }
  // Rotation by a random amount k
  cerr << "." << std::flush;
  for (long i=0; i<nTests; i++) {
    Ctxt c0 = cc[0];
    long k = rotationAmount(ea,publicKey,/*withMatrix=*/false);
    startFHEtimer("automorph");
    c0.smartAutomorph(k);
    c0.modDownToLevel(c0.findBaseLevel()); // mod-down if needed
    stopFHEtimer("automorph");
    ret += c0; // Just so the compiler doesn't optimize it away
  }

  // record the results
  td.addConst = getTime4func("addConstant") / getNumCalls4func("addConstant");
  td.add = getTime4func("add") / getNumCalls4func("add");
  td.multConst = getTime4func("multByConstant") / getNumCalls4func("multByConstant");
  td.mult = getTime4func("multiplyBy") / getNumCalls4func("multiplyBy");
  if (getNumCalls4func("multiplyBy2") > 0)
    td.multBy2 = getTime4func("multiplyBy2") / getNumCalls4func("multiplyBy2");
  else td.multBy2 = 0;
  td.autoNative = getTime4func("nativeAutomorph") / getNumCalls4func("nativeAutomorph");
  td.autoTypical = getTime4func("automorph") / getNumCalls4func("automorph");
  td.innerProd = getTime4func("innerProduct") / getNumCalls4func("innerProduct");
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

  virtual void get(RX& out, long i, long j) const {
    assert(i >= 0 && i < ea.size());
    assert(j >= 0 && j < ea.size());
    out = data[i][j];
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
  PlaintextMatrixBaseInterface *ptr = buildRandomMatrix(ea);
  Ctxt tmp = c[0];
  tmp.modDownToLevel(td.lvl);
  cerr << "." << std::flush;
  startFHEtimer("MatMul");
  ea.mat_mul(tmp, *ptr);      // multiply the ciphertext vector
  stopFHEtimer("MatMul");
  ret = tmp;
  delete ptr;

  for (long i=0; i<nTests; i++) {
    cerr << "." << std::flush;
    long nSlots = ea.size();
    long r = RandomBnd(nSlots);
    tmp = c[i % c.size()];
    tmp.modDownToLevel(td.lvl);
    // time rotation
    startFHEtimer("rotate");
    ea.rotate(tmp, r);
    stopFHEtimer("rotate");
    // time shift, amount between -nSlots/2 and +nSlots/2
    if (r>nSlots/2) r -= nSlots;
    startFHEtimer("shift");
    ea.shift(tmp, r);
    stopFHEtimer("shift");
    ret += tmp; // just so the compiler will not optimize it out
  }

  cerr << "." << std::flush;
  for (long i=0; i<nTests && i<ea.size(); i++) {
    tmp = c[i % c.size()];
    tmp.modDownToLevel(td.lvl);
    startFHEtimer("replicate");
    replicate(ea, tmp, i);
    stopFHEtimer("replicate");    
    ret += tmp; // just so the compiler will not optimize it out
  }

  cerr << "." << std::flush;
  ReplicateDummy handler;
  tmp = c[1];
  tmp.modDownToLevel(td.lvl);
  startFHEtimer("replicateAll");
  replicateAll(ea, tmp, &handler);
  stopFHEtimer("replicateAll");
  ret += tmp;

  cerr << "." << std::flush;
  Permut pi;
  randomPerm(pi, trees.getSize());
  tmp = c[2];
  tmp.modDownToLevel(td.lvl);

  PermNetwork net;
  startFHEtimer("permutation");
  net.buildNetwork(pi, trees);  // Build a permutation network for pi
  net.applyToCtxt(tmp);         // Apply permutation netwrok
  stopFHEtimer("permutation");
  ret +=tmp; // just so the compiler will not optimize it out

  // record the timing data
  td.rotate = getTime4func("rotate") / getNumCalls4func("rotate");
  td.shift = getTime4func("shift") / getNumCalls4func("shift");
  td.permute = getTime4func("permutation") / getNumCalls4func("permutation");
  td.matmul = getTime4func("MatMul") / getNumCalls4func("MatMul");
  td.replicate = getTime4func("replicate") / getNumCalls4func("replicate");
  td.replAll = getTime4func("replicateAll") / getNumCalls4func("replicateAll");
  resetAllTimers();
}


void  TimeIt(long m, long p, TimingData& data, bool high=false)
{
  resetAllTimers();
  long phim = phi_N(m);
  long L = (floor((7.2*phim)/(pSize* /*cc*/1.33* (110+/*k*/80))) -1)/2;
  if (L<2) L=2; // Make sure we have at least a few primes

  // Initialize a context with r=2,d=1
  startFHEtimer("init4");
  FHEcontext context(m, 2, 2);
  buildModChain(context, L, /*c=*/3);

  ZZX G; SetX(G); // G(X) = X
  EncryptedArray ea(context, G);
  startFHEtimer("init4");

  startFHEtimer("keyGen");
  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(64); // A Hamming-weight-64 secret key
  addSome1DMatrices(secretKey); // compute key-switching matrices
  stopFHEtimer("keyGen");

  // record timing results for init/geygen, etc.
  data.m = m;
  data.phim = phim;
  data.nSlots = ea.size();
  data.p = p;

  timeInit(m, p, /*r=*/1, /*d=*/1, L, /*nTests=*/4);
  timeInit(m, p, /*r=*/1, /*d=*/0, L, /*nTests=*/4);
  timeInit(m, p, /*r=*/2, /*d=*/1, L, /*nTests=*/4);
  timeInit(m, p, /*r=*/2, /*d=*/0, L, /*nTests=*/4);

  data.other.init2 = getTime4func("init2") / getNumCalls4func("init2");
  data.other.init4 = getTime4func("init4") / getNumCalls4func("init4");
  data.other.keyGen = getTime4func("keyGen") / getNumCalls4func("keyGen");
  data.other.encode2 = getTime4func("encode2") / getNumCalls4func("encode2");
  data.other.encode2d = getTime4func("encode2d") / getNumCalls4func("encode2d");
  data.other.encode4 = getTime4func("encode4") / getNumCalls4func("encode4");
  data.other.encode4d = getTime4func("encode4d") / getNumCalls4func("encode4d");
  data.other.encrypt = getTime4func("encrypt") / getNumCalls4func("encrypt");
  data.other.decrypt = getTime4func("decrypt") / getNumCalls4func("decrypt");
  data.other.decode2 = getTime4func("decode2") / getNumCalls4func("decode2");
  data.other.decode4 = getTime4func("decode4") / getNumCalls4func("decode4");
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
