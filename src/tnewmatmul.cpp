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

/**
 * @file matmul.h
 * @brief some matrix / linear algenra stuff
 */
#include "matmul.h"
#include "newmatmul.h"
#include <NTL/BasicThreadPool.h>




template<class type> class RandomMatrix_new : public  MatMul1D_derived<type> {
public:
  PA_INJECT(type) 

private:
  vector< vector< RX > > data;
  const EncryptedArray& ea;
  long dim;

public:
  virtual ~RandomMatrix_new() {}
  RandomMatrix_new(const EncryptedArray& _ea, long _dim): 
    ea(_ea), dim(_dim)
  {
    RBak bak; bak.save(); ea.getAlMod().restoreContext();
    long n = ea.size();
    long d = ea.getDegree();
    long D = ea.sizeOfDimension(dim);

    RandomStreamPush push;
    SetSeed(ZZ(123));

    data.resize(D);
    for (long i = 0; i < D; i++) {
      data[i].resize(D);
      for (long j = 0; j < D; j++) {
        random(data[i][j], d);
      }
    }
  }

  const EncryptedArray& getEA() const override { return ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }

  bool get(RX& out, long i, long j, long k) const override {
    long D = ea.sizeOfDimension(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }
};


static MatMul1D*
buildRandomMatrix_new(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomMatrix_new<PA_GF2>(ea, dim);
    }
    case PA_zz_p_tag: {
      return new RandomMatrix_new<PA_zz_p>(ea, dim);
    }
    default: return 0;
  }
}


template<class type> class RandomMatrix : public  MatMul<type> {
public:
  PA_INJECT(type) 

private:
  vector< vector< RX > > data;
  long dim;

public:
  virtual ~RandomMatrix() {}
  RandomMatrix(const EncryptedArray& _ea, long _dim, long g): 
    MatMul<type>(_ea,g), dim(_dim)
  {
    RBak bak; bak.save(); _ea.getAlMod().restoreContext();
    long n = _ea.size();
    long d = _ea.getDegree();
    long D = _ea.sizeOfDimension(dim);

    RandomStreamPush push;
    SetSeed(ZZ(123));

    data.resize(D);
    for (long i = 0; i < D; i++) {
      data[i].resize(D);
      for (long j = 0; j < D; j++) {
        random(data[i][j], d);
      }
    }
  }

  virtual bool get(RX& out, long i, long j) const {
    long D = this->getEA().sizeOfDimension(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }
};


static MatMulBase*
buildRandomMatrix(const EncryptedArray& ea, long dim, long giantStep)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomMatrix<PA_GF2>(ea, dim, giantStep);
    }
    case PA_zz_p_tag: {
      return new RandomMatrix<PA_zz_p>(ea, dim, giantStep);
    }
    default: return 0;
  }
}

template<class type> class RandomMultiMatrix_new : public  MatMul1D_derived<type> {
public:
  PA_INJECT(type) 

private:
  vector< vector< vector< RX > > > data;
  const EncryptedArray& ea;
  long dim;

public:
  virtual ~RandomMultiMatrix_new() {}
  RandomMultiMatrix_new(const EncryptedArray& _ea, long _dim): 
    ea(_ea), dim(_dim)
  {
    RBak bak; bak.save(); ea.getAlMod().restoreContext();
    long n = ea.size();
    long d = ea.getDegree();
    long D = ea.sizeOfDimension(dim);

    RandomStreamPush push;
    SetSeed(ZZ(123));

    data.resize(n/D);
    for (long k = 0; k < n/D; k++) {
      data[k].resize(D);
      for (long i = 0; i < D; i++) {
	data[k][i].resize(D);
	for (long j = 0; j < D; j++) {
	  random(data[k][i][j], d);
	}
      }
    }
  }

  const EncryptedArray& getEA() const override { return ea; }
  bool multipleTransforms() const override { return true; }
  long getDim() const override { return dim; }

  bool get(RX& out, long i, long j, long k) const override {
    long n = ea.size();
    long D = ea.sizeOfDimension(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    assert(k >= 0 && k < n/D);
    if (IsZero(data[k][i][j])) return true;
    out = data[k][i][j];
    return false;
  }
};


static MatMul1D*
buildRandomMultiMatrix_new(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomMultiMatrix_new<PA_GF2>(ea, dim);
    }
    case PA_zz_p_tag: {
      return new RandomMultiMatrix_new<PA_zz_p>(ea, dim);
    }
    default: return 0;
  }
}

template<class type> class RandomMultiMatrix : public  MatMul<type> {
public:
  PA_INJECT(type) 

private:
  vector< vector< vector< RX > > > data;
  long dim;

public:
  virtual ~RandomMultiMatrix() {}
  RandomMultiMatrix(const EncryptedArray& _ea, long _dim, long g)
    : MatMul<type>(_ea, g), dim(_dim)
  {
    RBak bak; bak.save(); _ea.getAlMod().restoreContext();
    long n = _ea.size();
    long d = _ea.getDegree();
    long D = _ea.sizeOfDimension(dim);

    RandomStreamPush push;
    SetSeed(ZZ(123));

    data.resize(n/D);
    for (long k = 0; k < n/D; k++) {
      data[k].resize(D);
      for (long i = 0; i < D; i++) {
	data[k][i].resize(D);
	for (long j = 0; j < D; j++) {
	  random(data[k][i][j], d);
	}
      }
    }
  }

  virtual bool multiGet(RX& out, long i, long j, long k) const
  {
    long D = this->getEA().sizeOfDimension(dim);
    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    assert(k >= 0 && k < this->getEA().size()/D);
    if (IsZero(data[k][i][j])) return true;
    out = data[k][i][j];
    return false;
  }
};

static MatMulBase*
buildRandomMultiMatrix(const EncryptedArray& ea, long dim, long giantStep)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomMultiMatrix<PA_GF2>(ea, dim, giantStep);
    }
    case PA_zz_p_tag: {
      return new RandomMultiMatrix<PA_zz_p>(ea, dim, giantStep);
    }
    default: return 0;
  }
}



//********************************

template<class type> 
class RandomBlockMatrix : public BlockMatMul<type> {
  PA_INJECT(type) 

  vector< vector< mat_R > > data;
  long dim;

public:
  ~RandomBlockMatrix() { /*cout << "destructor: random matrix\n";*/ }

  RandomBlockMatrix(const EncryptedArray& _ea, long _dim):
    BlockMatMul<type>(_ea), dim(_dim)
  {
    RBak bak; bak.save(); _ea.getAlMod().restoreContext();
    long n = _ea.size();
    long d = _ea.getDegree();
    long D = _ea.sizeOfDimension(dim);


    RandomStreamPush push;
    SetSeed(ZZ(123));

    data.resize(D);
    for (long i = 0; i < D; i++) {
      data[i].resize(D);
      for (long j = 0; j < D; j++) {
        data[i][j].SetDims(d, d);
        for (long u = 0; u < d; u++)
          for (long v = 0; v < d; v++) 
            random(data[i][j][u][v]);
      }
    }
  }

  virtual bool get(mat_R& out, long i, long j) const
  {
    const EncryptedArray& ea = this->getEA();
    long D = ea.sizeOfDimension(dim);
    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }
};

static MatMulBase*
buildRandomBlockMatrix(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomBlockMatrix<PA_GF2>(ea,dim);
    }
    case PA_zz_p_tag: {
      return new RandomBlockMatrix<PA_zz_p>(ea, dim);
    }
    default: return 0;
  }
}

template<class type> 
class RandomBlockMatrix_new : public BlockMatMul1D_derived<type> {
  PA_INJECT(type) 

  const EncryptedArray& ea;
  long dim;

  vector< vector< mat_R > > data;

public:

  RandomBlockMatrix_new(const EncryptedArray& _ea, long _dim):
    ea(_ea), dim(_dim)
  {
    RBak bak; bak.save(); _ea.getAlMod().restoreContext();
    long n = _ea.size();
    long d = _ea.getDegree();
    long D = _ea.sizeOfDimension(dim);


    RandomStreamPush push;
    SetSeed(ZZ(123));

    data.resize(D);
    for (long i = 0; i < D; i++) {
      data[i].resize(D);
      for (long j = 0; j < D; j++) {
        data[i][j].SetDims(d, d);
        for (long u = 0; u < d; u++)
          for (long v = 0; v < d; v++) 
            random(data[i][j][u][v]);
      }
    }
  }

  bool get(mat_R& out, long i, long j, long k) const override
  {
    long D = ea.sizeOfDimension(dim);
    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return ea; }
  long getDim() const override { return dim; }
  bool multipleTransforms() const override { return false; }
};

static BlockMatMul1D*
buildRandomBlockMatrix_new(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomBlockMatrix_new<PA_GF2>(ea,dim);
    }
    case PA_zz_p_tag: {
      return new RandomBlockMatrix_new<PA_zz_p>(ea, dim);
    }
    default: return 0;
  }
}


//********************************

template<class type> 
class RandomMultiBlockMatrix : public BlockMatMul<type> {
  PA_INJECT(type) 

  vector< vector< vector< mat_R > > > data;
  long dim;

public:
  virtual ~RandomMultiBlockMatrix() {}
  RandomMultiBlockMatrix(const EncryptedArray& _ea, long _dim):
    BlockMatMul<type>(_ea), dim(_dim)
  {
    RBak bak; bak.save(); _ea.getAlMod().restoreContext();
    long n = _ea.size();
    long d = _ea.getDegree();
    long D = _ea.sizeOfDimension(dim);

    RandomStreamPush push;
    SetSeed(ZZ(123));

    data.resize(n/D);
    for (long k = 0; k < n/D; k++) {
      data[k].resize(D);
      for (long i = 0; i < D; i++) {
	data[k][i].resize(D);
	for (long j = 0; j < D; j++) {
          data[k][i][j].SetDims(d, d);
	  for (long u = 0; u < d; u++)
	    for (long v = 0; v < d; v++) 
	      random(data[k][i][j][u][v]);
	}
      }
    }
  }

  virtual long size() const // how many transformations
  {
    const EncryptedArray& ea = this->getEA();
    long n = ea.size();
    long D = ea.sizeOfDimension(dim);
    return n/D;
  }

  virtual bool multiGet(mat_R& out, long i, long j, long k) const
  {
    const EncryptedArray& ea = this->getEA();
    long n = ea.size();
    long D = ea.sizeOfDimension(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    assert(k >= 0 && k < n/D);
    if (IsZero(data[k][i][j])) return true;
    out = data[k][i][j];
    return false;
  }
};

static MatMulBase*
buildRandomMultiBlockMatrix(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomMultiBlockMatrix<PA_GF2>(ea, dim);
    }
    case PA_zz_p_tag: {
      return new RandomMultiBlockMatrix<PA_zz_p>(ea, dim);
    }
    default: return 0;
  }
}


template<class type> 
class RandomMultiBlockMatrix_new : public BlockMatMul1D_derived<type> {
  PA_INJECT(type) 

  const EncryptedArray& ea;
  long dim;

  vector< vector< vector< mat_R > > > data;

public:

  RandomMultiBlockMatrix_new(const EncryptedArray& _ea, long _dim):
    ea(_ea), dim(_dim)
  {
    RBak bak; bak.save(); _ea.getAlMod().restoreContext();
    long n = _ea.size();
    long d = _ea.getDegree();
    long D = _ea.sizeOfDimension(dim);

    RandomStreamPush push;
    SetSeed(ZZ(123));

    data.resize(n/D);
    for (long k = 0; k < n/D; k++) {
      data[k].resize(D);
      for (long i = 0; i < D; i++) {
	data[k][i].resize(D);
	for (long j = 0; j < D; j++) {
          data[k][i][j].SetDims(d, d);
	  for (long u = 0; u < d; u++)
	    for (long v = 0; v < d; v++) 
	      random(data[k][i][j][u][v]);
	}
      }
    }
  }


  bool get(mat_R& out, long i, long j, long k) const override
  {
    long n = ea.size();
    long D = ea.sizeOfDimension(dim);

    assert(i >= 0 && i < D);
    assert(j >= 0 && j < D);
    assert(k >= 0 && k < n/D);
    if (IsZero(data[k][i][j])) return true;
    out = data[k][i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return ea; }
  long getDim() const override { return dim; }
  bool multipleTransforms() const override { return true; }
};

static BlockMatMul1D*
buildRandomMultiBlockMatrix_new(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomMultiBlockMatrix_new<PA_GF2>(ea, dim);
    }
    case PA_zz_p_tag: {
      return new RandomMultiBlockMatrix_new<PA_zz_p>(ea, dim);
    }
    default: return 0;
  }
}




template<class type> 
class RandomFullMatrix : public MatMulFull_derived<type> {
  PA_INJECT(type) 
  const EncryptedArray& ea;
  vector<vector<RX>> data;

public:
  RandomFullMatrix(const EncryptedArray& _ea): ea(_ea) {
    long n = ea.size();
    long d = ea.getDegree();
    long bnd = 2*n; // non-zero with probability 1/bnd

    RBak bak; bak.save(); ea.getContext().alMod.restoreContext();
    data.resize(n);
    for (long i: range(n)) {
      data[i].resize(n);
      for (long j: range(n)) {
        //bool zEntry = (RandomBnd(bnd) > 0);
        bool zEntry = false;
        if (zEntry) 
          clear(data[i][j]);
        else
          random(data[i][j], d);
      }
    }
  }

  bool get(RX& out, long i, long j) const override {
    assert(i >= 0 && i < ea.size());
    assert(j >= 0 && j < ea.size());
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return ea; }
};

static MatMulFull* buildRandomFullMatrix(EncryptedArray& ea)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: { return new RandomFullMatrix<PA_GF2>(ea); }
    case PA_zz_p_tag:{ return new RandomFullMatrix<PA_zz_p>(ea); }
    default: return nullptr;
  }
}



template<class type> 
class RandomFullBlockMatrix : public BlockMatMulFull_derived<type> {
  PA_INJECT(type) 
  const EncryptedArray& ea;
  vector<vector<mat_R>> data;

public:
  RandomFullBlockMatrix(const EncryptedArray& _ea): ea(_ea) {
    long n = ea.size();
    long d = ea.getDegree();
    long bnd = 2*n; // non-zero with probability 1/bnd

    RBak bak; bak.save(); ea.getContext().alMod.restoreContext();
    data.resize(n);
    for (long i: range(n)) {
      data[i].resize(n);
      for (long j: range(n)) {
        //bool zEntry = (RandomBnd(bnd) > 0);
        bool zEntry = false;

        data[i][j].SetDims(d, d);
        if (zEntry) 
          clear(data[i][j]);
        else {
	  for (long u = 0; u < d; u++)
	    for (long v = 0; v < d; v++) 
	      random(data[i][j][u][v]);
        }
      }
    }
  }

  bool get(mat_R& out, long i, long j) const override {
    assert(i >= 0 && i < ea.size());
    assert(j >= 0 && j < ea.size());
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return ea; }
};

static BlockMatMulFull* buildRandomFullBlockMatrix(EncryptedArray& ea)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: { return new RandomFullBlockMatrix<PA_GF2>(ea); }
    case PA_zz_p_tag:{ return new RandomFullBlockMatrix<PA_zz_p>(ea); }
    default: return nullptr;
  }
}

template<class Matrix>
void DoTest(const Matrix& mat, const EncryptedArray& ea, 
            const FHESecKey& secretKey)
{
  resetAllTimers();
  typename Matrix::ExecType mat_exec(mat, secretKey);
  mat_exec.upgrade();
  printAllTimers();

  // choose a random plaintext vector
  NewPlaintextArray v(ea);
  random(ea, v);

  // encrypt the random vector
  Ctxt ctxt(secretKey);
  ea.encrypt(ctxt, secretKey, v);
  Ctxt ctxt2 = ctxt;

  resetAllTimers();

  mat_exec.mul(ctxt);

  printAllTimers();


  mul(v, mat);     // multiply the plaintext vector

  NewPlaintextArray v1(ea);
  ea.decrypt(ctxt, secretKey, v1); // decrypt the ciphertext vector

  if (equals(ea, v, v1))        // check that we've got the right answer
    cout << "Nice!!\n";
  else
    cout << "Grrr@*\n";
}


int fhe_test_force_minimal = 0;


void  TestIt(FHEcontext& context, long dim, bool verbose, long full, long block)
{
  if (verbose) {
    context.zMStar.printout();
    cout << endl;
  }

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(/*w=*/64); // A Hamming-weight-w secret key


  if (fhe_test_force_minimal > 0) {
    addMinimal1DMatrices(secretKey);
    addMinimalFrbMatrices(secretKey);
  }
  else {
    addSome1DMatrices(secretKey);
    addSomeFrbMatrices(secretKey); 
  }

  // encrypted array with "full slots"
#if 1
  EncryptedArray ea(context, context.alMod);
#else
  long p = publicKey.getContext().zMStar.getP();
  ZZX G = makeIrredPoly(p, 10);
  EncryptedArray ea(context, G);
#endif


  if (full == 0 && block == 0) {
    std::unique_ptr< MatMul1D > ptr(buildRandomMatrix_new(ea,dim));
    DoTest(*ptr, ea, secretKey);
  }
  else if (full == 0 && block == 1) {
    std::unique_ptr< BlockMatMul1D > ptr(buildRandomBlockMatrix_new(ea,dim));
    DoTest(*ptr, ea, secretKey);
  }
  else if (full == 1 && block == 0) {
    std::unique_ptr< MatMulFull > ptr(buildRandomFullMatrix(ea));
    DoTest(*ptr, ea, secretKey);
  }
  else if (full == 1 && block == 1) {
    std::unique_ptr< BlockMatMulFull > ptr(buildRandomFullBlockMatrix(ea));
    DoTest(*ptr, ea, secretKey);
  }

}


int main(int argc, char *argv[]) 
{
  ArgMapping amap;

  long m=2047;
  amap.arg("m", m, "defines the cyclotomic polynomial Phi_m(X)");
  long p=2;
  amap.arg("p", p, "plaintext base");
  long r=1;
  amap.arg("r", r,  "lifting");
  long L=4;
  amap.arg("L", L, "# of levels in the modulus chain");
  long dim=0;
  amap.arg("dim", dim, "dimension along which to multiply");
  long verbose=0;
  amap.arg("verbose", verbose, "print timing and other info");
  long nt=1;
  amap.arg("nt", nt, "# threads");

  amap.arg("force_bsgs", fhe_test_force_bsgs, 
           "1 to force on, -1 to force off"); 
  amap.arg("force_hoist", fhe_test_force_hoist, 
           "-1 to force off"); 
  amap.arg("force_minimal", fhe_test_force_minimal, 
           "1 to force on"); 

  long full = 0; 
  amap.arg("full", full, "0: 1D, 1: block");

  long block = 0; 
  amap.arg("block", block, "0: 1D, 1: block");

  NTL::Vec<long> gens;
  amap.arg("gens", gens, "use specified vector of generators", NULL);
  amap.note("e.g., gens='[562 1871 751]'");
  NTL::Vec<long> ords;
  amap.arg("ords", ords, "use specified vector of orders", NULL);
  amap.note("e.g., ords='[4 2 -4]', negative means 'bad'");

  amap.parse(argc, argv);

  cout << "*** matmul1D: m=" << m
       << ", p=" << p
       << ", r=" << r
       << ", L=" << L
       << ", dim=" << dim
       << ", nt=" << nt
       // << ", gens=" << gens
       // << ", ords=" << ords
       << ", full=" << full
       << ", block=" << block
       << ", force_bsgs=" << fhe_test_force_bsgs
       << ", force_hoist=" << fhe_test_force_hoist
       << ", force_minimal=" << fhe_test_force_minimal
       << endl;

  vector<long> gens1, ords1;
  convert(gens1, gens);
  convert(ords1, ords);

  if (nt > 1) SetNumThreads(nt);

  setTimersOn();

  FHEcontext context(m, p, r, gens1, ords1);
  buildModChain(context, L, /*c=*/3);

  TestIt(context, dim, verbose, full, block);
  cout << endl;
  if (0 && verbose) {
    printAllTimers();
    cout << endl;
  }
}
