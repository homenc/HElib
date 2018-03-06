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
 * @file randomMatrices.h
 * @brief implementation of random matrices of various forms, used for testing
 */
#include "matmul.h"
#include <NTL/BasicThreadPool.h>

template<class type> class RandomMatrix : public  MatMul1D_derived<type> {
public:
  PA_INJECT(type) 

private:
  vector< vector< RX > > data;
  const EncryptedArray& ea;
  long dim;

public:
  virtual ~RandomMatrix() {}
  RandomMatrix(const EncryptedArray& _ea, long _dim): 
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
buildRandomMatrix(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomMatrix<PA_GF2>(ea, dim);
    }
    case PA_zz_p_tag: {
      return new RandomMatrix<PA_zz_p>(ea, dim);
    }
    default: return 0;
  }
}



template<class type> class RandomMultiMatrix : public  MatMul1D_derived<type> {
public:
  PA_INJECT(type) 

private:
  vector< vector< vector< RX > > > data;
  const EncryptedArray& ea;
  long dim;

public:
  virtual ~RandomMultiMatrix() {}
  RandomMultiMatrix(const EncryptedArray& _ea, long _dim): 
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
buildRandomMultiMatrix(const EncryptedArray& ea, long dim)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomMultiMatrix<PA_GF2>(ea, dim);
    }
    case PA_zz_p_tag: {
      return new RandomMultiMatrix<PA_zz_p>(ea, dim);
    }
    default: return 0;
  }
}



//********************************


template<class type> 
class RandomBlockMatrix : public BlockMatMul1D_derived<type> {
  PA_INJECT(type) 

  const EncryptedArray& ea;
  long dim;

  vector< vector< mat_R > > data;

public:

  RandomBlockMatrix(const EncryptedArray& _ea, long _dim):
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


//********************************



template<class type> 
class RandomMultiBlockMatrix : public BlockMatMul1D_derived<type> {
  PA_INJECT(type) 

  const EncryptedArray& ea;
  long dim;

  vector< vector< vector< mat_R > > > data;

public:

  RandomMultiBlockMatrix(const EncryptedArray& _ea, long _dim):
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

static MatMulFull* buildRandomFullMatrix(const EncryptedArray& ea)
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

static BlockMatMulFull* buildRandomFullBlockMatrix(const EncryptedArray& ea)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: { return new RandomFullBlockMatrix<PA_GF2>(ea); }
    case PA_zz_p_tag:{ return new RandomFullBlockMatrix<PA_zz_p>(ea); }
    default: return nullptr;
  }
}
