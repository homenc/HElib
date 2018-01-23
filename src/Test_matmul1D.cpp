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
#include "FHE.h"
#include "timing.h"
#include "EncryptedArray.h"
#include "matmul.h"

// Forward declerations
static MatMulBase* buildRandomMatrix(const EncryptedArray& ea,
                                     long dim, long giantStep);
static MatMulBase* buildRandomMultiMatrix(const EncryptedArray& ea,
                                          long dim, long giantStep);
static MatMulBase* buildRandomBlockMatrix(const EncryptedArray& ea,
                                          long dim);
static MatMulBase* buildRandomMultiBlockMatrix(const EncryptedArray& ea,
                                               long dim);


// The callback interface for the matrix-multiplication routines.

//! \cond FALSE (make doxygen ignore these classes)
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
//********************************

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
//! \endcond

void  TestIt(FHEcontext& context, long g, long dim, bool verbose)
{
  if (verbose) {
    context.zMStar.printout();
    cout << endl;
  }

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(/*w=*/64); // A Hamming-weight-w secret key

  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need
  EncryptedArray ea(context, context.alMod);

  // Test a "normal" matrix over the extension field
  {
    // choose a random plaintext square matrix
    std::unique_ptr< MatMulBase > ptr(buildRandomMatrix(ea, dim, g));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << " Multiplying 1D with MatMulBase... " << std::flush;
    matMul1D(v, *ptr, dim);
    matMul1D(ctxt2, *ptr, dim, cachezzX);
       // multiply ciphertext and build cache
    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";

    cout << " Multiplying 1D with MatMulBase+dcrt cache... " << std::flush;
    ctxt2 = ctxt;
    matMul1D(ctxt2, *ptr, dim, cacheDCRT); // upgrade cache and use in multiplication

    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }
  {
    // choose a random plaintext square matrix
    std::unique_ptr< MatMulBase > ptr(buildRandomMatrix(ea,dim,g));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << " Multiplying 1D with MatMulBase+zzx cache... " << std::flush;
    buildCache4MatMul1D(*ptr, dim, cachezzX);// build the cache
    matMul1D(ctxt, *ptr, dim);               // then use it
    matMul1D(v, *ptr, dim);     // multiply the plaintext vector

    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }

  // Test a "multi" matrix over the extension field
  {
    // choose a random plaintext square matrix
    std::unique_ptr< MatMulBase > ptr(buildRandomMultiMatrix(ea,dim,g));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << "\n Multiplying multi 1D with MatMulBase... " << std::flush;
    matMulti1D(v, *ptr, dim);
    matMulti1D(ctxt2, *ptr, dim, cachezzX); // multiply ciphertext and build cache
    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";

    cout <<" Multiplying multi 1D with MatMulBase+dcrt cache... "<< std::flush;
    ctxt2 = ctxt;
    matMulti1D(ctxt2, *ptr, dim, cacheDCRT); // upgrade cache and use in multiplication

    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))    // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }
  {
    // choose a random plaintext square matrix
    std::unique_ptr< MatMulBase > ptr(buildRandomMultiMatrix(ea,dim,g));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << " Multiplying multi 1D with MatMulBase+zzx cache... " << std::flush;
    buildCache4MatMulti1D(*ptr, dim, cachezzX);// build the cache
    matMulti1D(ctxt, *ptr, dim);               // then use it
    matMulti1D(v, *ptr, dim);     // multiply the plaintext vector

    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }

  // Test a "block matrix" over the base field
  {
    // choose a random plaintext square matrix
    shared_ptr<MatMulBase> ptr(buildRandomBlockMatrix(ea,dim));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << endl << " Multiplying 1D with BlockMatMul... " 
	 << std::flush;
    blockMatMul1D(ctxt2, *ptr, dim, cachezzX);
                                  // multiply ciphertext and build cache
    blockMatMul1D(v, *ptr, dim); // multiply the plaintext vector
    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";

    cout << " Multiplying 1D with BlockMatMul+dcrt cache... " << std::flush;
    ctxt2 = ctxt;
    blockMatMul1D(ctxt2, *ptr, dim, cacheDCRT);
                                  // upgrade cache and use in multiplication
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }
  {
    // choose a random plaintext square matrix
    shared_ptr<MatMulBase> ptr(buildRandomBlockMatrix(ea,dim));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << " Multiplying 1D with BlockMatMul+zzx cache... " << std::flush;
    buildCache4BlockMatMul1D(*ptr, dim, cachezzX);// build the cache
    blockMatMul1D(ctxt2, *ptr, dim);              // then use it
    blockMatMul1D(v, *ptr, dim); // multiply the plaintext vector
    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }
  // Test multiple "block" matrices over the base field
  {
    // choose a random plaintext square matrix
    shared_ptr<MatMulBase> ptr(buildRandomMultiBlockMatrix(ea,dim));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << " Multiplying multi 1D with BlockMatMul... " << std::flush;
    blockMatMulti1D(v, *ptr, dim);     // multiply the plaintext vector
    blockMatMulti1D(ctxt2, *ptr, dim); // multiply the ciphertext vector

    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";

  }
}


/* Testing the functionality of multiplying an encrypted vector by a
 * plaintext matrix, either over the extension- or the base-field/ring.
 *
 * Usage: Test_matmul1D [optional params]
 *
 *  m defines the cyclotomic polynomial Phi_m(X)
 *  p is the plaintext base [default=2]
 *  r is the lifting [default=1]
 *  g is the giant-step parameter [defauls=2]
 *  L is the # of primes in the modulus chain [default=4]
 *  dim is the dimension alng which we multiply [default=0]
 *  verbose print timing info [default=0]
 */
int main(int argc, char *argv[]) 
{
  ArgMapping amap;

  long m=2047;
  amap.arg("m", m, "defines the cyclotomic polynomial Phi_m(X)");
  long p=2;
  amap.arg("p", p, "plaintext base");
  long r=1;
  amap.arg("r", r,  "lifting");
  long g=2;
  amap.arg("g", g,  "giant-step parameter");
  long L=4;
  amap.arg("L", L, "# of levels in the modulus chain");
  long dim=0;
  amap.arg("dim", dim, "dimension along which to multiply");
  long verbose=0;
  amap.arg("verbose", verbose, "print timing and other info");

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
       << ", g=" << g
       << ", dim=" << dim
       // << ", gens=" << gens
       // << ", ords=" << ords
       << endl;

  vector<long> gens1, ords1;
  convert(gens1, gens);
  convert(ords1, ords);

  setTimersOn();

  FHEcontext context(m, p, r, gens1, ords1);
  buildModChain(context, L, /*c=*/3);

  TestIt(context, g, dim, verbose);
  cout << endl;
  if (verbose) {
    printAllTimers();
    cout << endl;
  }
}
