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
static MatMulBase* buildRandomMatrix(const EncryptedArray& ea, long dim);
static MatMulBase* buildRandomMultiMatrix(const EncryptedArray& ea, long dim);
static MatMulBase* buildRandomBlockMatrix(const EncryptedArray& ea, long dim);
static MatMulBase*
buildRandomMultiBlockMatrix(const EncryptedArray& ea, long dim);


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
  RandomMultiMatrix(const EncryptedArray& _ea, long _dim)
    : MatMul<type>(_ea), dim(_dim)
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

static MatMulBase* buildRandomMultiMatrix(const EncryptedArray& ea, long dim)
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
  RandomMatrix(const EncryptedArray& _ea, long _dim): 
    MatMul<type>(_ea), dim(_dim)
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

static MatMulBase* buildRandomMatrix(const EncryptedArray& ea, long dim)
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

static MatMulBase* buildRandomBlockMatrix(const EncryptedArray& ea, long dim)
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

void  TestIt(long m, long p, long r, long d, long L, long dim)
{
  cout << "*** TestIt: m=" << m
       << ", p=" << p
       << ", r=" << r
       << ", d=" << d
       << ", L=" << L
       << ", dim=" << dim
       << endl;

  FHEcontext context(m, p, r);
  buildModChain(context, L, /*c=*/3);

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
  addFrbMatrices(secretKey); // compute key-switching matrices that we need
  cout << "done\n";

  cout << "computing masks and tables for rotation...";
  EncryptedArray ea(context, G);
  cout << "done\n";

  // Test a "normal" matrix over the extension field
  {
    // choose a random plaintext square matrix
    std::unique_ptr< MatMulBase > ptr(buildRandomMatrix(ea, dim));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << " Multiplying 1D with MatMulBase... " << std::flush;
    matMul1D(v, *ptr, dim);
    matMul1D(ctxt2, *ptr, dim, cachezzX);// multiply ciphertext and build cache
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
    std::unique_ptr< MatMulBase > ptr(buildRandomMatrix(ea,dim));

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
    std::unique_ptr< MatMulBase > ptr(buildRandomMultiMatrix(ea,dim));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << "\n Multiplying multi 1D with MatMulBase... " << std::flush;
    matMulti1D(v, *ptr, dim);
    matMulti1D(ctxt2, *ptr, dim, cachezzX);// multiply ciphertext and build cache
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
    std::unique_ptr< MatMulBase > ptr(buildRandomMultiMatrix(ea,dim));

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
      cout << "Grrr...\n";

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
      cout << "Grrr...\n";
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


void usage(char *prog) 
{
  cout << "Usage: "<<prog<<" [ optional parameters ]...\n";
  cout << "  optional parameters have the form 'attr1=val1 attr2=val2 ...'\n";
  cout << "  e.g, 'm=2047 p=2 L=4'\n\n";
  cout << "  m defines the cyclotomic polynomial Phi_m(X)\n";
  cout << "  p is the plaintext base [default=2]" << endl;
  cout << "  r is the lifting [default=1]" << endl;
  cout << "  d is the degree of the field extension [default==1]\n";
  cout << "    (d == 0 => factors[0] defined the extension)\n";
  cout << "  L is the # of primes in the modulus chain [default=4]\n";
  exit(0);
}

/* Testing the functionality of multiplying an encrypted vector by a plaintext
 * matrix, either over the extension- or the base-field/ring.
 */
int main(int argc, char *argv[]) 
{
  argmap_t argmap;
  argmap["m"] = "2047";
  argmap["p"] = "2";
  argmap["r"] = "1";
  argmap["d"] = "0";
  argmap["L"] = "4";
  argmap["dim"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long m = atoi(argmap["m"]);
  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long d = atoi(argmap["d"]);
  long L = atoi(argmap["L"]);
  long dim = atoi(argmap["dim"]);

  setTimersOn();
  TestIt(m, p, r, d, L, dim);
  cout << endl;
}
