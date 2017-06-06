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

static MatMulBase* buildRandomMatrix(EncryptedArray& ea);
static MatMulBase* buildRandomBlockMatrix(const EncryptedArray& ea);


void  TestIt(long m, long p, long r, long d, long L, bool verbose)
{
  cout << "*** TestIt: m=" << m
       << ", p=" << p
       << ", r=" << r
       << ", d=" << d
       << ", L=" << L
       << endl;

  FHEcontext context(m, p, r);
  buildModChain(context, L, /*c=*/3);

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(/*w=*/64); // A Hamming-weight-w secret key

  ZZX G;
  if (d == 0)
    G = context.alMod.getFactorsOverZZ()[0];
  else
    G = makeIrredPoly(p, d); 

  if (verbose) {
    context.zMStar.printout();
    cout << endl;
    cout << "G = " << G << "\n";
  }

  addSome1DMatrices(secretKey); // compute key-switching matrices that we need
  addFrbMatrices(secretKey); // compute key-switching matrices that we need
  EncryptedArray ea(context, G);

  // Test a "dense" matrix over the extension field
  {
    // choose a random plaintext square matrix
    unique_ptr<MatMulBase> ptr(buildRandomMatrix(ea));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << " Multiplying with MatMulBase... " << std::flush;
    matMul(ctxt2, *ptr, cachezzX); // multiply ciphertext and build cache
    matMul(v, *ptr);     // multiply the plaintext vector

    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";

    cout << " Multiplying with MatMulBase+dcrt cache... " << std::flush;
    ctxt2 = ctxt;
    matMul(ctxt2, *ptr, cacheDCRT); // upgrade cache and use in multiplication

    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }
  {
    // choose a random plaintext square matrix
    unique_ptr<MatMulBase> ptr(buildRandomMatrix(ea));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    cout << " Multiplying with MatMulBase+zzx cache... " << std::flush;
    buildCache4MatMul(*ptr, cachezzX);// build the cache
    matMul(ctxt, *ptr);               // then use it
    matMul(v, *ptr);     // multiply the plaintext vector

    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }
  // Test a "diagonal sparse" matrix over the extension field
  {
    // choose a random plaintext square matrix
    unique_ptr<MatMulBase> ptr(buildRandomMatrix(ea));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << "\n Multiplying with Sparse MatMulBase... " << std::flush;
    matMul_sparse(ctxt2, *ptr, cachezzX); // multiply ciphertext and build cache
    matMul(v, *ptr);     // multiply the plaintext vector

    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";

    cout << " Multiplying with Sparse MatMulBase+dcrt cache... " << std::flush;
    ctxt2 = ctxt;
    matMul_sparse(ctxt2, *ptr, cacheDCRT); // upgrade the cache

    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }
  {
    // choose a random plaintext square matrix
    unique_ptr<MatMulBase> ptr(buildRandomMatrix(ea));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);

    cout << " Multiplying with Sparse MatMulBase+zzx cache... " << std::flush;
    buildCache4MatMul_sparse(*ptr, cachezzX); // build the cache
    matMul_sparse(ctxt, *ptr);                // then use it
    matMul(v, *ptr);                          // multiply plaintext vector

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
    std::unique_ptr<MatMulBase> ptr(buildRandomBlockMatrix(ea));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);
    Ctxt ctxt2 = ctxt;

    cout << "\n Multiplying with BlockMatMul... "  << std::flush;
    blockMatMul(ctxt2, *ptr, cachezzX); // multiply ciphertext and build cache
    blockMatMul(v, *ptr);      // multiply the plaintext vector

    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector

    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr...\n";

    cout << " Multiplying with BlockMatMul+dcrt cache... " << std::flush;
    ctxt2 = ctxt;
    blockMatMul(ctxt2, *ptr, cacheDCRT); // upgrade the cache

    ea.decrypt(ctxt2, secretKey, v1); // decrypt the ciphertext vector
    if (equals(ea, v, v1))        // check that we've got the right answer
      cout << "Nice!!\n";
    else
      cout << "Grrr@*\n";
  }
  {
    // choose a random plaintext square matrix
    std::unique_ptr<MatMulBase> ptr(buildRandomBlockMatrix(ea));

    // choose a random plaintext vector
    NewPlaintextArray v(ea);
    random(ea, v);

    // encrypt the random vector
    Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, v);

    cout << " Multiplying with BlockMatMul+zzx cache... " << std::flush;
    buildCache4BlockMatMul(*ptr, cachezzX); // build the cache
    blockMatMul(ctxt, *ptr);                // then use it
    blockMatMul(v, *ptr);     // multiply the plaintext vector

    NewPlaintextArray v1(ea);
    ea.decrypt(ctxt, secretKey, v1); // decrypt the ciphertext vector

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
  cout << "  verbose print timing info [default=0]\n";
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
  argmap["d"] = "1";
  argmap["L"] = "4";
  argmap["verbose"] = "0";

  // get parameters from the command line
  if (!parseArgs(argc, argv, argmap)) usage(argv[0]);

  long m = atoi(argmap["m"]);
  long p = atoi(argmap["p"]);
  long r = atoi(argmap["r"]);
  long d = atoi(argmap["d"]);
  long L = atoi(argmap["L"]);
  bool v = atoi(argmap["verbose"]);

  //  setTimersOn();
  setTimersOn();
  TestIt(m, p, r, d, L, v);
  cout << endl;
  if (v) {
    printAllTimers();
    cout << endl;
  }
}




template<class type> class RandomMatrix : public MatMul<type> {
  PA_INJECT(type) 
  vector< vector< RX > > data;

public:
  ~RandomMatrix() {/*cout << "destructor: random dense matrix\n";*/}
  RandomMatrix(const EncryptedArray& _ea): MatMul<type>(_ea) {
    long n = _ea.size();
    long d = _ea.getDegree();
    long bnd = 2*n; // non-zero with probability 1/bnd

    RBak bak; bak.save(); _ea.getContext().alMod.restoreContext();
    data.resize(n);
    for (long i = 0; i < n; i++) {
      data[i].resize(n);
      for (long j = 0; j < n; j++) {
        bool zEntry = (RandomBnd(bnd) > 0);
        if (zEntry) clear(data[i][j]);
        else        random(data[i][j], d);
      }
    }
  }

  virtual bool get(RX& out, long i, long j) const {
    assert(i >= 0 && i < this->getEA().size());
    assert(j >= 0 && j < this->getEA().size());
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }
};
static MatMulBase* buildRandomMatrix(EncryptedArray& ea)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: { return new RandomMatrix<PA_GF2>(ea); }
    case PA_zz_p_tag:{ return new RandomMatrix<PA_zz_p>(ea); }
    default: return nullptr;
  }
}


template<class type> class RandomBlockMatrix : public BlockMatMul<type> {
  PA_INJECT(type)

  std::vector< std::vector< mat_R > > data;

public:
  virtual ~RandomBlockMatrix() {}
  RandomBlockMatrix(const EncryptedArray& _ea): BlockMatMul<type>(_ea)
  { 
    RBak bak; bak.save(); _ea.getAlMod().restoreContext();
    long n = _ea.size();
    long d = _ea.getDegree();
    long bnd = 2*n; // non-zero with probability 1/bnd

    data.resize(n);
    for (long i = 0; i < n; i++) {
      data[i].resize(n);
      for (long j = 0; j < n; j++) {
        data[i][j].SetDims(d, d);

        bool zEntry = (RandomBnd(bnd) > 0);

        for (long u = 0; u < d; u++)
          for (long v = 0; v < d; v++) 
            if (zEntry) 
              clear(data[i][j][u][v]);
            else
              random(data[i][j][u][v]);
      }
    }
  }

  virtual bool get(mat_R& out, long i, long j) const {
    assert(i >= 0 && i < this->getEA().size());
    assert(j >= 0 && j < this->getEA().size());
    if (IsZero(data[i][j])) return true;
    out = data[i][j];
    return false;
  }

  const std::vector< std::vector< mat_R > >& getData() const {return data;}
};

static MatMulBase* buildRandomBlockMatrix(const EncryptedArray& ea)
{
  switch (ea.getTag()) {
    case PA_GF2_tag: {
      return new RandomBlockMatrix<PA_GF2>(ea);
    }

    case PA_zz_p_tag: {
      return new RandomBlockMatrix<PA_zz_p>(ea);
    }

    default: return 0;
  }
}
