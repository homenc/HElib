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
#include <NTL/BasicThreadPool.h>

#if (defined(__unix__) || defined(__unix) || defined(unix))
#include <sys/time.h>
#include <sys/resource.h>
#endif

// Implementation of the various random matrices is found here
#include "randomMatrices.h"
/*
 * Defined in this file are the following class templates:
 *
 *   class RandomMatrix: public MatMul1D_derived<type>
 *   class RandomMultiMatrix: public MatMul1D_derived<type>
 *   class RandomBlockMatrix: public BlockMatMul1D_derived<type>
 *   class RandomMultiBlockMatrix: public BlockMatMul1D_derived<type>
 *   class RandomFullMatrix: public MatMulFull_derived<type>
 *   class RandomFullBlockMatrix : public BlockMatMulFull_derived<type>
 *
 * Each of them has a corresponding build function, namely:
 *
 *   MatMul1D* buildRandomMatrix(const EncryptedArray& ea, long dim);
 *   MatMul1D* buildRandomMultiMatrix(const EncryptedArray& ea, long dim);
 *   BlockMatMul1D* buildRandomBlockMatrix(const EncryptedArray& ea, long dim);
 *   BlockMatMul1D* buildRandomMultiBlockMatrix(const EncryptedArray& ea, long dim);
 *   MatMulFull* buildRandomFullMatrix(const EncryptedArray& ea);
 *   BlockMatMulFull* buildRandomFullBlockMatrix(const EncryptedArray& ea);
 */

template<class Matrix>
bool DoTest(const Matrix& mat, const EncryptedArray& ea, 
            const FHESecKey& secretKey, bool minimal, bool verbose)
{
  FHE_NTIMER_START(EncodeMartix_MatMul);
  typename Matrix::ExecType mat_exec(mat, minimal);
  mat_exec.upgrade();
  FHE_NTIMER_STOP(EncodeMartix_MatMul);

  // choose a random plaintext vector and encrypt it
  NewPlaintextArray v(ea);
  random(ea, v);

  // encrypt the random vector
  Ctxt ctxt(secretKey);
  ea.encrypt(ctxt, secretKey, v);
  Ctxt ctxt2 = ctxt;

  mat_exec.mul(ctxt);

  mul(v, mat);     // multiply the plaintext vector

  NewPlaintextArray v1(ea);
  ea.decrypt(ctxt, secretKey, v1); // decrypt the ciphertext vector

  return equals(ea, v, v1);        // check that we've got the right answer
}

int ks_strategy = 0;
// 0 == default
// 1 == full
// 2 == BSGS
// 3 == minimal


void TestIt(FHEcontext& context, long dim, bool verbose, long full, long block)
{
  resetAllTimers();
  if (verbose)
    context.zMStar.printout();

  FHESecKey secretKey(context);
  const FHEPubKey& publicKey = secretKey;
  secretKey.GenSecKey(/*w=*/64); // A Hamming-weight-w secret key

  bool minimal = ks_strategy == 3;

  // we call addSomeFrbMatrices for all strategies except minimal

  switch (ks_strategy) {
  case 0: 
    addSome1DMatrices(secretKey);
    addSomeFrbMatrices(secretKey);
    break;
  case 1: 
    add1DMatrices(secretKey);
    addSomeFrbMatrices(secretKey);
    break;
  case 2: 
    addBSGS1DMatrices(secretKey);
    addSomeFrbMatrices(secretKey);
    break;
  case 3: 
    addMinimal1DMatrices(secretKey);
    addMinimalFrbMatrices(secretKey);
    break;

   default:
     Error("bad ks_strategy");
   }

  // encrypted array with "full slots"
  EncryptedArray ea(context, context.alMod);

  if (verbose) {
    {std::stringstream ss;
    ss << publicKey;
    cout << "\n  |pubKey|="<<ss.tellp()<<" bytes";}
    cout << ", security=" << context.securityLevel()<<endl;
    cout << "  #threads="<<NTL::AvailableThreads();
    if (block)
      cout << "block-size="<<(ea.getPAlgebra().getOrdP());
    cout << ", vector-dimension="
         << (full? ea.size() : ea.sizeOfDimension(dim))<<", ";
    if (!full)
      cout << (ea.size()/ea.sizeOfDimension(dim))<<" products in parallel, ";
  }

  bool okSoFar = true;
  for (long i=0; i<5; i++) {
    if (full == 0 && block == 0) {
      std::unique_ptr< MatMul1D > ptr(buildRandomMatrix(ea,dim));
      if (!DoTest(*ptr, ea, secretKey, minimal, verbose))
        okSoFar = false;
    }
    else if (full == 0 && block == 1) {
      std::unique_ptr< BlockMatMul1D > ptr(buildRandomBlockMatrix(ea,dim));
      if (!DoTest(*ptr, ea, secretKey, minimal, verbose))
        okSoFar = false;
    }
    else if (full == 1 && block == 0) {
      std::unique_ptr< MatMulFull > ptr(buildRandomFullMatrix(ea));
      if (!DoTest(*ptr, ea, secretKey, minimal, verbose))
        okSoFar = false;
    }
    else if (full == 1 && block == 1) {
      std::unique_ptr< BlockMatMulFull > ptr(buildRandomFullBlockMatrix(ea));
      if (!DoTest(*ptr, ea, secretKey, minimal, verbose))
        okSoFar = false;
    }
  }
  cout << (okSoFar? "Nice!!\n\n" : "Grrr@*\n\n");

  if (verbose) {
    printAllTimers(cout);
#if (defined(__unix__) || defined(__unix) || defined(unix))
      struct rusage rusage;
      getrusage( RUSAGE_SELF, &rusage );
      cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << endl;
#endif
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
  long L=3;
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
  amap.arg("ks_strategy", ks_strategy,
           "0: default, 1:full, 2:bsgs, 3:minimal"); 

  long full = 0; 
  amap.arg("full", full, "0: 1D, 1: full");

  long block = 0; 
  amap.arg("block", block, "0: normal, 1: block");

  NTL::Vec<long> gens;
  amap.arg("gens", gens, "use specified vector of generators", NULL);
  amap.note("e.g., gens='[420 1105 1425]'");
  NTL::Vec<long> ords;
  amap.arg("ords", ords, "use specified vector of orders", NULL);
  amap.note("e.g., ords='[11 8 2]', negative means 'bad'");

  amap.parse(argc, argv);

  if (verbose) {
    cout << "*** Test_MatMul: m=" << m
	 << ", p=" << p
	 << ", r=" << r
	 << ", L=" << L
	 << ", dim=" << dim
	 << ", nt=" << nt
	 << ", full=" << full
	 << ", block=" << block
	 << ", force_bsgs=" << fhe_test_force_bsgs
	 << ", force_hoist=" << fhe_test_force_hoist
	 << ", ks_strategy=" << ks_strategy
	 << endl;
   }

  vector<long> gens1, ords1;
  convert(gens1, gens);
  convert(ords1, ords);

  if (nt > 1) SetNumThreads(nt);

  setTimersOn();

  FHEcontext context(m, p, r, gens1, ords1);
  buildModChain(context, L, /*c=*/3);

  TestIt(context, dim, verbose, full, block);
}
