/* Copyright (C) 2020 IBM Corp.
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
#include <helib/matmul.h>
#include <helib/NumbTh.h>
#include <NTL/BasicThreadPool.h>
#include <helib/fhe_stats.h>

#if (defined(__unix__) || defined(__unix) || defined(unix))
#include <sys/time.h>
#include <sys/resource.h>
#endif

#include <helib/ArgMap.h>

NTL_CLIENT
using namespace helib;



int ks_strategy = 0;
// 0 == default
// 1 == full
// 2 == BSGS
// 3 == minimal


void TestIt(Context& context, bool verbose)
{
  resetAllTimers();
  if (verbose) {
    context.zMStar.printout();
    std::cout << "# small primes = " << context.smallPrimes.card() << "\n";
    std::cout << "# ctxt primes = " << context.ctxtPrimes.card() << "\n";
    std::cout << "# bits in ctxt primes = " 
	 << long(context.logOfProduct(context.ctxtPrimes)/log(2.0) + 0.5) << "\n";
    std::cout << "# special primes = " << context.specialPrimes.card() << "\n";
    std::cout << "# bits in special primes = " 
	 << long(context.logOfProduct(context.specialPrimes)/log(2.0) + 0.5) << "\n";

    fhe_stats = true;
  }

  SecKey secretKey(context);
  secretKey.GenSecKey(); // A Hamming-weight-w secret key



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

  long n = context.getNSlots();

  std::vector<double> v(n);
  for (long j: range(n)) v[j] = RandomReal();

  MatMul_CKKS_Complex mat(context, 
                [n](long i, long j) 
                { return ((i+j)% n)/double(n); } );
  // Note the use of a "lambda": this allows for quite
  // general ways to describe a matrix with minimal fuss

  //const PubKey publicKey = secretKey;
  const PubKey& publicKey = secretKey;

  Ctxt ctxt(publicKey);
  PtxtArray ptxt(context, v); 
  // initialze ptxt with the vector v

  ptxt.encrypt(ctxt);
  // Encrypt ptxt. We have to supply *some* upper bound on the magnitude
  // of the slots.

  // perform the linear transformation on encrypted data
  EncodedMatMul_CKKS emat(mat);
  //emat.upgrade();
  ctxt *= emat;

  // we can also just do it this way:
  // ctxt *= mat;

  ptxt *= mat;
  // perform the linear transformation on plaintext data

  PtxtArray ptxt1(context);
  ptxt1.decrypt(ctxt, secretKey);

#if 0
  std::vector<double> w1;
  ptxt1.store(w1);
  // w1 is the result of performing the transformation on
  // encrypted data
  
  std::vector<double> w;
  ptxt.store(w);
  // w1 is the result of performing the transformation on
  // plaintext data

  std::cout << w << "\n";
  std::cout << w1 << "\n";
#endif

  std::cout << Distance(ptxt, ptxt1) << "\n";
  //std::cout << Norm(ptxt) << "\n";
  //std::cout << Norm(ptxt1) << "\n";

  if (verbose) {
    printAllTimers(cout);
#if (defined(__unix__) || defined(__unix) || defined(unix))
      struct rusage rusage;
      getrusage( RUSAGE_SELF, &rusage );
      std::cout << "  rusage.ru_maxrss="<<rusage.ru_maxrss << std::endl;
#endif
    print_stats(cout);
  }
}


int main(int argc, char *argv[]) 
{
  helog.setLogToStderr();

  ArgMap amap;

  long m=16;
  amap.arg("m", m, "defines the cyclotomic polynomial Phi_m(X)");
  long r=10;
  amap.arg("r", r,  "precision");
  long L=200;
  amap.arg("L", L, "# of bits in the modulus chain");
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


  amap.parse(argc, argv);

  if (verbose) {
    std::cout << "*** Test_MatMul: m=" << m
	 << ", r=" << r
	 << ", L=" << L
	 << ", nt=" << nt
	 << ", force_bsgs=" << fhe_test_force_bsgs
	 << ", force_hoist=" << fhe_test_force_hoist
	 << ", ks_strategy=" << ks_strategy
	 << std::endl;

   }

  if (nt > 1) SetNumThreads(nt);

  setTimersOn();

  Context context(m, /*p=*/-1, r);
  buildModChain(context, L, /*c=*/3);

  TestIt(context, verbose);
}
