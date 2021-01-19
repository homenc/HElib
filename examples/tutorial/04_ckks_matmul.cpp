/* Copyright (C) 2020-2021 IBM Corp.
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

// In the CKKS encryption scheme, since a ciphertext encrypts a vector of
// slots, it makes sense to multiply that vector by a matrix.  HElib provides
// highly optimized routines for multiplying an encrypted vector by a plaintext
// matrix.

#include <helib/helib.h>

using namespace std;
using namespace helib;

// To use these routines, we need to include an extra file:
#include <helib/matmul.h>

// In this example, we will also make some performance measurements.  HElib
// provides convenient "timers" to measure running time.  We will also be
// measuring space. For this, we will use the getrusage function, if available:

#if defined(__unix__) || defined(__unix) || defined(unix)
#include <sys/time.h>
#include <sys/resource.h>
void printMemoryUsage()
{
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  cout << "  ru_maxrss=" << r.ru_maxrss << endl;
}
#else
void printMemoryUsage() {}
#endif

int main(int argc, char* argv[])
{
  Context context =
      ContextBuilder<CKKS>().m(16 * 1024).bits(119).precision(30).c(2).build();

  cout << "securityLevel=" << context.securityLevel() << "\n";

  long n = context.getNSlots();

  SecKey secretKey(context);
  secretKey.GenSecKey();

  addSome1DMatrices(secretKey);

  const PubKey& publicKey = secretKey;

  //===========================================================================

  // Let's encrypt something!
  vector<double> v(n);
  for (long i = 0; i < n; i++)
    v[i] = sin(2.0 * PI * i / n);
  PtxtArray p(context, v);
  Ctxt c(publicKey);
  p.encrypt(c);

  cout << "c.capacity=" << c.capacity() << " ";
  cout << "c.errorBound=" << c.errorBound() << "\n";

  //===========================================================================

  // We define an n x n plaintext matrix as follows:
  MatMul_CKKS mat(context,
                  [n](long i, long j) { return ((i + j) % n) / double(n); });

  // Note that the second parameter of the MatMul_CKKS constructor is of type
  // std::function<double(long, long)>, meaning that it should be a
  // function-like object that takes two long's and returns a double.  In this
  // example, the actual parameter is a C++ "lambda" object. In input (i, j),
  // this should return the value of the matrix in row i and column j.

  // We now multiply ciphertext c by this matrix:
  c *= mat;

  // Note that this computes c = c*mat, where the slots of c are viewed as a
  // row vector

  cout << "c.capacity=" << c.capacity() << " ";
  cout << "c.errorBound=" << c.errorBound() << "\n";

  // We can multiply the plaintext p by the same matrix:
  p *= mat;

  //===========================================================================

  // Let's decrypt and compare:
  PtxtArray pp(context);
  pp.decrypt(c, secretKey);

  double distance = Distance(p, pp);
  cout << "distance=" << distance << "\n";

  //===========================================================================

  // If a given matrix is going to be used many times, one can obtain better
  // performance by doing a one-time pre-computation. Let's begin by
  // performing the same ciphertext/matrix multiplication, but this
  // time, let's measure the running time.  HElib provides a convenient
  // mechanism for doing this:

  Ctxt c0 = c;
  HELIB_NTIMER_START(mul0); // starts a timer called "mul0"
  c0 *= mat;
  HELIB_NTIMER_STOP(mul0); // stops the time "mul0"
  printNamedTimer(cout, "mul0");
  // On my machine, this took about 4.6s

  // A pre-computation is performed by "encoding" the matrix, as follows:
  HELIB_NTIMER_START(encode);
  EncodedMatMul_CKKS emat(mat);
  HELIB_NTIMER_STOP(encode);
  printNamedTimer(cout, "encode");
  // On my machine, this took about 2.4s

  // We can apply the encoded matrix to a ciphertext as follows:
  Ctxt c1 = c;

  {
    HELIB_NTIMER_START(mul1);
    c1 *= emat;
  } // The timer "mul1" automatically gets stopped when control exits the block
  printNamedTimer(cout, "mul1");
  // On my machine, this took about 2.4s

  // We can perform even more precomputation, but it takes up more space.
  // First, let's see how much space we are currently using:
  printMemoryUsage();
  // On my machine, the memory footprint is now about 340MB

  // Now do more pre-computation:
  {
    HELIB_NTIMER_START(upgrade);
    emat.upgrade();
  }
  printNamedTimer(cout, "upgrade");
  // On my machine, this took about 1.7s

  // And let's see how much the space increased:
  printMemoryUsage();
  // On my machine, the memory footprint is now about 850MB,
  // and so the upgrade costs about 510MB of space.

  // Now we apply the upgraded encoded matrix to the ciphertext:
  Ctxt c2 = c;

  {
    HELIB_NTIMER_START(mul2);
    c2 *= emat;
  }
  printNamedTimer(cout, "mul2");
  // On my machine, this took about 1.4s. Compared the original time of 4.6s,
  // we see a roughly 3.3x speedup.

  return 0;
}
