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

// In the CKKS encryption scheme, ciphertexts can encrypt vectors of complex
// numbers.

#include <helib/helib.h>

#include <helib/matmul.h>
// This is only needed if you want to do matrix multiplication

using namespace std;
using namespace helib;

int main(int argc, char* argv[])
{
  Context context =
      ContextBuilder<CKKS>().m(32 * 1024).bits(358).precision(30).c(6).build();

  cout << "securityLevel=" << context.securityLevel() << "\n";

  long n = context.getNSlots();

  SecKey secretKey(context);
  secretKey.GenSecKey();

  addSome1DMatrices(secretKey);
  // This only needs to be done if you want to do matrix multiplication

  addSomeFrbMatrices(secretKey);
  // This only needs to be done if you want to do conjugation

  const PubKey& publicKey = secretKey;

  //===========================================================================

  // Let's encrypt something!
  vector<std::complex<double>> v0(n);
  for (long i = 0; i < n; i++)
    v0[i] = std::complex<double>(cos(2.0 * PI * i / n), sin(2.0 * PI * i / n));

  // A PtxtArray can be initialized with a vector of complex numbers
  PtxtArray p0(context, v0);

  // Encryption works the same as with real numbers
  Ctxt c0(publicKey);
  p0.encrypt(c0);

  //===========================================================================

  // We next create another ciphertext that encrypts random complex numbers:

  PtxtArray p1(context);
  p1.randomComplex();
  // this fills each entry of p1 with a random number in the complex
  // unit circle

  Ctxt c1(publicKey);
  p1.encrypt(c1);

  //===========================================================================

  // We can perform homomorphic computations in the same way as we did before:

  Ctxt c2 = c0;
  c2 *= 2.5;
  c2 += c1;

  // Note that there is no direct support for combining a ciphertext with a
  // complex scalar. This can be achieved, however, by first converting the
  // complex scaler to a PtxtArray.  For example:

  PtxtArray I(context, std::complex<double>(0.0, 1.0));
  // I has the imaginary unit in each slot

  c2 *= I;

  cout << "c2.capacity=" << c2.capacity() << " ";
  cout << "c2.errorBound=" << c2.errorBound() << "\n";

  // Data movement operations, like rotate and shift, work exactly as before.

  // There is also support for multiplying a ciphertext by a plaintext matrix
  // of complex numbers.  In 04_ckks_matmul.cpp, we showed how you could
  // specify an n x n matrix of real numbers using the class MatMul_CKKS. One
  // can specify an n x n matrix of complex numbers as follows:

  MatMul_CKKS_Complex mat(context, [n](long i, long j) {
    return std::complex<double>(i, j) / double(n);
  });

  c2 *= mat;

  cout << "c2.capacity=" << c2.capacity() << " ";
  cout << "c2.errorBound=" << c2.errorBound() << "\n";

  //===========================================================================

  // One can homomorphically compute the complex conjugate of each slot
  // of a ciphertext as follows:

  conjugate(c2);

  cout << "c2.capacity=" << c2.capacity() << " ";
  cout << "c2.errorBound=" << c2.errorBound() << "\n";

  //===========================================================================

  // Let's decrypt the results:

  PtxtArray pp2(context);
  pp2.decryptComplex(c2, secretKey);

  // Note that if one just writes pp2.decrypt(c2, secretKey) instead of
  // pp2.decryptComplex(c2, secretKey), the imaginary part will be discarded.

  // We can store pp2 in a standard vector, as usual:

  std::vector<std::complex<double>> v2;
  pp2.store(v2);

  //===========================================================================

  // We can also perform the computation on plaintexts and compare:

  PtxtArray p2 = p0;
  p2 *= 2.5;
  p2 += p1;
  p2 *= I;
  p2 *= mat;
  conjugate(p2);

  double distance = Distance(p2, pp2);
  cout << "distance=" << distance << "\n";
}
