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

// In the CKKS encryption scheme, besides SIMD operations that act on the slots
// of a ciphertext in parallel, it is also possible to move data around among
// the slots of a ciphertext.

#include <helib/helib.h>

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

  // To support data movement, we need to add some information to the public
  // key. This is done as follows:
  addSome1DMatrices(secretKey);

  // Recall that SecKey is a subclass of PubKey. The call to addSome1DMatrices
  // needs data stored in the secret key, but the information it computes is
  // stored in the public key.

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

  // We can rotate the data in the slots by any amount.

  rotate(c, 2);
  // rotate c right by 2:
  // (c[0], ..., c[n-1]) = (c[n-2], c[n-1], c[0], c[1], ..., c[n-3])

  cout << "c.capacity=" << c.capacity() << " ";
  cout << "c.errorBound=" << c.errorBound() << "\n";

  rotate(c, -1);
  // rotate c left by 1
  // (c[0], ..., c[n-1]) = (c[1], c[2], ..., c[n-1], c[0])

  cout << "c.capacity=" << c.capacity() << " ";
  cout << "c.errorBound=" << c.errorBound() << "\n";

  //===========================================================================

  // We can shift the data in the slots by any amount.

  shift(c, 2);
  // rotate c right by 2:
  // (c[0], ..., c[n-1]) = (0, 0, c[0], c[1], ..., c[n-3])

  cout << "c.capacity=" << c.capacity() << " ";
  cout << "c.errorBound=" << c.errorBound() << "\n";

  shift(c, -1);
  // rotate c left by 1
  // (c[0], ..., c[n-1]) = (c[1], c[2], ..., c[n-1], 0)

  cout << "c.capacity=" << c.capacity() << " ";
  cout << "c.errorBound=" << c.errorBound() << "\n";

  //===========================================================================

  // We can also sum all of slots, leaving the sum in each slot

  totalSums(c);
  // (c[0], ..., c[n-1]) = (S, ..., S), where S = sum_{i=0}^{n-1} c[i]

  cout << "c.capacity=" << c.capacity() << " ";
  cout << "c.errorBound=" << c.errorBound() << "\n";

  // There are a number of other data movement operations available.

  //===========================================================================

  // Let's perform the same computation on the plaintext:

  rotate(p, 2);
  rotate(p, -1);
  shift(p, 2);
  shift(p, -1);
  totalSums(p);

  //===========================================================================

  // Let's decrypt and compare:
  PtxtArray pp(context);
  pp.decrypt(c, secretKey);

  double distance = Distance(p, pp);
  cout << "distance=" << distance << "\n";

  // For debugging, you can also make "approximate" comparisons as follows:
  if (pp == Approx(p))
    cout << "GOOD\n";
  else
    cout << "BAD\n";

  // Here, p is the "correct value" and you want to test if pp is "close" to it.

  // NOTES: The Approx function (which is really a class constructor) takes two
  // optional arguments:
  //   double tolerance; // default is 0.01
  //   double floor;     // default is 1.0
  //
  // The expression
  //   a == Approx(b, tolerance, floor)
  // is true iff Distance(a,b) <= tolerance*max(Norm(b),floor), The idea is
  // that it checks if the relative error is at most tolerance, unless Norm(b)
  // itself is too small (as determined by floor). Here, Norm(b) is the max
  // absolute value of the slots, and Distance(a,b) = Norm(a-b).
  //
  // In addition to PtxtArray's, you can compare values of type double or
  // complex<double>, and vectors of type double or complex<double>.

  return 0;
}
