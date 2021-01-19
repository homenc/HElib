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

#include <helib/helib.h>

using namespace std;
using namespace helib;

// In the CKKS encryption scheme, ciphertexts have a certain amount of "noise".
// This noise increases with the depth of a homomorphic computation, where by
// "depth", we mean the depth of the arithmetic circuit representing the
// computation. Noise negatively impacts a homomorphic computation in two ways:
// as it grows, it reduces both the *capacity* and the *accuracy* of a
// ciphertext.
//
// The capacity of a ciphertext starts out as some number which is a little
// less than to the bits parameter specified when building a Context object,
// and it is reduced by some amount by each homomorphic computation.  When the
// capacity drops below 1, the ciphertext can no longer be decrypted.
//
// The accuracy can be measured in terms of the *absolute error* of a
// ciphertext c compared to the plaintext p that it should encrypt.  The
// absolute error of c is defined to be max{|c[i]-p[i]| : i=1..n-1}. Here, n is
// the number of slots. The absolute error of a freshly encrypted ciphertext
// should be no more than (about) 2^{-precision}, where precision is a
// parameter specified in building a Context object.  The absolute error will
// grow as a homomorphic computation proceeds.
//
// Given a ciphertext c, one can obtain its capacity by invoking c.capacity(),
// and one can obtain a bound on its absolute error by invoking c.errorBound().

#include <helib/helib.h>

using namespace std;
using namespace helib;

int main(int argc, char* argv[])
{
  Context context =
      ContextBuilder<CKKS>().m(32 * 1024).bits(358).precision(20).c(6).build();

  cout << "securityLevel=" << context.securityLevel() << "\n";

  long n = context.getNSlots();

  SecKey secretKey(context);
  secretKey.GenSecKey();
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

  // Let's square c a few times and see what happens

  c *= c;
  cout << "c.capacity=" << c.capacity() << " ";
  cout << "c.errorBound=" << c.errorBound() << "\n";

  c *= c;
  cout << "c.capacity=" << c.capacity() << " ";
  cout << "c.errorBound=" << c.errorBound() << "\n";

  c *= c;
  cout << "c.capacity=" << c.capacity() << " ";
  cout << "c.errorBound=" << c.errorBound() << "\n";

  c *= c;
  cout << "c.capacity=" << c.capacity() << " ";
  cout << "c.errorBound=" << c.errorBound() << "\n";

  c *= c;
  cout << "c.capacity=" << c.capacity() << " ";
  cout << "c.errorBound=" << c.errorBound() << "\n";

  //===========================================================================

  // Let's perform the same computation on the plaintext:
  p *= p;
  p *= p;
  p *= p;
  p *= p;
  p *= p;

  //===========================================================================

  // Let's decrypt and compare:
  PtxtArray pp(context);
  pp.decrypt(c, secretKey);

  double distance = Distance(p, pp);
  cout << "distance=" << distance << "\n";

  //===========================================================================

  // On my machine, I get the following output:
  //
  // c.capacity=328.497 c.errorBound=1.28242e-06
  // c.capacity=289.748 c.errorBound=2.69423e-06
  // c.capacity=252.063 c.errorBound=5.71405e-06
  // c.capacity=213.502 c.errorBound=1.1591e-05
  // c.capacity=176.579 c.errorBound=2.37053e-05
  // c.capacity=139.634 c.errorBound=4.79147e-05
  // distance=1.84256e-05
  //
  // So we see that we start out with capacity about 328 (which is somewhat
  // less than the value 358 of the bits parameter), and an errorBound of
  // 1.28242e-06, which is slightly larger than 2^{-20} = 2^{-precision}.
  // After each squaring, capacity decreases by 37-39, while errorBound
  // increases by about a factor of 2 (i.e., we lose one bit of precision).
  // Finally, when we decrypt, we see the actual error (1.84256e-05) is
  // somewhat smaller than errorBound (4.79147e-05).
  //
  // Note that the values returned by capacity() and errorBound() may vary from
  // one run of the program to another, even if all the parameters and
  // plaintext data are the same.  However, they should not change by much from
  // one run to another.

  return 0;
}
