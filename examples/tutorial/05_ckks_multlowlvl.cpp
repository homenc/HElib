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

// In the CKKS encryption scheme (as well as in BGV), ciphertext multiplication
// is a two-step process.  The operation ctxt1 *= ctxt2 is equivalent to the
// following:
//   ctxt1.multLowLvl(ctxt2);
//   ctxt1.reLinearize();
// The operation ctxt1.multLowLvl(ctxt2) multiplies ctxt1 by ctxt2, but it
// leaves ctxt1 in a non-canonical state. The operation ctxt1.reLinearize()
// puts ctxt1 back into a canonical state.  As it happens,
// ctxt1.multLowLvl(ctxt2) is a very fast operation, while ctxt1.reLinearize()
// is a much slower operation.  In addition, some operations, such as
// ciphertext addition, can be applied directly to ciphertexts in non-canonical
// states, yielding ciphertexts also in a non-canonical state. This behavior
// can sometimes be exploited to achieve significant speedups, as illustrated
// here.

#include <helib/helib.h>

using namespace std;
using namespace helib;

int main(int argc, char* argv[])
{
  Context context =
      ContextBuilder<CKKS>().m(16 * 1024).bits(119).precision(20).c(2).build();

  cout << "securityLevel=" << context.securityLevel() << "\n";

  long n = context.getNSlots();

  SecKey secretKey(context);
  secretKey.GenSecKey();

  const PubKey& publicKey = secretKey;

  //===========================================================================

  // Let's encrypt a bunch of random ciphertexts

  int len = 3;

  vector<PtxtArray> p, q;
  for (int i = 0; i < len; i++) {
    p.emplace_back(context);
    p[i].random();
    q.emplace_back(context);
    q[i].random();
  }

  // p[i] is a random PtxtArray for i = 0..len-1
  // q[i] is a random PtxtArray for i = 0..len-1

  vector<Ctxt> c, d;
  for (int i = 0; i < len; i++) {
    c.emplace_back(publicKey);
    p[i].encrypt(c[i]);
    d.emplace_back(publicKey);
    q[i].encrypt(d[i]);
  }

  // c[i] encrypts p[i] for i = 0..len-1
  // d[i] encrypts q[i] for i = 0..len-1

  //===========================================================================

  // Now let's compute the inner product e = sum_{i=0}^{len-1} c[i]*d[i] using
  // multLowLvl. NOTE: this is for illustration purposes only, as HElib already
  // provides a function innerProduct that does the same thing in essentially
  // the same way.

  Ctxt e(publicKey);
  // We use the fact that a freshly constructed ciphertext acts like an
  // encryption of 0

  for (int i = 0; i < len; i++) {
    Ctxt tmp = c[i];
    tmp.multLowLvl(d[i]);
    // tmp is now c[i]*d[i] but in a non-canonical state

    e += tmp;
    // e is now c[0]*d[0] + ... c[i]*d[i] but in a non-canonical state
  }

  e.reLinearize();
  // This puts e back into a canonical state. In this example, we do not really
  // have to do this, but if e were to be used in other computations, it would
  // likely be more efficient to put e into a canonical state once and for all
  // at this point.  The point is, if we had written the above loop with
  // tmp *= d[i] instead of tmp.multLowLvl(d[i]), we would have preformed len
  // expensive reLinearize operations, instead of just one.

  //===========================================================================

  // Let's do the same computation on plaintexts to check the results.

  PtxtArray r(context);
  // We use the fact that a freshly constructed plaintext is 0

  for (int i = 0; i < len; i++) {
    PtxtArray tmp = p[i];
    tmp *= q[i];
    r += tmp;
  }

  // Let's decrypt and compare:
  PtxtArray rr(context);
  rr.decrypt(e, secretKey);

  double distance = Distance(r, rr);
  cout << "distance=" << distance << "\n";

  return 0;
}
