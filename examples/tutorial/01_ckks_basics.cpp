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

// In the CKKS encryption scheme, plaintexts are vectors of real or complex
// numbers.  The length, n, of these vectors is determined by the choice of
// parameters, as discussed below.  We often refer to the components of these
// vectors as "slots", which are indexed 0, ..., n-1.  We can add, subtract, or
// multiply two ciphertexts, and the corresponding operations are carried out
// slot by slot.  This is sometimes referred to as a "SIMD operation" (SIMD
// means Single Instruction Multiple Data), since we can effectively perform
// the same scalar operation in parallel on all n slots.

#include <helib/helib.h>

using namespace std;
using namespace helib;

int main(int argc, char* argv[])
{
  // To get started, we need to choose some parameters.  This is done by
  // initializing a Context object.  Since there are a lot of parameters, many
  // of them optional, HElib provides a "builder pattern" than lets you provide
  // these parameters "by name".

  Context context =

      // initialize a Context object using the builder pattern
      ContextBuilder<CKKS>()

          .m(16 * 1024)
          // m is the "cyclotomic index". For CKKS, m must be a power of 2.  As
          // m increases, you get more security and more slots, but the
          // performance degrades and the size of a ciphertext increases. See
          // table below for more information.

          .bits(119)
          // bits specifies the number of bits in the "ciphertext modulus".  As
          // bits increases, you get less security, but you can perform deeper
          // homomorphic computations; in addition, the size of a ciphertext
          // increases.  See table below for more information. Also see
          // 02_depth.cpp for more information about how depth and bits are
          // related.

          .precision(20)
          // precision specifies the number of bits of precision when data is
          // encoded, encrypted, or decrypted.  More precisely, each of these
          // operations are designed to add an error term of at most
          // 2^{-precision} to each slot.  As precision increases, the allowed
          // depth of homomorphic computations decreases (but security and
          // performance are not affected).  It is not recommended to use
          // precision greater than about 40 or so.

          .c(2)
          // c specifies the number of columns in key-switching matrices.  Yes,
          // it sounds very technical, and it is.  However, all you have to know
          // about this parameter is that as c increases, you get a little more
          // security, but performance degrades and the memory requirement for
          // the public key increases. c must be at least 2 and it is not
          // recommended to set c higher than 8.  See table below for more
          // information.

          .build();
  // last step of the builder pattern

  // The following table lists settings of m, bits, and c that yield (at least)
  // 128-bit security.  It is highly recommended to only use settings from this
  // table.
  //
  //	m	bits	c
  //	16384	119	2
  //	32768	358	6
  //	32768	299	3
  //	32768	239	2
  //	65536	725	8
  //	65536	717	6
  //	65536	669	4
  //	65536	613	3
  //	65536	558	2
  //	131072	1445	8
  //	131072	1435	6
  //	131072	1387	5
  //	131072	1329	4
  //	131072	1255	3
  //	131072	1098	2
  //	262144	2940	8
  //	262144	2870	6
  //	262144	2763	5
  //	262144	2646	4
  //	262144	2511	3
  //	262144	2234	2

  // We can print out the estimated security level.
  // This estimate is based on the LWE security estimator.
  cout << "securityLevel=" << context.securityLevel() << "\n";

  // Get the number of slots, n.  Note that for CKKS, we always have n=m/4.
  long n = context.getNSlots();

  // Construct a secret key. A secret key must be associated with a specific
  // Context, which is passed (by reference) to the constructor.  Programming
  // note: to avoid dangling pointers, the given Context object must not be
  // destroyed while any objects associated with it are still in use.
  SecKey secretKey(context);

  // Constructing a secret key does not actually do very much.  To actually
  // build a full-fledged secret key, we have to invoke the GenSecKey method.
  secretKey.GenSecKey();

  // In HElib, the SecKey class is actually a subclass if the PubKey class.  So
  // one way to initialize a public key object is like this:
  const PubKey& publicKey = secretKey;

  // TECHNICAL NOTE: Note the "&" in the declaration of publicKey. Since the
  // SecKey class is a subclass of PubKey, this particular PubKey object is
  // ultimately a SecKey object, and through the magic of C++ polymorphism,
  // encryptions done via publicKey will actually use the secret key, which has
  // certain advantages.  If one left out the "&", then encryptions done via
  // publicKey will NOT use the secret key.

  //===========================================================================

  // Let's encrypt something!
  // HElib provides a number of idioms for encrypting and decrypting.  We focus
  // on one particular idiom here.

  // We start by declaring a vector of length n, and we fill it with some
  // arbitrary numbers. Note that PI is defined by HElib.
  vector<double> v0(n);
  for (long i = 0; i < n; i++)
    v0[i] = sin(2.0 * PI * i / n);

  // Next, we load the plaintext vector v0 into a special type of container,
  // called a PtxtArray.  Note that a PtxtArray is associated with a Context
  // object, which is passed (by reference) to the constructor.
  PtxtArray p0(context, v0);

  // Note that many types of vectors can be loaded into a PtxtArray object
  // (including, vectors of int, long, double, or even complex<double>).  Also
  // note that constructing p0 and loading v0 into could have been done in two
  // separate steps:
  //   PtxtArray p0(context); p0.load(v0);

  // Next, we construct a ciphertext c0. A ciphertext is associated with a
  // PubKey object, which is passed (by reference) to the constructor.
  // Programming note: to avoid dangling pointers, the given PubKey object must
  // not be destroyed while any objects associated with it are still in use.
  Ctxt c0(publicKey);

  // Finally, we can encrypt p0 and store it in c0:
  p0.encrypt(c0);
  // Note that since a ciphertext is always associated with a public key, there
  // is no need to pass a public key as a separate parameter to the encryption
  // method.

  //===========================================================================

  // We next create another ciphertext c1, in a slightly different way.
  // First, we construct another PtxtArray p1:
  PtxtArray p1(context);

  // Next, we fill all n slots of p1 with random numbers in the interval [0,1]:
  p1.random();

  // Finally, we encrypt p1 and store it in c1, as above:
  Ctxt c1(publicKey);
  p1.encrypt(c1);

  //===========================================================================

  // We next create a ciphertext c2, in the same as was we did c1:
  PtxtArray p2(context);
  p2.random();
  Ctxt c2(publicKey);
  p2.encrypt(c2);

  //===========================================================================

  // Now we homorphically compute c3 = c0*c1 + c2*1.5:
  Ctxt c3 = c0;
  c3 *= c1;
  Ctxt c4 = c2;
  c4 *= 1.5;
  c3 += c4;

  // When this is done, if we denote the i-th slot of a ciphertext c by c[i],
  // then we have c3[i] = c0[i]*c1[i] + c2[i]*1.5 for i = 0..n-1.

  // More generally, for a Ctxt c, one can perform c *= d, c += d, or c -= d,
  // where d can be (among other things) a long, a double, or even a PtxtArray.

  //===========================================================================

  // Next we decrypt c3.
  // First, we construct a new PtxtArray pp3.
  PtxtArray pp3(context);

  // Next, we decrypt c3, storing the plaintext in p3:
  pp3.decrypt(c3, secretKey);

  // Finally, we store the PtxtArray p3 into a standard vector v3:
  vector<double> v3;
  pp3.store(v3);

  //===========================================================================

  // If we like, we can test the accuracy of the computation.
  // First, we perform the same computation directly on plaintexts.
  // The PtxtArray class allows this to be done very easily:
  PtxtArray p3 = p0;
  p3 *= p1;
  PtxtArray p4 = p2;
  p4 *= 1.5;
  p3 += p4;

  // Then, we compute the distance between p3 (computed on plaintexts) and pp3
  // (computed homomorphically on ciphertexts). This is computed as
  // max{ |p3[i]-pp3[i]| : i = 0..n-1 }
  double distance = Distance(p3, pp3);

  cout << "distance=" << distance << "\n";

  return 0;
}
