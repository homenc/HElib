/* Copyright (C) 2012-2020 IBM Corp.
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

#include <NTL/ZZ.h>

#include <helib/helib.h>
#include <helib/permutations.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

struct Parameters
{

  long test;
  long m;
  long p;
  long r;
  long depth;
  long L;
  long ord1;
  long ord2;
  long ord3;
  long ord4;
  long good1;
  long good2;
  long good3;
  long good4;

  Parameters(long test,
             long m,
             long p,
             long r,
             long depth,
             long L,
             long ord1,
             long ord2,
             long ord3,
             long ord4,
             long good1,
             long good2,
             long good3,
             long good4) :
      test(test),
      m(m),
      p(p),
      r(r),
      depth(depth),
      L(L),
      ord1(ord1),
      ord2(ord2),
      ord3(ord3),
      ord4(ord4),
      good1(good1),
      good2(good2),
      good3(good3),
      good4(good4){};

  // Let googletest know how to print the Parameters
  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "test=" << params.test << ","
              << "p=" << params.p << ","
              << "r=" << params.r << ","
              << "m=" << params.m << ","
              << "depth=" << params.depth << ","
              << "L=" << params.L << ","
              << "ord1=" << params.ord1 << ","
              << "ord2=" << params.ord2 << ","
              << "ord3=" << params.ord3 << ","
              << "ord4=" << params.ord4 << ","
              << "good1=" << params.good1 << ","
              << "good2=" << params.good2 << ","
              << "good3=" << params.good3 << ","
              << "good4=" << params.good4 << "}";
  };
};

class GTestPermutations : public ::testing::TestWithParam<Parameters>
{
protected:
  GTestPermutations() :
      test(GetParam().test),
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      depth(GetParam().depth),
      L(GetParam().L),
      ord1(GetParam().ord1),
      ord2(GetParam().ord2),
      ord3(GetParam().ord3),
      ord4(GetParam().ord4),
      good1(GetParam().good1),
      good2(GetParam().good2),
      good3(GetParam().good3),
      good4(GetParam().good4){};

  long test;
  long m;
  long p;
  long r;
  long depth;
  long L;
  long ord1;
  long ord2;
  long ord3;
  long ord4;
  long good1;
  long good2;
  long good3;
  long good4;

  virtual void SetUp() override { helib::setDryRun(helib_test::dry); };

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }
};

void testCube(NTL::Vec<helib::GenDescriptor>& vec, long widthBound)
{
  helib::GeneratorTrees trees;
  long cost = trees.buildOptimalTrees(vec, widthBound);
  if (!helib_test::noPrint) {
    std::cout << "@TestCube: trees=" << trees << std::endl;
    std::cout << " cost =" << cost << std::endl;
  }
  NTL::Vec<long> dims;
  trees.getCubeDims(dims);
  helib::CubeSignature sig(dims);

  for (long cnt = 0; cnt < 3; cnt++) {
    helib::Permut pi;
    helib::randomPerm(pi, trees.getSize());

    helib::PermNetwork net;
    net.buildNetwork(pi, trees);

    helib::HyperCube<long> cube1(sig), cube2(sig);
    for (long i = 0; i < cube1.getSize(); i++)
      cube1[i] = i;
    helib::HyperCube<long> cube3 = cube1;
    helib::applyPermToVec(cube2.getData(),
                          cube1.getData(),
                          pi); // direct application
    net.applyToCube(cube3);    // applying permutation network

    const auto getErrorMessage = [&cube1, &cube2, &cube3]() {
      std::ostringstream os;
      if (cube1.getSize() < 100 && !helib_test::noPrint)
        os << "in=" << cube1.getData() << std::endl
           << "out1=" << cube2.getData() << ", out2=" << cube3.getData()
           << std::endl
           << std::endl;
      return os;
    };

    ASSERT_EQ(cube2, cube3) << getErrorMessage().str();
  }
}

void testCtxt(long m, long p, long widthBound = 0, long L = 0, long r = 1);

void testCtxt(long m, long p, long widthBound, long L, long r)
{
  if (!helib_test::noPrint)
    std::cout << "@testCtxt(m=" << m << ",p=" << p << ",depth=" << widthBound
              << ",r=" << r << ")";

  helib::Context context(m, p, r);
  helib::EncryptedArray ea(context); // Use G(X)=X for this ea object

  // Some arbitrary initial plaintext array
  std::vector<long> in(ea.size());
  for (long i = 0; i < ea.size(); i++)
    in[i] = i % p;

  // Setup generator-descriptors for the PAlgebra generators
  NTL::Vec<helib::GenDescriptor> vec(NTL::INIT_SIZE, ea.dimension());
  for (long i = 0; i < ea.dimension(); i++)
    vec[i] = helib::GenDescriptor(/*order=*/ea.sizeOfDimension(i),
                                  /*good=*/ea.nativeDimension(i),
                                  /*genIdx=*/i);

  // Some default for the width-bound, if not provided
  if (widthBound <= 0)
    widthBound = 1 + log2((double)ea.size());

  // Get the generator-tree structures and the corresponding hypercube
  helib::GeneratorTrees trees;
  long cost = trees.buildOptimalTrees(vec, widthBound);
  if (!helib_test::noPrint) {
    context.zMStar.printout();
    std::cout << ": trees=" << trees << std::endl;
    std::cout << " cost =" << cost << std::endl;
  }
  //  NTL::Vec<long> dims;
  //  trees.getCubeDims(dims);
  //  helib::CubeSignature sig(dims);

  // 1/2 prime per level should be more or less enough, here we use 1 per layer
  if (L <= 0)
    L = (1 + trees.numLayers()) * context.BPL();
  helib::buildModChain(context, /*nLevels=*/L, /*nDigits=*/3);
  if (!helib_test::noPrint)
    std::cout << "**Using " << L << " and " << context.ctxtPrimes.card()
              << " Ctxt-primes)\n";

  // Generate a sk/pk pair
  helib::SecKey secretKey(context);
  const helib::PubKey& publicKey = secretKey;
  secretKey.GenSecKey(); // A +-1/0 secret key
  helib::Ctxt ctxt(publicKey);

  for (long cnt = 0; cnt < 3; cnt++) {
    helib::resetAllTimers();
    // Choose a random permutation
    helib::Permut pi;
    helib::randomPerm(pi, trees.getSize());

    // Build a permutation network for pi
    helib::PermNetwork net;
    net.buildNetwork(pi, trees);

    // make sure we have the key-switching matrices needed for this network
    helib::addMatrices4Network(secretKey, net);

    // Apply the permutation pi to the plaintext
    std::vector<long> out1(ea.size());
    std::vector<long> out2(ea.size());
    helib::applyPermToVec(out1, in, pi); // direct application

    // Encrypt plaintext array, then apply permutation network to ciphertext
    ea.encrypt(ctxt, publicKey, in);
    if (!helib_test::noPrint)
      std::cout << "  ** applying permutation network to ciphertext... "
                << std::flush;
    double t = NTL::GetTime();
    net.applyToCtxt(ctxt, ea); // applying permutation network
    t = NTL::GetTime() - t;
    if (!helib_test::noPrint)
      std::cout << "done in " << t << " seconds" << std::endl;
    ea.decrypt(ctxt, secretKey, out2);

    ASSERT_EQ(out1, out2);
    // printAllTimers();
  }
}

/* m = 31, p = 2, phi(m) = 30
  ord(p)=5
  generator 6 has order (== Z_m^*) of 6
  T = [1 6 5 30 25 26 ]

  m = 61, p = 3, phi(m) = 60
  ord(p)=10
  generator 13 has order (== Z_m^*) of 3
  generator 2 has order (!= Z_m^*) of 2
  T = [1 2 13 26 47 33 ]

  m = 683, p = 2, phi(m) = 682
  ord(p)=22
  generator 3 has order (== Z_m^*) of 31

  m = 47127, p = 2, phi(m) = 30008
  ord(p)=22
  generator 5 has order (== Z_m^*) of 682
  generator 13661 has order (== Z_m^*) of 2
*/

TEST_P(GTestPermutations, ciphertextPermutations)
{
  if (test == 0 || helib_test::dry != 0) {
    NTL::Vec<helib::GenDescriptor> vec;
    long nGens;
    if (ord2 <= 1)
      nGens = 1;
    else if (ord3 <= 1)
      nGens = 2;
    else if (ord4 <= 1)
      nGens = 3;
    else
      nGens = 4;
    vec.SetLength(nGens);

    switch (nGens) {
    case 4:
      vec[3] = helib::GenDescriptor(ord4, good4, /*genIdx=*/3);
      // FALLTHROUGH
    case 3:
      vec[2] = helib::GenDescriptor(ord3, good3, /*genIdx=*/2);
      // FALLTHROUGH
    case 2:
      vec[1] = helib::GenDescriptor(ord2, good2, /*genIdx=*/1);
      // FALLTHROUGH
    default:
      vec[0] = helib::GenDescriptor(ord1, good1, /*genIdx=*/0);
    }
    if (!helib_test::noPrint) {
      std::cout << "***Testing ";
      if (helib::isDryRun())
        std::cout << "(dry run) ";
      for (long i = 0; i < vec.length(); i++)
        std::cout << "(" << vec[i].order << "," << vec[i].good << ")";
      std::cout << ", depth=" << depth << "\n";
    }
    ASSERT_NO_FATAL_FAILURE(testCube(vec, depth));
  } else {
    helib::setTimersOn();
    if (!helib_test::noPrint) {
      std::cout << "***Testing m=" << m << ", p=" << p << ", depth=" << depth
                << std::endl;
    }
    ASSERT_NO_FATAL_FAILURE(testCtxt(m, p, depth, L, r));
  }
}

INSTANTIATE_TEST_SUITE_P(
    defaultParameters,
    GTestPermutations,
    ::testing::Values(
        // FAST
        // Parameters(1, 91, 2, 1, 5, 0, 30, 0, 0, 0, 1, 1, 1, 1)
        // SLOW
        Parameters(1, 4369, 2, 1, 5, 0, 30, 0, 0, 0, 1, 1, 1, 1)));

} // namespace

#if 0
  std::cout << "***Testing m=31, p=2, width=3\n"; // (6 good)
  testCtxt(/*m=*/31, /*p=*/2, /*width=*/3);

  std::cout << "\n***Testing m=61, p=3, width=3\n"; // (3 good), (2, bad)
  testCtxt(/*m=*/61, /*p=*/3, /*width=*/3);

  std::cout << "\n***Testing m=683, p=2, width=5\n"; // (31, good)
  testCtxt(/*m=*/683, /*p=*/2, /*width=*/5);

  //  std::cout << "\n***Testing m=47127, p=2, width=11\n"; // (682,good),(2,good)
  //  testCtxt(/*m=*/47127, /*p=*/2, /*width=*/11);

  // Test 1: a single good small prime-order generator (3)
  {
  NTL::Vec<helib::GenDescriptor> vec(INIT_SIZE, 1);
  vec[0] = helib::GenDescriptor(/*order=*/3, /*good=*/true, /*genIdx=*/0);
  std::cout << "***Testing (3,good), width=1\n";
  testCube(vec, /*width=*/1);
  }

  // Test 2: a single bad larger prime-order generator (31)
  {
  NTL::Vec<helib::GenDescriptor> vec(INIT_SIZE, 1);
  vec[0] = helib::GenDescriptor(/*order=*/31, /*good=*/false, /*genIdx=*/0);
  std::cout << "\n***Testing (31,bad), width=5\n";
  testCube(vec, /*width=*/5);
  }

  // Test 3: two generators with small prime orders (2,3), both bad
  {
  NTL::Vec<helib::GenDescriptor> vec(INIT_SIZE, 2);
  vec[0] = helib::GenDescriptor(/*order=*/2, /*good=*/false, /*genIdx=*/0);
  vec[1] = helib::GenDescriptor(/*order=*/3, /*good=*/false, /*genIdx=*/1);
  std::cout << "\n***Testing [(2,bad),(3,bad)], width=3\n";
  testCube(vec, /*width=*/3);
  }

  // Test 4: two generators with small prime orders (2,3), one good
  {
  NTL::Vec<helib::GenDescriptor> vec(INIT_SIZE, 2);
  vec[0] = helib::GenDescriptor(/*order=*/3, /*good=*/true, /*genIdx=*/0);
  vec[1] = helib::GenDescriptor(/*order=*/2, /*good=*/false, /*genIdx=*/1);
  std::cout << "\n***Testing [(3,good),(2,bad)], width=3\n";
  testCube(vec, /*width=*/3);
  }

  // Test 5: a single good composite-order generator (6)
  {
  NTL::Vec<helib::GenDescriptor> vec(INIT_SIZE, 1);
  vec[0] = helib::GenDescriptor(/*order=*/6, /*good=*/true, /*genIdx=*/0);
  std::cout << "\n***Testing (6,good), width=3\n";
  testCube(vec, /*width=*/3);
  }

  // Test 6: (6,good),(2,bad)
  {
  NTL::Vec<helib::GenDescriptor> vec(INIT_SIZE, 2);
  vec[0] = helib::GenDescriptor(/*order=*/6,/*good=*/true, /*genIdx=*/0);
  vec[1] = helib::GenDescriptor(/*order=*/ 2, /*good=*/false,/*genIdx=*/1);
  std::cout << "\n**Testing [(6,good),(2,bad)], width=5\n";
  testCube(vec, /*width=*/5);
  }

  // Test 7: the "general case", (682,good),(2,bad)
  {
  NTL::Vec<helib::GenDescriptor> vec(INIT_SIZE, 2);
  vec[0] = helib::GenDescriptor(/*order=*/682,/*good=*/true, /*genIdx=*/0);
  vec[1] = helib::GenDescriptor(/*order=*/ 2, /*good=*/false,/*genIdx=*/1);
  std::cout << "\n**Testing [(682,good),(2,bad)], width=11\n";
  testCube(vec, /*width=*/11);
  }
#endif
