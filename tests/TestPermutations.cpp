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

/* TestPermutations.cpp - Applying plaintext permutation to encrypted vector
 */
#include <NTL/ZZ.h>

#include <helib/NumbTh.h>
#include <helib/timing.h>
#include <helib/permutations.h>
#include <helib/EncryptedArray.h>
#include <helib/ArgMap.h>
#include <helib/debugging.h>
#include "../src/macro.h"

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

struct Parameters
{
  Parameters(std::vector<long> orders, std::vector<long> good, long depth) :
      depth(depth)
  {
    gens.SetLength(orders.size());
    for (std::size_t i = 0; i < orders.size(); ++i) {
      gens[i] = helib::GenDescriptor(orders[i], good[i], /*genIdx=*/i);
    }
  }

  NTL::Vec<helib::GenDescriptor> gens;
  const long depth;

  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    os << "{gens = ";
    for (long i = 0; i < params.gens.length(); ++i) {
      os << "(" << params.gens[i].order << ", " << params.gens[i].good << "), ";
    }
    return os << "depth = " << params.depth << "}";
  }
};

struct BGVParameters
{
  BGVParameters(long m, long p, long r, long bits, long depth) :
      m(m), p(p), r(r), bits(bits), depth(depth){};

  const long m;
  const long p;
  const long r;
  const long bits;
  const long depth;

  friend std::ostream& operator<<(std::ostream& os, const BGVParameters& params)
  {
    return os << "{"
              << "m = " << params.m << ", "
              << "p = " << params.p << ", "
              << "r = " << params.r << ", "
              << "bits = " << params.bits << ", "
              << "depth = " << params.depth << "}";
  }
};

struct CKKSParameters : public BGVParameters
{
  CKKSParameters(long m, long precision, long bits, long depth) :
      BGVParameters(m, /*p=*/-1, /*r=*/precision, bits, depth){};
};

class TestPermutationsGeneral : public ::testing::TestWithParam<Parameters>
{
protected:
  const NTL::Vec<helib::GenDescriptor> gens;
  const long depth;

  TestPermutationsGeneral() : gens(GetParam().gens), depth(GetParam().depth) {}
};

class TestPermutationsBGV : public ::testing::TestWithParam<BGVParameters>
{
protected:
  const long m;
  const long p;
  const long r;
  const long bits;
  const long depth;

  helib::Context context;
  helib::SecKey secretKey;
  helib::PubKey publicKey;
  const helib::EncryptedArray& ea;

  TestPermutationsBGV() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      bits(GetParam().bits),
      depth(GetParam().depth),
      context(helib::ContextBuilder<helib::BGV>()
                  .m(m)
                  .p(p)
                  .r(r)
                  .bits(bits)
                  .build()),
      secretKey(context),
      publicKey(
          (secretKey.GenSecKey(), addSome1DMatrices(secretKey), secretKey)),
      ea(context.getEA())
  {}

  virtual void SetUp() override
  {
    helib::setupDebugGlobals(&secretKey, context.shareEA());
  }

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }
};

class TestPermutationsCKKS : public ::testing::TestWithParam<CKKSParameters>
{
protected:
  const long m;
  const long r;
  const long bits;
  const long depth;

  helib::Context context;
  helib::SecKey secretKey;
  helib::PubKey publicKey;
  const helib::EncryptedArrayCx& ea;

  TestPermutationsCKKS() :
      m(GetParam().m),
      r(GetParam().r),
      bits(GetParam().bits),
      depth(GetParam().depth),
      context(helib::ContextBuilder<helib::CKKS>()
                  .m(m)
                  .precision(r)
                  .bits(bits)
                  .build()),
      secretKey(context),
      publicKey(
          (secretKey.GenSecKey(), addSome1DMatrices(secretKey), secretKey)),
      ea(context.getEA().getCx())
  {}

  virtual void SetUp() override
  {
    helib::setupDebugGlobals(&secretKey, context.shareEA());
  }

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }
};

TEST_P(TestPermutationsBGV, ciphertextPermutations)
{
  std::vector<long> in(ea.size());
  // Arbitrary initial data
  for (std::size_t i = 0; i < in.size(); ++i) {
    in[i] = i % p;
  }

  // Set up generator-descriptors for the PAlgebra generators
  NTL::Vec<helib::GenDescriptor> gens;
  gens.SetLength(ea.dimension());
  for (long i = 0; i < ea.dimension(); ++i) {
    gens[i] = helib::GenDescriptor(/*order=*/ea.sizeOfDimension(i),
                                   /*good=*/ea.nativeDimension(i),
                                   /*genIdx=*/i);
  }

  // Get the gnerator-tree structures and the corresponding hypercube
  helib::GeneratorTrees trees;
  long cost = trees.buildOptimalTrees(gens, depth);

  if (helib_test::verbose) {
    std::cout << "trees = " << trees << "\n"
              << "cost = " << cost << std::endl;
  }

  for (long cnt = 0; cnt < 3; ++cnt) {
    // Choose a random permutation
    helib::Permut pi;
    helib::randomPerm(pi, trees.getSize());

    // Build a permutation network for pi
    helib::PermNetwork net;
    net.buildNetwork(pi, trees);

    // Make sure we have the key-switching matrices needed for this network
    addMatrices4Network(secretKey, net);

    // Apply the permutation pi to the plaintext
    std::vector<long> out1(ea.size());
    std::vector<long> out2(ea.size());
    helib::applyPermToVec(out1, in, pi); // Direct application

    // Encrypt plaintext array, then apply permutation network to the ciphertext
    helib::Ctxt ctxt(publicKey);
    ea.encrypt(ctxt, publicKey, in);
    if (helib_test::verbose) {
      std::cout << "  ** Applying permutation network to ciphertext... "
                << std::flush;
    }

    double t = NTL::GetTime();
    net.applyToCtxt(ctxt, ea); // Applying permutation network
    t = NTL::GetTime() - t;

    if (helib_test::verbose) {
      std::cout << "done in " << t << " seconds" << std::endl;
    }

    ea.decrypt(ctxt, secretKey, out2);

    EXPECT_EQ(out1, out2);
  }
}

TEST_P(TestPermutationsBGV, ciphertextPermutationsWithNewAPI)
{
  helib::PermIndepPrecomp pip(context, depth);

  if (helib_test::verbose) {
    std::cout << "depth = " << pip.getDepth() << "\n"
              << "cost = " << pip.getCost() << "\n";
  }

  helib::Permut pi;
  helib::randomPerm(pi, context.getNSlots());

  helib::PermPrecomp pp(pip, pi);
  helib::Ctxt ctxt(publicKey);
  helib::PtxtArray v(context);
  v.random();
  v.encrypt(ctxt);

  pp.apply(ctxt); // Ciphertext permutation
  pp.apply(v);    // Plaintext permutation

  helib::PtxtArray w(context);
  w.decrypt(ctxt, secretKey);

  EXPECT_EQ(w, v);
}

// This test is in TestPermutations for now as this is where
// this issue was discovered.
TEST(TestPermutationsCKKS, ckksFailIfRBitsTooLarge)
{
  EXPECT_THROW(
      helib::ContextBuilder<helib::CKKS>().precision(HELIB_SP_NBITS).build(),
      helib::InvalidArgument);
}

TEST_P(TestPermutationsCKKS, ciphertextPermutationsWithNewAPI)
{
  helib::PermIndepPrecomp pip(context, depth);

  if (helib_test::verbose) {
    std::cout << "depth = " << pip.getDepth() << "\n"
              << "cost = " << pip.getCost() << "\n";
  }

  helib::Permut pi;
  helib::randomPerm(pi, context.getNSlots());

  helib::PermPrecomp pp(pip, pi);
  helib::Ctxt ctxt(publicKey);
  helib::PtxtArray v(context);
  v.random();
  v.encrypt(ctxt, context.getNSlots()); // CKKS encryption

  pp.apply(ctxt); // Ciphertext permutation
  pp.apply(v);    // Plaintext permutation

  helib::PtxtArray w(context);
  w.decrypt(ctxt, secretKey);

  // TODO-FB investigate the use of EXPECT_NEAR with an error threshold
  EXPECT_TRUE(w == helib::Approx(v));
}

TEST_P(TestPermutationsGeneral, testCube)
{
  helib::GeneratorTrees trees;
  long cost = trees.buildOptimalTrees(gens, depth);

  if (helib_test::verbose) {
    std::cout << "trees = " << trees << "\n"
              << "cost  = " << cost << std::endl;
  }

  NTL::Vec<long> dims;
  trees.getCubeDims(dims);
  helib::CubeSignature sig(dims);

  for (long count = 0; count < 3; ++count) {
    helib::Permut pi;
    helib::randomPerm(pi, trees.getSize());

    helib::PermNetwork net;
    net.buildNetwork(pi, trees);

    helib::HyperCube<long> cube1(sig), cube2(sig);
    for (long i = 0; i < cube1.getSize(); ++i) {
      cube1[i] = i;
    }
    helib::HyperCube<long> cube3 = cube1;
    // Direct application
    helib::applyPermToVec(cube2.getData(), cube1.getData(), pi);
    // Applying permutation network
    net.applyToCube(cube3);

    EXPECT_EQ(cube2, cube3)
        << "input = " << cube1.getData() << "\noutput1 = " << cube2.getData()
        << "\noutput2 = " << cube2.getData();
  }
}

INSTANTIATE_TEST_SUITE_P(variousParameters,
                         TestPermutationsBGV,
                         ::testing::Values(BGVParameters(/*m=*/4369,
                                                         /*p=*/2,
                                                         /*r=*/1,
                                                         /*bits=*/1000,
                                                         /*depth=*/5)));

INSTANTIATE_TEST_SUITE_P(
    variousParameters,
    TestPermutationsCKKS,
    ::testing::Values(CKKSParameters(/*m=*/8192,
                                     /*precision=*/HELIB_SP_NBITS > 52
                                         ? 52
                                         : HELIB_SP_NBITS - 1,
                                     /*bits=*/1000,
                                     /*depth=*/5)));

INSTANTIATE_TEST_SUITE_P(
    variousParameters,
    TestPermutationsGeneral,
    ::testing::Values(Parameters(/*orders=*/{3}, /*good=*/{1}, /*depth=*/1),
                      Parameters({31}, {0}, 5),
                      Parameters({2, 3}, {0, 0}, 3),
                      Parameters({2, 3}, {1, 0}, 3),
                      Parameters({6}, {1}, 3),
                      Parameters({6, 2}, {1, 0}, 5),
                      Parameters({682, 2}, {1, 0}, 11)));

} // namespace
