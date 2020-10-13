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

/* TestPermutations.cpp - Applying plaintext permutation to encrypted vector
 */
#include <NTL/ZZ.h>

#include <helib/NumbTh.h>
#include <helib/timing.h>
#include <helib/permutations.h>
#include <helib/EncryptedArray.h>
#include <helib/ArgMap.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

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
  CKKSParameters(long m, long r, long bits, long depth) :
    BGVParameters(m, /*p=*/-1, r, bits, depth){};
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
      context(m, p, r),
      secretKey((buildModChain(context, bits), context)),
      publicKey((secretKey.GenSecKey(), addSome1DMatrices(secretKey), secretKey)),
      ea(*(context.ea))
    {}

  virtual void SetUp() override
  {
    helib::setupDebugGlobals(&secretKey, context.ea);
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
      context(m, /*p=*/-1, r),
      secretKey((buildModChain(context, bits), context)),
      publicKey((secretKey.GenSecKey(), addSome1DMatrices(secretKey), secretKey)),
      ea(context.ea->getCx())
    {}

  virtual void SetUp() override
  {
    helib::setupDebugGlobals(&secretKey, context.ea);
  }

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }
};

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

  // Ciphertext permutation
  pp.apply(ctxt);
  // Plaintext permutation
  pp.apply(v);

  helib::PtxtArray w(context);
  w.decrypt(ctxt, secretKey);

  EXPECT_EQ(w, v);
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

  // Ciphertext permutation
  pp.apply(ctxt);
  // Plaintext permutation
  pp.apply(v);

  helib::PtxtArray w(context);
  w.decrypt(ctxt, secretKey);

  EXPECT_TRUE(w == helib::Approx(v));
}

INSTANTIATE_TEST_SUITE_P(variousParameters,
                         TestPermutationsBGV,
                         ::testing::Values(BGVParameters(/*m=*/4369, /*p=*/2, /*r=*/1, /*bits=*/1000, /*depth=*/5)));

INSTANTIATE_TEST_SUITE_P(variousParameters,
                         TestPermutationsCKKS,
                         ::testing::Values(CKKSParameters(/*m=*/8192, /*r=*/20, /*bits=*/1000, /*depth=*/5)));

} // namespace
