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
// testPacking.cxx - testing unpack/repack functionality
#include <helib/intraSlot.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

struct Parameters
{
  long p;    // plaintext base
  long n;    // number of packed ciphertexts
  long r;    // lifting
  long L;    // # of levels in the modulus chain
  long m;    // use specified value as modulus
  long seed; // PRG seed

  Parameters(long p, long n, long r, long L, long m, long seed) :
      p(p), n(n), r(r), L(L), m(m), seed(seed){};

  // Let googletest know how to print the Parameters
  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "p=" << params.p << ","
              << "n=" << params.n << ","
              << "r=" << params.r << ","
              << "L=" << params.L << ","
              << "m=" << params.m << ","
              << "seed=" << params.seed << "}";
  };
};

class GTestIntraSlot : public ::testing::TestWithParam<Parameters>
{

  static helib::Context& setupContext(helib::Context& context, long L)
  {
    if (helib_test::verbose) {
      context.zMStar.printout();
    }
    buildModChain(context, L, 3);
    return context;
  };

protected:
  GTestIntraSlot() :
      p(GetParam().p),
      n(GetParam().n),
      r(GetParam().r),
      L(GetParam().L),
      m(GetParam().m),
      seed(GetParam().seed),
      context(m, p, r),
      secretKey(setupContext(context, L)),
      publicKey(secretKey){};

  long p;
  long n;
  long r;
  long L;
  long m;
  long seed;
  helib::Context context;
  helib::SecKey secretKey;
  const helib::PubKey& publicKey;

  void SetUp() override
  {
    SetSeed(NTL::ZZ(seed));
    secretKey.GenSecKey(); // A +-1/0 secret key
    helib::addSome1DMatrices(
        secretKey); // compute key-switching matrices that we need
    helib::addFrbMatrices(secretKey);

    helib::setupDebugGlobals(&secretKey, context.ea);
  };

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }
};

TEST_P(GTestIntraSlot, packingAndUnpackingWorks)
{
  NTL::ZZX G = context.alMod.getFactorsOverZZ()[0];
  helib::EncryptedArray ea(context, G);

  long d = ea.getDegree(); // size of each slot

  std::vector<helib::Ctxt> unpacked(d * n - 1, helib::Ctxt(publicKey));

  // generate (almost) d*n ciphertexts, with only integers in the slots
  std::vector<helib::PlaintextArray> p1(helib::lsize(unpacked),
                                        helib::PlaintextArray(ea));
  for (long i = 0; i < helib::lsize(unpacked); i++) {
    std::vector<long> slots;
    ea.random(slots);
    encode(ea, p1[i], slots);
    ea.encrypt(unpacked[i], publicKey, p1[i]);
  }

  // Pack (almost) d*n ciphertexts into only n of them
  std::vector<helib::Ctxt> ct(n, helib::Ctxt(publicKey));
  repack(helib::CtPtrs_vectorCt(ct), helib::CtPtrs_vectorCt(unpacked), ea);

  // Unpack them back
  std::vector<helib::zzX> unpackSlotEncoding;
  helib::buildUnpackSlotEncoding(unpackSlotEncoding, ea);
  unpack(helib::CtPtrs_vectorCt(unpacked),
         helib::CtPtrs_vectorCt(ct),
         ea,
         unpackSlotEncoding);

  helib::PlaintextArray p2(ea);
  for (long i = 0; i < helib::lsize(unpacked); i++) {
    ea.decrypt(unpacked[i], secretKey, p2);
    ASSERT_TRUE(equals(ea, p1[i], p2)) << "p2[" << i << "]=" << p2;
  }
}
INSTANTIATE_TEST_SUITE_P(someParameters,
                         GTestIntraSlot,
                         ::testing::Values(
                             // FAST
                             Parameters(2, 2, 1, 10, 91, 0)));
} // namespace
