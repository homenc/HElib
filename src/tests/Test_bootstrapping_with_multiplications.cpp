/* Copyright (C) 2012-2019 IBM Corp.
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
#include <helib.h>
#include "debugging.h"

#include "gtest/gtest.h"
#include "test_common.h"
#include <random>

namespace {
struct Parameters {
  Parameters(long m, long p, long r, long c, long bits,
             long t, int c_m, long n, std::vector<long> mvec,
             std::vector<long> gens, std::vector<long> ords) :
    m(m),
    p(p),
    r(r),
    c(c),
    bits(bits),
    t(t),
    c_m(c_m),
    n(n),
    mvec(convert<NTL::Vec<long>>(mvec)),
    gens(gens),
    ords(ords) {
  if (mvec.empty(), gens.empty() || ords.empty())
      throw helib::LogicError("mvec, gens, and ords must be non-empty");
    };

    const long m;
    const long p;
    const long r;
    const long c;
    const long bits;
    const long t;
    const int c_m;
    const long n;
    const NTL::Vec<long> mvec;
    const std::vector<long> gens;
    const std::vector<long> ords;

  friend std::ostream &operator<<(std::ostream &os, const Parameters &params) {
    return os << "{"
              << "m=" << params.m << ","
              << "p=" << params.p << ","
              << "r=" << params.r << ","
              << "c=" << params.c << ","
              << "bits=" << params.bits << ","
              << "skHwt=" << params.t << ","
              << "gens=" << params.gens << ","
              << "ords=" << params.ords << ","
              << "mvec=" << params.mvec << ","
              << "c_m=" << params.c_m << ","
              << "computation depth=" << params.n
              << "}";
  }

};

class Test_bootstrapping_with_multiplications : public ::testing::TestWithParam<Parameters> {
  protected:
    const long n; // Number of multiplications to perform

    FHEcontext context;
    FHESecKey secretKey;
    FHEPubKey publicKey;
    const EncryptedArray& ea;

    Test_bootstrapping_with_multiplications () :
      n(GetParam().n),
      context(GetParam().m, GetParam().p, GetParam().r,
              GetParam().gens, GetParam().ords),
      secretKey(postContextSetup(context, GetParam().c_m, GetParam().bits, GetParam().c,
                                 GetParam().t, GetParam().mvec)),
      publicKey(keySetup(secretKey)),
      ea(*(context.ea))
    {}

    static FHEcontext& postContextSetup(FHEcontext& context, int c_m, long bits, long c, long t, NTL::Vec<long> mvec) {
      context.zMStar.set_cM(c_m / 100);
      buildModChain(context, bits, c, true, t);
      context.makeBootstrappable(mvec, t, 0);
      return context;
    }

    static FHESecKey& keySetup(FHESecKey& secretKey) {
      secretKey.GenSecKey();
      addSome1DMatrices(secretKey);
      addFrbMatrices(secretKey);
      secretKey.genRecryptData();
      return secretKey;
    }

    virtual void SetUp() override
    {
      if(helib_test::verbose) {
        std::cout << "m=" << GetParam().m
                  << ", p=" << GetParam().p
                  << ", r=" << GetParam().r
                  << ", bits=" << GetParam().bits
                  << ", c=" << GetParam().c
                  << ", skHwt=" << GetParam().t
                  << ", c_m=" << GetParam().c_m 
                  << ", depth to compute=" << GetParam().n
                  << ", mvec=" << GetParam().mvec
                  << ", gens=" << GetParam().gens
                  << ", ords=" << GetParam().ords
                  << std::endl;
        ea.getPAlgebra().printout();
        std::cout << "ctxtPrimes="<<context.ctxtPrimes
                  << ", specialPrimes="<<context.specialPrimes<<std::endl<<std::endl;
      }
    }

    virtual void TearDown() override
    {
      if(helib_test::verbose) {
	printAllTimers();
      }
      cleanupGlobals();
    }
};

std::vector<long> generateRandomBinaryVector(long nslots)
{
  std::vector<long> ptxt(nslots);
  std::mt19937 gen(helib_test::random_seed);
  std::uniform_int_distribution<int> coinFlipDist(0,1);
  for(auto& num : ptxt)
    num = coinFlipDist(gen);
  return ptxt;
}

TEST_P(Test_bootstrapping_with_multiplications, correctly_performs_bootstrapping_with_no_multiplications) {
  const long nslots = ea.size();
  std::vector<long> ptxt(generateRandomBinaryVector(nslots)); // Random 0s and 1s

  std::vector<long> tmp_ptxt(ptxt);
  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, ptxt);
  if(!helib_test::noPrint) {
    CheckCtxt(ctxt, "Before recryption");
  }
  // Time the recryption step
  FHE_NTIMER_START(Bootstrap);
  // Recrypt/Bootstrap the ctxt
  publicKey.reCrypt(ctxt);
  FHE_NTIMER_STOP(Bootstrap);
  if(!helib_test::noPrint) {
    CheckCtxt(ctxt, "After recryption");
  }
  std::vector<long> decrypted(nslots);
  ea.decrypt(ctxt, secretKey, decrypted);

  EXPECT_EQ(decrypted, ptxt);
}

TEST_P(Test_bootstrapping_with_multiplications, correctly_performs_bootstrapping_with_multiplications) {
  const long nslots = ea.size();
  long p2r = context.alMod.getPPowR();
  std::vector<long> ptxt(generateRandomBinaryVector(nslots)); // Random 0s and 1s

  std::vector<long> tmp_ptxt(ptxt);
  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, ptxt);

  long depth = 0; // count to keep track of number of multiplications
  long round = 0; // count for number of bootstraps
  while (depth < n) {
    if(helib_test::verbose) {
      std::cout << "Round " << round << std::endl;
    }

    Ctxt tmp_ctxt(ctxt);
    if(!helib_test::noPrint) {
      CheckCtxt(ctxt, "Before multiplication");
    }
    FHE_NTIMER_START(RoundTotal);
    FHE_NTIMER_START(Multiplications);
    // Multiply the ciphertext with itself n times
    // until number of bits falls below threshold
    long count = 0; // count for number of multiplications this round
    while (ctxt.bitCapacity() >= 40 && depth + count < n) {
      if(helib_test::verbose) {
	std::cout << "multiplication " << count+1 << std::endl;
      }
      ctxt.multiplyBy(tmp_ctxt);
      count += 1;
    }
    FHE_NTIMER_STOP(Multiplications);
    if(!helib_test::noPrint) {
      CheckCtxt(ctxt, "Before recryption");
    }
    if (depth + count < n) {
      // Time the recryption step
      FHE_NTIMER_START(Bootstrap);
      // Recrypt/Bootstrap the ctxt
      publicKey.reCrypt(ctxt);
      FHE_NTIMER_STOP(Bootstrap);
    }
    FHE_NTIMER_STOP(RoundTotal);
    if(!helib_test::noPrint) {
      CheckCtxt(ctxt, "After recryption");
    }
    // Plaintext operation
    // Multiply with itself n times
    for (int j = 0; j < count; j++) {
      for (int i = 0; i < nslots; i++) {
	ptxt[i] *= tmp_ptxt[i];
	ptxt[i] = ptxt[i] % p2r;
      }
    }
    depth += count; // current depth of circuit computed
    round += 1;     // end of round
  }
  std::vector<long> decrypted(nslots);
  ea.decrypt(ctxt, secretKey, decrypted);

  EXPECT_EQ(decrypted, ptxt);

  if(!helib_test::noPrint) {
    printNamedTimer(std::cout << std::endl, "Setup");
    printNamedTimer(std::cout << std::endl, "Multiplications");
    printNamedTimer(std::cout << std::endl, "Bootstrap");
    printNamedTimer(std::cout << std::endl, "RoundTotal");
  }
}

INSTANTIATE_TEST_SUITE_P(typical_parameters, Test_bootstrapping_with_multiplications, ::testing::Values(
  //FAST
  Parameters(31*41, 2, 1, 2, 580, 64, 100, 50, {31, 41}, {1026, 249}, {30, -2}),
  Parameters(7*5*37, 17, 1, 3, 600, 64, 100, 50, {7, 5, 37}, {556, 1037}, {6, 4})
  //SLOW
  //Parameters(31775, 2, 1, 2, 580, 64, 100, 50, {41, 775}, {6976, 24806}, {40, 30}),
  //Parameters(35113, 2, 1, 2, 580, 64, 100, 50, {37, 949}, {16134, 8548}, {36, 24})
));

} // namespace
