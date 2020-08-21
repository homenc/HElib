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
#include <NTL/BasicThreadPool.h>
#include <helib/helib.h>
#include <NTL/lzz_pXFactoring.h>

#include <cstdio>

#include "gtest/gtest.h"
#include "test_common.h"

#include <helib/debugging.h>
#include <helib/replicate.h>
#include <helib/fhe_stats.h>

/**************

  1. c1.multiplyBy(c0)
  2. c0 += random constant
  3. c2 *= random constant
  4. tmp = c1
  5. ea.shift(tmp, random amount in [-nSlots/2, nSlots/2])
  6. c2 += tmp
  7. ea.rotate(c2, random amount in [1-nSlots, nSlots-1])
  8. c1.negate()
  9. c3.multiplyBy(c2)
  10. c0 -= c3

 **************/

namespace {

::testing::AssertionResult ciphertextMatches(const helib::EncryptedArray& ea,
                                             const helib::SecKey& sk,
                                             const helib::PlaintextArray& p,
                                             const helib::Ctxt& c)
{
  helib::PlaintextArray pp(ea);
  ea.decrypt(c, sk, pp);
  if (equals(ea, pp, p)) {
    return ::testing::AssertionSuccess();
  } else {
    return ::testing::AssertionFailure()
           << "Ciphertext does not match plaintext:" << std::endl
           << "p = " << p << std::endl
           << "pp = " << pp << std::endl;
  }
}

struct Parameters
{
  Parameters(long R,
             long p,
             long r,
             long d,
             long c,
             long k,
             long L,
             long s,
             long m,
             std::vector<long> mvec,
             std::vector<long> gens,
             std::vector<long> ords,
             long seed,
             long nt) :
      R(R),
      p(p),
      r(r),
      d(d),
      c(c),
      k(k),
      L(L),
      s(s),
      m(m),
      mvec(helib::convert<NTL::Vec<long>>(mvec)),
      gens(gens),
      ords(ords),
      seed(seed),
      nt(nt){};

  long R; // number of rounds
  long p; // plaintext base
  long r; // lifting
  long d; // degree of the field extension
  // Note: d == 0 => factors[0] defines extension
  long c;              // number of columns in the key-switching matrices
  long k;              // security parameter
  long L;              // # of bits in the modulus chain
  long s;              // minimum number of slots
  long m;              // use specified value as modulus
  NTL::Vec<long> mvec; // use product of the integers as modulus
  // Note: mvec takes priority over m
  std::vector<long> gens; // use specified vector of generators
  std::vector<long> ords; // use specified vector of orders
  // e.g., ords=[4 2 -4], negative means 'bad'
  long seed; // PRG seed
  long nt;   // num threads

  // Let googletest know how to print the Parameters
  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "R=" << params.R << ","
              << "p=" << params.p << ","
              << "r=" << params.r << ","
              << "d=" << params.d << ","
              << "c=" << params.c << ","
              << "k=" << params.k << ","
              << "L=" << params.L << ","
              << "s=" << params.s << ","
              << "m=" << params.m << ","
              << "mvec=" << params.mvec << ","
              << "gens=" << helib::vecToStr(params.gens) << ","
              << "ords=" << helib::vecToStr(params.ords) << ","
              << "seed=" << params.seed << ","
              << "nt=" << params.nt << "}";
  };
};

class GTestGeneral : public ::testing::TestWithParam<Parameters>
{
protected:
  long R;
  long p;
  long r;
  long d;
  long c;
  long k;
  long w;
  long L;
  long m;

  std::vector<long> gens;
  std::vector<long> ords;

  helib::Context context;

  helib::SecKey secretKey;

  const helib::PubKey& publicKey;

  GTestGeneral() :
      R(GetParam().R),
      p(GetParam().p),
      r(GetParam().r),
      d(GetParam().d),
      c(GetParam().c),
      k(GetParam().k),
      w(64),
      L(GetParam().L),
      m(helib::FindM(k,
                     L,
                     c,
                     p,
                     d,
                     GetParam().s,
                     GetParam().mvec.length() > 0
                         ? helib::computeProd(GetParam().mvec)
                         : GetParam().m,
                     true)),
      gens(GetParam().gens),
      ords(GetParam().ords),
      context(m, p, r, gens, ords),
      secretKey(
          (buildModChain(context,
                         L,
                         c,
                         /*willBeBootstrappable=*/false,
                         /*skHwt=*/0,
                         /*resolution=*/3,
                         /*bitsInSpecialPrimes=*/helib_test::special_bits),
           context)),
      publicKey(secretKey)
  {}

  virtual void SetUp() override
  {
    NTL::SetSeed(NTL::ZZ(GetParam().seed));
    NTL::SetNumThreads(GetParam().nt);
    secretKey.GenSecKey(w); // A Hamming-weight-w secret key
    helib::addSome1DMatrices(
        secretKey); // compute key-switching matrices that we need

    helib::setupDebugGlobals(&secretKey, context.ea);
    if (!helib_test::noPrint) {
      helib::fhe_stats = true;
    }
  };

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }
};

TEST_P(GTestGeneral, correctlyImplementsMixOfOperationsOverFourCiphertexts)
{
  char buffer[32];
  if (!helib_test::noPrint) {
    std::cout << "\n\n******** " << (helib_test::dry ? "(dry run):" : ":");
    std::cout << " R=" << R << ", p=" << p << ", r=" << r << ", d=" << d
              << ", c=" << c << ", k=" << k << ", w=" << w << ", L=" << L
              << ", m=" << m << ", gens=" << helib::vecToStr(gens)
              << ", ords=" << helib::vecToStr(ords) << std::endl;
  }

  NTL::ZZX G;
  if (d == 0)
    G = context.alMod.getFactorsOverZZ()[0];
  else
    G = helib::makeIrredPoly(p, d);

  if (!helib_test::noPrint) {
    context.zMStar.printout();
    std::cout << std::endl;

    std::cout << "security=" << context.securityLevel() << std::endl;
    std::cout << "# small primes = " << context.smallPrimes.card() << "\n";
    std::cout << "# ctxt primes = " << context.ctxtPrimes.card() << "\n";
    std::cout << "# bits in ctxt primes = "
              << long(context.logOfProduct(context.ctxtPrimes) / log(2.0) + 0.5)
              << "\n";
    std::cout << "# special primes = " << context.specialPrimes.card() << "\n";
    std::cout << "# bits in special primes = "
              << long(context.logOfProduct(context.specialPrimes) / log(2.0) +
                      0.5)
              << "\n";
    std::cout << "G = " << G << "\n";
  }
  std::shared_ptr<helib::EncryptedArray> ea_ptr =
      std::make_shared<helib::EncryptedArray>(context, G);
  helib::EncryptedArray& ea(*ea_ptr);
  long nslots = ea.size();

  // Debugging additions
  helib::setupDebugGlobals(&secretKey, ea_ptr);

  helib::PlaintextArray p0(ea);
  helib::PlaintextArray p1(ea);
  helib::PlaintextArray p2(ea);
  helib::PlaintextArray p3(ea);

  helib::random(ea, p0);
  helib::random(ea, p1);
  helib::random(ea, p2);
  helib::random(ea, p3);

  helib::Ctxt c0(publicKey), c1(publicKey), c2(publicKey), c3(publicKey);
  ea.encrypt(c0, publicKey, p0);
  // {ZZX ppp0; ea.encode(ppp0, p0); c0.DummyEncrypt(ppp0);} // dummy encryption
  ea.encrypt(c1, publicKey, p1); // real encryption
  ea.encrypt(c2, publicKey, p2); // real encryption
  ea.encrypt(c3, publicKey, p3); // real encryption

  helib::resetAllTimers();

  HELIB_NTIMER_START(Circuit);

  for (long i = 0; i < R; i++) {

    if (!helib_test::noPrint)
      std::cout << "*** round " << i << "..." << std::endl;

    long shamt = NTL::RandomBnd(2 * (nslots / 2) + 1) - (nslots / 2);
    // random number in [-nslots/2..nslots/2]
    long rotamt = NTL::RandomBnd(2 * nslots - 1) - (nslots - 1);
    // random number in [-(nslots-1)..nslots-1]

    // two random constants
    helib::PlaintextArray const1(ea);
    helib::PlaintextArray const2(ea);
    helib::random(ea, const1);
    helib::random(ea, const2);

    NTL::ZZX const1_poly, const2_poly;
    ea.encode(const1_poly, const1);
    ea.encode(const2_poly, const2);

    mul(ea, p1, p0); // c1.multiplyBy(c0)
    c1.multiplyBy(c0);
    if (!helib_test::noPrint)
      CheckCtxt(c1, "c1*=c0");
    EXPECT_TRUE(ciphertextMatches(ea, secretKey, p1, c1));

    add(ea, p0, const1); // c0 += random constant
    c0.addConstant(const1_poly);
    if (!helib_test::noPrint)
      CheckCtxt(c0, "c0+=k1");
    EXPECT_TRUE(ciphertextMatches(ea, secretKey, p0, c0));

    mul(ea, p2, const2); // c2 *= random constant
    c2.multByConstant(const2_poly);
    if (!helib_test::noPrint)
      CheckCtxt(c2, "c2*=k2");
    EXPECT_TRUE(ciphertextMatches(ea, secretKey, p2, c2));

    helib::PlaintextArray tmp_p(p1); // tmp = c1
    helib::Ctxt tmp(c1);
    sprintf(buffer, "tmp=c1>>=%d", (int)shamt);
    shift(ea,
          tmp_p,
          shamt); // ea.shift(tmp, random amount in [-nSlots/2,nSlots/2])
    ea.shift(tmp, shamt);
    if (!helib_test::noPrint)
      CheckCtxt(tmp, buffer);
    EXPECT_TRUE(ciphertextMatches(ea, secretKey, tmp_p, tmp));

    add(ea, p2, tmp_p); // c2 += tmp
    c2 += tmp;
    if (!helib_test::noPrint)
      CheckCtxt(c2, "c2+=tmp");
    EXPECT_TRUE(ciphertextMatches(ea, secretKey, p2, c2));

    sprintf(buffer, "c2>>>=%d", (int)rotamt);
    rotate(ea,
           p2,
           rotamt); // ea.rotate(c2, random amount in [1-nSlots, nSlots-1])
    ea.rotate(c2, rotamt);
    if (!helib_test::noPrint)
      CheckCtxt(c2, buffer);
    EXPECT_TRUE(ciphertextMatches(ea, secretKey, p2, c2));

    ::helib::negate(ea, p1); // c1.negate()
    c1.negate();
    if (!helib_test::noPrint)
      CheckCtxt(c1, "c1=-c1");
    EXPECT_TRUE(ciphertextMatches(ea, secretKey, p1, c1));

    mul(ea, p3, p2); // c3.multiplyBy(c2)
    c3.multiplyBy(c2);
    if (!helib_test::noPrint)
      CheckCtxt(c3, "c3*=c2");
    EXPECT_TRUE(ciphertextMatches(ea, secretKey, p3, c3));

    sub(ea, p0, p3); // c0 -= c3
    c0 -= c3;
    if (!helib_test::noPrint)
      CheckCtxt(c0, "c0=-c3");
    EXPECT_TRUE(ciphertextMatches(ea, secretKey, p0, c0));
  }

  c0.cleanUp();
  c1.cleanUp();
  c2.cleanUp();
  c3.cleanUp();

  HELIB_NTIMER_STOP(Circuit);

  if (!helib_test::noPrint) {
    std::cout << std::endl;
    helib::printAllTimers();
    std::cout << std::endl;
  }
  helib::resetAllTimers();
  HELIB_NTIMER_START(Check);

  EXPECT_TRUE(ciphertextMatches(ea, secretKey, p0, c0));
  EXPECT_TRUE(ciphertextMatches(ea, secretKey, p1, c1));
  EXPECT_TRUE(ciphertextMatches(ea, secretKey, p2, c2));
  EXPECT_TRUE(ciphertextMatches(ea, secretKey, p3, c3));

  HELIB_NTIMER_STOP(Check);

  std::cout << std::endl;
  if (!helib_test::noPrint) {
    helib::printAllTimers();
    std::cout << std::endl;
  }
}

// These helper functions and classes are for the test
// rotate1DWithBadDimensions
static bool check_replicate(const helib::Ctxt& c1,
                            const helib::Ctxt& c0,
                            long i,
                            const helib::SecKey& sKey,
                            const helib::EncryptedArray& ea)
{
  helib::PlaintextArray pa0(ea), pa1(ea);
  ea.decrypt(c0, sKey, pa0);
  ea.decrypt(c1, sKey, pa1);
  replicate(ea, pa0, i);

  return equals(ea, pa1, pa0); // returns true if replication succeeded
}

class StopReplicate
{};

// A class that handles the replicated ciphertexts one at a time
class ReplicateTester : public helib::ReplicateHandler
{
public:
  const helib::SecKey& sKey;
  const helib::EncryptedArray& ea;
  const helib::PlaintextArray& pa;
  long B;

  double t_last, t_total;
  long pos;
  bool error;

  ReplicateTester(const helib::SecKey& _sKey,
                  const helib::EncryptedArray& _ea,
                  const helib::PlaintextArray& _pa,
                  long _B) :
      sKey(_sKey), ea(_ea), pa(_pa), B(_B)
  {
    t_last = NTL::GetTime();
    t_total = 0.0;
    pos = 0;
    error = false;
  }

  // This method is called for every replicated ciphertext: in the i'th time
  // that it is called, the cipehrtext will have in all the slots the content
  // of the i'th input slot. In this test program we only decrypt and check
  // the result, in a real program it will do something with the cipehrtext.
  virtual void handle(const helib::Ctxt& ctxt)
  {

    double t_new = NTL::GetTime();
    double t_elapsed = t_new - t_last;
    t_total += t_elapsed;

    // Decrypt and check
    helib::PlaintextArray pa1 = pa;
    replicate(ea, pa1, pos);
    helib::PlaintextArray pa2(ea);

    ea.decrypt(ctxt, sKey, pa2);
    if (!equals(ea, pa1, pa2)) {
      error = true; // record the error, if any
    }
    t_last = NTL::GetTime();

    pos++;
    if (B > 0 && pos >= B)
      throw StopReplicate();
  }
};

TEST_P(GTestGeneral, rotate1D)
{
  if (!helib_test::noPrint) {
    std::cout << "\n\n******** TestIt"
              << (helib_test::dry ? "(dry run):" : ":");
    std::cout << " R=" << R << ", p=" << p << ", r=" << r << ", d=" << d
              << ", c=" << c << ", k=" << k << ", w=" << w << ", L=" << L
              << ", m=" << m << ", gens=" << helib::vecToStr(gens)
              << ", ords=" << helib::vecToStr(ords) << std::endl;
  }

  NTL::ZZX G;
  if (d == 0)
    G = context.alMod.getFactorsOverZZ()[0];
  else
    G = helib::makeIrredPoly(p, d);

  if (!helib_test::noPrint) {
    context.zMStar.printout();
    std::cout << std::endl;

    std::cout << "security=" << context.securityLevel() << std::endl;
    std::cout << "# small primes = " << context.smallPrimes.card() << "\n";
    std::cout << "# ctxt primes = " << context.ctxtPrimes.card() << "\n";
    std::cout << "# bits in ctxt primes = "
              << long(context.logOfProduct(context.ctxtPrimes) / log(2.0) + 0.5)
              << "\n";
    std::cout << "# special primes = " << context.specialPrimes.card() << "\n";
    std::cout << "# bits in special primes = "
              << long(context.logOfProduct(context.specialPrimes) / log(2.0) +
                      0.5)
              << "\n";
    std::cout << "G = " << G << "\n";
  }

  std::shared_ptr<helib::EncryptedArray> ea_ptr =
      std::make_shared<helib::EncryptedArray>(context, G);
  // Alias to avoid issues with previous code
  helib::EncryptedArray& ea = *ea_ptr;
  long nslots = ea.size();

  // Debugging additions
  helib::setupDebugGlobals(&secretKey, ea_ptr);

  helib::PlaintextArray p0(ea);
  helib::PlaintextArray p1(ea);
  helib::PlaintextArray p2(ea);

  random(ea, p0);
  random(ea, p1);
  random(ea, p2);

  long dim = context.zMStar.numOfGens() > 1 ? 1 : 0; // Dimension of rotation

  helib::Ctxt c0(publicKey), c1(publicKey), c2(publicKey);
  ea.encrypt(c0, publicKey, p0);
  ea.rotate1D(c0, dim, -1);
  ea.rotate1D(c0, dim, 1);

  EXPECT_TRUE(ciphertextMatches(ea, secretKey, p0, c0));

  long pos = 5; // Position of element to replicate

  ea.encrypt(c1, publicKey, p1);
  helib::Ctxt c1r = c1;
  replicate(ea, c1r, pos);
  EXPECT_TRUE(check_replicate(c1r, c1, pos, secretKey, ea));

  ea.encrypt(c2, publicKey, p2);
  ReplicateTester* handler = new ReplicateTester(secretKey, ea, p2, /*B=*/100);
  try {
    replicateAll(ea, c2, handler);
  } catch (StopReplicate) {
    std::cout << "We do not test replication more than nslots <= B (" << nslots
              << " <= " << handler->B << ")" << std::endl;
  }
  EXPECT_FALSE(handler->error);
  delete handler;
}

// clang-format off
INSTANTIATE_TEST_SUITE_P(variousParameters, GTestGeneral, ::testing::Values(
    //         R, p, r, d, c,  k,   L, s,  m,        mvec,        gens,     ords, seed, nt
    //DEEP
    Parameters(1, 2, 2, 1, 2, 10, 500, 0, 91,          {},          {},       {},    0,  1),
    Parameters(1, 2, 1, 2, 2, 10, 500, 0, 91,          {},          {},       {},    0,  1),
    Parameters(2, 7, 2, 1, 2, 10, 500, 0, 91,          {},          {},       {},    0,  1),
    // Algebra with a bad dimension
    Parameters(1, 2, 1, 1, 2, 80, 500, 0, 2761, {11, 251}, {1256, 254}, {10, -5},    0,  1)
    //FAST
    //Parameters(1, 2, 1, 1, 2, 80, 500, 0, 91,        {},          {},       {},    0,  1)
    ));
// clang-format on

// parameters to get an example where phi(m) is very
// close to m:
// m=18631 L=10 R=5

} // anonymous namespace
