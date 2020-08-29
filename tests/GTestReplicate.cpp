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

/* GTestReplicate.cpp - Testing the functionality of replicating one
 * slot from a vector access the whole vector (or replicating each slot
 * to a full ciphertext)
 */

#include <cassert>
#include <NTL/lzz_pXFactoring.h>

#include <helib/helib.h>
#include <helib/replicate.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {
struct Parameters
{

  long m; // cyclotomic ring
  long p; // plaintext base
  long r; // lifting
  long d; // degree of the field extension
  // d == 0 => factors[0] defines extension
  long L;   // heuristic
  long bnd; // recursion bound for replication
  long B;   // bound for # of replications

  Parameters(long m, long p, long r, long d, long L, long bnd, long B) :
      m(m), p(p), r(r), d(d), L(L), bnd(bnd), B(B){};

  friend std::ostream& operator<<(std::ostream& os, const Parameters params)
  {
    return os << "{"
              << "m=" << params.m << ","
              << "p=" << params.p << ","
              << "r=" << params.r << ","
              << "d=" << params.d << ","
              << "L=" << params.L << ","
              << "bnd=" << params.bnd << ","
              << "B=" << params.B << "}";
  };
};

class GTestReplicate : public ::testing::TestWithParam<Parameters>
{
protected:
  static void printContextAndG(const helib::Context& context, const NTL::ZZX& G)
  {
    if (!helib_test::noPrint) {
      context.zMStar.printout();
      std::cout << std::endl;
      std::cout << "G = " << G << "\n";
    }
  };

  static NTL::ZZX createG(const helib::Context& context, long p, long d)
  {
    return (d == 0) ? context.alMod.getFactorsOverZZ()[0]
                    : helib::makeIrredPoly(p, d);
  };

  GTestReplicate() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      d(GetParam().d),
      L(GetParam().L),
      bnd(GetParam().bnd),
      B(GetParam().B),
      context((helib::setDryRun(helib_test::dry), helib::setTimersOn(), m),
              p,
              r),
      secretKey((buildModChain(context, L, /*c=*/2), context)),
      G(createG(context, p, d)),
      publicKey((printContextAndG(context, G),
                 secretKey.GenSecKey(), // A +-1/0 secret key
                 addSome1DMatrices(secretKey),
                 secretKey)),
      ea(context, G),
      xp0(ea),
      xp1(ea),
      xc0(publicKey),
      xc1(publicKey){};

  virtual void SetUp() override
  {
    random(ea, xp0);
    random(ea, xp1);

    ea.encrypt(xc0, publicKey, xp0);
    ea.encode(poly_xp1, xp1);
    xc1 = xc0;

    helib::setupDebugGlobals(&secretKey, context.ea);
  };

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }

  const long m;
  const long p;
  const long r;
  const long d;
  const long L;
  const long bnd;
  const long B;
  helib::Context context;
  helib::SecKey secretKey;
  NTL::ZZX G;
  const helib::PubKey& publicKey;
  helib::EncryptedArray ea;
  helib::PlaintextArray xp0;
  helib::PlaintextArray xp1;
  helib::Ctxt xc0;
  helib::Ctxt xc1;
  NTL::ZZX poly_xp1; // To be default constructed
};

::testing::AssertionResult replicationSucceeds(const helib::Ctxt& c1,
                                               const helib::Ctxt& c0,
                                               long i,
                                               const helib::SecKey& sKey,
                                               const helib::EncryptedArray& ea)
{
  helib::PlaintextArray pa0(ea), pa1(ea);
  ea.decrypt(c0, sKey, pa0);
  ea.decrypt(c1, sKey, pa1);
  replicate(ea, pa0, i);

  if (equals(ea, pa1, pa0)) // replication succeeded
    return ::testing::AssertionSuccess();
  else
    return ::testing::AssertionFailure() << "Replication failed";
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
  // that it is called, the ciphertext will have in all the slots the content
  // of the i'th input slot. In this test program we only decrypt and check
  // the result, in a real program it will do something with the ciphertext.
  virtual void handle(const helib::Ctxt& ctxt)
  {

    double t_new = NTL::GetTime();
    double t_elapsed = t_new - t_last;
    t_total += t_elapsed;

    // Decrypt and check
    helib::PlaintextArray pa1 = pa;
    helib::replicate(ea, pa1, pos);
    helib::PlaintextArray pa2(ea);

    if (pos == 0 && !helib_test::noPrint)
      CheckCtxt(ctxt, "replicateAll");

    ea.decrypt(ctxt, sKey, pa2);
    if (!equals(ea, pa1, pa2))
      error = true; // record the error, if any
    t_last = NTL::GetTime();

    pos++;
    if (B > 0 && pos >= B)
      throw StopReplicate();
  }
};

TEST_P(GTestReplicate, replicateWorks)
{
  if (!helib_test::noPrint) {
    std::cout << "** Testing replicate():\n";
    CheckCtxt(xc1, "before replicate");
  }
  helib::replicate(ea, xc1, ea.size() / 2);
  EXPECT_TRUE(replicationSucceeds(xc1, xc0, ea.size() / 2, secretKey, ea));
  if (!helib_test::noPrint)
    CheckCtxt(xc1, "after replicate");
}

TEST_P(GTestReplicate, repeatedReplicationWorks)
{
  // Get some timing results
  for (long i = 0; i < 20 && i < ea.size(); i++) {
    xc1 = xc0;
    HELIB_NTIMER_START(replicate);
    helib::replicate(ea, xc1, i);
    EXPECT_TRUE(replicationSucceeds(xc1, xc0, i, secretKey, ea));
    HELIB_NTIMER_STOP(replicate);
  }
  if (!helib_test::noPrint) {
    helib::printAllTimers();
  }
}

TEST_P(GTestReplicate, replicateAllReplicatesAccurately)
{
  if (!helib_test::noPrint) {
    std::cout << "\n** Testing replicateAll()... " << std::flush;
  }
#ifdef HELIB_DEBUG
  helib::replicateVerboseFlag = true;
#else
  helib::replicateVerboseFlag = false;
#endif
  ReplicateTester handler(secretKey, ea, xp0, B);
  try {
    HELIB_NTIMER_START(replicateAll);
    helib::replicateAll(ea, xc0, &handler, bnd);
  } catch (StopReplicate) {
  }
  EXPECT_FALSE(handler.error);
  if (!helib_test::noPrint) {
    std::cout << "total time=" << handler.t_total << " ("
              << ((B > 0) ? B : ea.size()) << " vectors)\n";
  }
}

INSTANTIATE_TEST_SUITE_P(typicalParameters,
                         GTestReplicate,
                         ::testing::Values(
                             // SLOW
                             Parameters(1247, 2, 1, 1, 250, 64, 0)
                             // FAST
                             // Parameters(91, 2, 1, 1, 250, 64, 0)
                             ));

} // namespace
