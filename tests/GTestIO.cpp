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

#include <fstream>
#include <unistd.h>

#include <NTL/ZZX.h>
#include <NTL/vector.h>

#include <helib/helib.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {
struct Parameters
{
  long r;  // lifting
  long p;  // plaintext base
  long c;  // number of columns in the key-switching matrices
  long mm; // cyclotomic index

  Parameters(long r, long p, long c, long mm) : r(r), p(p), c(c), mm(mm){};

  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "r=" << params.r << ","
              << "p=" << params.p << ","
              << "c=" << params.c << ","
              << "mm=" << params.mm << "}";
  };
};

class GTestIO : public ::testing::TestWithParam<Parameters>
{
protected:
#define N_TESTS 3
  static constexpr long ms[N_TESTS][10] = {
      // nSlots  m   phi(m) ord(2)
      //  {   2,    7,    6,    3,   0,0,0,0,0,0},
      {6, 31, 30, 5, 0, 0, 0, 0, 0, 0},
      {6, 127, 126, 7, 127, 1, 108, 24, 6, -3},     // gens=108(6), 24(!3)
      {60, 1023, 600, 10, 11, 93, 838, 584, 10, 6}, // gens=838(10),584(6)
      //  {  378,  5461,  5292, 14}, // gens=3(126),509(3)
      //  {  630,  8191,  8190, 13}, // gens=39(630)
      //  {  600, 13981, 12000, 20}, // gens=10(30),23(10),3(!2)
      //  {  682, 15709, 15004, 22} // gens=5(682)
  };

  static std::string keyFilePath;
  static std::string getTestResourcePath()
  {
    std::string testResourcePath;
    std::string executablePath(helib_test::path_of_executable);
    std::size_t found;
    if ((found = executablePath.find_last_of("/")) == std::string::npos) {
      // Then there's no / in the filepath
      testResourcePath = "../test_resources";
    } else {
      testResourcePath = executablePath.substr(0, found) + "/../test_resources";
    }
    return testResourcePath;
  };

  GTestIO() :
      r(GetParam().r),
      p(GetParam().p),
      c(GetParam().c),
      w(64),
      L(100),
      mm(GetParam().mm),
      useTable(mm == 0 && p == 2),
      ptxtSpace(NTL::power_long(p, r)),
      numTests(useTable ? N_TESTS : 1),
      contexts(numTests),
      sKeys(numTests),
      ctxts(numTests),
      eas(numTests),
      ptxts(numTests){};

  const long r;
  const long p;
  const long c;
  const long w;
  const long L;
  const long mm;
  const bool useTable;
  const long ptxtSpace;
  const long numTests;

  std::vector<std::unique_ptr<helib::Context>> contexts;
  std::vector<std::unique_ptr<helib::SecKey>> sKeys;
  std::vector<std::unique_ptr<helib::Ctxt>> ctxts;
  std::vector<std::unique_ptr<helib::EncryptedArray>> eas;
  std::vector<std::vector<NTL::ZZX>> ptxts;

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }

public:
  static void SetUpTestCase()
  {
    keyFilePath = getTestResourcePath() + "/iotest.txt";
  };

  static void TearDownTestCase()
  {
    unlink(keyFilePath.c_str()); // clean up before exiting
  };
};

constexpr long GTestIO::ms[N_TESTS][10];
std::string GTestIO::keyFilePath;

// Testing the I/O of the important classes of the library
// (context, keys, ciphertexts).
TEST_P(GTestIO, importantClassesRemainConsistentUnderIo)
{
  // first loop: generate stuff and write it to file

  // open file for writing
  {
    std::fstream keyFile(keyFilePath, std::fstream::out | std::fstream::trunc);
    ASSERT_TRUE(keyFile.is_open());
    for (long i = 0; i < numTests; i++) {
      long m = (mm == 0) ? ms[i][1] : mm;

      if (!helib_test::noPrint) {
        std::cout << "Testing IO: m=" << m << ", p^r=" << p << "^" << r
                  << std::endl;
      }

      NTL::Vec<long> mvec(NTL::INIT_SIZE, 2);
      mvec[0] = ms[i][4];
      mvec[1] = ms[i][5];
      std::vector<long> gens(2);
      gens[0] = ms[i][6];
      gens[1] = ms[i][7];
      std::vector<long> ords(2);
      ords[0] = ms[i][8];
      ords[1] = ms[i][9];

      if (useTable && gens[0] > 0)
        contexts[i].reset(new helib::Context(m, p, r, gens, ords));
      else
        contexts[i].reset(new helib::Context(m, p, r));
      if (!helib_test::noPrint) {
        contexts[i]->zMStar.printout();
      }

      helib::buildModChain(*contexts[i], L, c); // Set the modulus chain
      if (mm == 0 && m == 1023)
        contexts[i]->makeBootstrappable(mvec);

      // Output the Context to file
      helib::writeContextBase(keyFile, *contexts[i]);
      if (!helib_test::noPrint) {
        helib::writeContextBase(std::cout, *contexts[i]);
        std::cout << std::endl;
      }
      keyFile << *contexts[i] << std::endl;

      sKeys[i].reset(new helib::SecKey(*contexts[i]));
      sKeys[i]->GenSecKey(0, ptxtSpace); // A +-1/0 secret key
      addSome1DMatrices(
          *sKeys[i]); // compute key-switching matrices that we need
      const helib::PubKey publicKey = *sKeys[i];
      eas[i].reset(new helib::EncryptedArray(*contexts[i]));

      long nslots = eas[i]->size();

      // Output the secret key to file, twice. Below we will have two copies
      // of most things.
      keyFile << *sKeys[i] << std::endl;
      ;
      keyFile << *sKeys[i] << std::endl;
      ;

      std::vector<NTL::ZZX> b;
      long p2r = eas[i]->getContext().alMod.getPPowR();
      NTL::ZZX poly = helib::RandPoly(
          0,
          NTL::to_ZZ(p2r)); // choose a random constant polynomial
      eas[i]->decode(ptxts[i], poly);

      ctxts[i].reset(new helib::Ctxt(publicKey));
      eas[i]->encrypt(*ctxts[i], publicKey, ptxts[i]);
      eas[i]->decrypt(*ctxts[i], *sKeys[i], b);
      ASSERT_EQ(ptxts[i].size(), b.size());
      for (long j = 0; j < nslots; j++) {
        ASSERT_EQ(ptxts[i][j], b[j]);
      }

      // output the plaintext
      keyFile << "[ ";
      for (long j = 0; j < nslots; j++)
        keyFile << ptxts[i][j] << " ";
      keyFile << "]\n";

      eas[i]->encode(poly, ptxts[i]);
      keyFile << poly << std::endl;

      // Output the ciphertext to file
      keyFile << *ctxts[i] << std::endl;
      keyFile << *ctxts[i] << std::endl;
      // std::cerr << "okay " << i << std::endl<< std::endl;
    }
    keyFile.close();
  }
  // std::cerr << "so far, so good\n\n";

  // second loop: read from input and repeat the computation

  // open file for read
  {
    std::fstream keyFile(keyFilePath, std::fstream::in);
    for (long i = 0; i < numTests; i++) {

      // Read context from file
      unsigned long m1, p1, r1;
      std::vector<long> gens, ords;
      helib::readContextBase(keyFile, m1, p1, r1, gens, ords);
      helib::Context tmpContext(m1, p1, r1, gens, ords);
      keyFile >> tmpContext;
      ASSERT_EQ(*contexts[i], tmpContext);
      // std::cerr << i << ": context matches input\n";

      // We define some things below wrt *contexts[i], not tmpContext.
      // This is because the various operator== methods check equality of
      // references, not equality of the referenced Context objects.
      helib::Context& context = *contexts[i];
      helib::SecKey secretKey(context);
      helib::SecKey secretKey2(tmpContext);
      const helib::PubKey& publicKey = secretKey;
      const helib::PubKey& publicKey2 = secretKey2;

      keyFile >> secretKey;
      keyFile >> secretKey2;
      ASSERT_EQ(secretKey, *sKeys[i]);
      // std::cerr << "   secret key matches input\n";

      helib::EncryptedArray ea(context);
      helib::EncryptedArray ea2(tmpContext);

      long nslots = ea.size();

      // Read the plaintext from file
      std::vector<NTL::ZZX> a;
      a.resize(nslots);
      ASSERT_EQ(nslots, (long)ptxts[i].size());
      helib::seekPastChar(keyFile, '['); // defined in NumbTh.cpp
      for (long j = 0; j < nslots; j++) {
        keyFile >> a[j];
        ASSERT_EQ(a[j], ptxts[i][j]);
      }
      helib::seekPastChar(keyFile, ']');
      // std::cerr << "   ptxt matches input\n";

      // Read the encoded plaintext from file
      NTL::ZZX poly1, poly2;
      keyFile >> poly1;
      eas[i]->encode(poly2, a);
      ASSERT_EQ(poly1, poly2);
      // std::cerr << "   eas[i].encode(a)==poly1 okay\n";

      ea.encode(poly2, a);
      ASSERT_EQ(poly1, poly2);
      // std::cerr << "   ea.encode(a)==poly1 okay\n";

      ea2.encode(poly2, a);
      ASSERT_EQ(poly1, poly2);
      // std::cerr << "   ea2.encode(a)==poly1 okay\n";

      eas[i]->decode(a, poly1);
      ASSERT_EQ(nslots, (long)a.size());
      for (long j = 0; j < nslots; j++) {
        ASSERT_EQ(a[j], ptxts[i][j]);
      }
      // std::cerr << "   eas[i].decode(poly1)==ptxts[i] okay\n";

      ea.decode(a, poly1);
      ASSERT_EQ(nslots, (long)a.size());
      for (long j = 0; j < nslots; j++) {
        ASSERT_EQ(a[j], ptxts[i][j]);
      }
      // std::cerr << "   ea.decode(poly1)==ptxts[i] okay\n";

      ea2.decode(a, poly1);
      ASSERT_EQ(nslots, (long)a.size());
      for (long j = 0; j < nslots; j++) {
        ASSERT_EQ(a[j], ptxts[i][j]);
      }
      // std::cerr << "   ea2.decode(poly1)==ptxts[i] okay\n";

      // Read ciphertext from file
      helib::Ctxt ctxt(publicKey);
      helib::Ctxt ctxt2(publicKey2);
      keyFile >> ctxt;
      keyFile >> ctxt2;
      ASSERT_TRUE(ctxts[i]->equalsTo(ctxt, /*comparePkeys=*/false));
      // std::cerr << "   ctxt matches input\n";

      sKeys[i]->Decrypt(poly2, *ctxts[i]);
      ASSERT_EQ(poly1, poly2);
      // std::cerr << "   sKeys[i]->decrypt(*ctxts[i]) == poly1 okay\n";

      secretKey.Decrypt(poly2, *ctxts[i]);
      ASSERT_EQ(poly1, poly2);
      // std::cerr << "   secretKey.decrypt(*ctxts[i]) == poly1 okay\n";

      secretKey.Decrypt(poly2, ctxt);
      ASSERT_EQ(poly1, poly2);
      // std::cerr << "   secretKey.decrypt(ctxt) == poly1 okay\n";

      secretKey2.Decrypt(poly2, ctxt2);
      ASSERT_EQ(poly1, poly2);
      // std::cerr << "   secretKey2.decrypt(ctxt2) == poly1 okay\n";

      eas[i]->decrypt(ctxt, *sKeys[i], a);
      ASSERT_EQ(nslots, (long)a.size());
      for (long j = 0; j < nslots; j++) {
        ASSERT_EQ(a[j], ptxts[i][j]);
      }
      // std::cerr << "   eas[i].decrypt(ctxt, *sKeys[i])==ptxts[i] okay\n";

      ea.decrypt(ctxt, secretKey, a);
      ASSERT_EQ(nslots, (long)a.size());
      for (long j = 0; j < nslots; j++) {
        ASSERT_EQ(a[j], ptxts[i][j]);
      }
      // std::cerr << "   ea.decrypt(ctxt, secretKey)==ptxts[i] okay\n";

      ea2.decrypt(ctxt2, secretKey2, a);
      ASSERT_EQ(nslots, (long)a.size());
      for (long j = 0; j < nslots; j++) {
        ASSERT_EQ(a[j], ptxts[i][j]);
      }
      // std::cerr << "   ea2.decrypt(ctxt2, secretKey2)==ptxts[i] okay\n";
    }
  }
}

INSTANTIATE_TEST_SUITE_P(someSmallParameters,
                         GTestIO,
                         ::testing::Values(
                             // FAST
                             Parameters(1, 2, 2, 91)));

} // namespace
