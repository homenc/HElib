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
#include <cstring>
#include <fstream>
#include <utility>
#include <unistd.h>

#include <NTL/ZZX.h>
#include <NTL/vector.h>

#include <helib/helib.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

bool isLittleEndian()
{
  int i = 1;
  return static_cast<bool>(*reinterpret_cast<char*>(&i));
}

void cleanupFiles(const char* file)
{
  if (unlink(file)) {
    std::cerr << "Delete of " << file << " failed." << std::endl;
  }
}

template <typename... Files>
void cleanupFiles(const char* file, Files... files)
{
  cleanupFiles(file);
  cleanupFiles(files...);
}

// Compare two binary files, expecting equality of length and content
::testing::AssertionResult filesAreEqual(std::string filename1,
                                         std::string filename2)
{
  std::ifstream file1(filename1);
  std::ifstream file2(filename2);

  if (!file1.is_open() && !file2.is_open()) {
    throw std::runtime_error("Could not open one of the following files: " +
                             filename1 + " and/or " + filename2 + ".");
  }

  std::fstream::pos_type file1size, file2size;

  file1size = file1.seekg(0, std::ifstream::end).tellg();
  file1.seekg(0, std::ifstream::beg);

  file2size = file2.seekg(0, std::ifstream::end).tellg();
  file2.seekg(0, std::ifstream::beg);

  // Quick compare sizes.
  if (file1size != file2size) {
    file1.close();
    file2.close();
    return ::testing::AssertionFailure()
           << "Files " << filename1 << " and " << filename2
           << " not the same size :(" << std::endl;
  }

  // Now compare byte blocks at a time.
  const size_t BLOCKSIZE = 4096; // 4 kB

  char buffer1[BLOCKSIZE];
  char buffer2[BLOCKSIZE];
  size_t curBlckSz = 0;

  for (size_t i = file1size, cnt = 0; i > 0; i -= curBlckSz, cnt++) {

    curBlckSz = std::min(BLOCKSIZE, i);

    file1.read(buffer1, curBlckSz);
    file2.read(buffer2, curBlckSz);

    if (memcmp(buffer1, buffer2, curBlckSz)) {
      return ::testing::AssertionFailure()
             << "Block " << cnt << " (block size: " << BLOCKSIZE << " bytes) "
             << cnt << " does not match :(" << std::endl;
    }
  }
  return ::testing::AssertionSuccess(); // Files are the same!
}

struct Parameters
{
  Parameters(long m,
             long r,
             long p,
             long c,
             long L,
             const std::string& sampleFilePrefix,
             bool cleanup) :
      m(m),
      r(r),
      p(p),
      c(c),
      L(L),
      sampleFilePrefix(sampleFilePrefix),
      cleanup(cleanup){};
  long m;
  long r;
  long p;
  long c;
  long L;
  std::string sampleFilePrefix;
  bool cleanup;

  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "m=" << params.m << ","
              << "r=" << params.r << ","
              << "p=" << params.p << ","
              << "c=" << params.c << ","
              << "L=" << params.L << ","
              << "sampleFilePrefix=" << params.sampleFilePrefix << ","
              << "cleanup=" << params.cleanup << "}";
  };
};

class GTestBinIO : public ::testing::TestWithParam<Parameters>
{
protected:
  const long m;
  const long r;
  const long p;
  const long c;
  const long w;
  const long L;
  const std::string sampleFilePrefix;
  const bool cleanup;

  static std::string testResourcePath;
  static std::string asciiFile1;
  static std::string asciiFile2;
  static std::string binFile1;
  static std::string otherEndianFileOut;

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

  GTestBinIO() :
      m(GetParam().m),
      r(GetParam().r),
      p(GetParam().p),
      c(GetParam().c),
      w(64),
      L(GetParam().L),
      sampleFilePrefix(GetParam().sampleFilePrefix),
      cleanup(GetParam().cleanup){};

  void SetUp() override{
      // Nothing clear to do here for now the way the test is written
  };

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }

public:
  static void SetUpTestCase()
  {
    testResourcePath = getTestResourcePath();
    asciiFile1 = testResourcePath + "/iotest_ascii1.txt";
    asciiFile2 = testResourcePath + "/iotest_ascii2.txt";
    binFile1 = testResourcePath + "/iotest_bin.bin";
    otherEndianFileOut = testResourcePath + "/iotest_ascii3.bin";
  };
};

std::string GTestBinIO::testResourcePath;
std::string GTestBinIO::asciiFile1;
std::string GTestBinIO::asciiFile2;
std::string GTestBinIO::binFile1;
std::string GTestBinIO::otherEndianFileOut;

TEST_P(GTestBinIO, implementsBinaryFileIoCorrectly)
{
  { // 1. Write ASCII and bin files
    std::ofstream asciiFile(asciiFile1);
    std::ofstream binFile(binFile1, std::ios::binary);
    ASSERT_TRUE(asciiFile.is_open());

    std::unique_ptr<helib::Context> context(new helib::Context(m, p, r));
    helib::buildModChain(*context, L, c); // Set the modulus chain

    if (!helib_test::noPrint) {
      std::cout << "Test to write out ASCII and Binary Files.\n";
      context->zMStar.printout(); // Printout context params
      std::cout << "\tSecurity Level: " << context->securityLevel()
                << std::endl;
    }
    std::unique_ptr<helib::SecKey> secKey(new helib::SecKey(*context));
    helib::PubKey* pubKey = (helib::PubKey*)secKey.get();
    secKey->GenSecKey(w);
    helib::addSome1DMatrices(*secKey);
    helib::addFrbMatrices(*secKey);

    helib::setupDebugGlobals(secKey.get(), context->ea);

    // ASCII
    if (!helib_test::noPrint) {
      std::cout << "\tWriting ASCII1 file " << asciiFile1 << std::endl;
    }
    helib::writeContextBase(asciiFile, *context);
    asciiFile << *context << std::endl << std::endl;
    asciiFile << *pubKey << std::endl << std::endl;
    asciiFile << *secKey << std::endl << std::endl;

    // Bin
    if (!helib_test::noPrint) {
      std::cout << "\tWriting Binary file " << binFile1 << std::endl;
    }
    helib::writeContextBaseBinary(binFile, *context);
    helib::writeContextBinary(binFile, *context);
    helib::writePubKeyBinary(binFile, *pubKey);
    helib::writeSecKeyBinary(binFile, *secKey);

    asciiFile.close();
    binFile.close();
  }
  { // 2. Read in bin files and write out ASCII.
    if (!helib_test::noPrint) {
      std::cout << "Test to read binary file and write it out as ASCII"
                << std::endl;
    }

    std::ifstream inFile(binFile1, std::ios::binary);
    std::ofstream outFile(asciiFile2);

    // Read in context,
    std::unique_ptr<helib::Context> context =
        helib::buildContextFromBinary(inFile);
    helib::readContextBinary(inFile, *context);

    // Read in SecKey and PubKey.
    std::unique_ptr<helib::SecKey> secKey(new helib::SecKey(*context));
    helib::PubKey* pubKey = (helib::PubKey*)secKey.get();

    helib::setupDebugGlobals(secKey.get(), context->ea);

    helib::readPubKeyBinary(inFile, *pubKey);
    helib::readSecKeyBinary(inFile, *secKey);

    // ASCII
    if (!helib_test::noPrint) {
      std::cout << "\tWriting ASCII2 file." << std::endl;
    }
    helib::writeContextBase(outFile, *context);
    outFile << *context << std::endl << std::endl;
    outFile << *pubKey << std::endl << std::endl;
    outFile << *secKey << std::endl << std::endl;

    inFile.close();
    outFile.close();
  }
  { // 3. Compare byte-wise the two ASCII files
    if (!helib_test::noPrint) {
      std::cout << "Comparing the two ASCII files\n";
    }

    ASSERT_TRUE(filesAreEqual(asciiFile1, asciiFile2));
  }
  { // 4. Read in binary and perform operation.
    if (!helib_test::noPrint) {
      std::cout << "Test reading in Binary files and performing an operation "
                   "between two ctxts\n";
    }

    std::ifstream inFile(binFile1, std::ios::binary);

    // Read in context,
    std::unique_ptr<helib::Context> context =
        helib::buildContextFromBinary(inFile);
    helib::readContextBinary(inFile, *context);

    // Read in PubKey.
    std::unique_ptr<helib::SecKey> secKey(new helib::SecKey(*context));
    helib::PubKey* pubKey = (helib::PubKey*)secKey.get();
    helib::readPubKeyBinary(inFile, *pubKey);
    helib::readSecKeyBinary(inFile, *secKey);
    inFile.close();

    helib::setupDebugGlobals(secKey.get(), context->ea);

    // Get the ea
    const helib::EncryptedArray& ea = *context->ea;

    // Setup some ptxts and ctxts.
    helib::Ctxt c1(*pubKey), c2(*pubKey);
    helib::PlaintextArray p1(ea), p2(ea);

    random(ea, p1);
    random(ea, p2);

    ea.encrypt(c1, *pubKey, p1);
    ea.encrypt(c2, *pubKey, p2);

    // Operation multiply and add.
    mul(ea, p1, p2);
    c1.multiplyBy(c2);
    // c1 *= c2;

    // Decrypt and Compare.
    helib::PlaintextArray pp1(ea);
    ea.decrypt(c1, *secKey, pp1);

    ASSERT_TRUE(equals(ea, p1, pp1)) << "BAD";

    if (cleanup) {
      if (!helib_test::noPrint)
        std::cout << "Clean up. Deleting created files." << std::endl;
      cleanupFiles(asciiFile1.c_str(), asciiFile2.c_str(), binFile1.c_str());
    }
  }
  { // 5. Read in binary from opposite little endian and print ASCII and compare
    bool littleEndian = isLittleEndian();

    std::string otherEndianFileIn = testResourcePath + sampleFilePrefix +
                                    (littleEndian ? "_BE.bin" : "_LE.bin");
    std::string otherEndianASCII = testResourcePath + sampleFilePrefix +
                                   (littleEndian ? "_BE.txt" : "_LE.txt");

    if (!helib_test::noPrint) {
      std::cout << "Test to read in" << (littleEndian ? " BE " : " LE ")
                << "binary file and write it out as ASCII" << std::endl;
    }

    if (sampleFilePrefix.empty()) {
      if (!helib_test::noPrint)
        std::cout << "\tSample prefix not provided, test not done."
                  << std::endl;
    } else {
      if (!helib_test::noPrint) {
        std::cout << "\tSample file used: " << otherEndianFileIn << std::endl;
      }

      std::ifstream inFile(otherEndianFileIn, std::ios::binary);

      ASSERT_TRUE(inFile.is_open())
          << "file " << otherEndianFileIn << " could not be opened.";
      std::ofstream outFile(otherEndianFileOut);

      // Read in context,
      std::unique_ptr<helib::Context> context =
          helib::buildContextFromBinary(inFile);
      helib::readContextBinary(inFile, *context);

      // Read in SecKey and PubKey.
      std::unique_ptr<helib::SecKey> secKey(new helib::SecKey(*context));
      helib::PubKey* pubKey = (helib::PubKey*)secKey.get();

      helib::readPubKeyBinary(inFile, *pubKey);
      helib::readSecKeyBinary(inFile, *secKey);
      inFile.close();

      // ASCII
      if (!helib_test::noPrint) {
        std::cout << "\tWriting other endian file." << std::endl;
      }
      helib::writeContextBase(outFile, *context);
      outFile << *context << std::endl << std::endl;
      outFile << *pubKey << std::endl << std::endl;
      outFile << *secKey << std::endl << std::endl;
      outFile.close();

      // Compare byte-wise the two ASCII files
      if (!helib_test::noPrint) {
        std::cout << "Comparing the two ASCII files\n";
      }

      ASSERT_TRUE(filesAreEqual(otherEndianASCII, otherEndianFileOut));

      if (cleanup) {
        if (!helib_test::noPrint) {
          std::cout << "Clean up. Deleting created files." << std::endl;
        }
        cleanupFiles(otherEndianFileOut.c_str());
      }
    }
  }
}

INSTANTIATE_TEST_SUITE_P(
    representativeParameters,
    GTestBinIO,
    ::testing::Values(
        // SLOW
        Parameters(127, 2, 2, 2, 300, std::string{}, true),
        Parameters(127, 1, 257, 2, 300, std::string{}, true)
        // FAST
        // Parameters(91, 1, 2, 2, 300, std::string{}, true)
        ));

} // namespace
