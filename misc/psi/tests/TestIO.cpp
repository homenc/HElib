/* Copyright (C) 2022 Intel Corporation
 * SPDX-License-Identifier: Apache-2.0
 */

#include <gtest/gtest.h>
#include <helib/helib.h>
#include <helib/debugging.h>
#include <helib/partialMatch.h>

#include "io.h"

namespace {

using Ptxt = helib::Ptxt<helib::BGV>;

struct BGVParameters
{
  BGVParameters(unsigned m, unsigned p, unsigned r, unsigned bits) :
      m(m), p(p), r(r), bits(bits){};

  const unsigned m;
  const unsigned p;
  const unsigned r;
  const unsigned bits;

  friend std::ostream& operator<<(std::ostream& os, const BGVParameters& params)
  {
    return os << "{"
              << "m = " << params.m << ", "
              << "p = " << params.p << ", "
              << "r = " << params.r << ", "
              << "bits = " << params.bits << "}";
  }
};

class TestIO : public ::testing::TestWithParam<BGVParameters>
{
protected:
  const unsigned long m;
  const unsigned long p;
  const unsigned long r;
  const unsigned long bits;
  helib::Context context;
  helib::SecKey secretKey;
  helib::PubKey publicKey;
  const helib::EncryptedArray& ea;

  TestIO() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      bits(GetParam().bits),
      context(helib::ContextBuilder<helib::BGV>()
                  .m(m)
                  .p(p)
                  .r(r)
                  .bits(bits)
                  .build()),
      secretKey(context),
      publicKey((secretKey.GenSecKey(),
                 addFrbMatrices(secretKey),
                 addSome1DMatrices(secretKey),
                 secretKey)),
      ea(context.getEA())
  {}

  virtual void SetUp() override
  {
    helib::setupDebugGlobals(&secretKey, context.shareEA());
  };

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }

  virtual ~TestIO() = default;
};

TEST_P(TestIO, readerThrowsForNonVectorPtxtQuery)
{
  std::vector<long> d1 = {1, 2, 3, 4};
  std::vector<long> d2 = {5, 6, 7, 8};
  std::vector<long> d3 = {9, 10, 11, 12};
  std::vector<long> d4 = {13, 14, 15, 16};

  Ptxt p1(context, d1);
  Ptxt p2(context, d2);
  Ptxt p3(context, d3);
  Ptxt p4(context, d4);
  Ptxt zero_ptxt(context);
  helib::Matrix<Ptxt> ptxt_query(zero_ptxt, 2, 2);
  ptxt_query(0, 0) = p1;
  ptxt_query(1, 0) = p2;
  ptxt_query(0, 1) = p3;
  ptxt_query(1, 1) = p4;

  std::string ptxtFile = "query.ptxt";
  writeResultsToFile(ptxtFile, ptxt_query);

  // This is just to allow us to call readQueryFromFile
  std::stringstream ss;
  ss << context;
  std::shared_ptr<helib::Context> contextp(helib::Context::readPtrFromJSON(ss));

  EXPECT_THROW(readQueryFromFile<Ptxt>(ptxtFile, publicKey),
               std::runtime_error);
  std::remove(ptxtFile.c_str());
}

TEST_P(TestIO, readerThrowsForNonVectorCtxtQuery)
{
  std::vector<long> d1 = {1, 2, 3, 4};
  std::vector<long> d2 = {5, 6, 7, 8};
  std::vector<long> d3 = {9, 10, 11, 12};
  std::vector<long> d4 = {13, 14, 15, 16};

  Ptxt p1(context, d1);
  Ptxt p2(context, d2);
  Ptxt p3(context, d3);
  Ptxt p4(context, d4);
  Ptxt zero_ptxt(context);
  helib::Matrix<Ptxt> ptxt_query(zero_ptxt, 2, 2);
  ptxt_query(0, 0) = p1;
  ptxt_query(1, 0) = p2;
  ptxt_query(0, 1) = p3;
  ptxt_query(1, 1) = p4;

  helib::Ctxt c1(publicKey);
  helib::Ctxt c2(publicKey);
  helib::Ctxt c3(publicKey);
  helib::Ctxt c4(publicKey);

  publicKey.Encrypt(c1, p1);
  publicKey.Encrypt(c2, p2);
  publicKey.Encrypt(c3, p3);
  publicKey.Encrypt(c4, p4);

  std::string ctxtFile = "query.ctxt";
  helib::Ctxt zero_ctxt(publicKey);
  helib::Matrix<helib::Ctxt> encrypted_query(zero_ctxt, 2, 2);
  encrypted_query(0, 0) = c1;
  encrypted_query(1, 0) = c2;
  encrypted_query(0, 1) = c3;
  encrypted_query(1, 1) = c4;
  writeResultsToFile(ctxtFile, encrypted_query);

  // This is just to allow us to call readQueryFromFile
  std::stringstream ss;
  ss << context;
  std::shared_ptr<helib::Context> contextp(helib::Context::readPtrFromJSON(ss));

  EXPECT_THROW(readQueryFromFile<Ptxt>(ctxtFile, publicKey),
               std::runtime_error);
  std::remove(ctxtFile.c_str());
}

TEST_P(TestIO, readingPtxtQueryWorksCorrectly)
{
  std::vector<long> d1 = {1, 2, 3, 4};
  std::vector<long> d2 = {5, 6, 7, 8};
  std::vector<long> d3 = {9, 10, 11, 12};
  std::vector<long> d4 = {13, 14, 15, 16};

  Ptxt p1(context, d1);
  Ptxt p2(context, d2);
  Ptxt p3(context, d3);
  Ptxt p4(context, d4);
  Ptxt zero_ptxt(context);
  helib::Matrix<Ptxt> ptxt_query(zero_ptxt, 1, 4);
  ptxt_query(0, 0) = p1;
  ptxt_query(0, 1) = p2;
  ptxt_query(0, 2) = p3;
  ptxt_query(0, 3) = p4;

  std::string ptxtFile = "query.ptxt";
  writeResultsToFile(ptxtFile, ptxt_query);

  // This is just to allow us to call readQueryFromFile
  std::stringstream ss;
  ss << context;
  std::shared_ptr<helib::Context> contextp(helib::Context::readPtrFromJSON(ss));

  auto deserialized_query = readQueryFromFile<Ptxt>(ptxtFile, publicKey);
  std::remove(ptxtFile.c_str());

  // Check dimensions
  EXPECT_EQ(ptxt_query.size(), deserialized_query.size());
  EXPECT_EQ(ptxt_query.dims(0), deserialized_query.dims(0));
  EXPECT_EQ(ptxt_query.dims(1), deserialized_query.dims(1));

  // Check elements
  for (std::size_t i = 0; i < deserialized_query.dims(0); ++i) {
    for (std::size_t j = 0; j < deserialized_query.dims(1); ++j) {
      EXPECT_EQ(ptxt_query(i, j), deserialized_query(i, j));
    }
  }
}

TEST_P(TestIO, readingCtxtQueryWorksCorrectly)
{
  std::vector<long> d1 = {1, 2, 3, 4};
  std::vector<long> d2 = {5, 6, 7, 8};
  std::vector<long> d3 = {9, 10, 11, 12};
  std::vector<long> d4 = {13, 14, 15, 16};

  Ptxt p1(context, d1);
  Ptxt p2(context, d2);
  Ptxt p3(context, d3);
  Ptxt p4(context, d4);
  Ptxt zero_ptxt(context);
  helib::Matrix<Ptxt> ptxt_query(zero_ptxt, 1, 4);
  ptxt_query(0, 0) = p1;
  ptxt_query(0, 1) = p2;
  ptxt_query(0, 2) = p3;
  ptxt_query(0, 3) = p4;

  helib::Ctxt c1(publicKey);
  helib::Ctxt c2(publicKey);
  helib::Ctxt c3(publicKey);
  helib::Ctxt c4(publicKey);

  publicKey.Encrypt(c1, p1);
  publicKey.Encrypt(c2, p2);
  publicKey.Encrypt(c3, p3);
  publicKey.Encrypt(c4, p4);

  std::string ctxtFile = "query.ctxt";
  helib::Ctxt zero_ctxt(publicKey);
  helib::Matrix<helib::Ctxt> encrypted_query(zero_ctxt, 1, 4);
  encrypted_query(0, 0) = c1;
  encrypted_query(0, 1) = c2;
  encrypted_query(0, 2) = c3;
  encrypted_query(0, 3) = c4;
  writeResultsToFile(ctxtFile, encrypted_query);

  // This is just to allow us to call readDbFromFile
  std::stringstream ss;
  ss << context;
  std::shared_ptr<helib::Context> contextp(helib::Context::readPtrFromJSON(ss));

  auto deserialized_query = readQueryFromFile<helib::Ctxt>(ctxtFile, publicKey);
  std::remove(ctxtFile.c_str()); // Remove file

  // Check dimensions
  EXPECT_EQ(encrypted_query.size(), deserialized_query.size());
  EXPECT_EQ(encrypted_query.dims(0), deserialized_query.dims(0));
  EXPECT_EQ(encrypted_query.dims(1), deserialized_query.dims(1));

  // Check elements
  for (std::size_t i = 0; i < encrypted_query.dims(0); ++i) {
    for (std::size_t j = 0; j < encrypted_query.dims(1); ++j) {
      Ptxt decrypted_entry(context);
      secretKey.Decrypt(decrypted_entry, deserialized_query(i, j));
      EXPECT_EQ(ptxt_query(i, j), decrypted_entry);
    }
  }
}

TEST_P(TestIO, readingPtxtDBWorksCorrectly)
{
  std::vector<long> d1 = {1, 2, 3, 4};
  std::vector<long> d2 = {5, 6, 7, 8};
  std::vector<long> d3 = {9, 10, 11, 12};
  std::vector<long> d4 = {13, 14, 15, 16};

  Ptxt p1(context, d1);
  Ptxt p2(context, d2);
  Ptxt p3(context, d3);
  Ptxt p4(context, d4);
  Ptxt zero_ptxt(context);
  helib::Matrix<Ptxt> ptxt_db(zero_ptxt, 2, 2);
  ptxt_db(0, 0) = p1;
  ptxt_db(1, 0) = p2;
  ptxt_db(0, 1) = p3;
  ptxt_db(1, 1) = p4;

  std::string ptxtFile = "db.ptxt";
  writeResultsToFile(ptxtFile, ptxt_db);

  // This is just to allow us to call readDbFromFile
  std::stringstream ss;
  ss << context;
  std::shared_ptr<helib::Context> contextp(helib::Context::readPtrFromJSON(ss));

  auto deserialized_db = readDbFromFile<Ptxt>(ptxtFile, contextp, publicKey);
  auto ptxt_matrix = deserialized_db.getData();
  std::remove(ptxtFile.c_str());

  // Check dimensions
  EXPECT_EQ(ptxt_db.size(), ptxt_matrix.size());
  EXPECT_EQ(ptxt_db.dims(0), ptxt_matrix.dims(0));
  EXPECT_EQ(ptxt_db.dims(1), ptxt_matrix.dims(1));

  // Check elements
  for (std::size_t i = 0; i < ptxt_matrix.dims(0); ++i) {
    for (std::size_t j = 0; j < ptxt_matrix.dims(1); ++j) {
      EXPECT_EQ(ptxt_db(i, j), ptxt_matrix(i, j));
    }
  }
}

TEST_P(TestIO, readingCtxtDBWorksCorrectly)
{
  std::vector<long> d1 = {1, 2, 3, 4};
  std::vector<long> d2 = {5, 6, 7, 8};
  std::vector<long> d3 = {9, 10, 11, 12};
  std::vector<long> d4 = {13, 14, 15, 16};

  Ptxt p1(context, d1);
  Ptxt p2(context, d2);
  Ptxt p3(context, d3);
  Ptxt p4(context, d4);
  Ptxt zero_ptxt(context);
  helib::Matrix<Ptxt> ptxt_db(zero_ptxt, 2, 2);
  ptxt_db(0, 0) = p1;
  ptxt_db(1, 0) = p2;
  ptxt_db(0, 1) = p3;
  ptxt_db(1, 1) = p4;

  helib::Ctxt c1(publicKey);
  helib::Ctxt c2(publicKey);
  helib::Ctxt c3(publicKey);
  helib::Ctxt c4(publicKey);

  publicKey.Encrypt(c1, p1);
  publicKey.Encrypt(c2, p2);
  publicKey.Encrypt(c3, p3);
  publicKey.Encrypt(c4, p4);

  std::string ctxtFile = "db.ctxt";
  helib::Ctxt zero_ctxt(publicKey);
  helib::Matrix<helib::Ctxt> encrypted_db(zero_ctxt, 2, 2);
  encrypted_db(0, 0) = c1;
  encrypted_db(1, 0) = c2;
  encrypted_db(0, 1) = c3;
  encrypted_db(1, 1) = c4;
  writeResultsToFile(ctxtFile, encrypted_db);

  // This is just to allow us to call readDbFromFile
  std::stringstream ss;
  ss << context;
  std::shared_ptr<helib::Context> contextp(helib::Context::readPtrFromJSON(ss));

  auto deserialized_db =
      readDbFromFile<helib::Ctxt>(ctxtFile, contextp, publicKey);
  auto encrypted_matrix = deserialized_db.getData();
  std::remove(ctxtFile.c_str()); // Remove file

  // Check dimensions
  EXPECT_EQ(encrypted_db.size(), encrypted_matrix.size());
  EXPECT_EQ(encrypted_db.dims(0), encrypted_matrix.dims(0));
  EXPECT_EQ(encrypted_db.dims(1), encrypted_matrix.dims(1));

  // Check elements
  for (std::size_t i = 0; i < encrypted_db.dims(0); ++i) {
    for (std::size_t j = 0; j < encrypted_db.dims(1); ++j) {
      Ptxt decrypted_entry(context);
      secretKey.Decrypt(decrypted_entry, encrypted_matrix(i, j));
      EXPECT_EQ(ptxt_db(i, j), decrypted_entry);
    }
  }
}

TEST_P(TestIO, ctxtDBDeserializationAgreesWithPtxtDeserialization)
{
  std::vector<long> d1 = {1, 2, 3, 4};
  std::vector<long> d2 = {5, 6, 7, 8};
  std::vector<long> d3 = {9, 10, 11, 12};
  std::vector<long> d4 = {13, 14, 15, 16};

  Ptxt p1(context, d1);
  Ptxt p2(context, d2);
  Ptxt p3(context, d3);
  Ptxt p4(context, d4);
  Ptxt zero_ptxt(context);
  helib::Matrix<Ptxt> ptxt_db(zero_ptxt, 2, 2);
  ptxt_db(0, 0) = p1;
  ptxt_db(1, 0) = p2;
  ptxt_db(0, 1) = p3;
  ptxt_db(1, 1) = p4;

  helib::Ctxt c1(publicKey);
  helib::Ctxt c2(publicKey);
  helib::Ctxt c3(publicKey);
  helib::Ctxt c4(publicKey);

  publicKey.Encrypt(c1, p1);
  publicKey.Encrypt(c2, p2);
  publicKey.Encrypt(c3, p3);
  publicKey.Encrypt(c4, p4);

  helib::Ctxt zero_ctxt(publicKey);
  helib::Matrix<helib::Ctxt> encrypted_db(zero_ctxt, 2, 2);
  encrypted_db(0, 0) = c1;
  encrypted_db(1, 0) = c2;
  encrypted_db(0, 1) = c3;
  encrypted_db(1, 1) = c4;

  std::string ctxtFile = "db.ctxt";
  writeResultsToFile(ctxtFile, encrypted_db);

  std::string ptxtFile = "db.ptxt";
  writeResultsToFile(ptxtFile, ptxt_db);

  // This is just to allow us to call readDbFromFile
  std::stringstream ss;
  ss << context;
  std::shared_ptr<helib::Context> contextp(helib::Context::readPtrFromJSON(ss));

  auto deserialized_ctxt_db =
      readDbFromFile<helib::Ctxt>(ctxtFile, contextp, publicKey);
  auto encrypted_matrix = deserialized_ctxt_db.getData();
  std::remove(ctxtFile.c_str()); // Remove file

  auto deserialized_ptxt_db =
      readDbFromFile<Ptxt>(ptxtFile, contextp, publicKey);
  auto ptxt_matrix = deserialized_ptxt_db.getData();
  std::remove(ptxtFile.c_str()); // Remove file

  // Check dimensions
  EXPECT_EQ(ptxt_matrix.size(), encrypted_matrix.size());
  EXPECT_EQ(ptxt_matrix.dims(0), encrypted_matrix.dims(0));
  EXPECT_EQ(ptxt_matrix.dims(1), encrypted_matrix.dims(1));

  // Check elements
  for (std::size_t i = 0; i < encrypted_db.dims(0); ++i) {
    for (std::size_t j = 0; j < encrypted_db.dims(1); ++j) {
      Ptxt decrypted_entry(context);
      secretKey.Decrypt(decrypted_entry, encrypted_matrix(i, j));
      EXPECT_EQ(ptxt_matrix(i, j), decrypted_entry);
    }
  }
}

INSTANTIATE_TEST_SUITE_P(variousParameters,
                         TestIO,
                         ::testing::Values(BGVParameters(24, 37, 1, 100)));

} // namespace
