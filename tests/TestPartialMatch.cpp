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

#include <bitset>

#include <helib/helib.h>
#include <helib/partialMatch.h>
#include <helib/debugging.h>

#include "gtest/gtest.h"
#include "test_common.h"

namespace {

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

class TestPartialMatch : public ::testing::TestWithParam<BGVParameters>
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

  TestPartialMatch() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      bits(GetParam().bits),
      context(m, p, r),
      secretKey((buildModChain(context, bits), context)),
      publicKey((secretKey.GenSecKey(),
                 addSome1DMatrices(secretKey),
                 addFrbMatrices(secretKey),
                 secretKey)),
      ea(*(context.ea))
  {}

  virtual void SetUp() override
  {
    if (helib_test::verbose) {
      ea.getPAlgebra().printout();
      std::cout << "r = " << context.alMod.getR() << std::endl;
      std::cout << "ctxtPrimes=" << context.ctxtPrimes
                << ", specialPrimes=" << context.specialPrimes << std::endl
                << std::endl;
    }

    helib::setupDebugGlobals(&secretKey, context.ea);
  }

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }
};

TEST(TestPartialMatch, calculateMasksThrowsWhenPassedEmptyQueryOrDatabase)
{
  helib::Context context(1024, 1087, 1);
  const helib::EncryptedArray& ea = *(context.ea);

  helib::Matrix<helib::Ptxt<helib::BGV>> empty(0l, 0l);
  helib::Matrix<helib::Ptxt<helib::BGV>> nonempty(1l, 1l);

  EXPECT_THROW(calculateMasks(ea, empty, nonempty), helib::InvalidArgument);
  EXPECT_THROW(calculateMasks(ea, nonempty, empty), helib::InvalidArgument);
}

TEST(TestPartialMatch, calculateMasksThrowsIfQueryAndDatabaseHaveUnequalWidth)
{
  helib::Context context(1024, 1087, 1);
  const helib::EncryptedArray& ea = *(context.ea);

  helib::Matrix<helib::Ptxt<helib::BGV>> query(1l, 2l);
  helib::Matrix<helib::Ptxt<helib::BGV>> database(2l, 3l);

  EXPECT_THROW(calculateMasks(ea, query, database), helib::InvalidArgument);
}

TEST(TestPartialMatch, calculateMasksReturnsZerosWhenNoMatchFound)
{
  helib::Context context(1024, 1087, 1);
  const helib::EncryptedArray& ea = *(context.ea);

  helib::Matrix<helib::Ptxt<helib::BGV>> query(
      helib::Ptxt<helib::BGV>(context, std::vector<long>(ea.size(), 0)),
      1l,
      5l);
  helib::Matrix<helib::Ptxt<helib::BGV>> database(
      helib::Ptxt<helib::BGV>(context, std::vector<long>(ea.size(), 1)),
      2l,
      5l);

  helib::Matrix<helib::Ptxt<helib::BGV>> result(
      calculateMasks(ea, query, database));
  helib::Matrix<helib::Ptxt<helib::BGV>> expected_result(
      helib::Ptxt<helib::BGV>(context, std::vector<long>(ea.size(), 0)),
      2l,
      5l);

  EXPECT_EQ(result.dims(0), database.dims(0));
  EXPECT_EQ(result.dims(1), database.dims(1));
  for (size_t i = 0; i < result.dims(0); ++i)
    for (size_t j = 0; j < result.dims(1); ++j) {
      EXPECT_EQ(result(i, j), expected_result(i, j));
    }
}

TEST(TestPartialMatch, calculateMasksReturnsOnesWhenEverythingMatches)
{
  helib::Context context(1024, 1087, 1);
  const helib::EncryptedArray& ea = *(context.ea);

  helib::Matrix<helib::Ptxt<helib::BGV>> query(
      helib::Ptxt<helib::BGV>(context, std::vector<long>(ea.size(), 1)),
      1l,
      5l);
  helib::Matrix<helib::Ptxt<helib::BGV>> database(
      helib::Ptxt<helib::BGV>(context, std::vector<long>(ea.size(), 1)),
      2l,
      5l);

  helib::Matrix<helib::Ptxt<helib::BGV>> result(
      calculateMasks(ea, query, database));
  helib::Matrix<helib::Ptxt<helib::BGV>> expected_result(
      helib::Ptxt<helib::BGV>(context, std::vector<long>(ea.size(), 1)),
      2l,
      5l);

  EXPECT_EQ(result.dims(0), database.dims(0));
  EXPECT_EQ(result.dims(1), database.dims(1));
  for (size_t i = 0; i < result.dims(0); ++i)
    for (size_t j = 0; j < result.dims(1); ++j) {
      EXPECT_EQ(result(i, j), expected_result(i, j));
    }
}

TEST(TestPartialMatch, calculateMasksWorksCorrectly)
{
  helib::Context context(1024, 1087, 1);
  const helib::EncryptedArray& ea = *(context.ea);

  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_database(1l, 5l);
  // columns/features
  std::vector<std::vector<long>> plaintext_database_numbers = {
      {2, 1, 3, 2, 2, 1, 4, 2, 3, 4, 1, 2},
      {2, 1, 3, 2, 2, 1, 4, 2, 3, 4, 1, 2},
      {5, 8, 1, 4, 7, 1, 7, 9, 5, 6, 3, 4},
      {9, 3, 7, 3, 1, 4, 9, 5, 1, 0, 1, 1},
      {1, 9, 3, 4, 5, 7, 5, 4, 5, 1, 8, 4}};
  for (int i = 0; i < 5; ++i)
    plaintext_database(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_numbers[i]);

  // columns/features
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_query(1l, 5l);
  std::vector<std::vector<long>> plaintext_query_numbers = {
      {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6},
      {7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7},
      {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
      {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
      {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}};
  for (int i = 0; i < 5; ++i)
    plaintext_query(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_query_numbers[i]);

  helib::Matrix<helib::Ptxt<helib::BGV>> expected_result(1l, 5l);
  std::vector<std::vector<long>> expected_result_numbers = {
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1},
      {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1}};
  for (int i = 0; i < 5; ++i) {
    expected_result(0, i) =
        helib::Ptxt<helib::BGV>(context, expected_result_numbers[i]);
    // Enters this loop if the algebra has more than 12 slots
    for (int j = 12; j < ea.size(); ++j)
      expected_result(0, i)[j] = 1L;
  }

  auto mask = calculateMasks(ea, plaintext_query, plaintext_database);

  EXPECT_EQ(mask.dims(0), plaintext_database.dims(0));
  EXPECT_EQ(mask.dims(1), plaintext_database.dims(1));
  for (size_t i = 0; i < mask.dims(0); ++i)
    for (size_t j = 0; j < mask.dims(1); ++j) {
      EXPECT_EQ(mask(i, j), expected_result(i, j));
    }
}

TEST(TestPartialMatch, calculateScoresThrowsIfIndexSetsAndOffsetsDifferInSize)
{
  helib::Context context(1024, 1087, 1);

  std::vector<std::vector<long>> index_sets(4, std::vector<long>{});
  std::vector<long> offsets(3, 0);
  std::vector<helib::Matrix<long>> weights(4, helib::Matrix<long>(2l, 2l));
  helib::Matrix<helib::Ptxt<helib::BGV>> mask(helib::Ptxt<helib::BGV>(context),
                                              4l,
                                              4l);

  EXPECT_THROW(calculateScores(index_sets, offsets, weights, mask),
               helib::InvalidArgument);
}

TEST(TestPartialMatch, calculateScoresThrowsIfWeightsAndOffsetsDifferInSize)
{
  helib::Context context(1024, 1087, 1);

  std::vector<std::vector<long>> index_sets(4, std::vector<long>{});
  std::vector<long> offsets(4, 0);
  std::vector<helib::Matrix<long>> weights(3, helib::Matrix<long>(2l, 2l));
  helib::Matrix<helib::Ptxt<helib::BGV>> mask(helib::Ptxt<helib::BGV>(context),
                                              4l,
                                              4l);

  EXPECT_THROW(calculateScores(index_sets, offsets, weights, mask),
               helib::InvalidArgument);
}

TEST(TestPartialMatch, calculateScoresThrowsIfWeightsAndIndexSetsDifferInSize)
{
  helib::Context context(1024, 1087, 1);

  std::vector<std::vector<long>> index_sets = {{2, 3}};
  std::vector<long> offsets = {5};
  // Weights is 4x1 so differs in size to index_sets.
  std::vector<helib::Matrix<long>> weights = {{{1}, {3}, {2}, {4}}};
  helib::Matrix<helib::Ptxt<helib::BGV>> mask(helib::Ptxt<helib::BGV>(context),
                                              4l,
                                              4l);

  EXPECT_THROW(calculateScores(index_sets, offsets, weights, mask),
               helib::InvalidArgument);
}

TEST(TestPartialMatch, calculateScoresThrowsIfMaskIsNotAColumnVector)
{
  helib::Context context(1024, 1087, 1);

  std::vector<std::vector<long>> index_sets = {{2, 3}};
  std::vector<long> offsets = {5};
  // Weights is now a 2x2 matrix but should be a column vector.
  std::vector<helib::Matrix<long>> weights = {{{1, 3}, {2, 4}}};
  helib::Matrix<helib::Ptxt<helib::BGV>> mask(helib::Ptxt<helib::BGV>(context),
                                              4l,
                                              4l);

  EXPECT_THROW(calculateScores(index_sets, offsets, weights, mask),
               helib::InvalidArgument);
}

TEST(TestPartialMatch, calculateScoresWorksCorrectly)
{
  helib::Context context(1024, 1087, 1);

  helib::Matrix<helib::Ptxt<helib::BGV>> mask(2l, 5l);
  std::vector<std::vector<long>> mask_numbers = {
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1},
      {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1}};
  for (int i = 0; i < 5; ++i) {
    // Repeating the same mask to ensure the logic is correct.
    mask(0, i) = helib::Ptxt<helib::BGV>(context, mask_numbers[i]);
    mask(1, i) = helib::Ptxt<helib::BGV>(context, mask_numbers[i]);
  }

  // Cook up some arbitrary Fs, taus, and mus
  std::vector<std::vector<long>> Fs = {{0, 3, 4}, {0, 1, 2}, {2, 3}, {2, 4}};
  std::vector<long> mus = {0, 3, 1, 2};
  std::vector<helib::Matrix<long>> taus = {{{1}, {2}, {3}},
                                           {{3}, {1}, {2}},
                                           {{1}, {1}},
                                           {{2}, {1}}};

  auto scores = calculateScores(Fs, mus, taus, mask);

  helib::Ptxt<helib::BGV> expected_result(context);
  expected_result[0] = 0;
  expected_result[1] = 0;
  expected_result[2] = 0;
  expected_result[3] = 27;
  expected_result[4] = 24;
  expected_result[5] = 0;
  expected_result[6] = 0;
  expected_result[7] = 27;
  expected_result[8] = 24;
  expected_result[9] = 0;
  expected_result[10] = 24;
  expected_result[11] = 90;

  EXPECT_EQ(scores.dims(0), mask.dims(0));
  EXPECT_EQ(scores.dims(1), 1l);
  for (size_t i = 0; i < scores.dims(0); ++i)
    for (size_t j = 0; j < scores.dims(1); ++j) {
      EXPECT_EQ(scores(i, j), expected_result);
    }
}

TEST_P(TestPartialMatch, calculateScoresForCtxtAndPtxtAreEquivalent)
{
  helib::Matrix<helib::Ptxt<helib::BGV>> mask(2l, 5l);
  std::vector<std::vector<long>> mask_numbers = {
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1},
      {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1}};
  for (int i = 0; i < 5; ++i) {
    // Repeating the same mask to ensure the logic is correct.
    mask(0, i) = helib::Ptxt<helib::BGV>(context, mask_numbers[i]);
    mask(1, i) = helib::Ptxt<helib::BGV>(context, mask_numbers[i]);
  }

  helib::Matrix<helib::Ctxt> encrypted_mask(helib::Ctxt(publicKey), 2l, 5l);
  for (std::size_t i = 0; i < mask.dims(0); ++i)
    for (std::size_t j = 0; j < mask.dims(1); ++j)
      publicKey.Encrypt(encrypted_mask(i, j), mask(i, j));

  // Cook up some arbitrary Fs, taus, and mus
  std::vector<std::vector<long>> Fs = {{0, 3, 4}, {0, 1, 2}, {2, 3}, {2, 4}};
  std::vector<long> mus = {0, 3, 1, 2};
  std::vector<helib::Matrix<long>> taus = {{{1}, {2}, {3}},
                                           {{3}, {1}, {2}},
                                           {{1}, {1}},
                                           {{2}, {1}}};

  auto plaintext_scores = calculateScores(Fs, mus, taus, mask);
  auto encrypted_scores = calculateScores(Fs, mus, taus, encrypted_mask);

  helib::Matrix<helib::Ptxt<helib::BGV>> results(
      helib::Ptxt<helib::BGV>(context),
      encrypted_scores.dims(0),
      encrypted_scores.dims(1));
  results.entrywiseOperation<helib::Ctxt>(
      encrypted_scores,
      [&](auto& ptxt, const auto& ctxt) -> decltype(auto) {
        secretKey.Decrypt(ptxt, ctxt);
        return ptxt;
      });

  EXPECT_EQ(plaintext_scores.dims(0), results.dims(0));
  EXPECT_EQ(plaintext_scores.dims(1), results.dims(1));
  for (size_t i = 0; i < plaintext_scores.dims(0); ++i)
    for (size_t j = 0; j < plaintext_scores.dims(1); ++j) {
      EXPECT_EQ(plaintext_scores(i, j), results(i, j));
    }
}

TEST_P(TestPartialMatch, calculateMasksWorksCorrectlyForCtxtAndPtxt)
{
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_database(2l, 5l);
  // columns/features
  std::vector<std::vector<long>> plaintext_database_numbers = {
      {2, 1, 3, 2, 2, 1, 4, 2, 3, 4, 1, 2},
      {2, 1, 3, 2, 2, 1, 4, 2, 3, 4, 1, 2},
      {5, 8, 1, 4, 7, 1, 7, 9, 5, 6, 3, 4},
      {9, 3, 7, 3, 1, 4, 9, 5, 1, 0, 1, 1},
      {1, 9, 3, 4, 5, 7, 5, 4, 5, 1, 8, 4}};
  for (int i = 0; i < 5; ++i) {
    plaintext_database(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_numbers[i]);
    plaintext_database(1, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_numbers[i]);
  }

  // columns/features
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_query(1l, 5l);
  std::vector<std::vector<long>> plaintext_query_numbers = {
      {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6},
      {7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7},
      {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
      {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
      {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}};
  for (int i = 0; i < 5; ++i)
    plaintext_query(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_query_numbers[i]);

  helib::Matrix<helib::Ctxt> encrypted_query(helib::Ctxt(publicKey), 1l, 5l);
  for (std::size_t i = 0; i < plaintext_query.dims(0); ++i)
    for (std::size_t j = 0; j < plaintext_query.dims(1); ++j)
      publicKey.Encrypt(encrypted_query(i, j), plaintext_query(i, j));

  auto plaintext_masks =
      calculateMasks(ea, plaintext_query, plaintext_database);
  auto encrypted_masks =
      calculateMasks(ea, encrypted_query, plaintext_database);

  helib::Matrix<helib::Ptxt<helib::BGV>> results(
      helib::Ptxt<helib::BGV>(context),
      encrypted_masks.dims(0),
      encrypted_masks.dims(1));
  results.entrywiseOperation<helib::Ctxt>(
      encrypted_masks,
      [&](auto& ptxt, const auto& ctxt) -> decltype(auto) {
        secretKey.Decrypt(ptxt, ctxt);
        return ptxt;
      });

  EXPECT_EQ(plaintext_masks.dims(0), results.dims(0));
  EXPECT_EQ(plaintext_masks.dims(1), results.dims(1));
  for (size_t i = 0; i < plaintext_masks.dims(0); ++i)
    for (size_t j = 0; j < plaintext_masks.dims(1); ++j) {
      EXPECT_EQ(plaintext_masks(i, j), results(i, j));
    }
}

TEST_P(TestPartialMatch,
       calculateMasksAndCalculateScoresWorksCorrectlyForCtxtAndPtxt)
{
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_database(2l, 5l);
  // columns/features
  std::vector<std::vector<long>> plaintext_database_numbers = {
      {2, 1, 3, 2, 2, 1, 4, 2, 3, 4, 1, 2},
      {2, 1, 3, 2, 2, 1, 4, 2, 3, 4, 1, 2},
      {5, 8, 1, 4, 7, 1, 7, 9, 5, 6, 3, 4},
      {9, 3, 7, 3, 1, 4, 9, 5, 1, 0, 1, 1},
      {1, 9, 3, 4, 5, 7, 5, 4, 5, 1, 8, 4}};
  for (int i = 0; i < 5; ++i) {
    plaintext_database(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_numbers[i]);
    plaintext_database(1, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_numbers[i]);
  }

  // columns/features
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_query(1l, 5l);
  std::vector<std::vector<long>> plaintext_query_numbers = {
      {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6},
      {7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7},
      {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
      {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
      {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}};
  for (int i = 0; i < 5; ++i)
    plaintext_query(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_query_numbers[i]);

  helib::Matrix<helib::Ctxt> encrypted_query(helib::Ctxt(publicKey), 1l, 5l);
  for (std::size_t i = 0; i < plaintext_query.dims(0); ++i)
    for (std::size_t j = 0; j < plaintext_query.dims(1); ++j)
      publicKey.Encrypt(encrypted_query(i, j), plaintext_query(i, j));

  auto plaintext_mask = calculateMasks(ea, plaintext_query, plaintext_database);
  auto encrypted_mask = calculateMasks(ea, encrypted_query, plaintext_database);

  // Cook up some arbitrary Fs, taus, and mus
  std::vector<std::vector<long>> Fs = {{0, 3, 4}, {0, 1, 2}, {2, 3}, {2, 4}};
  std::vector<long> mus = {0, 3, 1, 2};
  std::vector<helib::Matrix<long>> taus = {{{1}, {2}, {3}},
                                           {{3}, {1}, {2}},
                                           {{1}, {1}},
                                           {{2}, {1}}};

  auto plaintext_scores = calculateScores(Fs, mus, taus, plaintext_mask);
  auto encrypted_scores = calculateScores(Fs, mus, taus, encrypted_mask);

  helib::Matrix<helib::Ptxt<helib::BGV>> results(
      helib::Ptxt<helib::BGV>(context),
      encrypted_scores.dims(0),
      encrypted_scores.dims(1));
  results.entrywiseOperation<helib::Ctxt>(
      encrypted_scores,
      [&](auto& ptxt, const auto& ctxt) -> decltype(auto) {
        secretKey.Decrypt(ptxt, ctxt);
        return ptxt;
      });

  EXPECT_EQ(plaintext_scores.dims(0), results.dims(0));
  EXPECT_EQ(plaintext_scores.dims(1), results.dims(1));
  for (size_t i = 0; i < plaintext_scores.dims(0); ++i)
    for (size_t j = 0; j < plaintext_scores.dims(1); ++j) {
      EXPECT_EQ(plaintext_scores(i, j), results(i, j));
    }
}

TEST_P(TestPartialMatch, partialMatchEncodeEncodesIntegersCorrectly)
{
  for (uint32_t input = 0; input < 1Lu << 30; input = 2 * (input + 1)) {
    helib::PolyMod output = partialMatchEncode(input, context);
    auto vec = static_cast<std::vector<long>>(output);
    uint32_t result = 0;
    ASSERT_EQ(vec.size(), context.zMStar.getOrdP());
    for (long i = 0, multiplier = 1; i < long(vec.size()); multiplier *= p, ++i)
      result += multiplier * vec[i];
    EXPECT_EQ(input, result);
  }
}

TEST_P(TestPartialMatch, databaseLookupWorksCorrectly)
{
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_database(2l, 5l);
  // columns/features
  std::vector<std::vector<long>> plaintext_database_numbers = {
      {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6},
      {7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7},
      {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
      {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
      {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}};
  for (int i = 0; i < 5; ++i) {
    plaintext_database(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_numbers[i]);
    plaintext_database(1, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_numbers[i]);
  }

  // columns/features
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_query(1l, 5l);
  std::vector<std::vector<long>> plaintext_query_numbers = {
      {6, 6, 6, 2, 6, 6, 4, 9, 6, 2, 1, 4},
      {7, 8, 1, 3, 9, 4, 3, 8, 6, 2, 1, 4},
      {2, 2, 2, 5, 9, 1, 2, 2, 6, 2, 1, 4},
      {1, 1, 5, 1, 1, 6, 1, 8, 6, 2, 1, 4},
      {4, 4, 9, 4, 2, 4, 6, 4, 6, 2, 1, 4}};
  for (int i = 0; i < 5; ++i)
    plaintext_query(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_query_numbers[i]);

  auto mask = calculateMasks(ea, plaintext_query, plaintext_database);

  // Fs, taus, and mus to perform the query, q = (c0 & c2) ^ (c3 & c4)
  std::vector<std::vector<long>> Fs = {{0, 1, 2, 3, 4},
                                       {0, 1, 2, 3, 4},
                                       {0, 1, 2, 3, 4},
                                       {0, 1, 2, 3, 4}};
  std::vector<helib::Matrix<long>> taus = {
      // The taus are the cartesian product of the above query q
      {{1}, {0}, {0}, {1}, {0}},        // Either 0-th or 3rd column matches
      {{1}, {0}, {0}, {0}, {1}},        // Either 0-th or 4th column matches
      {{0}, {0}, {1}, {1}, {0}},        // Either 2nd or 3rd column matches
      {{0}, {0}, {1}, {0}, {1}}};       // Either 2nd or 4th column matches
  std::vector<long> mus = {0, 0, 0, 0}; // Offset

  auto scores = calculateScores(Fs, mus, taus, mask);

  // FLT on scores
  scores.apply([&](auto& ptxt) {
    ptxt.power(context.alMod.getPPowR() - 1);
    return ptxt;
  });

  // Initialize all as one so no need to worry about padding matching
  helib::Ptxt<helib::BGV> expected_result(context, NTL::ZZX(1l));
  expected_result[0] = 1;  // First row identical      -> q = true
  expected_result[1] = 1;  // Only c1 differs          -> q = true
  expected_result[2] = 1;  // Only (c0 and c2) matches -> q = true
  expected_result[3] = 1;  // Only (c3 and c4) matches -> q = true
  expected_result[4] = 0;  // Only c0 and c3 match     -> q = false
  expected_result[5] = 0;  // Only c0 and c4 match     -> q = false
  expected_result[6] = 0;  // Only c2 and c3 match     -> q = false
  expected_result[7] = 0;  // Only c2 and c4 match     -> q = false
  expected_result[8] = 0;  // Only c0 matches          -> q = false
  expected_result[9] = 0;  // Only c2 matches          -> q = false
  expected_result[10] = 0; // Only c3 matches          -> q = false
  expected_result[11] = 0; // Only c4 matches          -> q = false

  EXPECT_EQ(scores.dims(0), mask.dims(0));
  EXPECT_EQ(scores.dims(1), 1l);
  for (size_t i = 0; i < scores.dims(0); ++i)
    for (size_t j = 0; j < scores.dims(1); ++j) {
      EXPECT_EQ(scores(i, j), expected_result);
    }
}

TEST_P(TestPartialMatch, databaseLookupWorksCorrectlyForCtxtAndPtxt)
{
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_database(2l, 5l);
  // columns/features
  std::vector<std::vector<long>> plaintext_database_numbers = {
      {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6},
      {7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7},
      {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
      {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
      {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}};
  for (int i = 0; i < 5; ++i) {
    plaintext_database(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_numbers[i]);
    plaintext_database(1, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_numbers[i]);
  }

  // columns/features
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_query(1l, 5l);

  // TODO - nicer way of doing this
  std::vector<std::vector<long>> plaintext_query_numbers = {
      {6, 6, 6, 6, 6, 1, 6, 9, 6, 6, 6, 6},
      {4, 8, 1, 6, 9, 4, 3, 8, 2, 9, 2, 5},
      {2, 3, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2},
      {1, 1, 1, 1, 8, 1, 1, 1, 5, 2, 1, 9},
      {4, 4, 4, 4, 4, 4, 1, 4, 4, 2, 4, 1}};
  for (int i = 0; i < 5; ++i)
    plaintext_query(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_query_numbers[i]);

  // Encrypt the query
  helib::Matrix<helib::Ctxt> encrypted_query(helib::Ctxt(publicKey), 1l, 5l);
  for (std::size_t i = 0; i < plaintext_query.dims(0); ++i)
    for (std::size_t j = 0; j < plaintext_query.dims(1); ++j)
      publicKey.Encrypt(encrypted_query(i, j), plaintext_query(i, j));

  auto plaintext_mask = calculateMasks(ea, plaintext_query, plaintext_database);
  auto encrypted_mask = calculateMasks(ea, encrypted_query, plaintext_database);

  // Fs, taus, and mus to perform the query, q = c0 & c2 & (c3 | c4)
  std::vector<std::vector<long>> Fs = {{0, 1, 2, 3, 4},
                                       {0, 1, 2, 3, 4},
                                       {0, 1, 2, 3, 4},
                                       {0, 1, 2, 3, 4}};
  std::vector<helib::Matrix<long>> taus = {
      {{1}, {0}, {0}, {0}, {0}}, // Only 0-th column matches
      {{0}, {0}, {1}, {0}, {0}}, // Only 2nd column matches
      {{0}, {0}, {0}, {1}, {1}}, // Either 3rd or 4th column matches
      {{0}, {0}, {0}, {0}, {0}}};
  std::vector<long> mus = {0, 0, 0, 1}; // Offset to cancel out the tau_{F3} = 0

  auto plaintext_scores = calculateScores(Fs, mus, taus, plaintext_mask);
  auto encrypted_scores = calculateScores(Fs, mus, taus, encrypted_mask);

  // FLT on scores
  plaintext_scores.apply([&](auto& ptxt) {
    ptxt.power(context.alMod.getPPowR() - 1);
    return ptxt;
  });
  encrypted_scores.apply([&](auto& ctxt) {
    ctxt.power(context.alMod.getPPowR() - 1);
    return ctxt;
  });

  // Decrypt the result
  helib::Matrix<helib::Ptxt<helib::BGV>> results(
      helib::Ptxt<helib::BGV>(context),
      encrypted_scores.dims(0),
      encrypted_scores.dims(1));
  results.entrywiseOperation<helib::Ctxt>(
      encrypted_scores,
      [&](auto& ptxt, const auto& ctxt) -> decltype(auto) {
        secretKey.Decrypt(ptxt, ctxt);
        return ptxt;
      });

  EXPECT_EQ(plaintext_scores.dims(0), results.dims(0));
  EXPECT_EQ(plaintext_scores.dims(1), results.dims(1));
  for (size_t i = 0; i < plaintext_scores.dims(0); ++i)
    for (size_t j = 0; j < plaintext_scores.dims(1); ++j) {
      EXPECT_EQ(plaintext_scores(i, j), results(i, j));
    }
}

TEST(TestPartialMatch, databaseLookupQueryAPIGeneratesPostFix)
{
  const helib::QueryExpr& name = helib::makeQueryExpr(0);
  const helib::QueryExpr& age = helib::makeQueryExpr(1);
  const helib::QueryExpr& height = helib::makeQueryExpr(2);
  const helib::QueryExpr& weight = helib::makeQueryExpr(3);

  helib::QueryExpr res = name && (age || (height && weight));

  EXPECT_EQ("0 1 2 3 && || &&", res->eval());
}

TEST(TestPartialMatch, databaseLookupQueryAPIGeneratesMusAndTaus)
{
  const helib::QueryExpr& name = helib::makeQueryExpr(0);
  const helib::QueryExpr& age = helib::makeQueryExpr(1);
  const helib::QueryExpr& height = helib::makeQueryExpr(2);
  const helib::QueryExpr& weight = helib::makeQueryExpr(3);

  helib::QueryBuilder qb(name && (age || (height && weight)));
  // 0 && (1 || (2 && 3)) = 0 && (1 || 2) && (1 && 3)

  long columns = 5;

  helib::Query_t query = qb.build(columns);

  std::vector<std::vector<long>> expected_Fs = {{0, 1, 2, 3, 4},
                                                {0, 1, 2, 3, 4},
                                                {0, 1, 2, 3, 4}};

  std::vector<helib::Matrix<long>> expected_taus = {
      {{1}, {0}, {0}, {0}, {0}},  // Only 0-th column
      {{0}, {1}, {1}, {0}, {0}},  // Either 1st or 2nd column
      {{0}, {1}, {0}, {1}, {0}}}; // Either 1st or 3rd column

  std::vector<long> expected_mus = {0, 0, 0};

  EXPECT_EQ(expected_Fs.size(), query.Fs.size());
  EXPECT_EQ(expected_mus.size(), query.mus.size());
  EXPECT_EQ(expected_taus.size(), query.taus.size());
  EXPECT_EQ(columns, query.taus[0].size());
  for (size_t i = 0; i < expected_Fs.size(); ++i)
    EXPECT_EQ(expected_Fs[i], query.Fs[i]);
  for (size_t i = 0; i < expected_mus.size(); ++i)
    EXPECT_EQ(expected_mus[i], query.mus[i]) << "*** i= " << i;
  for (size_t i = 0; i < expected_taus.size(); ++i)
    EXPECT_TRUE(expected_taus[i] == query.taus[i]) << "*** i= " << i;
}

TEST_P(TestPartialMatch, databaseLookupWorksWithQueryAPI)
{
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_database(2l, 5l);
  // columns/features
  // TODO - nicer way of doing this
  std::vector<std::vector<long>> plaintext_database_numbers = {
      {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6},
      {7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7},
      {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
      {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
      {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}};
  for (int i = 0; i < 5; ++i) {
    plaintext_database(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_numbers[i]);
    plaintext_database(1, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_numbers[i]);
  }
  // Have to pass in a no-op deleter because this context is handled by the
  // fixture
  helib::Database<helib::Ptxt<helib::BGV>> database(plaintext_database,
                                                    context);

  // columns/features
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_query_data(1l, 5l);
  std::vector<std::vector<long>> plaintext_query_numbers = {
      {6, 6, 6, 6, 6, 1, 6, 9, 6, 6, 6, 6},
      {4, 8, 1, 6, 9, 4, 3, 8, 2, 9, 2, 5},
      {2, 3, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2},
      {1, 1, 1, 1, 8, 1, 1, 1, 5, 2, 1, 9},
      {4, 4, 4, 4, 4, 4, 1, 4, 4, 2, 4, 1}};
  for (int i = 0; i < 5; ++i)
    plaintext_query_data(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_query_numbers[i]);

  // Encrypt the query
  helib::Matrix<helib::Ctxt> encrypted_query_data(helib::Ctxt(publicKey),
                                                  1l,
                                                  5l);
  for (std::size_t i = 0; i < plaintext_query_data.dims(0); ++i)
    for (std::size_t j = 0; j < plaintext_query_data.dims(1); ++j)
      publicKey.Encrypt(encrypted_query_data(i, j), plaintext_query_data(i, j));

  // Build the Mus and Taus from our query
  const helib::QueryExpr& name = helib::makeQueryExpr(0);
  const helib::QueryExpr& age = helib::makeQueryExpr(1);
  const helib::QueryExpr& height = helib::makeQueryExpr(2);
  const helib::QueryExpr& weight = helib::makeQueryExpr(3);

  helib::QueryBuilder qb(name && (age || (height && weight)));
  // 0 && (1 || (2 && 3)) = 0 && (1 || 2) && (1 && 3)

  long columns = plaintext_database.dims(1);

  helib::Query_t lookup_query(qb.build(columns));

  auto plaintext_result = database.contains(lookup_query, plaintext_query_data);
  auto encrypted_result = database.contains(lookup_query, encrypted_query_data);

  // Decrypt the result
  helib::Matrix<helib::Ptxt<helib::BGV>> results(
      helib::Ptxt<helib::BGV>(context),
      encrypted_result.dims(0),
      encrypted_result.dims(1));
  results.entrywiseOperation<helib::Ctxt>(
      encrypted_result,
      [&](auto& ptxt, const auto& ctxt) -> decltype(auto) {
        secretKey.Decrypt(ptxt, ctxt);
        return ptxt;
      });

  EXPECT_EQ(plaintext_result, results);
}

TEST_P(TestPartialMatch, scoringWorksWithDatabaseAndQueryAPIs)
{
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_database(2l, 5l);
  // columns/features
  std::vector<std::vector<long>> plaintext_database_numbers = {
      {2, 1, 3, 2, 2, 1, 4, 2, 3, 4, 1, 2},
      {2, 1, 3, 2, 2, 1, 4, 2, 3, 4, 1, 2},
      {5, 8, 1, 4, 7, 1, 7, 9, 5, 6, 3, 4},
      {9, 3, 7, 3, 1, 4, 9, 5, 1, 0, 1, 1},
      {1, 9, 3, 4, 5, 7, 5, 4, 5, 1, 8, 4}};
  for (int i = 0; i < 5; ++i) {
    plaintext_database(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_numbers[i]);
    plaintext_database(1, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_numbers[i]);
  }
  // Have to pass in a no-op deleter because this context is handled by the
  // fixture
  helib::Database<helib::Ptxt<helib::BGV>> database(plaintext_database,
                                                    context);

  // columns/features
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_query(1l, 5l);
  std::vector<std::vector<long>> plaintext_query_numbers = {
      {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6},
      {7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7},
      {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
      {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
      {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}};
  for (int i = 0; i < 5; ++i)
    plaintext_query(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_query_numbers[i]);

  // Encrypt the query
  helib::Matrix<helib::Ctxt> encrypted_query(helib::Ctxt(publicKey), 1l, 5l);
  for (std::size_t i = 0; i < plaintext_query.dims(0); ++i)
    for (std::size_t j = 0; j < plaintext_query.dims(1); ++j)
      publicKey.Encrypt(encrypted_query(i, j), plaintext_query(i, j));

  // Cook up some arbitrary Fs, taus, and mus
  std::vector<std::vector<long>> Fs = {{0, 3, 4}, {0, 1, 2}, {2, 3}, {2, 4}};
  std::vector<long> mus = {0, 3, 1, 2};
  std::vector<helib::Matrix<long>> taus = {{{1}, {2}, {3}},
                                           {{3}, {1}, {2}},
                                           {{1}, {1}},
                                           {{2}, {1}}};

  // Build the query
  helib::Query_t weighted_query(Fs, mus, taus, false);

  auto plaintext_scores = database.getScore(weighted_query, plaintext_query);
  auto encrypted_scores = database.getScore(weighted_query, encrypted_query);

  helib::Matrix<helib::Ptxt<helib::BGV>> results(
      helib::Ptxt<helib::BGV>(context),
      encrypted_scores.dims(0),
      encrypted_scores.dims(1));
  results.entrywiseOperation<helib::Ctxt>(
      encrypted_scores,
      [&](auto& ptxt, const auto& ctxt) -> decltype(auto) {
        secretKey.Decrypt(ptxt, ctxt);
        return ptxt;
      });

  EXPECT_EQ(plaintext_scores.dims(0), results.dims(0));
  EXPECT_EQ(plaintext_scores.dims(1), results.dims(1));
  for (size_t i = 0; i < plaintext_scores.dims(0); ++i)
    for (size_t j = 0; j < plaintext_scores.dims(1); ++j) {
      EXPECT_EQ(plaintext_scores(i, j), results(i, j));
    }
}

TEST_P(TestPartialMatch, scoringWorksWithEncryptedDatabaseAndQueryAPIs)
{
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_database_data(2l, 5l);
  // columns/features
  std::vector<std::vector<long>> plaintext_database_data_numbers = {
      {2, 1, 3, 2, 2, 1, 4, 2, 3, 4, 1, 2},
      {2, 1, 3, 2, 2, 1, 4, 2, 3, 4, 1, 2},
      {5, 8, 1, 4, 7, 1, 7, 9, 5, 6, 3, 4},
      {9, 3, 7, 3, 1, 4, 9, 5, 1, 0, 1, 1},
      {1, 9, 3, 4, 5, 7, 5, 4, 5, 1, 8, 4}};
  for (int i = 0; i < 5; ++i) {
    plaintext_database_data(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_data_numbers[i]);
    plaintext_database_data(1, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_database_data_numbers[i]);
  }
  // Encrypt the database_data
  helib::Matrix<helib::Ctxt> encrypted_database_data(helib::Ctxt(publicKey),
                                                     2l,
                                                     5l);
  for (std::size_t i = 0; i < plaintext_database_data.dims(0); ++i)
    for (std::size_t j = 0; j < plaintext_database_data.dims(1); ++j)
      publicKey.Encrypt(encrypted_database_data(i, j),
                        plaintext_database_data(i, j));

  // Have to pass in a no-op deleter because this context is handled by the
  // fixture
  helib::Database<helib::Ptxt<helib::BGV>> plaintext_database(
      plaintext_database_data,
      context);

  // Have to pass in a no-op deleter because this context is handled by the
  // fixture
  helib::Database<helib::Ctxt> encrypted_database(encrypted_database_data,
                                                  context);

  // columns/features
  helib::Matrix<helib::Ptxt<helib::BGV>> plaintext_query(1l, 5l);
  std::vector<std::vector<long>> plaintext_query_numbers = {
      {6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6},
      {7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7},
      {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
      {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
      {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4}};
  for (int i = 0; i < 5; ++i)
    plaintext_query(0, i) =
        helib::Ptxt<helib::BGV>(context, plaintext_query_numbers[i]);

  // Encrypt the query
  helib::Matrix<helib::Ctxt> encrypted_query(helib::Ctxt(publicKey), 1l, 5l);
  for (std::size_t i = 0; i < plaintext_query.dims(0); ++i)
    for (std::size_t j = 0; j < plaintext_query.dims(1); ++j)
      publicKey.Encrypt(encrypted_query(i, j), plaintext_query(i, j));

  // Cook up some arbitrary Fs, taus, and mus
  std::vector<std::vector<long>> Fs = {{0, 3, 4}, {0, 1, 2}, {2, 3}, {2, 4}};
  std::vector<long> mus = {0, 3, 1, 2};
  std::vector<helib::Matrix<long>> taus = {{{1}, {2}, {3}},
                                           {{3}, {1}, {2}},
                                           {{1}, {1}},
                                           {{2}, {1}}};

  // Build the query
  helib::Query_t weighted_query(Fs, mus, taus, false);

  auto plaintext_scores =
      plaintext_database.getScore(weighted_query, plaintext_query);
  auto encrypted_scores =
      encrypted_database.getScore(weighted_query, encrypted_query);

  // Decrypt result
  helib::Matrix<helib::Ptxt<helib::BGV>> results(
      helib::Ptxt<helib::BGV>(context),
      encrypted_scores.dims(0),
      encrypted_scores.dims(1));
  results.entrywiseOperation<helib::Ctxt>(
      encrypted_scores,
      [&](auto& ptxt, const auto& ctxt) -> decltype(auto) {
        secretKey.Decrypt(ptxt, ctxt);
        return ptxt;
      });

  // Compare decrypted and plaintext results
  EXPECT_EQ(plaintext_scores.dims(0), results.dims(0));
  EXPECT_EQ(plaintext_scores.dims(1), results.dims(1));
  for (size_t i = 0; i < plaintext_scores.dims(0); ++i)
    for (size_t j = 0; j < plaintext_scores.dims(1); ++j) {
      EXPECT_EQ(plaintext_scores(i, j), results(i, j));
    }
}

INSTANTIATE_TEST_SUITE_P(variousParameters,
                         TestPartialMatch,
                         ::testing::Values(BGVParameters(1024, 1087, 1, 700)));

} // namespace
