/* Copyright (C) 2019-2020 IBM Corp.
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

#include <gtest/gtest.h>
#include <NTL/ZZX.h>
#include <sstream>
#include <helib/PolyMod.h>
#include <helib/exceptions.h>

namespace {
TEST(TestPolyMod, canBeDefaultConstructed) { helib::PolyMod poly; }

TEST(TestPolyMod, canBeCopyConstructed)
{
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, -1);
  auto ring = std::make_shared<helib::PolyModRing>(13, 1, G);
  helib::PolyMod poly(ring);
  helib::PolyMod copied_poly(poly);

  EXPECT_EQ(copied_poly, poly);
}

TEST(TestPolyMod, isValidCorrectlyReportsValidity)
{
  NTL::ZZX G;
  NTL::SetX(G);
  helib::PolyMod poly;
  auto ring = std::make_shared<helib::PolyModRing>(13, 1, G);
  helib::PolyMod valid_poly(ring);
  EXPECT_FALSE(poly.isValid());
  EXPECT_TRUE(valid_poly.isValid());
}

TEST(TestPolyMod, getGAndGetp2rReturnCorrectValues)
{
  NTL::ZZX G;
  NTL::SetX(G);
  long p = 13;
  auto ring = std::make_shared<helib::PolyModRing>(p, 1, G);
  helib::PolyMod poly(ring);
  EXPECT_EQ(poly.getG(), G);
  EXPECT_EQ(poly.getp2r(), p);
}

TEST(TestPolyMod, throwsIfUsedWithoutSetup)
{
  NTL::ZZX G;
  NTL::SetX(G);
  helib::PolyMod poly;
  auto ring = std::make_shared<helib::PolyModRing>(13, 1, G);
  helib::PolyMod valid_poly(ring);

  EXPECT_THROW(poly += NTL::ZZX(1), helib::LogicError);
  EXPECT_THROW(poly *= NTL::ZZX(1), helib::LogicError);
  EXPECT_THROW(poly -= NTL::ZZX(1), helib::LogicError);
  // Cast to avoid clang unused warning
  // EXPECT_THROW((void)(poly == NTL::ZZX(1)), helib::LogicError);
  EXPECT_THROW(poly += valid_poly, helib::LogicError);
  EXPECT_THROW(poly *= valid_poly, helib::LogicError);
  EXPECT_THROW(poly -= valid_poly, helib::LogicError);
  // Cast to avoid clang unused warning
  // EXPECT_THROW((void)(poly == valid_poly), helib::LogicError);
  EXPECT_THROW(valid_poly += poly, helib::LogicError);
  EXPECT_THROW(valid_poly *= poly, helib::LogicError);
  EXPECT_THROW(valid_poly -= poly, helib::LogicError);
  // Cast to avoid clang unused warning
  // EXPECT_THROW((void)(valid_poly == poly), helib::LogicError);
}

TEST(TestPolyMod, defaultConstructedPolyModDiffersFromZeroPolyMod)
{
  const long p2r = 13;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  // G = x^2

  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod poly(ring);
  helib::PolyMod default_poly;

  EXPECT_FALSE(default_poly == poly);
  EXPECT_FALSE(default_poly == NTL::ZZX(0));
  EXPECT_FALSE(default_poly == 0l);
  EXPECT_FALSE(poly == default_poly);
  EXPECT_TRUE(default_poly != poly);
  EXPECT_TRUE(default_poly != NTL::ZZX(0));
  EXPECT_TRUE(default_poly != 0l);
  EXPECT_TRUE(poly != default_poly);
}

TEST(TestPolyMod, PolyModsWithTheSameDataAreEqual)
{
  const long p2r = 13;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  // G = x^2
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod poly1(ring);
  helib::PolyMod poly2(ring);

  EXPECT_EQ(poly1, poly2);

  poly1 = {1, 2, 3};
  poly2 = {1, 2, 3};

  EXPECT_EQ(poly1, poly2);

  helib::PolyMod default1;
  helib::PolyMod default2;

  EXPECT_EQ(default1, default2);
}

TEST(TestPolyMod, throwsIfGOrP2rDiffer)
{
  const long p2r = 13;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, -1);
  auto ring1 = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  auto ring2 = std::make_shared<helib::PolyModRing>(11, 1, G); // p2r differs
  auto ring3 = std::make_shared<helib::PolyModRing>(p2r, 1, G + 1); // G differs
  helib::PolyMod poly1(ring1);
  helib::PolyMod poly2(ring2); // p2r differs
  helib::PolyMod poly3(ring3); // G differs

  EXPECT_THROW(poly1 += poly2, helib::LogicError);
  EXPECT_THROW(poly1 *= poly3, helib::LogicError);
}

TEST(TestPolyMod, canBeConstructedWithPolyModRing)
{
  const long p2r = 13;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, -1);
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod poly(ring);
}

TEST(TestPolyMod, canBeConstructedWithALong)
{
  const long p2r = 13;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, -1);
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod poly(1L, ring);
}

TEST(TestPolyMod, canBeAssignedFromALong)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, 1);
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  long data = 6;
  helib::PolyMod poly(ring);
  poly = data;
}

TEST(TestPolyMod, canBeAssignedFromAVectorOfLong)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, 1);
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  std::vector<long> data(NTL::deg(G), 6);
  helib::PolyMod poly(data, ring);
}

TEST(TestPolyMod, canBeAssignedFromStdVectorOfLong)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, 1);
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  std::vector<long> data(5, 6);
  helib::PolyMod poly(ring);
  poly = data;
}

TEST(TestPolyMod, canBeAssignedFromBraceEnclosedListOfLongs)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, 1);
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod poly(ring);
  poly = {5L, 6L};
}

TEST(TestPolyMod, canBeAssignedFromOtherPolyMod)
{
  const long p2r = 13;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, -1);
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod poly(ring);
  helib::PolyMod poly2 = poly;
}

TEST(TestPolyMod, defaultsToZeroWhenSetUpCorrectly)
{
  const long p2r = 13;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, -1);
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod poly(ring);
  EXPECT_EQ(poly, 0);
}

TEST(TestPolyMod, defaultsToZZXZeroWhenSetUpCorrectly)
{
  const long p2r = 13;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, -1);
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod poly(ring);
  EXPECT_EQ(poly, NTL::ZZX(0));
}

TEST(TestPolyMod, preservesPassedInZZX)
{
  const long p2r = 13;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX input;
  NTL::SetCoeff(input, 3, 4);
  NTL::SetCoeff(input, 2, 3);
  NTL::SetCoeff(input, 1, 2);
  NTL::SetCoeff(input, 0, 1);
  helib::PolyMod poly(input, ring);

  EXPECT_EQ(poly, input);
}

TEST(TestPolyMod, reducesInputModuloP2rAndG)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  // G = x^2
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX input;
  NTL::SetCoeff(input, 3, 4);
  NTL::SetCoeff(input, 2, 3);
  NTL::SetCoeff(input, 1, 2);
  NTL::SetCoeff(input, 0, 1);
  helib::PolyMod poly(input, ring);

  NTL::ZZX expected_result;
  NTL::SetX(expected_result);
  expected_result *= 2;
  expected_result += 1;

  EXPECT_EQ(poly, expected_result);
}

TEST(TestPolyMod, canBeAssignedFromZZX)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod poly(ring);
  NTL::ZZX input;
  NTL::SetCoeff(input, 3, 4);
  NTL::SetCoeff(input, 2, 3);
  NTL::SetCoeff(input, 1, 2);
  NTL::SetCoeff(input, 0, 1);
  poly = input;

  EXPECT_EQ(poly, input);
}

TEST(TestPolyMod, canBeConvertedToZZX)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod poly(ring);
  NTL::ZZX input;
  NTL::SetCoeff(input, 0, 1);
  NTL::SetCoeff(input, 2, 1);
  poly = input;
  NTL::ZZX converted = static_cast<NTL::ZZX>(poly);

  EXPECT_EQ(converted, input);
}

TEST(TestPolyMod, canBeConvertedToLong)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod poly(ring);
  poly = 5l;
  long converted = static_cast<long>(poly);

  EXPECT_EQ(converted, 5l);
}

TEST(TestPolyMod, canBeConvertedToCoeffVector)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod poly(ring);

  std::vector<long> coeffs{0, 1, 2, 0};
  poly = coeffs;
  std::vector<long> result = static_cast<std::vector<long>>(poly);
  ASSERT_EQ(coeffs.size(), result.size());
  ASSERT_EQ(NTL::deg(G), result.size());
  for (std::size_t i = 0; i < coeffs.size(); ++i) {
    EXPECT_EQ(coeffs[i], result[i]);
  }
}

TEST(TestPolyMod, canBeSerializedAndDeserialized)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  // G = x^2
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX input;
  NTL::SetCoeff(input, 0, 1);
  NTL::SetCoeff(input, 1, 2);
  helib::PolyMod pre_serialized_poly(input, ring);

  std::stringstream ss;
  ss << pre_serialized_poly;
  std::string serialized = "[1, 2]";
  EXPECT_EQ(ss.str(), serialized);
  helib::PolyMod deserialized(ring);
  ss >> deserialized;

  EXPECT_EQ(pre_serialized_poly, deserialized);
}

TEST(TestPolyMod, deserializeFunctionThrowsIfMoreElementsThanDegree)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  // G = x^2
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod dest(ring);

  std::stringstream ss;
  ss << "[1, 2, 3, 4]";

  EXPECT_THROW(deserialize(ss, dest), helib::IOError);
}

TEST(TestPolyMod, rightShiftOperatorThrowsIfMoreElementsThanDegree)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  // G = x^2
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  helib::PolyMod dest(ring);

  std::stringstream ss;
  ss << "[1, 2, 3, 4]";

  EXPECT_THROW(ss >> dest, helib::IOError);
}

TEST(TestPolyMod, serializeFunctionSerializesCorrectly)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G; // G^2
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX data;
  NTL::SetCoeff(data, 0, 3);
  NTL::SetCoeff(data, 1, 1);
  helib::PolyMod poly(data, ring);

  std::stringstream ss;
  serialize(ss, poly);
  std::string expected_string = "[3, 1]";

  EXPECT_EQ(ss.str(), expected_string);
}

TEST(TestPolyMod, deserializeFunctionDeserializesCorrectly)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G; // G^2
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX data;
  NTL::SetCoeff(data, 0, 3);
  NTL::SetCoeff(data, 1, 1);
  helib::PolyMod expected_result(data, ring);

  std::string string_poly = "[3, 1]";
  helib::PolyMod deserialized_poly(ring);
  std::stringstream str(string_poly);
  deserialize(str, deserialized_poly);

  EXPECT_EQ(deserialized_poly, expected_result);
}

TEST(TestPolyMod, minusOperatorCorrectlyNegates)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  // G = x^2
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX input;
  NTL::SetCoeff(input, 1, 2);
  NTL::SetCoeff(input, 0, -1);
  helib::PolyMod poly(input, ring);

  poly = -poly;

  NTL::ZZX expected_result;
  NTL::SetCoeff(expected_result, 1, 3);
  NTL::SetCoeff(expected_result, 0, 1);

  EXPECT_EQ(poly, expected_result);
}

TEST(TestPolyMod, plusOperatorWithOtherPolyModWorks)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX augend_input;
  NTL::SetCoeff(augend_input, 1, 2);
  NTL::SetCoeff(augend_input, 0, 1);
  NTL::ZZX addend_input;
  NTL::SetCoeff(addend_input, 2, 1);
  NTL::SetCoeff(addend_input, 0, 2);
  helib::PolyMod augend(augend_input, ring);
  helib::PolyMod addend(addend_input, ring);
  helib::PolyMod summand(ring);

  summand = augend + addend; // (2x + 1) + (x^2 + 2)

  NTL::ZZX expected_result;
  expected_result = augend_input + addend_input;

  EXPECT_EQ(summand, expected_result);
}

TEST(TestPolyMod, plusOperatorWithZZXWorks)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX augend_input;
  NTL::SetCoeff(augend_input, 1, 2);
  NTL::SetCoeff(augend_input, 0, 1);
  NTL::ZZX addend;
  NTL::SetCoeff(addend, 2, 1);
  NTL::SetCoeff(addend, 0, 2);
  helib::PolyMod augend(augend_input, ring);
  helib::PolyMod summand(ring);

  summand = augend + addend; // (2x + 1) + (x^2 + 2)

  NTL::ZZX expected_result;
  expected_result = augend_input + addend;

  EXPECT_EQ(summand, expected_result);
}

TEST(TestPolyMod, plusOperatorWithLongWorks)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX augend_input;
  NTL::SetCoeff(augend_input, 1, 2);
  NTL::SetCoeff(augend_input, 0, 1);
  long addend = 4;
  helib::PolyMod augend(augend_input, ring);
  helib::PolyMod summand(ring);

  summand = augend + addend; // (2x + 1) + 4

  NTL::ZZX expected_result;
  expected_result = augend_input + addend;

  EXPECT_EQ(summand, expected_result);
}

TEST(TestPolyMod, minusOperatorWithOtherPolyModWorks)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX minuend_input;
  NTL::SetCoeff(minuend_input, 1, 2);
  NTL::SetCoeff(minuend_input, 0, 4);
  NTL::ZZX subtrahend_input;
  NTL::SetCoeff(subtrahend_input, 2, 1);
  NTL::SetCoeff(subtrahend_input, 0, 2);
  helib::PolyMod minuend(minuend_input, ring);
  helib::PolyMod subtrahend(subtrahend_input, ring);
  helib::PolyMod difference(ring);

  difference = minuend - subtrahend; // (2x + 4) - (x^2 + 2)

  NTL::ZZX expected_result;
  expected_result = minuend_input - subtrahend_input;

  EXPECT_EQ(difference, expected_result);
}

TEST(TestPolyMod, minusOperatorWithZZXWorks)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX minuend_input;
  NTL::SetCoeff(minuend_input, 1, 2);
  NTL::SetCoeff(minuend_input, 0, 4);
  NTL::ZZX subtrahend;
  NTL::SetCoeff(subtrahend, 2, 1);
  NTL::SetCoeff(subtrahend, 0, 2);
  helib::PolyMod minuend(minuend_input, ring);
  helib::PolyMod difference(ring);

  difference = minuend - subtrahend; // (2x + 4) - (x^2 + 2)

  NTL::ZZX expected_result;
  expected_result = minuend_input - subtrahend;

  EXPECT_EQ(difference, expected_result);
}

TEST(TestPolyMod, minusOperatorWithLongWorks)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX minuend_input;
  NTL::SetCoeff(minuend_input, 1, 2);
  NTL::SetCoeff(minuend_input, 0, 4);
  long subtrahend = 4;
  helib::PolyMod minuend(minuend_input, ring);
  helib::PolyMod difference(ring);

  difference = minuend - subtrahend; // (2x + 4) - 4

  NTL::ZZX expected_result;
  expected_result = minuend_input - subtrahend;

  EXPECT_EQ(difference, expected_result);
}

TEST(TestPolyMod, negateCorrectlyNegates)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  // G = x^2
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX input;
  NTL::SetCoeff(input, 1, 2);
  NTL::SetCoeff(input, 0, -1);
  helib::PolyMod poly(input, ring);

  poly.negate();

  NTL::ZZX expected_result;
  NTL::SetCoeff(expected_result, 1, 3);
  NTL::SetCoeff(expected_result, 0, 1);

  EXPECT_EQ(poly, expected_result);
}

TEST(TestPolyMod, timesOperatorWithOtherPolyModWorks)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX multiplicand_input;
  NTL::SetCoeff(multiplicand_input, 1, 2);
  NTL::SetCoeff(multiplicand_input, 0, 1);
  NTL::ZZX multiplier_input;
  NTL::SetCoeff(multiplier_input, 2, 1);
  NTL::SetCoeff(multiplier_input, 0, 2);
  helib::PolyMod multiplicand(multiplicand_input, ring);
  helib::PolyMod multiplier(multiplier_input, ring);
  helib::PolyMod product(ring);

  product = multiplicand * multiplier; // (2x + 1)(x^2 + 2)

  NTL::ZZX expected_result;
  expected_result = multiplicand_input * multiplier_input;

  EXPECT_EQ(product, expected_result);
}

TEST(TestPolyMod, timesOperatorWithZZXWorks)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX multiplicand_input;
  NTL::SetCoeff(multiplicand_input, 1, 2);
  NTL::SetCoeff(multiplicand_input, 0, 1);
  NTL::ZZX multiplier;
  NTL::SetCoeff(multiplier, 2, 1);
  NTL::SetCoeff(multiplier, 0, 2);
  helib::PolyMod multiplicand(multiplicand_input, ring);
  helib::PolyMod product(ring);

  product = multiplicand * multiplier; // (2x + 1)(x^2 + 2)

  NTL::ZZX expected_result;
  expected_result = multiplicand_input * multiplier;

  EXPECT_EQ(product, expected_result);
}

TEST(TestPolyMod, timesOperatorWithLongWorks)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX multiplicand_input;
  NTL::SetCoeff(multiplicand_input, 1, 2);
  NTL::SetCoeff(multiplicand_input, 0, 1);
  long multiplier = 4;
  helib::PolyMod multiplicand(multiplicand_input, ring);
  helib::PolyMod product(ring);

  product = multiplicand * multiplier; // (2x + 1) * 4

  NTL::ZZX expected_result;
  expected_result = multiplicand_input * multiplier;

  EXPECT_EQ(product, expected_result);
}

TEST(TestPolyMod, timesEqualsWithZZXWorks)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX multiplicand_input;
  NTL::SetCoeff(multiplicand_input, 1, 2);
  NTL::SetCoeff(multiplicand_input, 0, 1);
  NTL::ZZX multiplier;
  NTL::SetCoeff(multiplier, 2, 1);
  NTL::SetCoeff(multiplier, 0, 2);
  helib::PolyMod product(multiplicand_input, ring);

  product *= multiplier; // (2x + 1)(x^2 + 2)

  NTL::ZZX expected_result;
  expected_result = multiplicand_input * multiplier;

  EXPECT_EQ(product, expected_result);
}

TEST(TestPolyMod, timesEqualsWithOtherPolyModWorks)
{
  const long p2r = 5;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  G *= G;
  // G = x^4
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX multiplicand_input;
  NTL::SetCoeff(multiplicand_input, 1, 2);
  NTL::SetCoeff(multiplicand_input, 0, 1);
  NTL::ZZX multiplier_input;
  NTL::SetCoeff(multiplier_input, 2, 1);
  NTL::SetCoeff(multiplier_input, 0, 2);
  helib::PolyMod product(multiplicand_input, ring);
  helib::PolyMod multiplier(multiplier_input, ring);

  product *= multiplier; // (2x + 1)(x^2 + 2)

  NTL::ZZX expected_result;
  expected_result = multiplicand_input * multiplier_input;

  EXPECT_EQ(product, expected_result);
}

TEST(TestPolyMod, plusEqualsWithZZXWorks)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  // G = x^2
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX augend_input;
  NTL::SetCoeff(augend_input, 1, 2);
  NTL::SetCoeff(augend_input, 0, 3);
  NTL::ZZX addend;
  NTL::SetCoeff(addend, 1, 4);
  NTL::SetCoeff(addend, 0, 2);
  helib::PolyMod sum(augend_input, ring);

  sum += addend; // (2x + 3) + (4x + 2)

  NTL::ZZX expected_result;
  expected_result = augend_input + addend;

  EXPECT_EQ(sum, expected_result);
}

TEST(TestPolyMod, plusEqualsWithOtherPolyModWorks)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  // G = x^2
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX augend_input;
  NTL::SetCoeff(augend_input, 1, 2);
  NTL::SetCoeff(augend_input, 0, 3);
  NTL::ZZX addend_input;
  NTL::SetCoeff(addend_input, 1, 4);
  NTL::SetCoeff(addend_input, 0, 2);
  helib::PolyMod sum(augend_input, ring);
  helib::PolyMod addend(addend_input, ring);

  sum += addend; // (2x + 3) + (4x + 2)

  NTL::ZZX expected_result;
  expected_result = augend_input + addend_input;

  EXPECT_EQ(sum, expected_result);
}

TEST(TestPolyMod, minusEqualsWithZZXWorks)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  // G = x^2
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX minuend_input;
  NTL::SetCoeff(minuend_input, 1, 4);
  NTL::SetCoeff(minuend_input, 0, 3);
  NTL::ZZX subtrahend;
  NTL::SetCoeff(subtrahend, 1, 1);
  NTL::SetCoeff(subtrahend, 0, -3);
  helib::PolyMod difference(minuend_input, ring);

  difference -= subtrahend; // (4x + 3) - (x - 3)

  NTL::ZZX expected_result;
  expected_result = minuend_input - subtrahend;

  EXPECT_EQ(difference, expected_result);
}

TEST(TestPolyMod, minusEqualsWithOtherPolyModWorks)
{
  const long p2r = 7;
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  // G = x^2
  auto ring = std::make_shared<helib::PolyModRing>(p2r, 1, G);
  NTL::ZZX minuend_input;
  NTL::SetCoeff(minuend_input, 1, 4);
  NTL::SetCoeff(minuend_input, 0, 3);
  NTL::ZZX subtrahend_input;
  NTL::SetCoeff(subtrahend_input, 1, 1);
  NTL::SetCoeff(subtrahend_input, 0, -3);
  helib::PolyMod difference(minuend_input, ring);
  helib::PolyMod subtrahend(subtrahend_input, ring);

  difference -= subtrahend; // (4x + 3) - (x - 3)

  NTL::ZZX expected_result;
  expected_result = minuend_input - subtrahend_input;

  EXPECT_EQ(difference, expected_result);
}

} // namespace
