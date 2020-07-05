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
#include <helib/Context.h>
#include <helib/PolyModRing.h>
#include <helib/EncryptedArray.h>
#include <helib/exceptions.h>

namespace {

TEST(TestPolyModRing, canBeConstructed)
{
  const long p = 5;
  const long r = 3;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, -1);
  helib::PolyModRing poly(p, r, G);
  EXPECT_EQ(poly.p, p);
  EXPECT_EQ(poly.r, r);
  EXPECT_EQ(poly.G, G);
  EXPECT_EQ(poly.p2r, pow(p, r));
}

TEST(TestPolyModRing, canBeCopyConstructed)
{
  const long p = 5;
  const long r = 3;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, -1);
  helib::PolyModRing poly(p, r, G);
  helib::PolyModRing copy(poly);
  EXPECT_EQ(poly, copy);
}

TEST(TestPolyModRing, canBeMoveConstructed)
{
  const long p = 5;
  const long r = 3;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, -1);
  helib::PolyModRing poly(p, r, G);
  helib::PolyModRing copied(poly);
  helib::PolyModRing moved(std::move(poly));
  EXPECT_EQ(copied, moved);
}

TEST(TestPolyModRing, equalsAndNotEqualsAreCorrect)
{
  const long p = 5;
  const long r = 3;
  NTL::ZZX G;
  NTL::SetCoeff(G, 1, 1);
  NTL::SetCoeff(G, 0, -1);
  helib::PolyModRing poly1(p, r, G);
  helib::PolyModRing poly2(p + 2, r, G);
  helib::PolyModRing poly3(p, r + 2, G);
  helib::PolyModRing poly4(p, r, G + 2);
  helib::PolyModRing polyA(p, r, G);

  EXPECT_EQ(poly1, polyA);
  EXPECT_FALSE(poly1 != polyA);

  EXPECT_NE(poly1, poly2);
  EXPECT_NE(poly1, poly3);
  EXPECT_NE(poly1, poly4);

  EXPECT_FALSE(poly1 == poly2);
  EXPECT_FALSE(poly1 == poly3);
  EXPECT_FALSE(poly1 == poly4);
}

} // namespace
