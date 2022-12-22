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

/* Copyright (C) 2022 Intel Corporation
* SPDX-License-Identifier: Apache-2.0
*
* Added tests for separated SK, PK and Key switching matrices
*/

/* Note this file only tests JSON (de)serialization.*/
/* If you are searching for binary (de)serialization go to TestBinIO.cpp*/

#include <tuple>

#include <helib/helib.h>
#include <helib/debugging.h>
#include <io.h>

#include "test_common.h"
#include "gtest/gtest.h"

namespace {

// Returns the data vector and the expected string for BGV
static std::pair<std::vector<NTL::ZZX>, std::string> createData(long size,
                                                                long p2r,
                                                                long deg)
{
  std::vector<NTL::ZZX> data(size);
  std::stringstream ss;
  ss << "[";
  for (long i = 0; i < size; ++i) {
    NTL::ZZX input;
    NTL::SetCoeff(input, 0, i % p2r);
    ss << "[" << (i % p2r);
    long num = (i + 2) % p2r;
    if (deg != 1 && num != 0) {
      NTL::SetCoeff(input, 1, num);
      ss << "," << num;
    }
    data[i] = input;
    ss << "]" << (i != size - 1 ? "," : "");
  }
  ss << "]";

  return {std::move(data), ss.str()};
}

// Returns the data vector and the expected string for CKKS
static std::vector<std::complex<double>> createData(long size)
{
  std::vector<std::complex<double>> data(size);
  for (long i = 0; i < size; ++i) {
    data[i] = {(i + 2) / 10.0, (2 * i) / 2.5};
  }

  return data;
}

static std::string addHeader(const std::string& data,
                             const std::string& scheme,
                             const std::string& type = "Ptxt")
{
  return "{\"HElibVersion\":\"" + std::string(helib::version::asString) +
         "\",\"content\":{\"scheme\":\"" + scheme + "\",\"slots\":" + data +
         "},\"serializationVersion\":\"0.0.1\",\"type\":\"" + type + "\"}";
}

struct BGVParameters
{
  const long m;
  const long p;
  const long r;
  const long bits;
  const std::vector<long> gens;
  const std::vector<long> ords;
  const std::vector<long> mvec;

  BGVParameters(long m,
                long p,
                long r,
                long bits,
                std::vector<long> gens,
                std::vector<long> ords,
                std::vector<long> mvec) :
      m(m), p(p), r(r), bits(bits), gens(gens), ords(ords), mvec(mvec){};

  friend std::ostream& operator<<(std::ostream& os, const BGVParameters& params)
  {
    return os << "{"
              << "m = " << params.m << ", "
              << "p = " << params.p << ", "
              << "r = " << params.r << ", "
              << "bits = " << params.bits << ", "
              << "gens = " << helib::vecToStr(params.gens) << ", "
              << "ords = " << helib::vecToStr(params.ords) << ", "
              << "mvec = " << helib::vecToStr(params.mvec) << "}";
  }
};

struct CKKSParameters
{
  const long m;
  const long precision;
  const long bits;

  CKKSParameters(long m, long precision, long bits) :
      m(m), precision(precision), bits(bits){};

  friend std::ostream& operator<<(std::ostream& os,
                                  const CKKSParameters& params)
  {
    return os << "{"
              << "m = " << params.m << ", "
              << "precision = " << params.precision << ", "
              << "bits = " << params.bits << "}";
  }
};

class TestIO_BGV : public ::testing::TestWithParam<BGVParameters>
{
protected:
  const long m;
  const long p;
  const long r;
  const long bits;
  const std::vector<long> gens;
  const std::vector<long> ords;
  const std::vector<long> mvec;

  helib::Context context;
  helib::SecKey secretKey;
  helib::PubKey publicKey;
  const helib::EncryptedArray& ea;

  TestIO_BGV() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      bits(GetParam().bits),
      gens(GetParam().gens),
      ords(GetParam().ords),
      mvec(GetParam().mvec),
      context(helib::ContextBuilder<helib::BGV>()
                  .m(m)
                  .p(p)
                  .r(r)
                  .bits(bits)
                  .gens(gens)
                  .ords(ords)
                  .mvec(mvec)
                  .build()),
      secretKey(context),
      publicKey(
          (secretKey.GenSecKey(), addSome1DMatrices(secretKey), secretKey)),
      ea(context.getEA())
  {}

  virtual void SetUp() override
  {
    helib::setupDebugGlobals(&secretKey, context.shareEA());
  }

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }
};

class TestIO_CKKS : public ::testing::TestWithParam<CKKSParameters>
{
protected:
  const long m;
  const long precision;
  const long bits;

  helib::Context context;
  helib::SecKey secretKey;
  helib::PubKey publicKey;
  const helib::EncryptedArray& ea;

  TestIO_CKKS() :
      m(GetParam().m),
      precision(GetParam().precision),
      bits(GetParam().bits),
      context(helib::ContextBuilder<helib::CKKS>()
                  .m(m)
                  .precision(precision)
                  .bits(bits)
                  .build()),
      secretKey(context),
      publicKey(
          (secretKey.GenSecKey(), addSome1DMatrices(secretKey), secretKey)),
      ea(context.getEA())
  {}

  virtual void SetUp() override
  {
    helib::setupDebugGlobals(&secretKey, context.shareEA());
  }

  virtual void TearDown() override { helib::cleanupDebugGlobals(); }
};

static const json makeExpectedJSON(const helib::Context& context)
{
  // Generators
  long numGens = context.getZMStar().numOfGens();
  std::vector<long> gens(numGens);
  for (long i = 0; i < numGens; ++i) {
    gens[i] = context.getZMStar().ZmStarGen(i);
  }

  // Orders of the generators
  std::vector<long> ords(numGens);
  for (long i = 0; i < numGens; ++i) {
    // Check for good/bad dimension
    if (context.getZMStar().SameOrd(i)) {
      ords[i] = context.getZMStar().OrderOf(i);
    } else { // Bad dimensions are negated
      ords[i] = -context.getZMStar().OrderOf(i);
    }
  }

  // Prime chain
  std::vector<long> qs(context.numPrimes());
  for (std::size_t i = 0; i < qs.size(); ++i) {
    qs[i] = context.ithPrime(i);
  }

  json j = {{"m", context.getM()},
            {"p", context.getP()},
            {"r", context.getR()},
            {"gens", gens},
            {"ords", ords},
            {"stdev", context.getStdev()},
            {"scale", context.getScale()},
            {"smallPrimes", unwrap(context.getSmallPrimes().writeToJSON())},
            {"specialPrimes", unwrap(context.getSpecialPrimes().writeToJSON())},
            {"qs", qs},
            {"digits", writeVectorToJSON(context.getDigits())},
            {"hwt_param", context.getHwt()},
            {"e_param", context.getE()},
            {"ePrime_param", context.getEPrime()},
            {"mvec", context.getRcData().mvec},
            {"build_cache", context.getRcData().build_cache},
            {"alsoThick", context.getRcData().alsoThick}};
  return helib::toTypedJson<helib::Context>(j);
}

TEST(TestIO, TypedJSONLabelsAreTypeNames)
{
  EXPECT_EQ(helib::ContextBuilder<helib::BGV>::typeName, "ContextBuilder");
  EXPECT_EQ(helib::Context::typeName, "Context");
  EXPECT_EQ(helib::Ctxt::typeName, "Ctxt");
  EXPECT_EQ(helib::PubKey::typeName, "PubKey");
  EXPECT_EQ(helib::SecKey::typeName, "SecKey");
  EXPECT_EQ(helib::KeySwitch::typeName, "KeySwitch");
  EXPECT_EQ(helib::Ptxt<helib::BGV>::typeName, "Ptxt");
  EXPECT_EQ(helib::PtxtArray::typeName, "PtxtArray");
}

TEST(TestIO, toTypedJSONWorks)
{
  json jcont = 42;

  EXPECT_EQ(
      helib::toTypedJson<helib::ContextBuilder<helib::BGV>>(jcont).at("type"),
      helib::ContextBuilder<helib::BGV>::typeName);
  EXPECT_EQ(helib::toTypedJson<helib::Context>(jcont).at("type"),
            helib::Context::typeName);
  EXPECT_EQ(helib::toTypedJson<helib::Ctxt>(jcont).at("type"),
            helib::Ctxt::typeName);
  EXPECT_EQ(helib::toTypedJson<helib::PubKey>(jcont).at("type"),
            helib::PubKey::typeName);
  EXPECT_EQ(helib::toTypedJson<helib::SecKey>(jcont).at("type"),
            helib::SecKey::typeName);
  EXPECT_EQ(helib::toTypedJson<helib::KeySwitch>(jcont).at("type"),
            helib::KeySwitch::typeName);
  EXPECT_EQ(helib::toTypedJson<helib::Ptxt<helib::BGV>>(jcont).at("type"),
            helib::Ptxt<helib::BGV>::typeName);
  EXPECT_EQ(helib::toTypedJson<helib::PtxtArray>(jcont).at("type"),
            helib::PtxtArray::typeName);

  EXPECT_EQ(
      helib::toTypedJson<helib::Context>(jcont).at("serializationVersion"),
      helib::jsonSerializationVersion);

  EXPECT_EQ(helib::toTypedJson<helib::Context>(jcont).at("HElibVersion"),
            helib::version::asString);

  EXPECT_EQ(helib::toTypedJson<helib::KeySwitch>(jcont).at("content"), jcont);
}

TEST(TestIO, fromTypedJSONWorks)
{
  json jcont = 42;

  EXPECT_EQ(helib::fromTypedJson<helib::ContextBuilder<helib::BGV>>(
                helib::toTypedJson<helib::ContextBuilder<helib::BGV>>(jcont)),
            jcont);
  EXPECT_EQ(helib::fromTypedJson<helib::Context>(
                helib::toTypedJson<helib::Context>(jcont)),
            jcont);
  EXPECT_EQ(
      helib::fromTypedJson<helib::Ctxt>(helib::toTypedJson<helib::Ctxt>(jcont)),
      jcont);
  EXPECT_EQ(helib::fromTypedJson<helib::PubKey>(
                helib::toTypedJson<helib::PubKey>(jcont)),
            jcont);
  EXPECT_EQ(helib::fromTypedJson<helib::SecKey>(
                helib::toTypedJson<helib::SecKey>(jcont)),
            jcont);
  EXPECT_EQ(helib::fromTypedJson<helib::KeySwitch>(
                helib::toTypedJson<helib::KeySwitch>(jcont)),
            jcont);
  EXPECT_EQ(helib::fromTypedJson<helib::Ptxt<helib::BGV>>(
                helib::toTypedJson<helib::Ptxt<helib::BGV>>(jcont)),
            jcont);
  EXPECT_EQ(helib::fromTypedJson<helib::PtxtArray>(
                helib::toTypedJson<helib::PtxtArray>(jcont)),
            jcont);

  EXPECT_EQ(helib::fromTypedJson<helib::Context>(
                helib::toTypedJson<helib::Context>(jcont)),
            jcont);

  EXPECT_EQ(helib::fromTypedJson<helib::KeySwitch>(
                helib::toTypedJson<helib::KeySwitch>(jcont)),
            jcont);
}

TEST(TestIO, fromTypedJSONThrowsWhenMetadataIsWrong)
{
  json jcont = 42;
  json tj_orig = helib::toTypedJson<helib::Context>(jcont);
  json tj = tj_orig;
  tj.at("type") = "wrong";

  EXPECT_THROW(helib::fromTypedJson<helib::ContextBuilder<helib::BGV>>(tj),
               helib::IOError);
  EXPECT_THROW(helib::fromTypedJson<helib::Context>(tj), helib::IOError);
  EXPECT_THROW(helib::fromTypedJson<helib::Ctxt>(tj), helib::IOError);
  EXPECT_THROW(helib::fromTypedJson<helib::PubKey>(tj), helib::IOError);
  EXPECT_THROW(helib::fromTypedJson<helib::SecKey>(tj), helib::IOError);
  EXPECT_THROW(helib::fromTypedJson<helib::KeySwitch>(tj), helib::IOError);
  EXPECT_THROW(helib::fromTypedJson<helib::Ptxt<helib::BGV>>(tj),
               helib::IOError);
  EXPECT_THROW(helib::fromTypedJson<helib::PtxtArray>(tj), helib::IOError);

  tj = tj_orig;
  std::stringstream ss;
  ss << helib::jsonSerializationVersion << ".wrong";
  tj.at("serializationVersion") = ss.str();
  EXPECT_THROW(helib::fromTypedJson<helib::Context>(tj), helib::IOError);

  ss.str("");
  ss.clear();
  tj = tj_orig;
  ss << helib::version::asString << ".wrong";
  tj.at("HElibVersion") = ss.str();
  EXPECT_THROW(helib::fromTypedJson<helib::Context>(tj), helib::IOError);

  tj = tj_orig;
  EXPECT_EQ(helib::fromTypedJson<helib::Context>(tj), jcont);
}

TEST(TestIO, serializeComplexNumbers)
{
  std::complex<double> cd1 = 0;
  json j1 = cd1;
  EXPECT_EQ(j1.dump(), "[0.0,0.0]");

  std::complex<double> cd2 = 42.5;
  json j2 = cd2;
  EXPECT_EQ(j2.dump(), "[42.5,0.0]");

  std::complex<double> cd3 = {0, 12.2};
  json j3 = cd3;
  EXPECT_EQ(j3.dump(), "[0.0,12.2]");

  std::complex<double> cd4 = {10.9, 9.1};
  json j4 = cd4;
  EXPECT_EQ(j4.dump(), "[10.9,9.1]");
}

TEST(TestIO, deserializeComplexNumbers)
{
  std::stringstream ss1("[0.0,0.0]");
  std::complex<double> cd1 = 0;
  json j1;
  ss1 >> j1;
  EXPECT_EQ(j1.get<std::complex<double>>(), cd1);

  std::stringstream ss2("[42.5,0.0]");
  std::complex<double> cd2 = 42.5;
  json j2;
  ss2 >> j2;
  EXPECT_EQ(j2.get<std::complex<double>>(), cd2);

  std::stringstream ss3("[0.0,12.2]");
  std::complex<double> cd3 = {0, 12.2};
  json j3;
  ss3 >> j3;
  EXPECT_EQ(j3.get<std::complex<double>>(), cd3);

  std::stringstream ss4("[10.9,9.1]");
  std::complex<double> cd4 = {10.9, 9.1};
  json j4;
  ss4 >> j4;
  EXPECT_EQ(j4.get<std::complex<double>>(), cd4);

  std::stringstream ss5("[0]");
  std::complex<double> cd5 = 0;
  json j5;
  ss5 >> j5;
  EXPECT_EQ(j5.get<std::complex<double>>(), cd5);

  std::stringstream ss6("[42.5]");
  std::complex<double> cd6 = 42.5;
  json j6;
  ss6 >> j6;
  EXPECT_EQ(j6.get<std::complex<double>>(), cd6);

  std::stringstream ss7("[42.5]");
  std::complex<double> cd7 = {0, 42.5};
  json j7;
  ss7 >> j7;
  EXPECT_NE(j7.get<std::complex<double>>(), cd7);

  std::stringstream ss8("42.5");
  std::complex<double> cd8 = 42.5;
  json j8;
  ss8 >> j8;
  EXPECT_EQ(j8.get<std::complex<double>>(), cd8);

  std::stringstream ss9("[]");
  std::complex<double> cd9 = 0;
  json j9;
  ss9 >> j9;
  EXPECT_EQ(j9.get<std::complex<double>>(), cd9);

  std::stringstream ss10("[1,2,3]");
  std::complex<double> cd10;
  json j10;
  ss10 >> j10;
  EXPECT_THROW(std::from_json(j10, cd10), helib::IOError);
}

TEST_P(TestIO_BGV, polyModWriteToJSONWorks)
{
  std::stringstream ss1, ss2, ss3;
  helib::PolyMod s1(context.getSlotRing()), s2(context.getSlotRing()),
      s3(context.getSlotRing());
  s1 = 0;
  s2 = 42;
  s3 = {1, 2, 3};
  ss1 << s1;
  ss2 << s2;
  ss3 << s3;

  EXPECT_EQ(ss1.str(), "[0]");
  EXPECT_EQ(ss2.str(), "[0]");
  EXPECT_EQ(ss3.str(), "[1,0,1]");

  helib::PolyMod invalid;
  EXPECT_THROW(ss1 << invalid, helib::LogicError);
}

TEST_P(TestIO_BGV, polyModReadFromJSONWorks)
{
  auto ring = context.getSlotRing();
  std::stringstream ss1("[0]"), ss2("[0]"), ss3("[1,0,1]");
  helib::PolyMod s1(ring), s2(ring), s3(ring), s4(ring);
  s1 = 0;
  s2 = 42;
  s3 = {1, 2, 3};
  EXPECT_EQ(helib::PolyMod::readFromJSON(ss1, ring), s1);
  EXPECT_EQ(helib::PolyMod::readFromJSON(ss2, ring), s2);
  EXPECT_EQ(helib::PolyMod::readFromJSON(ss3, ring), s3);

  std::stringstream ss4;
  ss4 << "[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, "
         "20, 21, 22, 23, 24, 25]";
  EXPECT_THROW(helib::PolyMod::readFromJSON(ss4, ring), helib::IOError);
}

TEST_P(TestIO_BGV, polyModReadJSONPolyModObjects)
{
  auto ring = context.getSlotRing();
  std::stringstream ss1("[0]"), ss2("[0]"), ss3("[1,0,1]");
  helib::PolyMod s1(ring), s2(ring), s3(ring), s4(ring), ds1(ring), ds2(ring),
      ds3(ring), ds4(ring), invalid;
  s1 = 0;
  s2 = 42;
  s3 = {1, 2, 3};
  ds1.readJSON(ss1);
  ds2.readJSON(ss2);
  ds3.readJSON(ss3);
  EXPECT_EQ(ds1, s1);
  EXPECT_EQ(ds2, s2);
  EXPECT_EQ(ds3, s3);

  std::stringstream ss4;
  ss4 << "[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, "
         "20, 21, 22, 23, 24, 25]";
  EXPECT_THROW(ds4.readJSON(ss4), helib::IOError);
  ss4 << "[1]";
  EXPECT_THROW(invalid.readJSON(ss4), helib::LogicError);
}

TEST_P(TestIO_BGV,
       polyModReadFromJSONFunctionThrowsIfMoreEqualsElementsThanDegree)
{
  auto ring = context.getSlotRing();

  std::stringstream ss;
  ss << "[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, "
        "20, 21, 22, 23, 24, 25]";

  EXPECT_THROW(helib::PolyMod::readFromJSON(ss, ring), helib::IOError);

  ss.str("");
  ss.clear();

  json j;
  to_json(j, ring->G);
  ss << j;

  EXPECT_THROW(helib::PolyMod::readFromJSON(ss, ring), helib::IOError);
}

TEST_P(TestIO_BGV, polyModReadFromJSONFunctionThrowsIfFloatingPointCoefficient)
{
  // G = x^2
  auto ring = context.getSlotRing();

  std::stringstream ss;
  ss << "[1.5, 2.2]";

  EXPECT_THROW(helib::PolyMod::readFromJSON(ss, ring), helib::IOError);
}

TEST_P(TestIO_BGV, rpolyModRightShiftOperatorThrowsIfMoreElementsThanDegree)
{
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G;
  // G = x^2
  auto ring = context.getSlotRing();
  helib::PolyMod dest(ring);

  std::stringstream ss;
  ss << "[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, "
        "20, 21, 22, 23, 24, 25]";

  EXPECT_THROW(ss >> dest, helib::IOError);
}

TEST_P(TestIO_BGV, polyModWriteToJSONSerializesCorrectly)
{
  NTL::ZZX G;
  NTL::SetX(G);
  G *= G; // G^2
  auto ring = context.getSlotRing();
  NTL::ZZX data;
  NTL::SetCoeff(data, 0, 2);
  NTL::SetCoeff(data, 1, 1);
  helib::PolyMod poly(data, ring);

  std::stringstream ss;
  poly.writeToJSON(ss);
  std::string expected_string = "[0,1]";

  EXPECT_EQ(ss.str(), expected_string);
}

TEST_P(TestIO_BGV, polyModReadFromJSONDeserializesCorrectly)
{
  auto ring = context.getSlotRing();
  NTL::ZZX data;
  NTL::SetCoeff(data, 0, 3);
  NTL::SetCoeff(data, 1, 1);
  helib::PolyMod expected_result(data, ring);

  std::string string_poly = "[3, 1]";
  std::stringstream str(string_poly);
  helib::PolyMod deserialized_poly = helib::PolyMod::readFromJSON(str, ring);

  EXPECT_EQ(deserialized_poly, expected_result);
}

TEST_P(TestIO_BGV, serializeContextWithStreamOperator)
{
  std::stringstream strm;
  const json expected_json = makeExpectedJSON(context);

  EXPECT_NO_THROW(strm << context);
  EXPECT_EQ(strm.str(), expected_json.dump());
}

TEST_P(TestIO_BGV, deserializingContextDoesNotThrow)
{
  std::stringstream str;

  context.writeToJSON(str);

  EXPECT_NO_THROW(helib::Context::readFromJSON(str));
}

TEST_P(TestIO_BGV, handlingContextDeserializationWithMissingField)
{
  std::stringstream ss;

  ss << context;

  json j;
  ss >> j;

  // Clears the stream
  ss.str("");
  ss.clear();
  j.at("content").erase("m");

  ss << j;
  EXPECT_THROW(helib::Context::readFromJSON(ss), helib::RuntimeError);
}

TEST_P(TestIO_BGV, readContextFromDeserializeCorrectly)
{
  std::stringstream ss;

  ss << context;

  helib::Context deserialized_context = helib::Context::readFromJSON(ss);

  EXPECT_EQ(deserialized_context, context);
}

TEST_P(TestIO_BGV, canPerformOperationWithDeserializedContext)
{
  std::stringstream ss;

  ss << context;

  helib::Context deserialized_context = helib::Context::readFromJSON(ss);

  EXPECT_NO_THROW(helib::PubKey{deserialized_context});
  EXPECT_NO_THROW(helib::SecKey{deserialized_context});

  helib::Ctxt ctxt{helib::PubKey(deserialized_context)};

  EXPECT_NO_THROW(ctxt.square());
  EXPECT_NO_THROW(ctxt += ctxt);
  EXPECT_NO_THROW(ctxt.reLinearize());
  EXPECT_NO_THROW(deserialized_context.getEA().rotate(ctxt, 1));
}

TEST_P(TestIO_BGV, serializeKeysWithStreamOperator)
{
  std::stringstream str;

  EXPECT_NO_THROW(str << publicKey);
  EXPECT_NO_THROW(str << secretKey);
}

TEST_P(TestIO_BGV, deserializeKeysWithStreamOperator)
{
  std::stringstream str;

  str << publicKey;

  EXPECT_NO_THROW(str >> publicKey);

  str.str("");
  str.clear();

  str << secretKey;

  EXPECT_NO_THROW(str >> secretKey);
}

TEST_P(TestIO_BGV, readKeysFromDeserializeCorrectly)
{
  std::stringstream ss;

  ss << publicKey;

  helib::PubKey deserialized_pk(context);
  ss >> deserialized_pk;

  EXPECT_EQ(publicKey, deserialized_pk);

  ss.str("");
  ss.clear();

  ss << secretKey;

  helib::SecKey deserialized_sk = helib::SecKey::readFromJSON(ss, context);

  EXPECT_EQ(secretKey, deserialized_sk);
}

TEST_P(TestIO_BGV, secretKeyOnlyThrowsOnMismatchContext)
{
  std::stringstream str;
  secretKey.writeToJSON(str, /*sk_only=*/true);
  helib::Context new_context(helib::ContextBuilder<helib::BGV>()
                                 .m(41)
                                 .p(p)
                                 .r(r)
                                 .bits(bits)
                                 .gens(gens)
                                 .ords(ords)
                                 .mvec(mvec)
                                 .build());
  helib::SecKey deserialized_sk(new_context);
  EXPECT_THROW(deserialized_sk.readFromJSON(str, new_context, /*sk_only=*/true),
               helib::LogicError);
}

TEST_P(TestIO_BGV, handlingPublicKeyDeserializationWithMissingField)
{
  std::stringstream ss;

  ss << publicKey;

  json j;
  ss >> j;

  // Clears the stream
  ss.str("");
  ss.clear();
  j.at("content").erase("KS_strategy");

  ss << j;
  EXPECT_THROW(helib::PubKey::readFromJSON(ss, context), helib::RuntimeError);
}

TEST_P(TestIO_BGV, handlingSecretKeyDeserializationWithMissingField)
{
  std::stringstream ss;

  ss << secretKey;

  json j;
  ss >> j;

  // Clears the stream
  ss.str("");
  ss.clear();
  j.at("content").erase("PubKey");

  ss << j;
  EXPECT_THROW(helib::SecKey::readFromJSON(ss, context), helib::RuntimeError);
}

TEST_P(TestIO_BGV, canEncryptWithDeserializedPublicKey)
{
  std::stringstream ss;

  ss << publicKey;

  helib::PubKey deserialized_pk = helib::PubKey::readFromJSON(ss, context);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(deserialized_pk);

  EXPECT_NO_THROW(ptxt.encrypt(ctxt));

  decrypted_result.decrypt(ctxt, secretKey);

  EXPECT_EQ(ptxt, decrypted_result);
}

TEST_P(TestIO_BGV, canEncryptWithDeserializedSecretKeyOnly)
{
  std::stringstream ss;

  secretKey.writeToJSON(ss, /*sk_only=*/true);
  helib::SecKey deserialized_sk =
      helib::SecKey::readFromJSON(ss, context, /*sk_only=*/true);
  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(deserialized_sk);

  EXPECT_NO_THROW(ptxt.encrypt(ctxt));

  decrypted_result.decrypt(ctxt, secretKey);

  EXPECT_EQ(ptxt, decrypted_result);
}

TEST_P(TestIO_BGV, canDecryptWithDeserializedSecretKey)
{
  std::stringstream ss;

  ss << secretKey;

  helib::SecKey deserialized_sk = helib::SecKey::readFromJSON(ss, context);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();

  helib::Ctxt ctxt(publicKey);
  ptxt.encrypt(ctxt);

  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt, deserialized_sk));
  EXPECT_EQ(ptxt, decrypted_result);
}

TEST_P(TestIO_BGV, canDecryptWithDeserializedSecretKeyOnly)
{
  std::stringstream ss;

  secretKey.writeToJSON(ss, /*sk_only=*/true);
  helib::SecKey deserialized_sk =
      helib::SecKey::readFromJSON(ss, context, /*sk_only=*/true);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();

  helib::Ctxt ctxt(publicKey);
  ptxt.encrypt(ctxt);

  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt, deserialized_sk));
  EXPECT_EQ(ptxt, decrypted_result);
}

TEST_P(TestIO_BGV, serializeCiphertextWithStreamOperator)
{
  std::stringstream str;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(str << ctxt);
}

TEST_P(TestIO_BGV, deserializeCiphertextWithStreamOperator)
{
  std::stringstream str;
  helib::Ctxt ctxt(publicKey);

  str << ctxt;

  EXPECT_NO_THROW(str >> ctxt);
}

TEST_P(TestIO_BGV, readCiphertextFromDeserializeCorrectly)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  ss << ctxt;

  helib::Ctxt deserialized_ctxt(publicKey);
  ss >> deserialized_ctxt;

  EXPECT_EQ(ctxt, deserialized_ctxt);
}

TEST_P(TestIO_BGV, handlingCiphertextDeserializationWithMissingField)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  ss << ctxt;

  json j;
  ss >> j;

  // Clears the stream
  ss.str("");
  ss.clear();
  j.at("content").erase("noiseBound");

  ss << j;
  EXPECT_THROW(helib::Ctxt::readFromJSON(ss, publicKey), helib::RuntimeError);
}

TEST_P(TestIO_BGV, canPerformOperationsOnDeserializedCiphertext)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  ss << ctxt;
  helib::Ctxt deserialized_ctxt = helib::Ctxt::readFromJSON(ss, publicKey);
  helib::PtxtArray ptxt(ea);

  EXPECT_NO_THROW(deserialized_ctxt *= ctxt);
  EXPECT_NO_THROW(deserialized_ctxt += ctxt);
  EXPECT_NO_THROW(deserialized_ctxt.reLinearize());
  EXPECT_NO_THROW(ea.rotate(deserialized_ctxt, 1));
  EXPECT_NO_THROW(ptxt.decrypt(deserialized_ctxt, secretKey));
}

TEST_P(TestIO_BGV, ptxtWritesDataCorrectlyToOstream)
{
  const long p2r = context.getSlotRing()->p2r;
  const long d = context.getOrdP();
  helib::PolyMod poly(context.getSlotRing());
  std::vector<helib::PolyMod> data(context.getEA().size(), poly);
  std::stringstream ss;
  ss << "[";
  for (long i = 0; i < helib::lsize(data); ++i) {
    NTL::ZZX input;
    NTL::SetCoeff(input, 0, i % p2r);
    if (d != 1) {
      NTL::SetCoeff(input, 1, (i + 2) % p2r);
    }
    data[i] = input;
    ss << data[i] << (i != helib::lsize(data) - 1 ? "," : "");
  }
  ss << "]";
  helib::Ptxt<helib::BGV> ptxt(context, data);
  std::string expected =
      "{\"HElibVersion\":\"" + std::string(helib::version::asString) +
      "\",\"content\":{\"scheme\":\"BGV\","
      "\"slots\":" +
      ss.str() + "},\"serializationVersion\":\"0.0.1\",\"type\":\"Ptxt\"}";
  std::ostringstream os;
  os << ptxt;

  EXPECT_EQ(os.str(), expected);
}

TEST_P(TestIO_BGV, ptxtReadsDataCorrectlyFromIstream)
{
  helib::PolyMod poly(context.getSlotRing());
  std::vector<helib::PolyMod> data(context.getEA().size(), poly);
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {i, i + 2};
  }
  helib::Ptxt<helib::BGV> ptxt(context);
  std::stringstream ss;
  ss << "[";
  for (auto it = data.begin(); it != data.end(); it++) {
    ss << *it;
    if (it != data.end() - 1) {
      ss << ", ";
    }
  }
  ss << "]";

  std::string expected = addHeader(ss.str(), "BGV");
  ss.str(expected);
  ss >> ptxt;

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }
}

TEST_P(TestIO_BGV, ptxtReadsJSONVectorFromIstream)
{
  helib::PolyMod poly(context.getSlotRing());
  std::vector<helib::PolyMod> data(context.getEA().size(), poly);
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {i, i + 2};
  }
  helib::Ptxt<helib::BGV> ptxt(context);
  std::stringstream ss;
  ss << "[";
  for (auto it = data.begin(); it != data.end(); it++) {
    ss << *it;
    if (it != data.end() - 1) {
      ss << ", ";
    }
  }
  ss << "]";

  ss >> ptxt;

  for (std::size_t i = 0; i < ptxt.size(); ++i) {
    EXPECT_EQ(ptxt[i], data[i]);
  }
}

TEST_P(TestIO_BGV, ptxtWriteToJSONFunctionSerializesPtxtCorrectly)
{
  helib::PolyMod poly(context.getSlotRing());
  std::vector<helib::PolyMod> data(context.getEA().size(), poly);
  std::stringstream ptxt_string_stream;
  ptxt_string_stream << "[";
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = 2 * i;
    ptxt_string_stream << "[" << helib::mcMod(2 * i, context.getPPowR()) << "]";
    if (i < helib::lsize(data) - 1)
      ptxt_string_stream << ",";
  }
  ptxt_string_stream << "]";
  helib::Ptxt<helib::BGV> ptxt(context, data);

  std::stringstream ss;
  ptxt.writeToJSON(ss);

  EXPECT_EQ(ss.str(), addHeader(ptxt_string_stream.str(), "BGV"));
}

TEST_P(TestIO_BGV, ptxtReadFromJSONFunctionDeserializesPtxtCorrectly)
{
  helib::PolyMod poly(context.getSlotRing());
  std::vector<helib::PolyMod> data(context.getEA().size(), poly);
  std::stringstream ptxt_string_stream;
  ptxt_string_stream << "[";
  for (long i = 0; i < helib::lsize(data); ++i) {
    NTL::ZZX tmp;
    ptxt_string_stream << "[";
    for (long j = 0; j < context.getOrdP(); ++j) {
      NTL::SetCoeff(tmp, j, j * j);
      ptxt_string_stream << j * j;
      if (j < context.getOrdP() - 1)
        ptxt_string_stream << ",";
    }
    data[i] = tmp;
    ptxt_string_stream << "]";
    if (i < helib::lsize(data) - 1)
      ptxt_string_stream << ",";
  }
  ptxt_string_stream << "]";
  std::string expected = addHeader(ptxt_string_stream.str(), "BGV");
  ptxt_string_stream.str(expected);
  helib::Ptxt<helib::BGV> ptxt(context, data);

  helib::Ptxt<helib::BGV> deserialized_ptxt =
      helib::Ptxt<helib::BGV>::readFromJSON(ptxt_string_stream, context);

  EXPECT_EQ(ptxt, deserialized_ptxt);
}

TEST_P(TestIO_BGV, ptxtReadFromJSONFunctionThrowsIfMoreElementsThanSlots)
{
  helib::PolyMod poly(context.getSlotRing());
  std::vector<helib::PolyMod> data(context.getEA().size() + 1, poly);
  std::stringstream ptxt_string_stream;
  ptxt_string_stream << "[";
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {i, i * i};
    ptxt_string_stream << "[" << i << ", " << i * i << "]";
    if (i < helib::lsize(data) - 1)
      ptxt_string_stream << ",";
  }
  ptxt_string_stream << "]";
  std::string ptxt_string = addHeader(ptxt_string_stream.str(), "BGV");
  ptxt_string_stream.str(ptxt_string);

  EXPECT_THROW(
      helib::Ptxt<helib::BGV>::readFromJSON(ptxt_string_stream, context),
      helib::IOError);
}

TEST_P(TestIO_BGV, ptxtRightShiftOperatorThrowsIfMoreElementsThanSlots)
{
  helib::PolyMod poly(context.getSlotRing());
  std::vector<helib::PolyMod> data(context.getEA().size() + 1, poly);
  std::stringstream ptxt_string_stream;
  ptxt_string_stream << "[";
  for (long i = 0; i < helib::lsize(data); ++i) {
    data[i] = {i, i * i};
    ptxt_string_stream << "[" << i << ", " << i * i << "]";
    if (i < helib::lsize(data) - 1)
      ptxt_string_stream << ", ";
  }
  ptxt_string_stream << "]";
  std::string s = addHeader(ptxt_string_stream.str(), "BGV");
  ptxt_string_stream.str(s);

  helib::Ptxt<helib::BGV> deserialized_ptxt(context);
  EXPECT_THROW(ptxt_string_stream >> deserialized_ptxt, helib::IOError);
}

TEST_P(TestIO_BGV, ptxtWriteToJSONCorrectlyWritesMetdata)
{
  helib::Ptxt<helib::BGV> ptxt(context);
  std::stringstream ss;
  ss << ptxt;
  json j;
  ss >> j;

  EXPECT_TRUE(j.contains("type"));
  EXPECT_EQ(j.at("type").get<std::string>(), helib::Ptxt<helib::BGV>::typeName);

  EXPECT_TRUE(j.contains("serializationVersion"));
  EXPECT_EQ(j.at("serializationVersion").get<std::string>(),
            helib::jsonSerializationVersion);

  EXPECT_TRUE(j.contains("HElibVersion"));
  EXPECT_EQ(j.at("HElibVersion").get<std::string>(), helib::version::asString);

  j = j.at("content");
  EXPECT_TRUE(j.contains("scheme"));
  EXPECT_EQ(j.at("scheme").get<std::string>(), helib::BGV::schemeName);

  EXPECT_TRUE(j.contains("slots"));
  EXPECT_EQ(j.at("slots").size(), ea.size());
}

TEST_P(TestIO_BGV, ptxtReadFromJSONFailsWhenBadMetadata)
{
  helib::PolyMod poly(context.getSlotRing());
  helib::Ptxt<helib::BGV> original_ptxt{context};
  std::stringstream ptxt_string_stream;
  for (long i = 0; i < original_ptxt.lsize(); ++i) {
    original_ptxt[i] = {i, i * i};
  }

  helib::JsonWrapper jw = original_ptxt.writeToJSON();

  helib::Ptxt<helib::BGV> dest_ptxt;

  EXPECT_THROW(dest_ptxt.readJSON(jw), helib::RuntimeError);

  dest_ptxt = helib::Ptxt<helib::BGV>(context);

  json jmod = unwrap(jw);
  jmod.at("type") = "wrong";
  EXPECT_THROW(dest_ptxt.readJSON(helib::wrap(jmod)), helib::IOError);

  jmod = unwrap(jw);
  jmod.at("serializationVersion") = "wrong";
  EXPECT_THROW(dest_ptxt.readJSON(helib::wrap(jmod)), helib::IOError);

  jmod = unwrap(jw);
  jmod.at("HElibVersion") = "wrong";
  EXPECT_THROW(dest_ptxt.readJSON(helib::wrap(jmod)), helib::IOError);

  jmod = unwrap(jw);
  jmod.at("content").at("scheme") = "wrong";
  EXPECT_THROW(dest_ptxt.readJSON(helib::wrap(jmod)), helib::IOError);
}

TEST_P(TestIO_BGV, ptxtReadFromJSONCorrectlyPadsData)
{
  std::stringstream ss(addHeader("[]", "BGV"));
  helib::Ptxt<helib::BGV> deserialized_ptxt(context);
  ss >> deserialized_ptxt;

  EXPECT_EQ(deserialized_ptxt.size(), ea.size());
  for (long i = 0; i < ea.size(); ++i) {
    EXPECT_EQ(deserialized_ptxt[i], 0);
  }

  ss.str("");
  ss.clear();

  std::vector<long> data(context.getEA().size() / 2, 1);
  json j = data;
  ss.str(addHeader(j.dump(), "BGV"));
  ss >> deserialized_ptxt;

  EXPECT_EQ(deserialized_ptxt.size(), ea.size());
  for (long i = 0; i < ea.size(); ++i) {
    if (i < static_cast<long>(data.size())) {
      EXPECT_EQ(deserialized_ptxt[i], 1);
    } else {
      EXPECT_EQ(deserialized_ptxt[i], 0);
    }
  }
}

TEST_P(TestIO_BGV, ptxtReadsManyPtxtsFromStream)
{
  helib::PolyMod poly(context.getSlotRing());
  std::vector<helib::PolyMod> data1(context.getEA().size(), poly);
  std::vector<helib::PolyMod> data2(context.getEA().size(), poly);
  std::vector<helib::PolyMod> data3(context.getEA().size(), poly);
  for (long i = 0; i < helib::lsize(data1); ++i) {
    data1[i] = {i, i + 2};
    data2[i] = {2 * i, 2 * (i + 2)};
    data3[i] = {3 * i, 3 * (i + 2)};
  }
  helib::Ptxt<helib::BGV> ptxt1(context, data1);
  helib::Ptxt<helib::BGV> ptxt2(context, data2);
  helib::Ptxt<helib::BGV> ptxt3(context, data3);

  std::stringstream ss;
  ss << ptxt1 << std::endl;
  ss << ptxt2 << std::endl;
  ss << ptxt3 << std::endl;

  helib::Ptxt<helib::BGV> deserialized1(context);
  helib::Ptxt<helib::BGV> deserialized2(context);
  helib::Ptxt<helib::BGV> deserialized3(context);
  ss >> deserialized1;
  ss >> deserialized2;
  ss >> deserialized3;

  EXPECT_EQ(ptxt1, deserialized1);
  EXPECT_EQ(ptxt2, deserialized2);
  EXPECT_EQ(ptxt3, deserialized3);
}

// This test should be moved to a set of PtxtArray tests
TEST(TestPtxtArray, ptxtArrayLoadsDataWithLessSlotsThanMax)
{
  // nslots = 18
  helib::Context context =
      helib::ContextBuilder<helib::BGV>().m(127).p(2).build();

  std::vector<long> data(10, 1);
  helib::PtxtArray pa(context);

  EXPECT_NO_THROW(pa.store(data));
}

TEST_P(TestIO_BGV, ptxtArrayWritesDataCorrectlyToOstream)
{
  const long p2r = context.getSlotRing()->p2r;
  const long d = context.getOrdP();

  auto [data, expectedString] = createData(context.getNSlots(), p2r, d);

  std::string expected = addHeader(expectedString, "BGV", "PtxtArray");

  helib::PtxtArray pa(context, data);
  std::ostringstream os;
  os << pa;

  EXPECT_EQ(os.str(), expected);
}

TEST_P(TestIO_BGV, ptxtArrayReadsDataCorrectlyFromIstream)
{
  const long p2r = context.getSlotRing()->p2r;
  const long d = context.getOrdP();

  auto [data, expectedString] = createData(context.getNSlots(), p2r, d);

  helib::PtxtArray pa(context);
  std::string expected = addHeader(expectedString, "BGV", "PtxtArray");
  std::stringstream ss(expectedString);
  ss >> pa;

  EXPECT_EQ(pa, helib::PtxtArray(context, data));
}

TEST_P(TestIO_BGV, ptxtArrayReadsJSONVectorFromIstream)
{
  const long p2r = context.getSlotRing()->p2r;
  const long d = context.getOrdP();

  auto [data, expectedString] = createData(context.getNSlots(), p2r, d);

  helib::PtxtArray pa(context);
  std::stringstream ss(expectedString);
  ss >> pa;

  EXPECT_EQ(pa, helib::PtxtArray(context, data));
}

TEST_P(TestIO_BGV, ptxtArrayWriteToJSONFunctionSerializesPtxtCorrectly)
{
  const long p2r = context.getSlotRing()->p2r;
  const long d = context.getOrdP();

  auto [data, expectedString] = createData(context.getNSlots(), p2r, d);

  helib::PtxtArray pa(context, data);

  std::stringstream ss;
  pa.writeToJSON(ss);

  EXPECT_EQ(ss.str(), addHeader(expectedString, "BGV", "PtxtArray"));
}

TEST_P(TestIO_BGV, ptxtArrayReadFromJSONFunctionDeserializesPtxtCorrectly)
{
  const long p2r = context.getSlotRing()->p2r;
  const long d = context.getOrdP();

  auto [data, expectedString] = createData(context.getNSlots(), p2r, d);

  std::string expected = addHeader(expectedString, "BGV", "PtxtArray");
  std::stringstream ss(expectedString);

  helib::PtxtArray pa(context, data);
  helib::PtxtArray deserialized_pa =
      helib::PtxtArray::readFromJSON(ss, context);

  EXPECT_EQ(pa, deserialized_pa);
}

TEST_P(TestIO_BGV, ptxtArrayReadFromJSONFunctionThrowsIfMoreElementsThanSlots)
{
  const long p2r = context.getSlotRing()->p2r;
  const long d = context.getOrdP();
  std::string expectedString;

  std::tie(std::ignore, expectedString) =
      createData(context.getNSlots() + 1, p2r, d);

  std::stringstream ss(addHeader(expectedString, "BGV", "PtxtArray"));

  EXPECT_THROW(helib::PtxtArray::readFromJSON(ss, context), helib::IOError);
}

TEST_P(TestIO_BGV, ptxtArrayRightShiftOperatorThrowsIfMoreElementsThanSlots)
{
  const long p2r = context.getSlotRing()->p2r;
  const long d = context.getOrdP();
  std::string expectedString;

  std::tie(std::ignore, expectedString) =
      createData(context.getNSlots() + 1, p2r, d);

  std::stringstream ss(addHeader(expectedString, "BGV", "PtxtArray"));
  helib::PtxtArray pa(context);

  EXPECT_THROW(ss >> pa, helib::IOError);
}

TEST_P(TestIO_BGV, ptxtArrayWriteToJSONCorrectlyWritesMetdata)
{
  helib::PtxtArray pa(context);
  std::stringstream ss;
  ss << pa;
  json j;
  ss >> j;

  EXPECT_TRUE(j.contains("type"));
  EXPECT_EQ(j.at("type").get<std::string>(), helib::PtxtArray::typeName);

  EXPECT_TRUE(j.contains("serializationVersion"));
  EXPECT_EQ(j.at("serializationVersion").get<std::string>(),
            helib::jsonSerializationVersion);

  EXPECT_TRUE(j.contains("HElibVersion"));
  EXPECT_EQ(j.at("HElibVersion").get<std::string>(), helib::version::asString);

  j = j.at("content");
  EXPECT_TRUE(j.contains("scheme"));
  EXPECT_EQ(j.at("scheme").get<std::string>(), helib::BGV::schemeName);

  EXPECT_TRUE(j.contains("slots"));
  EXPECT_EQ(j.at("slots").size(), ea.size());
}

TEST_P(TestIO_BGV, ptxtArrayReadFromJSONFailsWhenBadMetadata)
{
  const long p2r = context.getSlotRing()->p2r;
  const long d = context.getOrdP();
  std::vector<NTL::ZZX> data;

  std::tie(data, std::ignore) = createData(context.getNSlots(), p2r, d);

  helib::PtxtArray pa(context, data);
  helib::JsonWrapper jw = pa.writeToJSON();

  helib::PtxtArray destPa(context);
  json jmod = unwrap(jw);
  jmod.at("type") = "wrong";
  EXPECT_THROW(destPa.readJSON(helib::wrap(jmod)), helib::IOError);

  jmod = unwrap(jw);
  jmod.at("serializationVersion") = "wrong";
  EXPECT_THROW(destPa.readJSON(helib::wrap(jmod)), helib::IOError);

  jmod = unwrap(jw);
  jmod.at("HElibVersion") = "wrong";
  EXPECT_THROW(destPa.readJSON(helib::wrap(jmod)), helib::IOError);

  jmod = unwrap(jw);
  jmod.at("content") = "wrong";
  EXPECT_THROW(destPa.readJSON(helib::wrap(jmod)), helib::IOError);
}

TEST_P(TestIO_BGV, ptxtArrayReadFromJSONCorrectlyPadsData)
{
  std::stringstream ss(addHeader("[]", "BGV", "PtxtArray"));
  helib::PtxtArray pa(context);
  ss >> pa;

  EXPECT_EQ(pa.size(), ea.size());

  ss.str("");
  ss.clear();

  std::vector<long> data(context.getNSlots() / 2, 1);
  json j = data;
  ss.str(addHeader(j.dump(), "BGV", "PtxtArray"));
  ss >> pa;

  EXPECT_EQ(pa.size(), ea.size());
  std::vector<long> deserialized;
  pa.store(deserialized);

  EXPECT_EQ(pa.size(), deserialized.size());

  for (long i = 0; i < pa.size(); ++i) {
    if (i < static_cast<long>(data.size())) {
      EXPECT_EQ(deserialized[i], 1);
    } else {
      EXPECT_EQ(deserialized[i], {});
    }
  }
}

TEST_P(TestIO_BGV, ptxtArrayReadsManyPtxtsFromStream)
{
  helib::PtxtArray pa1(context);
  helib::PtxtArray pa2(context);
  helib::PtxtArray pa3(context);
  pa1.random();
  pa2.random();
  pa3.random();

  std::stringstream ss;
  ss << pa1 << std::endl;
  ss << pa2 << std::endl;
  ss << pa3 << std::endl;

  helib::PtxtArray deserialized1(context);
  helib::PtxtArray deserialized2(context);
  helib::PtxtArray deserialized3(context);
  ss >> deserialized1;
  ss >> deserialized2;
  ss >> deserialized3;

  EXPECT_EQ(pa1, deserialized1);
  EXPECT_EQ(pa2, deserialized2);
  EXPECT_EQ(pa3, deserialized3);
}

TEST_P(TestIO_CKKS, serializeContextWithStreamOperator)
{
  std::stringstream strm;
  const json expected_json = makeExpectedJSON(context);

  EXPECT_NO_THROW(strm << context);
  EXPECT_EQ(strm.str(), expected_json.dump());
}

TEST_P(TestIO_CKKS, deserializingContextDoesNotThrow)
{
  std::stringstream str;

  context.writeToJSON(str);

  EXPECT_NO_THROW(helib::Context::readFromJSON(str));
}

TEST_P(TestIO_CKKS, handlingContextDeserializationWithMissingField)
{
  std::stringstream ss;

  ss << context;

  json j;
  ss >> j;

  // Clears the stream
  ss.str("");
  ss.clear();
  j.at("content").erase("m");

  ss << j;
  EXPECT_THROW(helib::Context::readFromJSON(ss), helib::RuntimeError);
}

TEST_P(TestIO_CKKS, canPerformOperationWithDeserializedContext)
{
  std::stringstream ss;

  ss << context;

  helib::Context deserialized_context = helib::Context::readFromJSON(ss);

  EXPECT_NO_THROW(helib::PubKey{deserialized_context});
  EXPECT_NO_THROW(helib::SecKey{deserialized_context});

  helib::Ctxt ctxt{helib::PubKey(deserialized_context)};

  EXPECT_NO_THROW(ctxt.square());
  EXPECT_NO_THROW(ctxt += ctxt);
  EXPECT_NO_THROW(ctxt.reLinearize());
  EXPECT_NO_THROW(deserialized_context.getEA().rotate(ctxt, 1));
}

TEST_P(TestIO_CKKS, readContextFromDeserializeCorrectly)
{
  std::stringstream ss;

  ss << context;

  helib::Context deserialized_context = helib::Context::readFromJSON(ss);

  EXPECT_EQ(deserialized_context, context);
}

TEST_P(TestIO_CKKS, serializeKeysWithStreamOperator)
{
  std::stringstream str;

  EXPECT_NO_THROW(str << publicKey);
  EXPECT_NO_THROW(str << secretKey);
}

TEST_P(TestIO_CKKS, deserializeKeysWithStreamOperator)
{
  std::stringstream str;

  str << publicKey;

  EXPECT_NO_THROW(str >> publicKey);

  str.str("");
  str.clear();

  str << secretKey;

  EXPECT_NO_THROW(str >> secretKey);
}

TEST_P(TestIO_CKKS, readKeysFromDeserializeCorrectly)
{
  std::stringstream ss;

  ss << publicKey;

  helib::PubKey deserialized_pk(context);
  ss >> deserialized_pk;

  EXPECT_EQ(publicKey, deserialized_pk);

  ss.str("");
  ss.clear();

  ss << secretKey;

  helib::SecKey deserialized_sk = helib::SecKey::readFromJSON(ss, context);

  EXPECT_EQ(secretKey, deserialized_sk);
}

TEST_P(TestIO_CKKS, secretKeyOnlyThrowsOnMismatchContext)
{
  std::stringstream str;
  secretKey.writeToJSON(str, /*sk_only=*/true);
  helib::Context new_context(helib::ContextBuilder<helib::CKKS>()
                                 .m(32)
                                 .precision(precision)
                                 .bits(bits)
                                 .build());
  helib::SecKey deserialized_sk(new_context);
  EXPECT_THROW(deserialized_sk.readFromJSON(str, new_context, /*sk_only=*/true),
               helib::LogicError);
}

TEST_P(TestIO_CKKS, handlingPublicKeyDeserializationWithMissingField)
{
  std::stringstream ss;

  ss << publicKey;

  json j;
  ss >> j;

  // Clears the stream
  ss.str("");
  ss.clear();
  j.at("content").erase("KS_strategy");

  ss << j;
  EXPECT_THROW(helib::PubKey::readFromJSON(ss, context), helib::RuntimeError);
}

TEST_P(TestIO_CKKS, handlingSecretKeyDeserializationWithMissingField)
{
  std::stringstream ss;

  ss << secretKey;

  json j;
  ss >> j;

  // Clears the stream
  ss.str("");
  ss.clear();
  j.at("content").erase("PubKey");

  ss << j;
  EXPECT_THROW(helib::SecKey::readFromJSON(ss, context), helib::RuntimeError);
}

TEST_P(TestIO_CKKS, canEncryptWithDeserializedPublicKey)
{
  std::stringstream ss;

  ss << publicKey;

  helib::PubKey deserialized_pk = helib::PubKey::readFromJSON(ss, context);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();

  helib::Ctxt ctxt(deserialized_pk);

  EXPECT_NO_THROW(ptxt.encrypt(ctxt));

  decrypted_result.decrypt(ctxt, secretKey);

  EXPECT_EQ(ptxt, helib::Approx(decrypted_result));
}

TEST_P(TestIO_CKKS, canEncryptWithDeserializedSecretKeyOnly)
{
  std::stringstream ss;

  secretKey.writeToJSON(ss, /*sk_only=*/true);
  helib::SecKey deserialized_sk =
      helib::SecKey::readFromJSON(ss, context, /*sk_only=*/true);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();

  helib::Ctxt ctxt(deserialized_sk);
  EXPECT_NO_THROW(ptxt.encrypt(ctxt));

  decrypted_result.decrypt(ctxt, secretKey);
  EXPECT_EQ(ptxt, helib::Approx(decrypted_result));
}

TEST_P(TestIO_CKKS, canDecryptWithDeserializedSecretKey)
{
  std::stringstream ss;

  ss << secretKey;

  helib::SecKey deserialized_sk = helib::SecKey::readFromJSON(ss, context);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();

  helib::Ctxt ctxt(publicKey);
  ptxt.encrypt(ctxt);

  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt, deserialized_sk));
  EXPECT_EQ(ptxt, helib::Approx(decrypted_result));
}

TEST_P(TestIO_CKKS, canDecryptWithDeserializedSecretKeyOnly)
{
  std::stringstream ss;

  secretKey.writeToJSON(ss, /*sk_only=*/true);
  helib::SecKey deserialized_sk =
      helib::SecKey::readFromJSON(ss, context, /*sk_only=*/true);
  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();

  helib::Ctxt ctxt(publicKey);
  ptxt.encrypt(ctxt);

  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt, deserialized_sk));
  EXPECT_EQ(ptxt, helib::Approx(decrypted_result));
}

TEST_P(TestIO_CKKS, serializeCiphertextWithStreamOperator)
{
  std::stringstream str;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(str << ctxt);
}

TEST_P(TestIO_CKKS, deserializeCiphertextWithStreamOperator)
{
  std::stringstream str;
  helib::Ctxt ctxt(publicKey);

  str << ctxt;

  EXPECT_NO_THROW(str >> ctxt);
}

TEST_P(TestIO_CKKS, readCiphertextFromDeserializeCorrectly)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  ss << ctxt;

  helib::Ctxt deserialized_ctxt(publicKey);
  ss >> deserialized_ctxt;

  EXPECT_EQ(ctxt, deserialized_ctxt);
}

TEST_P(TestIO_CKKS, handlingCiphertextDeserializationWithMissingField)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  ss << ctxt;

  json j;
  ss >> j;

  // Clears the stream
  ss.str("");
  ss.clear();
  j.at("content").erase("noiseBound");

  ss << j;
  EXPECT_THROW(helib::Ctxt::readFromJSON(ss, publicKey), helib::RuntimeError);
}

TEST_P(TestIO_CKKS, canPerformOperationsOnDeserializedCiphertext)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  ss << ctxt;
  helib::Ctxt deserialized_ctxt = helib::Ctxt::readFromJSON(ss, publicKey);
  helib::PtxtArray ptxt(ea);

  EXPECT_NO_THROW(deserialized_ctxt *= ctxt);
  EXPECT_NO_THROW(deserialized_ctxt += ctxt);
  EXPECT_NO_THROW(deserialized_ctxt.reLinearize());
  EXPECT_NO_THROW(ea.rotate(deserialized_ctxt, 1));
  EXPECT_NO_THROW(ptxt.decrypt(deserialized_ctxt, secretKey));
}

TEST_P(TestIO_CKKS, ptxtWritesDataCorrectlyToOstream)
{
  std::vector<std::complex<double>> data(context.getEA().size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  json j = data;
  std::string expected = addHeader(j.dump(), "CKKS");
  std::ostringstream os;
  os << ptxt;

  EXPECT_EQ(os.str(), expected);
}

TEST_P(TestIO_CKKS, ptxtReadsDataCorrectlyFromIstream)
{
  std::vector<std::complex<double>> data(context.getEA().size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context);
  std::stringstream ss;
  json j = data;
  ss << addHeader(j.dump(), "CKKS");
  std::istringstream is(ss.str());
  is >> ptxt;

  COMPARE_CXDOUBLE_VECS(ptxt, data);
}

TEST_P(TestIO_CKKS, ptxtReadsSquareBracketsDataCorrectly)
{
  std::vector<std::complex<double>> data(context.getEA().size());
  std::stringstream ss;

  ss << "[" << std::setprecision(std::numeric_limits<double>::max_digits10);
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
    ss << "[" << data[i].real() << "," << data[i].imag() << "]"
       << (i != data.size() - 1 ? "," : "");
  }
  ss << "]";

  helib::Ptxt<helib::CKKS> ptxt(context);
  std::istringstream is(ss.str());
  is >> ptxt;

  COMPARE_CXDOUBLE_VECS(ptxt, data);
}

TEST_P(TestIO_CKKS, ptxtWriteToJSONSerializesPtxtCorrectly)
{
  std::vector<std::complex<double>> data(context.getEA().size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context, data);
  json j = data;
  std::string expected = addHeader(j.dump(), "CKKS");
  std::stringstream serialized_ptxt;
  ptxt.writeToJSON(serialized_ptxt);

  EXPECT_EQ(serialized_ptxt.str(), expected);
}

TEST_P(TestIO_CKKS, ptxtRightShiftDeserializeCorrectly)
{
  std::vector<std::complex<double>> data(context.getEA().size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context);
  json j = data;
  std::string expected = addHeader(j.dump(), "CKKS");
  std::istringstream is(expected);
  is >> ptxt;

  COMPARE_CXDOUBLE_VECS(ptxt, data);
}

TEST_P(TestIO_CKKS, ptxtReadsJSONVectorFromIstream)
{
  std::vector<std::complex<double>> data(context.getEA().size());
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context);
  json j = data;
  std::istringstream is(j.dump());
  is >> ptxt;

  COMPARE_CXDOUBLE_VECS(ptxt, data);
}

TEST_P(TestIO_CKKS, ptxtReadFromJSONThrowsIfMoreElementsThanSlots)
{
  std::vector<std::complex<double>> data(context.getEA().size() + 1);
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  std::stringstream ss;
  ss << "[";
  ss << std::setprecision(std::numeric_limits<double>::digits10);
  for (auto it = data.begin(); it != data.end(); it++) {
    ss << "[" << it->real() << ", " << it->imag() << "]";
    if (it != data.end() - 1) {
      ss << ", ";
    }
  }
  ss << "]";
  std::istringstream is(ss.str());
  EXPECT_THROW(helib::Ptxt<helib::CKKS>::readFromJSON(is, context),
               helib::IOError);
}

TEST_P(TestIO_CKKS, ptxtRightShiftOperatorThrowsIfMoreElementsThanSlots)
{
  std::vector<std::complex<double>> data(context.getEA().size() + 1);
  for (std::size_t i = 0; i < data.size(); ++i) {
    data[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
  }
  helib::Ptxt<helib::CKKS> ptxt(context);
  std::stringstream ss;
  ss << "[";
  ss << std::setprecision(std::numeric_limits<double>::digits10);
  for (auto it = data.begin(); it != data.end(); it++) {
    ss << "[" << it->real() << ", " << it->imag() << "]";
    if (it != data.end() - 1) {
      ss << ", ";
    }
  }
  ss << "]";
  std::istringstream is(addHeader(ss.str(), "CKKS"));
  EXPECT_THROW(is >> ptxt, helib::IOError);
}

TEST_P(TestIO_CKKS, ptxtReadsManyPtxtsFromStream)
{
  std::vector<std::complex<double>> data1(context.getEA().size());
  std::vector<std::complex<double>> data2(context.getEA().size());
  std::vector<std::complex<double>> data3(context.getEA().size());
  for (std::size_t i = 0; i < data1.size(); ++i) {
    data1[i] = {(i * i) / 10.0, (i * i * i) / 7.5};
    data2[i] = data1[i] * 2.0;
    data3[i] = data1[i] * 3.5;
  }
  helib::Ptxt<helib::CKKS> ptxt1(context, data1);
  helib::Ptxt<helib::CKKS> ptxt2(context, data2);
  helib::Ptxt<helib::CKKS> ptxt3(context, data3);

  std::stringstream ss;
  ss << ptxt1 << std::endl;
  ss << ptxt2 << std::endl;
  ss << ptxt3 << std::endl;

  helib::Ptxt<helib::CKKS> deserialized1(context);
  helib::Ptxt<helib::CKKS> deserialized2(context);
  helib::Ptxt<helib::CKKS> deserialized3(context);
  ss >> deserialized1;
  ss >> deserialized2;
  ss >> deserialized3;

  COMPARE_CXDOUBLE_VECS(ptxt1, deserialized1);
  COMPARE_CXDOUBLE_VECS(ptxt2, deserialized2);
  COMPARE_CXDOUBLE_VECS(ptxt3, deserialized3);
}

TEST_P(TestIO_CKKS, ptxtArrayWritesDataCorrectlyToOstream)
{
  auto data = createData(context.getNSlots());
  json j = data;
  std::string expected = addHeader(j.dump(), "CKKS", "PtxtArray");

  helib::PtxtArray pa(context, data);
  std::ostringstream os;
  os << pa;

  EXPECT_EQ(os.str(), expected);
}

TEST_P(TestIO_CKKS, ptxtArrayReadsDataCorrectlyFromIstream)
{
  auto data = createData(context.getNSlots());
  helib::PtxtArray pa(context);
  json j = data;
  std::stringstream ss;
  ss << addHeader(j.dump(), "CKKS", "PtxtArray");
  std::istringstream is(ss.str());
  is >> pa;

  EXPECT_EQ(pa, helib::PtxtArray(context, data));
}

TEST_P(TestIO_CKKS, ptxtArrayReadsSquareBracketsDataCorrectly)
{
  auto data = createData(context.getNSlots());

  std::stringstream ss;
  ss << "[" << std::setprecision(std::numeric_limits<double>::digits10);
  for (auto it = data.begin(); it != data.end(); ++it) {
    ss << "[" << it->real() << "," << it->imag() << "]"
       << (it != data.end() - 1 ? "," : "");
  }
  ss << "]";

  helib::PtxtArray pa(context);
  std::istringstream is(ss.str());
  is >> pa;

  EXPECT_EQ(pa, helib::PtxtArray(context, data));
}

TEST_P(TestIO_CKKS, ptxtArrayWriteToJSONSerializesPtxtCorrectly)
{
  auto data = createData(context.getNSlots());
  helib::PtxtArray pa(context, data);
  json j = data;
  std::string expected = addHeader(j.dump(), "CKKS", "PtxtArray");
  std::stringstream ss;
  pa.writeToJSON(ss);

  EXPECT_EQ(ss.str(), expected);
}

TEST_P(TestIO_CKKS, ptxtArrayRightShiftDeserializeCorrectly)
{
  auto data = createData(context.getNSlots());
  helib::PtxtArray pa(context);
  json j = data;
  std::istringstream is(addHeader(j.dump(), "CKKS", "PtxtArray"));
  is >> pa;

  EXPECT_EQ(pa, helib::PtxtArray(context, data));
}

TEST_P(TestIO_CKKS, ptxtArrayReadsJSONVectorFromIstream)
{
  auto data = createData(context.getNSlots());
  helib::PtxtArray pa(context);
  json j = data;
  std::istringstream is(j.dump());
  is >> pa;

  EXPECT_EQ(pa, helib::PtxtArray(context, data));
}

TEST_P(TestIO_CKKS, ptxtArrayReadFromJSONThrowsIfMoreElementsThanSlots)
{
  auto data = createData(context.getNSlots() + 1);
  json j = data;
  std::istringstream is(j.dump());

  EXPECT_THROW(helib::PtxtArray::readFromJSON(is, context), helib::IOError);
}

TEST_P(TestIO_CKKS, ptxtArrayRightShiftOperatorThrowsIfMoreElementsThanSlots)
{
  auto data = createData(context.getNSlots() + 1);
  helib::PtxtArray pa(context);
  json j = data;
  std::istringstream is(j.dump());

  EXPECT_THROW(is >> pa, helib::IOError);
}

TEST_P(TestIO_CKKS, ptxtArrayReadsManyPtxtsFromStream)
{
  helib::PtxtArray pa1(context);
  helib::PtxtArray pa2(context);
  helib::PtxtArray pa3(context);
  pa1.random();
  pa2.random();
  pa3.random();

  std::stringstream ss;
  ss << pa1 << std::endl;
  ss << pa2 << std::endl;
  ss << pa3 << std::endl;

  helib::PtxtArray deserialized1(context);
  helib::PtxtArray deserialized2(context);
  helib::PtxtArray deserialized3(context);
  ss >> deserialized1;
  ss >> deserialized2;
  ss >> deserialized3;

  EXPECT_EQ(pa1, deserialized1);
  EXPECT_EQ(pa2, deserialized2);
  EXPECT_EQ(pa3, deserialized3);
}

TEST_P(TestIO_BGV, contextBuilderLogsCorrectly)
{
  long c = 3;
  std::vector<long> gens = {3};
  std::vector<long> ords = {-2};
  std::vector<long> mvec = {3};
  bool buildModChainFlag = false;
  long bits = 2;
  long skHwt = 64;
  long resolution = 1;
  long bitsInSpecialPrimes = 15;
  bool bootstrappableFlag = true;
  bool buildCacheFlag = true;
  bool thickFlag = true;

  // clang-format off
  auto cb = helib::ContextBuilder<helib::BGV>()
                          .m(m)
                          .p(p)
                          .r(r)
                          .c(c)
                          .gens(gens)
                          .ords(ords)
                          .buildModChain(buildModChainFlag)
                          .bits(bits)
                          .skHwt(skHwt)
                          .resolution(resolution)
                          .bitsInSpecialPrimes(bitsInSpecialPrimes)
                          .bootstrappable(bootstrappableFlag)
                          .mvec(mvec)
                          .buildCache(buildCacheFlag)
                          .thickboot();
  // clang-format off

  std::stringstream ss;
  json actual_json;
  ss << cb;
  ss >> actual_json; 
  
  // Content only
  const json expected_json = { 
                         { "scheme", "bgv" },
                         { "m", m },
                         { "p", p },
                         { "r", r },
                         { "c", c },
                         { "gens", gens },
                         { "ords", ords },
                         { "buildModChainFlag", buildModChainFlag },
                         { "bits", bits },
                         { "skHwt", skHwt },
                         { "resolution", resolution },
                         { "bitsInSpecialPrimes", bitsInSpecialPrimes },
                         { "bootstrappableFlag", bootstrappableFlag },
                         { "mvec", mvec },
                         { "buildCacheFlag", buildCacheFlag },
                         { "thickFlag", thickFlag }
                      };

  EXPECT_EQ(actual_json.at("content"), expected_json);
}

TEST_P(TestIO_CKKS, contextBuilderLogsCorrectly)
{
  long c = 3;
  std::vector<long> gens = {3};
  std::vector<long> ords = {-2};
  bool buildModChainFlag = false;
  long bits = 2;
  long skHwt = 64;
  long resolution = 1;
  long bitsInSpecialPrimes = 15;

  // clang-format off
  auto cb = helib::ContextBuilder<helib::CKKS>()
                          .m(m)
                          .precision(precision)
                          .c(c)
                          .gens(gens)
                          .ords(ords)
                          .buildModChain(buildModChainFlag)
                          .bits(bits)
                          .skHwt(skHwt)
                          .resolution(resolution)
                          .bitsInSpecialPrimes(bitsInSpecialPrimes);
  // clang-format off

  std::stringstream ss;
  json actual_json;
  ss << cb;
  ss >> actual_json; 

  // Content only
  const json expected_json = {
                          { "scheme", "ckks"},
                          { "m", m },
                          { "precision", precision },
                          { "c", c },
                          { "gens", gens },
                          { "ords", ords },
                          { "buildModChainFlag", buildModChainFlag },
                          { "bits", bits },
                          { "skHwt", skHwt },
                          { "resolution", resolution },
                          { "bitsInSpecialPrimes", bitsInSpecialPrimes } 
                       };

  EXPECT_EQ(actual_json.at("content"), expected_json);
}

INSTANTIATE_TEST_SUITE_P(Parameters,
                         TestIO_BGV,
                         ::testing::Values(BGVParameters(/*m=*/45,
                                                         /*p=*/2,
                                                         /*r=*/1,
                                                         /*bits=*/30,
                                                         /*gens=*/{},
                                                         /*ords=*/{},
                                                         /*mvec=*/{})));

INSTANTIATE_TEST_SUITE_P(Parameters,
                         TestIO_CKKS,
                         ::testing::Values(CKKSParameters(/*m=*/64,
                                                          /*precision=*/30,
                                                          /*bits=*/30)));

} // namespace
