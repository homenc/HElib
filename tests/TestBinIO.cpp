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

#include <cmath> // isinf
#include <cstring>
#include <sstream>
#include <helib/helib.h>
#include <helib/debugging.h>
#include <binio.h>

#include "test_common.h"
#include "gtest/gtest.h"
// Several of these tests relate to writing objects to stream and then reading
// them back without using test suite parameters, particularly using a
// "deserialized context". This causes challenges when validating against test
// suite parameters and variables, as the underlying contexts have different
// addresses. To accomodate this, in some places we will validate by writing
// test suite parameters and variables to stream and then reading them back
// using variables local to the specific test.
namespace {

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

class TestBinIO_BGV : public ::testing::TestWithParam<BGVParameters>
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

  TestBinIO_BGV() :
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

class TestBinIO_CKKS : public ::testing::TestWithParam<CKKSParameters>
{
protected:
  const long m;
  const long precision;
  const long bits;

  helib::Context context;
  helib::SecKey secretKey;
  helib::PubKey publicKey;
  const helib::EncryptedArray& ea;

  TestBinIO_CKKS() :
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

// Helper function for converting std::array<char> to string
std::string eyeCatcherToStr(
    const std::array<char, helib::EyeCatcher::SIZE>& eyeCatcher)
{
  return std::string(eyeCatcher.begin(), eyeCatcher.end());
}

TEST(TestBinIO, headerSizeIs24bytes)
{
  EXPECT_EQ(sizeof(helib::SerializeHeader<void>), 24);
}

TEST(TestBinIO, headerForContext)
{
  helib::SerializeHeader<helib::Context> header;

  EXPECT_TRUE(header.beginCatcher == helib::EyeCatcher::HEADER_BEGIN);
  EXPECT_TRUE(header.endCatcher == helib::EyeCatcher::HEADER_END);
  EXPECT_TRUE(header.version == helib::Binio::VERSION_0_0_1_0);
  EXPECT_EQ(header.structId, 5);
}

TEST(TestBinIO, headerEquals)
{
  helib::SerializeHeader<helib::Context> header1;
  helib::SerializeHeader<helib::Context> header2;

  EXPECT_TRUE(!memcmp(&header1, &header2, sizeof(header1)));
}

TEST(TestBinIO, headerSerializationDeserialization)
{
  helib::SerializeHeader<helib::Context> header;

  std::stringstream ss;

  header.writeTo(ss);

  auto DeserialisedHeader =
      helib::SerializeHeader<helib::Context>::readFrom(ss);

  EXPECT_TRUE(!memcmp(&header, &DeserialisedHeader, sizeof(header)));
}

TEST_P(TestBinIO_BGV, singleFunctionSerialization)
{
  std::stringstream str;

  EXPECT_NO_THROW(context.writeTo(str));
}

TEST_P(TestBinIO_BGV, singleFunctionDeserialization)
{
  std::stringstream str;

  context.writeTo(str);

  EXPECT_NO_THROW(context.readFrom(str));
}

TEST_P(TestBinIO_BGV, throwsWhenPreContextEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(context.writeTo(ss));

  // Delete pre-context eye catcher
  std::string s = ss.str();
  std::size_t pos = s.find(eyeCatcherToStr(helib::EyeCatcher::CONTEXT_BEGIN));
  s.erase(pos, pos + helib::EyeCatcher::SIZE);
  ss.str(s);

  EXPECT_THROW(context.readFrom(ss), helib::IOError);
}

TEST_P(TestBinIO_BGV, throwsWhenPostContextEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(context.writeTo(ss));

  // Delete post-context eye catcher
  std::string s = ss.str();
  s.erase(s.find(eyeCatcherToStr(helib::EyeCatcher::CONTEXT_END)),
          s.size() - 1);
  ss.str(s);

  EXPECT_THROW(context.readFrom(ss), helib::IOError);
}

TEST_P(TestBinIO_BGV, readContextFromDeserializeCorrectly)
{
  std::stringstream str;

  EXPECT_NO_THROW(context.writeTo(str));

  helib::Context deserialized_context = helib::Context::readFrom(str);

  EXPECT_EQ(context, deserialized_context);
}

TEST_P(TestBinIO_BGV, readContextPtrFromDeserializeCorrectly)
{
  std::stringstream str;

  EXPECT_NO_THROW(context.writeTo(str));

  helib::Context* deserialized_contextp = helib::Context::readPtrFrom(str);

  EXPECT_EQ(context, *deserialized_contextp);
}

TEST(TestBinIO_BGV, readContextFromDeserializeCorrectlyBootstrappable)
{
  // clang-format off
  helib::Context context = helib::ContextBuilder<helib::BGV>()
      .m(1271)
      .p(2)
      .r(1)
      .gens({1026, 249})
      .ords({30, -2})
      .bits(30)
      .bootstrappable(true)
      .mvec(helib::convert<NTL::Vec<long>>(std::vector<long>({31, 41})))
      .build();
  // clang-format on

  std::stringstream str;

  EXPECT_NO_THROW(context.writeTo(str));

  helib::Context deserialized_context = helib::Context::readFrom(str);

  EXPECT_EQ(context, deserialized_context);
  EXPECT_TRUE(deserialized_context.isBootstrappable());
}

TEST_P(TestBinIO_BGV, canPerformOperationWithDeserializedContext)
{
  std::stringstream ss;

  context.writeTo(ss);

  helib::Context deserialized_context = helib::Context::readFrom(ss);

  EXPECT_NO_THROW(helib::PubKey{deserialized_context});
  EXPECT_NO_THROW(helib::SecKey{deserialized_context});

  helib::Ctxt ctxt{helib::PubKey(deserialized_context)};

  EXPECT_NO_THROW(ctxt.square());
  EXPECT_NO_THROW(ctxt += ctxt);
  EXPECT_NO_THROW(ctxt.reLinearize());
  EXPECT_NO_THROW(deserialized_context.getEA().rotate(ctxt, 1));
}

TEST_P(TestBinIO_BGV, singleFunctionSerializationOfKeys)
{
  std::stringstream str;

  EXPECT_NO_THROW(publicKey.writeTo(str));
  EXPECT_NO_THROW(secretKey.writeTo(str));
}

TEST_P(TestBinIO_BGV, singleFunctionDeserializationOfKeys)
{
  std::stringstream str;

  publicKey.writeTo(str);

  EXPECT_NO_THROW(publicKey.readFrom(str, context));

  str.str("");
  str.clear();

  secretKey.writeTo(str);

  EXPECT_NO_THROW(secretKey.readFrom(str, context));
}

TEST_P(TestBinIO_BGV, throwsWhenPrePublicKeyEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(publicKey.writeTo(ss));

  // Delete pre-publicKey eye catcher
  std::string s = ss.str();
  std::size_t pos = s.find(eyeCatcherToStr(helib::EyeCatcher::PK_BEGIN));
  s.erase(pos, pos + helib::EyeCatcher::SIZE);
  ss.str(s);

  EXPECT_THROW(publicKey.readFrom(ss, context), helib::IOError);
}

TEST_P(TestBinIO_BGV, throwsWhenPostPublicKeyEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(publicKey.writeTo(ss));

  // Delete post-publicKey eye catcher
  std::string s = ss.str();
  s.erase(s.find(eyeCatcherToStr(helib::EyeCatcher::PK_END)), s.size() - 1);
  ss.str(s);

  EXPECT_THROW(publicKey.readFrom(ss, context), helib::IOError);
}

TEST_P(TestBinIO_BGV, throwsWhenPreSecretKeyEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(secretKey.writeTo(ss));

  // Delete pre-secretKey eye catcher
  std::string s = ss.str();
  std::size_t pos = s.find(eyeCatcherToStr(helib::EyeCatcher::SK_BEGIN));
  s.erase(pos, pos + helib::EyeCatcher::SIZE);
  ss.str(s);

  EXPECT_THROW(secretKey.readFrom(ss, context), helib::IOError);
}

TEST_P(TestBinIO_BGV, throwsWhenPostSecretKeyEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(secretKey.writeTo(ss));

  // Delete post-secretKey eye catcher
  std::string s = ss.str();
  s.erase(s.find(eyeCatcherToStr(helib::EyeCatcher::SK_END)), s.size() - 1);
  ss.str(s);

  EXPECT_THROW(secretKey.readFrom(ss, context), helib::IOError);
}

TEST_P(TestBinIO_BGV, throwsWhenPreOnlySecretKeyEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(secretKey.writeTo(ss, /*sk_only=*/true));

  // Delete pre-secretKey eye catcher
  std::string s = ss.str();
  std::size_t pos = s.find(eyeCatcherToStr(helib::EyeCatcher::SK_BEGIN));
  s.erase(pos, pos + helib::EyeCatcher::SIZE);
  ss.str(s);

  EXPECT_THROW(secretKey.readFrom(ss, context, /*sk_only=*/true),
               helib::IOError);
}

TEST_P(TestBinIO_BGV, throwsWhenPostOnlySecretKeyEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(secretKey.writeTo(ss, /*sk_only=*/true));

  // Delete post-secretKey eye catcher
  std::string s = ss.str();
  s.erase(s.find(eyeCatcherToStr(helib::EyeCatcher::SK_END)), s.size() - 1);
  ss.str(s);

  EXPECT_THROW(secretKey.readFrom(ss, context, /*sk_only=*/true),
               helib::IOError);
}

TEST_P(TestBinIO_BGV, readOnlySecretKeyThrowsWhenMismatchContext)
{
  std::stringstream str;
  secretKey.writeTo(str, /*sk_only=*/true);
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
  EXPECT_THROW(deserialized_sk.readFrom(str, new_context, /*sk_only=*/true),
               helib::LogicError);
}

TEST_P(TestBinIO_BGV, readKeysFromDeserializeCorrectly)
{
  std::stringstream str;

  EXPECT_NO_THROW(publicKey.writeTo(str));

  helib::PubKey deserialized_pk = helib::PubKey::readFrom(str, context);

  EXPECT_EQ(publicKey, deserialized_pk);

  str.str("");
  str.clear();

  EXPECT_NO_THROW(secretKey.writeTo(str));

  helib::SecKey deserialized_sk = helib::SecKey::readFrom(str, context);

  EXPECT_EQ(secretKey, deserialized_sk);
}

TEST_P(TestBinIO_BGV, readKeyPtrsFromDeserializeCorrectly)
{
  std::stringstream str;

  EXPECT_NO_THROW(publicKey.writeTo(str));

  std::shared_ptr<helib::PubKey> deserialized_pkp =
      std::make_shared<helib::PubKey>(helib::PubKey::readFrom(str, context));

  EXPECT_EQ(publicKey, *deserialized_pkp);

  str.str("");
  str.clear();

  EXPECT_NO_THROW(secretKey.writeTo(str));

  std::shared_ptr<helib::SecKey> deserialized_skp =
      std::make_shared<helib::SecKey>(helib::SecKey::readFrom(str, context));

  EXPECT_EQ(secretKey, *deserialized_skp);
}

TEST_P(TestBinIO_BGV, canEncryptWithDeserializedPublicKey)
{
  std::stringstream ss;

  publicKey.writeTo(ss);

  helib::PubKey deserialized_pk = helib::PubKey::readFrom(ss, context);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(deserialized_pk);

  EXPECT_NO_THROW(ptxt.encrypt(ctxt));

  decrypted_result.decrypt(ctxt, secretKey);

  EXPECT_EQ(ptxt, decrypted_result);
}

TEST_P(TestBinIO_BGV, canEncryptWithDeserializedPublicEncryptionKeyAndContext)
{
  std::stringstream ss;

  helib::PubKey& encPublicKey = secretKey;

  context.writeTo(ss);
  encPublicKey.writeTo(ss);

  helib::Context deserialized_context = helib::Context::readFrom(ss);
  helib::PubKey deserialized_enc_pk =
      helib::PubKey::readFrom(ss, deserialized_context);

  helib::PtxtArray ptxt(deserialized_context),
      decrypted_result(deserialized_context);
  ptxt.random();
  helib::Ctxt ctxt(deserialized_enc_pk);

  EXPECT_NO_THROW(ptxt.encrypt(ctxt));

  secretKey.writeTo(ss, /*sk_only=*/true);
  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss, deserialized_context, /*sk_only=*/true);
  decrypted_result.decrypt(ctxt, deserialized_sk);

  EXPECT_EQ(ptxt, decrypted_result);
}

TEST_P(TestBinIO_BGV, canEncryptWithDeserializedSecretKey)
{
  std::stringstream ss;

  secretKey.writeTo(ss);

  helib::SecKey deserialized_sk = helib::SecKey::readFrom(ss, context);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(deserialized_sk);

  EXPECT_NO_THROW(ptxt.encrypt(ctxt));

  decrypted_result.decrypt(ctxt, secretKey);

  EXPECT_EQ(ptxt, decrypted_result);
}

TEST_P(TestBinIO_BGV, canEncryptWithDeserializedSecretKeyOnly)
{
  std::stringstream ss;
  secretKey.writeTo(ss, /*sk_only=*/true);

  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss, context, /*sk_only=*/true);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(deserialized_sk);
  EXPECT_NO_THROW(ptxt.encrypt(ctxt));
  decrypted_result.decrypt(ctxt, secretKey);
  EXPECT_EQ(ptxt, decrypted_result);
}

TEST_P(TestBinIO_BGV, canDecryptWithDeserializedSecretKey)
{
  std::stringstream ss;
  secretKey.writeTo(ss);
  helib::SecKey deserialized_sk = helib::SecKey::readFrom(ss, context);
  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(publicKey);

  ptxt.encrypt(ctxt);

  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt, deserialized_sk));
  EXPECT_EQ(ptxt, decrypted_result);
}

TEST_P(TestBinIO_BGV, canDecryptWithDeserializedSecretKeyOnly)
{
  std::stringstream ss;

  secretKey.writeTo(ss, /*sk_only=*/true);

  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss, context, /*sk_only=*/true);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(publicKey);

  ptxt.encrypt(ctxt);

  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt, deserialized_sk));
  EXPECT_EQ(ptxt, decrypted_result);
}

TEST_P(TestBinIO_BGV, canDecryptWithDeserializedSecretKeyOnlyAndContext)
{
  std::stringstream ss;
  context.writeTo(ss);
  secretKey.writeTo(ss, /*sk_only=*/true);
  helib::Context deserialized_context = helib::Context::readFrom(ss);
  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss, deserialized_context, /*sk_only=*/true);
  helib::PtxtArray ptxt(deserialized_context),
      decrypted_result(deserialized_context);
  ptxt.random();
  helib::Ctxt ctxt(deserialized_sk);

  ptxt.encrypt(ctxt);

  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt, deserialized_sk));
  EXPECT_EQ(ptxt, decrypted_result);
}

TEST_P(TestBinIO_BGV, singleFunctionSerializationOfCiphertext)
{
  std::stringstream str;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(str));
}

TEST_P(TestBinIO_BGV, singleFunctionDeserializationOfCiphertext)
{
  std::stringstream str;
  helib::Ctxt ctxt(publicKey);
  ctxt.writeTo(str);

  EXPECT_NO_THROW(helib::Ctxt::readFrom(str, publicKey));
}

TEST_P(TestBinIO_BGV, singleFunctionDeserializationOfCiphertextInPlace)
{
  std::stringstream str;
  helib::Ctxt ctxt(publicKey);
  ctxt.writeTo(str);

  EXPECT_NO_THROW(ctxt.read(str));
}

TEST_P(TestBinIO_BGV, throwsWhenPreCiphertextEyeCatcherNotFound)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(ss));

  // Delete pre-ciphertext eye catcher
  std::string s = ss.str();
  std::size_t pos = s.find(eyeCatcherToStr(helib::EyeCatcher::CTXT_BEGIN));
  s.erase(pos, pos + helib::EyeCatcher::SIZE);
  ss.str(s);

  EXPECT_THROW(ctxt.readFrom(ss, publicKey), helib::IOError);
}

TEST_P(TestBinIO_BGV, throwsWhenPostCiphertextEyeCatcherNotFound)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(ss));

  // Delete post-ciphertext eye catcher
  std::string s = ss.str();
  s.erase(s.find(eyeCatcherToStr(helib::EyeCatcher::CTXT_END)), s.size() - 1);
  ss.str(s);

  EXPECT_THROW(ctxt.readFrom(ss, publicKey), helib::IOError);
}

TEST_P(TestBinIO_BGV, throwsWhenPreCiphertextEyeCatcherNotFoundInPlace)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(ss));

  // Delete pre-ciphertext eye catcher
  std::string s = ss.str();
  std::size_t pos = s.find(eyeCatcherToStr(helib::EyeCatcher::CTXT_BEGIN));
  s.erase(pos, pos + helib::EyeCatcher::SIZE);
  ss.str(s);

  EXPECT_THROW(ctxt.read(ss), helib::IOError);
}

TEST_P(TestBinIO_BGV, throwsWhenPostCiphertextEyeCatcherNotFoundInPlace)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(ss));

  // Delete post-ciphertext eye catcher
  std::string s = ss.str();
  s.erase(s.find(eyeCatcherToStr(helib::EyeCatcher::CTXT_END)), s.size() - 1);
  ss.str(s);

  EXPECT_THROW(ctxt.read(ss), helib::IOError);
}

TEST_P(TestBinIO_BGV, readCiphertextFromDeserializeCorrectly)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  ss << ctxt;

  helib::Ctxt deserialized_ctxt(publicKey);
  ss >> deserialized_ctxt;

  EXPECT_EQ(ctxt, deserialized_ctxt);
}

TEST_P(TestBinIO_BGV, readCiphertextInPlaceFromDeserializeCorrectly)
{
  std::stringstream str;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(str));

  helib::Ctxt deserialized_ctxt(publicKey);
  deserialized_ctxt.read(str);

  EXPECT_EQ(ctxt, deserialized_ctxt);
}

TEST_P(TestBinIO_BGV, readCiphertextAndReadCiphertextInPlaceAreEquivalent)
{
  std::stringstream str1;
  std::stringstream str2;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(str1));
  EXPECT_NO_THROW(ctxt.writeTo(str2));

  helib::Ctxt deserialized_ctxt = helib::Ctxt::readFrom(str1, publicKey);
  helib::Ctxt inplace_ctxt(publicKey);
  inplace_ctxt.read(str2);

  EXPECT_EQ(inplace_ctxt, deserialized_ctxt);
}

TEST_P(TestBinIO_BGV, canPerformOperationsOnDeserializedCiphertext)
{
  std::stringstream ss1, ss2;
  helib::Ctxt ctxt(publicKey);

  ctxt.writeTo(ss1);
  ctxt.writeTo(ss2);

  helib::Ctxt deserialized_ctxt = helib::Ctxt::readFrom(ss1, publicKey);
  helib::PtxtArray ptxt1(ea), ptxt2(ea);

  EXPECT_NO_THROW(deserialized_ctxt *= ctxt);
  EXPECT_NO_THROW(deserialized_ctxt += ctxt);
  EXPECT_NO_THROW(deserialized_ctxt.reLinearize());
  EXPECT_NO_THROW(ea.rotate(deserialized_ctxt, 1));
  EXPECT_NO_THROW(ptxt1.decrypt(deserialized_ctxt, secretKey));

  helib::Ctxt inplace_ctxt(publicKey);
  inplace_ctxt.read(ss2);

  EXPECT_NO_THROW(inplace_ctxt *= ctxt);
  EXPECT_NO_THROW(inplace_ctxt += ctxt);
  EXPECT_NO_THROW(inplace_ctxt.reLinearize());
  EXPECT_NO_THROW(ea.rotate(inplace_ctxt, 1));
  EXPECT_NO_THROW(ptxt2.decrypt(inplace_ctxt, secretKey));

  EXPECT_EQ(ptxt1, ptxt2);
}

TEST_P(
    TestBinIO_BGV,
    canPerformOperationsOnDeserializedCiphertextWithDeserializedContextAndEvalKey)
{
  std::stringstream ss1, ss2;
  helib::PtxtArray ptxt(context);
  ptxt.random();
  helib::Ctxt ctxt(publicKey);
  ptxt.encrypt(ctxt);

  ctxt.writeTo(ss1);
  context.writeTo(ss2);
  publicKey.writeTo(ss2);

  helib::Context deserialized_context = helib::Context::readFrom(ss2);
  helib::PubKey deserialized_eval_pk =
      helib::PubKey::readFrom(ss2, deserialized_context);

  helib::Ctxt deserialized_ctxt =
      helib::Ctxt::readFrom(ss1, deserialized_eval_pk);
  helib::Ctxt deserialized_ctxt_copy(deserialized_ctxt);
  EXPECT_NO_THROW(deserialized_ctxt *= deserialized_ctxt_copy);
  EXPECT_NO_THROW(deserialized_ctxt += deserialized_ctxt_copy);
  EXPECT_NO_THROW(deserialized_context.getEA().rotate(deserialized_ctxt, 1));

  helib::PtxtArray ptxt_copy(ptxt);
  ptxt *= ptxt_copy;
  ptxt += ptxt_copy;
  rotate(ptxt, 1);
  helib::PtxtArray decrypted_result(deserialized_context);
  std::stringstream ss3;
  secretKey.writeTo(ss3, /*sk_only=*/true);
  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss3, deserialized_context, /*sk_only=*/true);
  // to have consistent contexts, read ptxt to stream and then read with this
  // context
  ptxt.writeToJSON(ss3);
  helib::PtxtArray ptxt_result =
      helib::PtxtArray::readFromJSON(ss3, deserialized_context);
  decrypted_result.decrypt(deserialized_ctxt, deserialized_sk);
  EXPECT_EQ(decrypted_result, ptxt_result);
}

TEST_P(TestBinIO_BGV, decryptWithDeserializedSecretKeyOnlyAfterComputation)
{
  std::stringstream ss;

  secretKey.writeTo(ss, /*sk_only=*/true);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(publicKey);
  ptxt.encrypt(ctxt);

  EXPECT_NO_THROW(ptxt *= ptxt);
  EXPECT_NO_THROW(ptxt += ptxt);
  EXPECT_NO_THROW(rotate(ptxt, 1));

  EXPECT_NO_THROW(ctxt *= ctxt);
  EXPECT_NO_THROW(ctxt += ctxt);
  EXPECT_NO_THROW(ctxt.reLinearize());
  EXPECT_NO_THROW(ea.rotate(ctxt, 1));

  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss, context, /*sk_only=*/true);
  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt, deserialized_sk));

  EXPECT_EQ(ptxt, decrypted_result);
}

TEST_P(TestBinIO_BGV, decryptWithDeserializedSecretKeyOnlyAfterMultLowLvl)
{
  std::stringstream ss;

  secretKey.writeTo(ss, /*sk_only=*/true);

  helib::PtxtArray ptxt1(ea), ptxt2(ea), decrypted_result(ea);
  ptxt1.random();
  ptxt2.random();
  helib::Ctxt ctxt1(publicKey);
  ptxt1.encrypt(ctxt1);
  helib::Ctxt ctxt2(publicKey);
  ptxt2.encrypt(ctxt2);

  EXPECT_NO_THROW(ptxt1 *= ptxt2);

  EXPECT_NO_THROW(ctxt1.multLowLvl(ctxt2));

  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss, context, /*sk_only=*/true);
  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt1, deserialized_sk));

  EXPECT_EQ(ptxt1, decrypted_result);
}

TEST_P(TestBinIO_BGV,
       decryptDeserializedCiphertextWithDeserializedSecretKeyOnly)
{
  std::stringstream ss1, ss2;
  helib::Ctxt ctxt(publicKey);
  helib::PtxtArray ptxt(context);
  ptxt.decrypt(ctxt, secretKey);

  ctxt.writeTo(ss1);
  context.writeTo(ss2);
  secretKey.writeTo(ss2, /*sk_only=*/true);

  helib::Context deserialized_context = helib::Context::readFrom(ss2);
  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss2, deserialized_context, /*sk_only=*/true);
  helib::Ctxt deserialized_ctxt = helib::Ctxt::readFrom(ss1, deserialized_sk);
  helib::PtxtArray decrypted_result(deserialized_context);
  decrypted_result.decrypt(deserialized_ctxt, deserialized_sk);

  // for consistent contexts, write ptxt to stream and read with context
  std::stringstream ss3;
  ptxt.writeToJSON(ss3);
  helib::PtxtArray ptxt_result =
      helib::PtxtArray::readFromJSON(ss3, deserialized_context);

  EXPECT_EQ(ptxt_result, decrypted_result);
}

TEST_P(TestBinIO_CKKS, singleFunctionSerialization)
{
  std::stringstream str;

  EXPECT_NO_THROW(context.writeTo(str));
}

TEST_P(TestBinIO_CKKS, singleFunctionDeserialization)
{
  std::stringstream str;

  context.writeTo(str);

  EXPECT_NO_THROW(context.readFrom(str));
}

TEST_P(TestBinIO_CKKS, throwsWhenPreContextEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(context.writeTo(ss));

  // Delete pre-context eye catcher
  std::string s = ss.str();
  std::size_t pos = s.find(eyeCatcherToStr(helib::EyeCatcher::CONTEXT_BEGIN));
  s.erase(pos, pos + helib::EyeCatcher::SIZE);
  ss.str(s);

  EXPECT_THROW(context.readFrom(ss), helib::IOError);
}

TEST_P(TestBinIO_CKKS, throwsWhenPostContextEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(context.writeTo(ss));

  // Delete post-context eye catcher
  std::string s = ss.str();
  s.erase(s.find(eyeCatcherToStr(helib::EyeCatcher::CONTEXT_END)),
          s.size() - 1);
  ss.str(s);

  EXPECT_THROW(context.readFrom(ss), helib::IOError);
}

TEST_P(TestBinIO_CKKS, readContextFromDeserializeCorrectly)
{
  std::stringstream str;

  EXPECT_NO_THROW(context.writeTo(str));

  helib::Context deserialized_context = helib::Context::readFrom(str);

  EXPECT_EQ(context, deserialized_context);
}

TEST_P(TestBinIO_CKKS, readContextPtrFromDeserializeCorrectly)
{
  std::stringstream str;

  EXPECT_NO_THROW(context.writeTo(str));

  helib::Context* deserialized_contextp = helib::Context::readPtrFrom(str);

  EXPECT_EQ(context, *deserialized_contextp);
}

TEST_P(TestBinIO_CKKS, canPerformOperationWithDeserializedContext)
{
  std::stringstream ss;

  context.writeTo(ss);

  helib::Context deserialized_context = helib::Context::readFrom(ss);

  EXPECT_NO_THROW(helib::PubKey{deserialized_context});
  EXPECT_NO_THROW(helib::SecKey{deserialized_context});

  helib::Ctxt ctxt{helib::PubKey(deserialized_context)};

  EXPECT_NO_THROW(ctxt.square());
  EXPECT_NO_THROW(ctxt += ctxt);
  EXPECT_NO_THROW(ctxt.reLinearize());
  EXPECT_NO_THROW(deserialized_context.getEA().rotate(ctxt, 1));
}

TEST_P(TestBinIO_CKKS, singleFunctionSerializationOfKeys)
{
  std::stringstream str;

  EXPECT_NO_THROW(publicKey.writeTo(str));
  EXPECT_NO_THROW(secretKey.writeTo(str));
}

TEST_P(TestBinIO_CKKS, singleFunctionDeserializationOfKeys)
{
  std::stringstream str;

  publicKey.writeTo(str);

  EXPECT_NO_THROW(publicKey.readFrom(str, context));

  str.str("");
  str.clear();

  secretKey.writeTo(str);

  EXPECT_NO_THROW(secretKey.readFrom(str, context));
}

TEST_P(TestBinIO_CKKS, throwsWhenPrePublicKeyEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(publicKey.writeTo(ss));

  // Delete pre-publicKey eye catcher
  std::string s = ss.str();
  std::size_t pos = s.find(eyeCatcherToStr(helib::EyeCatcher::PK_BEGIN));
  s.erase(pos, pos + helib::EyeCatcher::SIZE);
  ss.str(s);

  EXPECT_THROW(publicKey.readFrom(ss, context), helib::IOError);
}

TEST_P(TestBinIO_CKKS, throwsWhenPostPublicKeyEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(publicKey.writeTo(ss));

  // Delete post-publicKey eye catcher
  std::string s = ss.str();
  s.erase(s.find(eyeCatcherToStr(helib::EyeCatcher::PK_END)), s.size() - 1);
  ss.str(s);

  EXPECT_THROW(publicKey.readFrom(ss, context), helib::IOError);
}

TEST_P(TestBinIO_CKKS, throwsWhenPreSecretKeyEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(secretKey.writeTo(ss));

  // Delete pre-secretKey eye catcher
  std::string s = ss.str();
  std::size_t pos = s.find(eyeCatcherToStr(helib::EyeCatcher::SK_BEGIN));
  s.erase(pos, pos + helib::EyeCatcher::SIZE);
  ss.str(s);

  EXPECT_THROW(secretKey.readFrom(ss, context), helib::IOError);
}

TEST_P(TestBinIO_CKKS, throwsWhenPostSecretKeyEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(secretKey.writeTo(ss));

  // Delete post-secretKey eye catcher
  std::string s = ss.str();
  s.erase(s.find(eyeCatcherToStr(helib::EyeCatcher::SK_END)), s.size() - 1);
  ss.str(s);

  EXPECT_THROW(secretKey.readFrom(ss, context), helib::IOError);
}

TEST_P(TestBinIO_CKKS, throwsWhenPreOnlySecretKeyEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(secretKey.writeTo(ss, /*sk_only=*/true));

  // Delete pre-secretKey eye catcher
  std::string s = ss.str();
  std::size_t pos = s.find(eyeCatcherToStr(helib::EyeCatcher::SK_BEGIN));
  s.erase(pos, pos + helib::EyeCatcher::SIZE);
  ss.str(s);

  EXPECT_THROW(secretKey.readFrom(ss, context, /*sk_only=*/true),
               helib::IOError);
}

TEST_P(TestBinIO_CKKS, throwsWhenPostOnlySecretKeyEyeCatcherNotFound)
{
  std::stringstream ss;

  EXPECT_NO_THROW(secretKey.writeTo(ss, /*sk_only=*/true));

  // Delete post-secretKey eye catcher
  std::string s = ss.str();
  s.erase(s.find(eyeCatcherToStr(helib::EyeCatcher::SK_END)), s.size() - 1);
  ss.str(s);

  EXPECT_THROW(secretKey.readFrom(ss, context, /*sk_only=*/true),
               helib::IOError);
}

TEST_P(TestBinIO_CKKS, readOnlySecretKeyThrowsWhenMismatchContext)
{
  std::stringstream str;
  secretKey.writeTo(str, /*sk_only=*/true);
  helib::Context new_context(helib::ContextBuilder<helib::CKKS>()
                                 .m(32)
                                 .precision(precision)
                                 .bits(bits)
                                 .build());
  helib::SecKey deserialized_sk(new_context);
  EXPECT_THROW(deserialized_sk.readFrom(str, new_context, /*sk_only=*/true),
               helib::LogicError);
}

TEST_P(TestBinIO_CKKS, readKeysFromDeserializeCorrectly)
{
  std::stringstream str;

  EXPECT_NO_THROW(publicKey.writeTo(str));

  helib::PubKey deserialized_pk = helib::PubKey::readFrom(str, context);

  EXPECT_EQ(publicKey, deserialized_pk);

  str.str("");
  str.clear();

  EXPECT_NO_THROW(secretKey.writeTo(str));

  helib::SecKey deserialized_sk = helib::SecKey::readFrom(str, context);

  EXPECT_EQ(secretKey, deserialized_sk);
}

TEST_P(TestBinIO_CKKS, readKeyPtrsFromDeserializeCorrectly)
{
  std::stringstream str;

  EXPECT_NO_THROW(publicKey.writeTo(str));

  std::shared_ptr<helib::PubKey> deserialized_pkp =
      std::make_shared<helib::PubKey>(helib::PubKey::readFrom(str, context));

  EXPECT_EQ(publicKey, *deserialized_pkp);

  str.str("");
  str.clear();

  EXPECT_NO_THROW(secretKey.writeTo(str));

  std::shared_ptr<helib::SecKey> deserialized_skp =
      std::make_shared<helib::SecKey>(helib::SecKey::readFrom(str, context));

  EXPECT_EQ(secretKey, *deserialized_skp);
}

TEST_P(TestBinIO_CKKS, canEncryptWithDeserializedPublicKey)
{
  std::stringstream ss;
  publicKey.writeTo(ss);
  helib::PubKey deserialized_pk = helib::PubKey::readFrom(ss, context);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(deserialized_pk);

  EXPECT_NO_THROW(ptxt.encrypt(ctxt));
  decrypted_result.decrypt(ctxt, secretKey);
  EXPECT_EQ(ptxt, helib::Approx(decrypted_result));
}

TEST_P(TestBinIO_CKKS, canEncryptWithDeserializedPublicEncryptionKeyAndContext)
{
  std::stringstream ss;

  helib::PubKey& encPublicKey = secretKey;
  context.writeTo(ss);
  encPublicKey.writeTo(ss);
  helib::Context deserialized_context = helib::Context::readFrom(ss);
  helib::PubKey deserialized_enc_pk =
      helib::PubKey::readFrom(ss, deserialized_context);
  helib::PtxtArray ptxt(deserialized_context),
      decrypted_result(deserialized_context);
  ptxt.random();
  helib::Ctxt ctxt(deserialized_enc_pk);
  EXPECT_NO_THROW(ptxt.encrypt(ctxt));

  // now make a secret key with the exact same context
  secretKey.writeTo(ss, /*sk_only=*/true);
  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss, deserialized_context, /*sk_only=*/true);
  decrypted_result.decrypt(ctxt, deserialized_sk);
  EXPECT_EQ(ptxt, helib::Approx(decrypted_result));
}

TEST_P(TestBinIO_CKKS, canEncryptWithDeserializedSecretKey)
{
  std::stringstream ss;
  secretKey.writeTo(ss);
  helib::SecKey deserialized_sk = helib::SecKey::readFrom(ss, context);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(deserialized_sk);

  EXPECT_NO_THROW(ptxt.encrypt(ctxt));
  decrypted_result.decrypt(ctxt, secretKey);
  EXPECT_EQ(ptxt, helib::Approx(decrypted_result));
}

TEST_P(TestBinIO_CKKS, canEncryptWithDeserializedSecretKeyOnly)
{
  std::stringstream ss;
  secretKey.writeTo(ss, /*sk_only=*/true);
  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss, context, /*sk_only=*/true);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(deserialized_sk);

  EXPECT_NO_THROW(ptxt.encrypt(ctxt));
  decrypted_result.decrypt(ctxt, secretKey);
  EXPECT_EQ(ptxt, helib::Approx(decrypted_result));
}

TEST_P(TestBinIO_CKKS, canDecryptWithDeserializedSecretKey)
{
  std::stringstream ss;
  secretKey.writeTo(ss);
  helib::SecKey deserialized_sk = helib::SecKey::readFrom(ss, context);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(publicKey);

  ptxt.encrypt(ctxt);

  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt, deserialized_sk));
  EXPECT_EQ(ptxt, helib::Approx(decrypted_result));
}

TEST_P(TestBinIO_CKKS, canDecryptWithDeserializedSecretKeyOnly)
{
  std::stringstream ss;
  secretKey.writeTo(ss, /*sk_only=*/true);

  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss, context, /*sk_only=*/true);
  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(publicKey);

  ptxt.encrypt(ctxt);

  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt, deserialized_sk));
  EXPECT_EQ(ptxt, helib::Approx(decrypted_result));
}

TEST_P(TestBinIO_CKKS, canDecryptWithDeserializedSecretKeyOnlyAndContext)
{
  std::stringstream ss;
  context.writeTo(ss);
  secretKey.writeTo(ss, /*sk_only=*/true);
  helib::Context deserialized_context = helib::Context::readFrom(ss);
  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss, deserialized_context, /*sk_only=*/true);
  helib::PtxtArray ptxt(deserialized_context),
      decrypted_result(deserialized_context);
  ptxt.random();
  helib::Ctxt ctxt(deserialized_sk);

  ptxt.encrypt(ctxt);

  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt, deserialized_sk));
  EXPECT_EQ(ptxt, helib::Approx(decrypted_result));
}

TEST_P(TestBinIO_CKKS, singleFunctionSerializationOfCiphertext)
{
  std::stringstream str;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(str));
}

TEST_P(TestBinIO_CKKS, singleFunctionDeserializationOfCiphertext)
{
  std::stringstream str;
  helib::Ctxt ctxt(publicKey);
  ctxt.writeTo(str);

  EXPECT_NO_THROW(helib::Ctxt::readFrom(str, publicKey));
}

TEST_P(TestBinIO_CKKS, singleFunctionDeserializationOfCiphertextInPlace)
{
  std::stringstream str;
  helib::Ctxt ctxt(publicKey);
  ctxt.writeTo(str);

  EXPECT_NO_THROW(ctxt.read(str));
}

TEST_P(TestBinIO_CKKS, throwsWhenPreCiphertextEyeCatcherNotFound)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(ss));

  // Delete pre-ciphertext eye catcher
  std::string s = ss.str();
  std::size_t pos = s.find(eyeCatcherToStr(helib::EyeCatcher::CTXT_BEGIN));
  s.erase(pos, pos + helib::EyeCatcher::SIZE);
  ss.str(s);

  EXPECT_THROW(ctxt.readFrom(ss, publicKey), helib::IOError);
}

TEST_P(TestBinIO_CKKS, throwsWhenPostCiphertextEyeCatcherNotFound)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(ss));

  // Delete post-ciphertext eye catcher
  std::string s = ss.str();
  s.erase(s.find(eyeCatcherToStr(helib::EyeCatcher::CTXT_END)), s.size() - 1);
  ss.str(s);

  EXPECT_THROW(ctxt.readFrom(ss, publicKey), helib::IOError);
}

TEST_P(TestBinIO_CKKS, throwsWhenPreCiphertextEyeCatcherNotFoundInPlace)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(ss));

  // Delete pre-ciphertext eye catcher
  std::string s = ss.str();
  std::size_t pos = s.find(eyeCatcherToStr(helib::EyeCatcher::CTXT_BEGIN));
  s.erase(pos, pos + helib::EyeCatcher::SIZE);
  ss.str(s);

  EXPECT_THROW(ctxt.read(ss), helib::IOError);
}

TEST_P(TestBinIO_CKKS, throwsWhenPostCiphertextEyeCatcherNotFoundInPlace)
{
  std::stringstream ss;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(ss));

  // Delete post-ciphertext eye catcher
  std::string s = ss.str();
  s.erase(s.find(eyeCatcherToStr(helib::EyeCatcher::CTXT_END)), s.size() - 1);
  ss.str(s);

  EXPECT_THROW(ctxt.read(ss), helib::IOError);
}

TEST_P(TestBinIO_CKKS, readCiphertextFromDeserializeCorrectly)
{
  std::stringstream str;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(str));

  helib::Ctxt deserialized_ctxt = helib::Ctxt::readFrom(str, publicKey);

  EXPECT_EQ(ctxt, deserialized_ctxt);
}

TEST_P(TestBinIO_CKKS, readCiphertextInPlaceFromDeserializeCorrectly)
{
  std::stringstream str;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(str));

  helib::Ctxt deserialized_ctxt(publicKey);
  deserialized_ctxt.read(str);

  EXPECT_EQ(ctxt, deserialized_ctxt);
}

TEST_P(TestBinIO_CKKS, readCiphertextAndReadCiphertextInPlaceAreEquivalent)
{
  std::stringstream str1;
  std::stringstream str2;
  helib::Ctxt ctxt(publicKey);

  EXPECT_NO_THROW(ctxt.writeTo(str1));
  EXPECT_NO_THROW(ctxt.writeTo(str2));

  helib::Ctxt deserialized_ctxt = helib::Ctxt::readFrom(str1, publicKey);
  helib::Ctxt inplace_ctxt(publicKey);
  inplace_ctxt.read(str2);

  EXPECT_EQ(inplace_ctxt, deserialized_ctxt);
}

TEST_P(TestBinIO_CKKS, canPerformOperationsOnDeserializedCiphertext)
{
  std::stringstream ss1, ss2;
  helib::Ctxt ctxt(publicKey);
  ctxt.writeTo(ss1);
  ctxt.writeTo(ss2);

  helib::Ctxt deserialized_ctxt = helib::Ctxt::readFrom(ss1, publicKey);
  helib::PtxtArray ptxt1(ea), ptxt2(ea);

  EXPECT_NO_THROW(deserialized_ctxt *= ctxt);
  EXPECT_NO_THROW(deserialized_ctxt += ctxt);
  EXPECT_NO_THROW(deserialized_ctxt.reLinearize());
  EXPECT_NO_THROW(ea.rotate(deserialized_ctxt, 1));
  EXPECT_NO_THROW(ptxt1.decrypt(deserialized_ctxt, secretKey));

  helib::Ctxt inplace_ctxt(publicKey);
  inplace_ctxt.read(ss2);

  EXPECT_NO_THROW(inplace_ctxt *= ctxt);
  EXPECT_NO_THROW(inplace_ctxt += ctxt);
  EXPECT_NO_THROW(inplace_ctxt.reLinearize());
  EXPECT_NO_THROW(ea.rotate(inplace_ctxt, 1));
  EXPECT_NO_THROW(ptxt2.decrypt(inplace_ctxt, secretKey));

  EXPECT_EQ(ptxt1, ptxt2);
}

TEST_P(
    TestBinIO_CKKS,
    canPerformOperationsOnDeserializedCiphertextWithDeserializedContextAndEvalKey)
{
  std::stringstream ss1, ss2;
  helib::PtxtArray ptxt(context);
  ptxt.random();
  helib::Ctxt ctxt(publicKey);
  ptxt.encrypt(ctxt);

  ctxt.writeTo(ss1);
  context.writeTo(ss2);
  publicKey.writeTo(ss2);

  helib::Context deserialized_context = helib::Context::readFrom(ss2);
  helib::PubKey deserialized_eval_pk =
      helib::PubKey::readFrom(ss2, deserialized_context);

  helib::Ctxt deserialized_ctxt =
      helib::Ctxt::readFrom(ss1, deserialized_eval_pk);
  helib::Ctxt deserialized_ctxt_copy(deserialized_ctxt);
  EXPECT_NO_THROW(deserialized_ctxt *= deserialized_ctxt_copy);
  EXPECT_NO_THROW(deserialized_ctxt += deserialized_ctxt_copy);
  EXPECT_NO_THROW(deserialized_context.getEA().rotate(deserialized_ctxt, 1));

  helib::PtxtArray ptxt_copy(ptxt);
  ptxt *= ptxt_copy;
  ptxt += ptxt_copy;
  rotate(ptxt, 1);
  helib::PtxtArray decrypted_result(deserialized_context);
  std::stringstream ss3;
  secretKey.writeTo(ss3, /*sk_only=*/true);
  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss3, deserialized_context, /*sk_only=*/true);
  // to have consistent contexts, write ptxt to stream and then read with this
  // context
  ptxt.writeToJSON(ss3);
  helib::PtxtArray ptxt_result =
      helib::PtxtArray::readFromJSON(ss3, deserialized_context);
  decrypted_result.decrypt(deserialized_ctxt, deserialized_sk);
  EXPECT_EQ(decrypted_result, helib::Approx(ptxt_result));
}

TEST_P(TestBinIO_CKKS, decryptWithDeserializedSecretKeyOnlyAfterComputation)
{
  std::stringstream ss;

  secretKey.writeTo(ss, /*sk_only=*/true);

  helib::PtxtArray ptxt(ea), decrypted_result(ea);
  ptxt.random();
  helib::Ctxt ctxt(publicKey);
  ptxt.encrypt(ctxt);

  EXPECT_NO_THROW(ptxt *= ptxt);
  EXPECT_NO_THROW(ptxt += ptxt);
  EXPECT_NO_THROW(rotate(ptxt, 1));

  EXPECT_NO_THROW(ctxt *= ctxt);
  EXPECT_NO_THROW(ctxt += ctxt);
  EXPECT_NO_THROW(ctxt.reLinearize());
  EXPECT_NO_THROW(ea.rotate(ctxt, 1));

  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss, context, /*sk_only=*/true);
  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt, deserialized_sk));

  EXPECT_EQ(ptxt, helib::Approx(decrypted_result));
}

TEST_P(TestBinIO_CKKS, decryptWithDeserializedSecretKeyOnlyAfterMultLowLvl)
{
  std::stringstream ss;

  secretKey.writeTo(ss, /*sk_only=*/true);

  helib::PtxtArray ptxt1(ea), ptxt2(ea), decrypted_result(ea);
  ptxt1.random();
  ptxt2.random();
  helib::Ctxt ctxt1(publicKey);
  ptxt1.encrypt(ctxt1);
  helib::Ctxt ctxt2(publicKey);
  ptxt2.encrypt(ctxt2);

  EXPECT_NO_THROW(ptxt1 *= ptxt2);

  EXPECT_NO_THROW(ctxt1.multLowLvl(ctxt2));

  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss, context, /*sk_only=*/true);
  EXPECT_NO_THROW(decrypted_result.decrypt(ctxt1, deserialized_sk));

  EXPECT_EQ(ptxt1, helib::Approx(decrypted_result));
}

TEST_P(TestBinIO_CKKS,
       decryptDeserializedCiphertextWithDeserializedSecretKeyOnly)
{
  std::stringstream ss1, ss2;
  helib::Ctxt ctxt(publicKey);
  helib::PtxtArray ptxt(context);
  ptxt.decrypt(ctxt, secretKey);

  ctxt.writeTo(ss1);
  context.writeTo(ss2);
  secretKey.writeTo(ss2, /*sk_only=*/true);

  helib::Context deserialized_context = helib::Context::readFrom(ss2);
  helib::SecKey deserialized_sk =
      helib::SecKey::readFrom(ss2, deserialized_context, /*sk_only=*/true);
  helib::Ctxt deserialized_ctxt = helib::Ctxt::readFrom(ss1, deserialized_sk);
  helib::PtxtArray decrypted_result(deserialized_context);
  decrypted_result.decrypt(deserialized_ctxt, deserialized_sk);

  // for consistent contexts, write ptxt to stream and read with context
  std::stringstream ss3;
  ptxt.writeToJSON(ss3);
  helib::PtxtArray ptxt_result =
      helib::PtxtArray::readFromJSON(ss3, deserialized_context);

  EXPECT_EQ(ptxt_result, decrypted_result);
}

INSTANTIATE_TEST_SUITE_P(Parameters,
                         TestBinIO_BGV,
                         ::testing::Values(BGVParameters(/*m=*/45,
                                                         /*p=*/2,
                                                         /*r=*/1,
                                                         /*bits=*/30,
                                                         /*gens=*/{},
                                                         /*ords=*/{},
                                                         /*mvec=*/{})));

INSTANTIATE_TEST_SUITE_P(Parameters,
                         TestBinIO_CKKS,
                         ::testing::Values(CKKSParameters(/*m=*/64,
                                                          /*precision=*/30,
                                                          /*bits=*/60)));

} // namespace
