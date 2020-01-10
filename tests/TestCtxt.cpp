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

// This test file does not fully cover Ctxt, and is intended as a starting
// point for testing new functionality of Ctxt going forward.
// The older tests with more extensive coverage can be found in the files
// with names matching "GTest*".

#include <helib/helib.h>

#include "test_common.h"
#include "gtest/gtest.h"

namespace {

struct BGVParameters
{
  BGVParameters(unsigned m, unsigned p, unsigned r, unsigned bits) :
      m(m),
      p(p),
      r(r),
      bits(bits){};

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

class TestCtxt : public ::testing::TestWithParam<BGVParameters>
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

  TestCtxt() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      bits(GetParam().bits),
      context(m, p, r),
      secretKey((buildModChain(context, bits), context)),
      publicKey((secretKey.GenSecKey(), secretKey)),
      ea(*(context.ea))
  {}
};

TEST_P(TestCtxt, timesEqualsWithLongWorks)
{
  helib::Ptxt<helib::BGV> ptxt(context, std::vector<long>(ea.size(), 5));
  helib::Ctxt ctxt(publicKey);
  publicKey.Encrypt(ctxt, ptxt);
  ctxt *= 2l;

  helib::Ptxt<helib::BGV> expected_result(context,
                                          std::vector<long>(ea.size(), 10));
  helib::Ptxt<helib::BGV> decrypted_result(context);
  secretKey.Decrypt(decrypted_result, ctxt);

  EXPECT_EQ(decrypted_result, expected_result);
}

TEST_P(TestCtxt, timesEqualsWithZZXWorks)
{
  helib::Ptxt<helib::BGV> ptxt(context, std::vector<long>(ea.size(), 5));
  helib::Ctxt ctxt(publicKey);
  publicKey.Encrypt(ctxt, ptxt);
  ctxt *= NTL::ZZX(2l);

  helib::Ptxt<helib::BGV> expected_result(context,
                                          std::vector<long>(ea.size(), 10));
  helib::Ptxt<helib::BGV> decrypted_result(context);
  secretKey.Decrypt(decrypted_result, ctxt);

  EXPECT_EQ(decrypted_result, expected_result);
}

INSTANTIATE_TEST_SUITE_P(variousParameters,
                         TestCtxt,
                         ::testing::Values(BGVParameters(17, 2, 1, 100)));

} // namespace
