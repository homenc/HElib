/* Copyright (C) 2019 IBM Corp.
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
#include <cmath> // isinf
#include "helib.h"

#include "test_common.h"
#include "gtest/gtest.h"

// TODO - Currently does not cover well the Context object.

namespace {

struct BGVParameters
{
  BGVParameters(unsigned m, unsigned p, unsigned r) : m(m), p(p), r(r){};

  const unsigned m;
  const unsigned p;
  const unsigned r;

  friend std::ostream& operator<<(std::ostream& os, const BGVParameters& params)
  {
    return os << "{"
              << "m = " << params.m << ", "
              << "p = " << params.p << ", "
              << "r = " << params.r << "}";
  }
};

class Test_with_BGV : public ::testing::TestWithParam<BGVParameters>
{
protected:
  Test_with_BGV() :
      m(GetParam().m),
      p(GetParam().p),
      r(GetParam().r),
      context(std::make_shared<helib::FHEcontext>(m, p, r))
  {}

  const unsigned long m;
  const unsigned long p;
  const unsigned long r;

  const std::shared_ptr<helib::FHEcontext> context;
};

TEST_P(Test_with_BGV,
       Context_throw_exception_when_calculating_security_before_modchain_built)
{
  EXPECT_THROW(context->securityLevel(), helib::LogicError);
}

TEST_P(Test_with_BGV, Context_calculating_security_after_modchain_built)
{
  buildModChain(*context, /*bits=*/100, /*c=*/2);
  double result = context->securityLevel();
  EXPECT_FALSE(std::isinf(result));
}

INSTANTIATE_TEST_SUITE_P(various_parameters,
                         Test_with_BGV,
                         ::testing::Values(BGVParameters(17, 2, 1)));

} // namespace
