/* Copyright (C) 2012-2017 IBM Corp.
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
#include <cassert>
#include <string>
#include <sstream>
#include <NTL/ZZ.h>
#include "NumbTh.h"
#include "FHEContext.h"
#include "debugging.h"

#include "gtest/gtest.h"
#include "test_common.h"

namespace {
struct Parameters {
    Parameters(long m, long p, long r, std::vector<long> gens, std::vector<long> ords) :
        m(m), // Cyclotomic index
              // e.g., m=1024, m=2047
        p(p), // plaintext base
              // use p=-1 for the complex field (CKKS)
        r(r), // lifting
        gens(gens), // use specified vector of generators
        ords(ords) // use specified vector of orders
    {};

    const long m;
    const long p;
    const long r;
    const std::vector<long> gens;
    const std::vector<long> ords;
    
    // Let googletest know how to print the Parameters
    friend std::ostream& operator<<(std::ostream& os, const Parameters& params) {
        return os << "{"
            << "m=" << params.m << ","
            << "p=" << params.p << ","
            << "r=" << params.r << ","
            << "gens=" << params.gens << ","
            << "ords=" << params.ords
            << "}";
    };
};

class GTest_PAlgebra : public ::testing::TestWithParam<Parameters> {
    protected:
        GTest_PAlgebra() :
            m(GetParam().m),
            p(GetParam().p),
            r(GetParam().r),
            gens(GetParam().gens),
            ords(GetParam().ords),
            context(m, p, r, gens, ords)
    {};

        const long m;
        const long p;
        const long r;
        const std::vector<long> gens;
        const std::vector<long> ords;
        FHEcontext context;

        static void printPrimeFactors(long m, const FHEcontext& context)
        {
            std::vector<long> f;
            factorize(f,m);
            std::cout << "factoring "<<m<<" gives [";
            for (const auto& factor : f)
                std::cout << factor << " ";
            std::cout << "]" << std::endl;
            context.zMStar.printout();
            std::cout << std::endl;
        }

        virtual void SetUp() override {
            buildModChain(context, 5, 2);
            if(!helib_test::noPrint) {
                printPrimeFactors(m, context);
            }
        };

        virtual void TearDown() override
        {
            cleanupGlobals();
        }
};

TEST_P(GTest_PAlgebra, reads_and_writes_contexts_as_strings)
{
  std::stringstream s1;
  writeContextBase(s1, context);
  s1 << context;

  std::string s2 = s1.str();

  if(!helib_test::noPrint) {
      std::cout << s2;
  }

  std::stringstream s3(s2);

  unsigned long m1, p1, r1;
  std::vector<long> gens, ords;
  readContextBase(s3, m1, p1, r1, gens, ords);

  FHEcontext c1(m1, p1, r1, gens, ords);
  s3 >> c1;

  EXPECT_EQ(context, c1);
};

INSTANTIATE_TEST_SUITE_P(small_parameters, GTest_PAlgebra, ::testing::Values(
            //FAST
            Parameters(91, 2, 1, std::vector<long>{}, std::vector<long>{})
            ));

} // namespace
