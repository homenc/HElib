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

#ifndef TEST_COMMON_H
#define TEST_COMMON_H
#include <helib/ArgMap.h>

namespace helib_test {

// Compare std::vector<std::complex<double>>
#define COMPARE_CXDOUBLE_VECS(vec1, vec2)                                      \
  EXPECT_EQ(vec1.size(), vec2.size());                                         \
  for (std::size_t i = 0; i < vec2.size(); ++i) {                              \
    EXPECT_DOUBLE_EQ(vec1[i].real(), vec2[i].real()) << "Slot " << i;          \
    EXPECT_DOUBLE_EQ(vec1[i].imag(), vec2[i].imag()) << "Slot " << i;          \
  }

// Compare std::vector<std::complex<float>>
#define COMPARE_CXFLOAT_VECS(vec1, vec2)                                       \
  EXPECT_EQ(vec1.size(), vec2.size());                                         \
  for (std::size_t i = 0; i < vec2.size(); ++i) {                              \
    EXPECT_FLOAT_EQ(vec1[i].real(), vec2[i].real()) << "Slot " << i;           \
    EXPECT_FLOAT_EQ(vec1[i].imag(), vec2[i].imag()) << "Slot " << i;           \
  }

extern char* path_of_executable;
extern bool noPrint;
extern bool verbose;
extern bool dry;
extern bool reset;
extern unsigned int random_seed;
extern long special_bits;

void parse_common_args(int argc, char* argv[]);

std::vector<std::pair<long, long>> getBadDimensionParams(long min_m,
                                                         long max_m,
                                                         long min_p,
                                                         long max_p,
                                                         long m_sparseness = 1,
                                                         long p_sparseness = 1);
std::vector<std::pair<long, long>> getGoodDimensionParams(
    long min_m,
    long max_m,
    long min_p,
    long max_p,
    long m_sparseness = 1,
    long p_sparseness = 1);
} // namespace helib_test

#endif /* ifndef TEST_COMMON_H */
