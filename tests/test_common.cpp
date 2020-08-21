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

#include "test_common.h"
#include <helib/Context.h>
#include <helib/EncryptedArray.h>
#include <iostream>
#include <random>

#include <NTL/BasicThreadPool.h>

namespace helib_test {

char* path_of_executable = nullptr;
bool noPrint = false;
bool verbose = false;
bool dry = false;
unsigned int random_seed = 0U;
long special_bits = 0L;

void parse_common_args(int argc, char* argv[])
{
  helib::ArgMap amap;
  path_of_executable = argv[0];
  amap.arg("dry", dry, "dry=1 for a dry-run");
  amap.arg("noPrint", noPrint, "suppress printouts");
  amap.arg("special_bits", special_bits, "# of bits in special primes");
  amap.arg("verbose", verbose, "print more information");
  amap.arg("seed", random_seed, "specify random seed for test data");
  amap.parse(argc, argv);
  if (random_seed == 0U) // Not specified: use random seed
    random_seed = std::random_device{}();
  // TODO: change this printout so that the random_seed is simply another
  // parameter of parameterised tests.  This will only be possible once
  // we parse the command-line args before gtest does.
  std::cout << "random seed: " << random_seed << std::endl;
  // TODO: Add in number of threads as an argument.
}

static inline long howManyThreads(long max = 0)
{
  if (max < 0)
    throw std::logic_error(std::string(__func__) + ": max number of threads (" +
                           std::to_string(max) +
                           ") cannot be less than zero."
                           " Zero means use the maximum available on system.");

  // Get system threads and clip to max
  long threads = std::thread::hardware_concurrency();
  if (max != 0)
    threads = (threads > max) ? max : threads;

  return threads ? threads : 1;
}

// Setting threads globally
long setThreads = []() -> long {
  NTL::SetNumThreads(howManyThreads(2));
  return NTL::AvailableThreads();
}();

// TODO: Should be a member of EncryptedArray?
bool hasBadDimension(const helib::Context& context)
{
  for (int i = 0; i < context.zMStar.numOfGens(); ++i)
    if (!context.ea->nativeDimension(i))
      return true;
  return false;
}

bool isPrime(const long num)
{
  for (long i = 2; i <= std::sqrt(num); ++i)
    if (num % i == 0)
      return false;
  return true;
}

std::vector<std::pair<long, long>> getParams(bool good,
                                             long min_m,
                                             long max_m,
                                             long min_p,
                                             long max_p,
                                             long m_sparseness,
                                             long p_sparseness)
{
  std::vector<std::pair<long, long>> params;
  std::vector<long> p_vals;
  std::vector<long> m_vals;

  for (long i = min_p; i < max_p; ++i)
    if (isPrime(i))
      p_vals.push_back(i);

  std::size_t write_head = 0;
  for (std::size_t j = 0; j < p_vals.size(); j += p_sparseness)
    p_vals[write_head++] = p_vals[j];
  p_vals.resize(write_head);

  for (long i = min_m; i < max_m; i += m_sparseness)
    m_vals.push_back(i);

  for (auto p : p_vals)
    for (auto m : m_vals) {
      if (m % p != 0) {
        helib::Context context(m, p, 1L);
        if (good ^ hasBadDimension(context)) {
          params.emplace_back(m, p);
        }
      }
    }
  return params;
}

std::vector<std::pair<long, long>> getBadDimensionParams(long min_m,
                                                         long max_m,
                                                         long min_p,
                                                         long max_p,
                                                         long m_sparseness,
                                                         long p_sparseness)
{
  return getParams(false,
                   min_m,
                   max_m,
                   min_p,
                   max_p,
                   m_sparseness,
                   p_sparseness);
}

std::vector<std::pair<long, long>> getGoodDimensionParams(long min_m,
                                                          long max_m,
                                                          long min_p,
                                                          long max_p,
                                                          long m_sparseness,
                                                          long p_sparseness)
{
  return getParams(true,
                   min_m,
                   max_m,
                   min_p,
                   max_p,
                   m_sparseness,
                   p_sparseness);
}
} // namespace helib_test
