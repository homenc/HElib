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

#ifndef HELIB_SET_H
#define HELIB_SET_H

#include <helib/SumRegister.h>
#include <NTL/BasicThreadPool.h>

namespace helib {

// Binary tree Summation of Ctxts.
// Destructive but more efficient algorithm.
// TODO: Can write a generic binSum to work with pointers
// and make this a vector of unique pointers
/**
 * @brief Performs a binary summation of a vector of elements.
 * @tparam TXT type of the elements of which to sum.
 * @param ctxtArray The array on which to perform the binary sum.
 * @note This function is destructive on the array.
 **/
template <typename TXT>
inline void binSumReduction(std::vector<TXT>& ctxtArray)
{
  int cnt = 0;
  int end = ctxtArray.size();

  while (end > 1) {
    ++cnt;
    int odd = end & 1;
    int comps = end >> 1;

    NTL_EXEC_RANGE(comps, first, last)
    for (unsigned i = first; i < last; i++) {
      ctxtArray.at(i) += ctxtArray.at(comps + i + odd);
    }
    NTL_EXEC_RANGE_END

    // Free the end of the vector.
    ctxtArray.erase(ctxtArray.begin() + comps + odd, ctxtArray.end());
    ctxtArray.shrink_to_fit();

    // Update end.
    end = (end + odd) >> 1;
  }
}

/**
 * @brief Given two sets, calculates and returns the set intersection.
 * @tparam TXT type of the query set. Must be a `Ptxt` or `Ctxt`.
 * @param query The query set of type `TXT` where the elements of the set are
 * held in the slots.
 * @param server_set The server set. A vector of integer polynomials.
 * @return A set of the same size as `query` holding the elements in the
 * intersecting set.
 **/
template <typename TXT>
inline TXT calculateSetIntersection(const TXT& query,
                                    const std::vector<NTL::ZZX>& server_set)
{
  long availableThreads =
      std::min(NTL::AvailableThreads(), long(server_set.size()));
  std::vector<TXT> interResult(availableThreads, query);

  NTL::PartitionInfo pinfo(server_set.size());

  NTL_EXEC_INDEX(availableThreads, index)
  long first, last;
  pinfo.interval(first, last, index);
  SumRegister<TXT> sumRegister(last - first);

  for (long i = first; i < last; ++i) {
    auto lquery = std::make_unique<TXT>(query);
    Ptxt<BGV> entry(query.getContext(), server_set[i]);
    *lquery -= entry;
    mapTo01(*query.getContext().ea, *lquery);
    lquery->negate();
    lquery->addConstant(NTL::ZZX(1L));
    sumRegister.add(lquery);
  }

  sumRegister.flush();

  assertTrue(sumRegister.hasResult(), "Sum Register did not have a result.");

  interResult.at(index) = *(sumRegister.getResult());
  NTL_EXEC_INDEX_END

  // Final binary sum to add the results of the sum registers
  binSumReduction<TXT>(interResult);
  return interResult.at(0) *= query;
}

} // namespace helib

#endif
