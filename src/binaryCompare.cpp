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
/**
 * @file binaryCompare.cpp
 * @brief Implementing integer comparison in binary representation.
 */
#include <algorithm>
#include <cassert>
#include <stdexcept>
// #include <numeric>
// #include <climits>
// #include <map>
// #include <atomic>
// #include <mutex>          // std::mutex, std::unique_lock

#include <NTL/BasicThreadPool.h>
#include "binaryArith.h"

#ifdef DEBUG_PRINTOUT
#include "debugging.h"
#endif

// returns new v[i] = \sum_{j=0}^i old v[i]
void runningSums(std::vector<Ctxt>& v)
{
  for (long i=1; i<lsize(v); i++) v[i] += v[i-1];
}

// Compares two integers in binary a,b.
// Returns max(a,b), min(a,b) and indicator bits mu=(a>b) and ni=(a<b)
void compareTwoNumbers(CtPtrs& max, CtPtrs& min, Ctxt& mu, Ctxt& ni,
                       const CtPtrs& aa, const CtPtrs& bb,
                       std::vector<zzX>* unpackSlotEncoding)
{
  FHE_TIMER_START;
  // make sure that lsize(b) >= lsize(a)
  const CtPtrs& a = (lsize(bb)>=lsize(aa))? aa : bb;
  const CtPtrs& b = (lsize(bb)>=lsize(aa))? bb : aa;
  long aSize = lsize(a);
  long bSize = lsize(b);
  if (aSize<1) { // a is empty
    mu.clear();
    ni.clear();
    ni.addConstant(ZZ(1L));
    vecCopy(max, b);
    setLengthZero(min);
    return;
  }

  // Begin by checking that we have enough levels
  if (findMinLevel({&a,&b}) < NTL::NumBits(long(bSize))+2)
    packedRecrypt(a,b,unpackSlotEncoding);
  if (findMinLevel({&a,&b}) < NTL::NumBits(long(bSize))+1)
    throw std::logic_error("not enough levels for comparison");

  // NOTE: this procedure minimizes the number of multiplications,
  //       but it uses 1-2 level too many. Can we optimize it?

  /* We first compute for each position i the value
   *    e_i = (a[i]==b[i]) = (a[i]+b[i]+1),
   * NOTE: e must be a std::vector in reverse order, so we
   *       can later use the incrementalProduct function of HElib.
   */

  const Ctxt zeroCtxt(ZeroCtxtLike, *(b.ptr2nonNull()));
  DoubleCRT one(zeroCtxt.getContext()); one += 1L;

  std::vector<Ctxt> e(bSize, zeroCtxt);
  for (long i=0; i<bSize; i++) {
    e[bSize-i-1] = *b[i]; // e = rev(b)
    e[bSize-i-1].addConstant(one, 1.0);
    if (i<aSize)
      e[bSize-i-1] += *a[i];
  }
  // Compute e*_i = prod_{j=i}^{bSize-1} e_j
  incrementalProduct(e); // e[bSize-i-1] = e*_i

  // Now compute the bits ag[bSize-i-1] = (ai>bi, and a==b upto bit i+1)
  std::vector<Ctxt> ag = e;
  //  for (long i=0; i<bSize; i++) {
  NTL_EXEC_RANGE(bSize, first, last)
  for (long i=first; i<last; i++) {
    if (i<aSize) {
      ag[bSize-i-1] = *a[i];
      ag[bSize-i-1].multiplyBy(*b[i]);
      ag[bSize-i-1] -= *a[i]; // now ag[bSize-i-1] = ai(bi-1) = (ai>bi)
      if (i<bSize-1)          // multiply by e*_{i+1}
        ag[bSize-i-1].multiplyBy(e[bSize-i-2]);
    }
    else ag[bSize-i-1].clear();
  }
  NTL_EXEC_RANGE_END
  runningSums(ag); // now ag[bSize-i-1] = (a>b upto bit i)

  // We are now ready to compute the bits of the result.

  mu = ag[bSize-1];  // a>b
  ni = ag[bSize-1];
  ni.addConstant(one, 1.0); // a <= b
  ni += e[bSize-1];         // a < b

  resize(max, bSize, zeroCtxt); // ensure enough space
  resize(min, aSize, zeroCtxt); // ensure enough space

  //  for (long i=0; i<bSize; i++) {
  NTL_EXEC_RANGE(bSize, first, last)
  for (long i=first; i<last; i++) {
    *max[i] = *b[i];
    if (i<aSize) {
      max[i]->multiplyBy(ag[bSize-i-1]);
      *max[i] -= *b[i];      // b[i] * (b>=a upto bit i)
      Ctxt tmp = ag[bSize-i-1];
      tmp.multiplyBy(*a[i]); // a[i] * (a>b upto bit i)
      *max[i] += tmp;

      *min[i] = *a[i];
      min[i]->multiplyBy(ag[bSize-i-1]);
      *min[i] -= *a[i];      // a[i] * (b>=a upto bit i)
      tmp = ag[bSize-i-1];
      tmp.multiplyBy(*b[i]); // b[i] * (a>b upto bit i)
      *min[i] += tmp;
    }
  }
  NTL_EXEC_RANGE_END
}
