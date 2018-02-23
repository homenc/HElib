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

// returns new v[i] = \sum_{j>=i} old v[i]
void runningSums(CtPtrs& v)
{
  FHE_TIMER_START;
  for (long i=lsize(v)-1; i>0; i--) *v[i-1] += *v[i];
}


// a recursive function that computes
//      e*[i] = prod_{j>=i} e[i]  and  g*[i] = e*[i+1] \cdot g[i]
// This function is optimized, so that instead of all the e*[i]'s
// it only computes e*[0] and the e*[i]'s that are used in the
// computation of the g*[i]'d
static void compProducts(const CtPtrs_slice& e, const CtPtrs_slice& g)
{
  long n = lsize(e);
  if (n <= 1) return; // nothing to do
#ifdef DEBUG_PRINTOUT
  cout << "compProducts(g["<<g.start<<".."<<(g.start+g.sz-1)<<"],e["
       << e.start<<".."<<(e.start+e.sz-1)<<"])\n";
#endif

  // split the array in two, second part has size the largest 2^l < n,
  // and first part is the rest

  long ell = NTL::NumBits(n-1) -1; // n/2 <= 2^l < n
  long n1 = n - (1UL<<ell);        // n1 \in [1, n/2]

  // Call the recursive procedure separately on the first and second parts
  compProducts(CtPtrs_slice(e,0,n1), CtPtrs_slice(g,0,n1));      // first half
  compProducts(CtPtrs_slice(e,n1,n-n1), CtPtrs_slice(g,n1,n-n1));// second half

  // Multiply the first product in the 2nd part into every product in the 1st
  NTL_EXEC_RANGE(1+n1, first, last)
  for (long i=first; i<last; i++) {
    if (i==0)               e[0]->multiplyBy(*e[n1]);
    else if (i-1<g.size()) g[i-1]->multiplyBy(*e[n1]);
  }
  NTL_EXEC_RANGE_END
#ifdef DEBUG_PRINTOUT
  cout << " g["<<g.start<<".."<<(g.start+g.sz-1)<<"], "
       << " e["<<e.start<<".."<<(e.start+e.sz-1)<<"]:\n";
  for (long i=0; i<g.size(); i++)
    decryptAndPrint((cout<<"   g["<<(i+g.start)<<"] ("<<((void*)g[i])<<"): "),
                    *g[i], *dbgKey, *dbgEa, FLAG_PRINT_POLY);
  for (long i=0; i<e.size(); i++)
    decryptAndPrint((cout<<"   e["<<(i+e.start)<<"] ("<<((void*)e[i])<<"): "),
                    *e[i], *dbgKey, *dbgEa, FLAG_PRINT_POLY);

  cout << endl;
#endif
}

// Compute aeqb[i] = (a==b upto bit i), agtb[i] = (aeqb[i+1] and ai>bi)
// We assume that b.size()>a.size()
static void
compEqGt(CtPtrs& aeqb, CtPtrs& agtb, const CtPtrs& a, const CtPtrs& b)
{
  FHE_TIMER_START;
  const Ctxt zeroCtxt(ZeroCtxtLike, *(b.ptr2nonNull()));
  DoubleCRT one(zeroCtxt.getContext()); one += 1L;
  
  resize(aeqb, lsize(b), zeroCtxt);
  resize(agtb, lsize(a), zeroCtxt);

  // First compute the local bits e[i]=(a[i]==b[i]), gt[i]=(a[i]>b[i])
  FHE_NTIMER_START(compEqGt1);
  long aSize = lsize(a);
  NTL_EXEC_RANGE(aSize, first, last)
  for (long i=first; i<last; i++) {
    *aeqb[i] = *b[i];               // b
    aeqb[i]->addConstant(one, 1.0); // b+1
    *agtb[i] = *aeqb[i];            // b+1
    *aeqb[i] += *a[i];              // a+b+1
    agtb[i]->multiplyBy(*a[i]);     // a(b+1)
  }
  NTL_EXEC_RANGE_END
  FHE_NTIMER_STOP(compEqGt1);

  // NOTE: Usually there isn't much gain in multi-threading the loop below,
  //    but computing b[i] can be expensive in some implementations of CtPtrs
  FHE_NTIMER_START(compEqGt2);
  if (lsize(b)-aSize >1) {
    NTL_EXEC_RANGE(lsize(b)-aSize, first, last)
    for (long i=first; i<last; i++) {
      *aeqb[i+aSize] = *b[i+aSize];         // b
      aeqb[i+aSize]->addConstant(one, 1.0); // b+1
    }
    NTL_EXEC_RANGE_END
  }
  else if (lsize(b)-aSize == 1) {
    *aeqb[aSize] = *b[aSize];         // b
    aeqb[aSize]->addConstant(one, 1.0); // b+1
  }
  FHE_NTIMER_STOP(compEqGt2);

#ifdef DEBUG_PRINTOUT
  for (long i=0; i<lsize(b); i++)
    decryptAndPrint((cout<<" e["<<i<<"]: "), *aeqb[i], *dbgKey, *dbgEa, FLAG_PRINT_POLY);
  for (long i=0; i<lsize(a); i++)
    decryptAndPrint((cout<<" ag["<<i<<"]: "), *agtb[i], *dbgKey, *dbgEa, FLAG_PRINT_POLY);
  cout << endl;
#endif

  // Call a recursive function to compute:
  // e*_i = \prod_{j>=i} aeqb_i, g*_i = aeqb*_{i+1} \cdot agtb_i
  FHE_NTIMER_START(compEqGt3);
  compProducts(CtPtrs_slice(aeqb,0), CtPtrs_slice(agtb,0));
  runningSums(agtb); // now ag[i] = (a>b upto bit i)
  FHE_NTIMER_STOP(compEqGt3);
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

  // Check that we have enough levels, try to bootstrap otherwise
  if (findMinLevel({&a,&b}) < NTL::NumBits(bSize+1)+2)
    packedRecrypt(a,b,unpackSlotEncoding);
  if (findMinLevel({&a,&b}) < NTL::NumBits(bSize)+1) // the bear minimum
    throw std::logic_error("not enough levels for comparison");

  // NOTE: this procedure minimizes the number of multiplications,
  //       but it may use one level too many. Can we optimize it?

  /* We first compute for each position i the values
   *   e[i] = (a==b upto position i)
   *   ag[i] = (a>b upto position i)
   */

  // We use max, min to hold the intermediate values e, ag
  CtPtrs& e = max;
  CtPtrs& ag = min;
  compEqGt(e, ag, a, b);

  // We are now ready to compute the bits of the result.

  FHE_NTIMER_START(compResults);
  mu = *ag[0];             // a > b
  ni = *ag[0];
  ni.addConstant(ZZ(1L));  // a <= b
  ni += *e[0];             // a < b

  NTL_EXEC_RANGE(aSize, first, last)
  for (long i=first; i<last; i++) {
    *max[i] = *a[i];
    *max[i] -= *b[i];
    max[i]->multiplyBy(*ag[i]);

    *min[i] = *max[i];
    *max[i] += *b[i];
    *min[i] -= *a[i];
  }
  NTL_EXEC_RANGE_END
  for (long i=aSize; i<bSize; i++)
    *max[i] = *b[i];
  FHE_NTIMER_STOP(compResults);
}
