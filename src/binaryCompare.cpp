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

  /* We first compute for each position i the local min and max:
   *    e_i = (a[i]==b[i]) = (a[i]+b[i]+1),
   *    f_i = (a[i]>b[i])  = a[i]*(b[i]+1), and
   *    g_i = (a[i]<b[i])  = b[i]*(a[i]+1).
   * NOTE: e must be a std::vector in reverse order, so we
   *       can use the incrementalProduct function of HElib.
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
  // NOTE: this procedure minimizes the number of multiplications,
  //       but it uses 1-2 level too many. Can we optimize it?
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
  mu = ag[bSize-1];  // a>b
  ni = ag[bSize-1];
  ni.addConstant(one, 1.0); // a <= b
  ni += e[bSize-1];         // a < b
}

// unused code below
#if 0
void basicMaxStep(NTL::Vec<Ctxt>& a, NTL::Vec<Ctxt>& data,
                  const vector< vector< GF2X > >& maskTable,
                  int dim, int ord, int i)
{
  const FHEcontext& context = a[0].getContext();
  const EncryptedArray& ea = *(context.ea);
  int aSize = lsize(a);
  int dataSize = lsize(data);

  // Build the masks M_i which have 1's in slots of the range [2^i, 2^{i+1}-1]
  int mIndex1 = pow(2,i);
  int mIndex2 = (ord>mIndex1*2)? mIndex1*2 : ord;
  const GF2X mt1 = maskTable[dim][mIndex1] - maskTable[dim][mIndex2];
  DoubleCRT m1(convert<zzX,GF2X>(mt1), context, a[0].getPrimeSet());

  cout << "basicMaxStep(dim="<<dim<<", amt="<<mIndex1<<")\n";

  // Begin by checking that we have enough levels
  int lvl = 100;
  for (long i=0; i<dataSize; i++) {
    int dlvl = data[i].findBaseLevel();
    if (dlvl < lvl) lvl = dlvl;
  }
  if (lvl < NTL::NumBits(long(aSize))+3) {
    packedRecrypt(a,data);
    if (a[0].findBaseLevel() < NTL::NumBits(long(aSize))+1)
    throw std::logic_error("not enough levels for basicMaxStep");
  }

  // Imitate power-of-two rotations in a non-power-of-two dimension
  NTL::Vec<Ctxt> a2 = a;       // copy of a
  for (int j=0; j<aSize; j++) {
    Ctxt tmp = a2[j];
    a2[j].multByConstant(m1);
    ea.rotate1D(a2[j], dim, -mIndex1); // keeps slots >2^i, shift left
    ea.rotate1D(tmp, dim, mIndex1);    // keeps slots <2^i, shift right
    tmp.multByConstant(m1);
    a2[j] += tmp;
  }
#if DEBUGGING
  vector<long> nums1, nums2;
  if (verbose) {
    decryptBinaryNums(nums1, a);
    decryptBinaryNums(nums2, a2);
    cout << " nums1="<<nums1<<endl;
    cout << " nums2="<<nums2<<endl;
  }
#endif

  // Compare the numbers to find the max and min
  Ctxt ni(ZeroCtxtLike, a[0]);
  {
    Ctxt mu(ZeroCtxtLike, a[0]);
    NTL::Vec<Ctxt> max(INIT_SIZE, aSize, ni);
    NTL::Vec<Ctxt> min(INIT_SIZE, aSize, ni);

    compTwoNumbers(max, min, mu, ni, a, a2);
#if DEBUGGING
  if (verbose) {
    decryptBinaryNums(nums1, max);
    decryptBinaryNums(nums2, min);
    cout << "   max="<<nums1<<endl;
    cout << "   min="<<nums2<<endl;
  }
#endif
    // Moves the larger numbers to slots closer to slot 0
    for (int j=0; j<aSize; j++) {
      a[j] = max[j];             // max
      max[j] -= min[j];          // max - min
      max[j].multByConstant(m1); // (max-min)*m1
      a[j] += max[j];            // max +(max-min)*m1 = max(1-m1) +min*m1
    }
#if DEBUGGING
    if (verbose) {
      decryptBinaryNums(nums1, a);
      cout << " merge="<<nums1<<endl;
    }
#endif
  } // release memory for min, max, mu

  a2 = data;       // copy of data
  for (int j=0; j<dataSize; j++) {
    Ctxt tmp = a2[j];
    a2[j].multByConstant(m1);
    ea.rotate1D(a2[j], dim, -mIndex1);
    ea.rotate1D(tmp, dim, mIndex1);
    tmp.multByConstant(m1);
    a2[j] += tmp;
  }

  // reassemble the associated data

  // replicate the bit ni across last dimension
  replicate1D(ea, ni, /*dim=*/ea.dimension()-1, /*pos=*/0);
  ni.addConstant(ZZ(1L));
  ni.multByConstant(m1); // (ni+1)*m1
  Ctxt tmp = ni;
  ea.rotate1D(tmp, dim, -mIndex1);
  ni += tmp;             // ni' = (ni+1)*m1 + ((ni+1)*m1)>>2^i
  
  for (int i=0; i<dataSize; i++) {
    a2[i] -= data[i];      // a2 - data
    a2[i].multiplyBy(ni);  // (a2-data)*ni'
    data[i] += a2[i];      // data +(a2-data)*ni' = data*(1-ni') +a2*ni'
  }
  //  decryptAndPrint(cout<<"   after assembling, data[0]=", data[0], *dbgKey, ea);
  //  testEncIndex(data[0], 104);
}

// Computes find-max-with-associated-data.
// Takes a packed ctxt of N binary numbers, finds the largest
// integer and moves it to the first ctxt slot.
// It moves the corresponding data to the same slots as its
// associated integer after the max is found.
// NOTE: This procedure only works when all values of our vector
//       are non-negative
// NOTE: I have only coded it for dimensions with order not
//       equal to a power of two
void findMax(NTL::Vec<Ctxt>& a, NTL::Vec<Ctxt>& data, long n)
{
  const FHEcontext& context = a[0].getContext();
  const EncryptedArray& ea = *(context.ea);
  int ord0 = ea.sizeOfDimension(0); // Order of generator 0
  int ell0 = NumBits(ord0-1);

  const EncryptedArrayDerived<PA_GF2>& ea2 = ea.getDerived(PA_GF2());
  const vector< vector< GF2X > >& maskTable = ea2.getTab().getMaskTable();

  if (n>2*ord0) // If more elements than twice the order of dimension 0
    basicMaxStep(a, data, maskTable, 1, ea.sizeOfDimension(1), 1);

  if (n>ord0)   // If more elements than the order of dimension 0
    basicMaxStep(a, data, maskTable, 1, ea.sizeOfDimension(1), 0);

  for (int i=ell0-1; i>=0; i--) // max along dimension 0
    basicMaxStep(a, data, maskTable, 0, ord0, i);
}

// Find the k largest integers from a set of n integers in binary
// representation using repeated calls of the findMax procedure
void findLargestK(NTL::Vec<NTL::Vec<Ctxt> >& result, NTL::Vec<Ctxt>& a,
                  NTL::Vec<Ctxt>& data, long n, long k)
{
  // Make sure we are not looking for more numbers than we have
  assert (k<=n);

  const FHEcontext& context = a[0].getContext();
  const EncryptedArray& ea = *(context.ea);

  // find the max number then zero out the result
  // and find the next largest integer until k
  // have been found

  ZZX zzMask;
  ea.encodeUnitSelector(zzMask, 0); // mask with 1 in slot 0
  const DoubleCRT mask(zzMask, context);

  // Make k find max calls
  result.SetLength(0);
#if DEBUGGING
  cout << " b4 findLargestK ";
  testEncIndex(data[0], n);
#endif
  for (int i=0; i<k; i++) {
#if DEBUGGING
    auto a2 = a;
    auto data2 = data;
    vector<long> nums;
    decryptBinaryNums(nums, a);
#endif
    findMax(a, data, n); // find largest element in a, rotate data accordingly
#if DEBUGGING
    vector<long> nums2;
    decryptBinaryNums(nums2, a);
    // cout << "  nums ="<<nums<<endl;
    // cout << "  nums2="<<nums2<<endl;
    cout << "- findMax #"<<i<<" ";

    long maxNum = *std::max_element(nums.begin(), nums.end());
    if (maxNum != nums2[0]) {
      cout << "  OOPS, max should be "<<maxNum<<" but got "<<nums2[0]<<endl;
      verbose = true;
      findMax(a2, data2, n);
      exit(0);
    }
    else
      cout << "  yay, max is indeed "<<maxNum<<endl;
    testEncIndex(data[0], n);
#endif
    result.append(data); // j=0 is the column of the max element
    // {
    //   ofstream ofs("./index_32767_"+to_string(i)+".dat",
    //                std::ofstream::out | std::ofstream::trunc);
    //   if (ofs.is_open()) ofs<<data[0]<<endl;
    // }
    for (int j=0; j<lsize(a); j++) { // zero out the largest element in a
      Ctxt tmp = a[j];           // copy of a
      tmp.multByConstant(mask);  // keep just slot 0
      a[j] -= tmp;               // all slots except slot 0
    }
  }
}
#endif // unused code
