/* Copyright (C) 2012-2020 IBM Corp.
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
 * @file tableLookup.cpp
 * @brief Code for homomorphic table lookup and fixed-point functions
 */
#include <limits>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <NTL/BasicThreadPool.h>
#include <helib/intraSlot.h>
#include <helib/tableLookup.h>

#ifdef HELIB_DEBUG
#include <helib/debugging.h>
#endif

namespace helib {

static void recursiveProducts(const CtPtrs& products,
                              const CtPtrs_slice& array);
static double pow2_double(long n); // compute 2^n as double

// For an n-size array, compute the 2^n products
//     products[j] = \prod_{i s.t. j_i=1} array[i]
//                   \times \prod_{i s.t. j_i=0}(a-array[i])
void computeAllProducts(/*Output*/ CtPtrs& products,
                        /*Index*/ const CtPtrs& array,
                        std::vector<zzX>* unpackSlotEncoding)
{
  HELIB_TIMER_START;
  long nBits = array.size();
  if (lsize(products) > 0) {
    long nBits2 = NTL::NumBits(lsize(products) - 1); // ceil(log_2(size))
    if (nBits > nBits2)
      nBits = nBits2; // ignore extra bits in 'array'
  }
  if (nBits < 1)
    return; // do nothing
  // Output cannot be bigger than 2^16
  assertTrue(nBits <= 16, "Output cannot be bigger than 2^16");

  if (lsize(products) == 0) // try to set the output size
    products.resize(1L << nBits, &array);
  for (long i = 0; i < lsize(products); i++)
    products[i]->clear();

  // Check that we have enough levels, try to bootstrap otherwise
  assertNotNull(array.ptr2nonNull(),
                "Invalid array (could not find non-null Ctxt)");
  long bpl = array.ptr2nonNull()->getContext().BPL();
  if (findMinBitCapacity(array) < (NTL::NumBits(nBits) + 1) * bpl) {
    const Ctxt* ct = array.ptr2nonNull(); // find some non-null Ctxt
    assertNotNull(unpackSlotEncoding,
                  "unpackSlotEncoding must not be null when bootstrapping");
    assertTrue(ct->getPubKey().isBootstrappable(),
               "Cannot bootstrap with non-bootstrappable public key");
    packedRecrypt(array,
                  *unpackSlotEncoding,
                  *(ct->getContext().ea),
                  /*belowLevel=*/nBits + 3);
  }
  if (findMinBitCapacity(array) < (NTL::NumBits(nBits) + 1) * bpl)
    throw LogicError("not enough levels for table lookup");

  // Call the recursive function that computes the products
  recursiveProducts(products, CtPtrs_slice(array, 0, nBits));
}

// The input is a plaintext table T[] and an array of encrypted bits
// I[], holding the binary representation of an index i into T.
// The output is the encrypted value T[i].
void tableLookup(Ctxt& out,
                 const std::vector<zzX>& table,
                 const CtPtrs& idx,
                 std::vector<zzX>* unpackSlotEncoding)
{
  HELIB_TIMER_START;
  out.clear();
  std::vector<Ctxt> products(lsize(table),
                             out); // to hold subset products of idx
  CtPtrs_vectorCt pWrap(products); // A wrapper

  // Compute all products of encrypted bits =: b_i
  computeAllProducts(pWrap, idx, unpackSlotEncoding);

  // Compute the sum b_i * T[i]
  NTL_EXEC_RANGE(lsize(table), first, last)
  for (long i = first; i < last; i++)
    products[i].multByConstant(table[i]); // p[i] = p[i]*T[i]
  NTL_EXEC_RANGE_END
  for (long i = 0; i < lsize(table); i++)
    out += products[i];
}

// A counterpart of tableLookup. The input is an encrypted table T[]
// and an array of encrypted bits I[], holding the binary representation
// of an index i into T.  This function increments by one the entry T[i].
void tableWriteIn(const CtPtrs& table,
                  const CtPtrs& idx,
                  std::vector<zzX>* unpackSlotEncoding)
{
  HELIB_TIMER_START;
  const Ctxt* ct = table.ptr2nonNull(); // find some non-null Ctxt
  long size = lsize(table);
  if (size == 0)
    return;
  std::vector<Ctxt> products(size, Ctxt(ZeroCtxtLike, *ct));
  CtPtrs_vectorCt pWrap(products); // A wrapper

  // Compute all products of encrypted bits =: b_i
  computeAllProducts(pWrap, idx, unpackSlotEncoding);

  // increment each entry of T[i] by products[i]
  NTL_EXEC_RANGE(lsize(table), first, last)
  for (long i = first; i < last; i++)
    *table[i] += products[i];
  NTL_EXEC_RANGE_END
}

// The function buildLookupTable is documented in tableLookup.h.
// The output is returned in T, size of T will be 2^{nbits_in}.
// For every signed integer x with bit-size 'nbits_in', we will have
//     T[x] = f(x * 2^{scale_in}) * 2^{-scale_out}),
// rounded to the nearest integer and truncated to 'nbits_out' bits.
// The bits are packed inside the slots, so it is assumed that each
// slot has enough room to fit these many bits. (Otherwise we only
// keep as many low-order bits as fit in a slot.)
void buildLookupTable(std::vector<zzX>& T, // encoded result is returned in T
                      std::function<double(double)> f,
                      long nbits_in, // number of precision bits
                      long scale_in, // scaling factor
                      long sign_in,  // 1: 2's complement signed, 0: unsigned
                      long nbits_out,
                      long scale_out,
                      long sign_out,
                      const EncryptedArray& ea)
{
  HELIB_TIMER_START;
  // tables of size > 2^{16} are not supported
  assertTrue(nbits_in <= 16, "tables of size > 2^{16} are not supported");
  long sz = 1L << nbits_in;
  T.resize(sz);

  double pow2_scale_in = pow2_double(scale_in);        // 2^{nbits_in}
  double pow2_neg_scale_out = pow2_double(-scale_out); // 2^{-nbits_out}

  // Compute the largest and smallest values that can be in T
  long largest_value, smallest_value;
  if (sign_out) { // values in T are encoded in 2's complement
    largest_value = (1L << (nbits_out - 1)) - 1;
    smallest_value = -(1L << (nbits_out - 1));
  } else { // values in T are all non-negative
    largest_value = (1L << nbits_out) - 1;
    smallest_value = 0;
  }

  for (long i = 0; i < sz; i++) { // Compute the entries of T
    long x;
    if (sign_in) { // indexes into T are in 2's complement
      long sign_bit = (1L << (nbits_in - 1)) & i;
      x = i - 2 * sign_bit;
    } else
      x = i; // indexes into T are all non-negative

    // Compute the value that should go in the table as rounded double
    double scaled_x = double(x) * pow2_scale_in;
    double y = round(f(scaled_x) * pow2_neg_scale_out);

    // saturated arithmetic (set to smallest or largest values)
    // NOTE: this should work fine even if y is an infinity
    long value;
    if (std::isnan(y))
      value = 0;
    else if (y > largest_value)
      value = largest_value;
    else if (y < smallest_value)
      value = smallest_value;
    else
      value = y;

    // convert to unsigned and mask to nbits_out bits
    unsigned long uvalue = value;
    uvalue &= ((1UL << nbits_out) - 1UL); // keep only bottom nbits_out bits

    packConstant(T[i], uvalue, nbits_out, ea);
  }
}

// A recursive function to compute, for an n-size array, the 2^n products
//     products[j] = \prod_{i s.t. j_i=1} array[i]
//                   \times \prod_{i s.t. j_i=0}(a-array[i])
// It is assume that 'products' size <= 2^n, else only 1st 2^n entries are set
static void recursiveProducts(const CtPtrs& products, const CtPtrs_slice& array)
{
  long nBits = lsize(array);
  long N = lsize(products);
  if (nBits == 0 || N == 0)
    return; // nothing to do

  if (N > (1L << nBits))
    N = (1L << nBits);
  else if (N < (1L << (nBits - 1)))
    nBits = NTL::NumBits(N - 1); // Ensure nBits <= ceil(log2(N))

  if (N <= 2) { // edge condition
    *products[0] = *array[0];
    products[0]->negate();
    products[0]->addConstant(NTL::ZZ(1)); // out[0] = 1-in
    if (N > 1)
      *products[1] = *array[0]; // out[1] = in
  }
  // optimization for n=2: a single multiplication instead of 4
  else if (N <= 4) {
    *products[0] = *array[1];           // x1
    products[0]->multiplyBy(*array[0]); // x1 x0

    *products[1] = *array[0];     // x0
    *products[1] -= *products[0]; // x0 - x1 x0 = (1-x1)x0

    *products[2] = *array[1];     // x1
    *products[2] -= *products[0]; // x1 - x1 x0 = x1(1-x0)

    if (N > 3)
      *products[3] = *products[0]; // x1 x0

    products[0]->addConstant(NTL::ZZ(1)); // 1 +x1 x0
    *products[0] -= *array[1];            // 1 +x1 x0 -x1
    *products[0] -= *array[0];            // 1 +x1 x0 -x1 -x0 = (1-x1)(1-x0)
  } else {                                // split the array into two parts;
    // first part is highest pow(2) < n, second part is what is left

    long n1 = 1L << (NTL::NumBits(nBits) - 1); // largest power of two <= n
    if (nBits <= n1)
      n1 = n1 / 2;               // largest power of two < n
    long k = 1L << n1;           // size of first part
    long l = 1L << (nBits - n1); // size of second part

    const Ctxt* ct = array.ptr2nonNull(); // find some non-null Ctxt
    std::vector<Ctxt> products1(k, Ctxt(ZeroCtxtLike, *ct));
    std::vector<Ctxt> products2(l, Ctxt(ZeroCtxtLike, *ct));

    // compute first part of the array
    recursiveProducts(CtPtrs_vectorCt(products1), CtPtrs_slice(array, 0, n1));

    // recursive call on second part of array
    recursiveProducts(CtPtrs_vectorCt(products2),
                      CtPtrs_slice(array, n1, nBits - n1));

    // multiplication to get all subset products
    NTL_EXEC_RANGE(lsize(products), first, last)
    for (long ii = first; ii < last; ii++) {
      long j = ii / k;
      long i = ii - j * k;
      *products[ii] = products1[i];
      products[ii]->multiplyBy(products2[j]);
    }
    NTL_EXEC_RANGE_END
  }
}

static double pow2_double(long n) // compute 2^n as double
{
  double res = 1;
  long abs_n = std::labs(n);

  for (long i = 0; i < abs_n; i++)
    res *= 2;
  if (n < 0)
    res = 1 / res;
  return res;
}

} // namespace helib
