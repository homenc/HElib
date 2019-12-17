/* Copyright (C) 2012-2019 IBM Corp.
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
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <NTL/BasicThreadPool.h>

#include <helib/helib.h>

#include <helib/intraSlot.h>
#include <helib/binaryArith.h>

#include "gtest/gtest.h"
#include "test_common.h"

#ifdef DEBUG_PRINTOUT
#include <helib/debugging.h>
#endif
// define flags FLAG_PRINT_ZZX, FLAG_PRINT_POLY, FLAG_PRINT_VEC, functions
//        decryptAndPrint(ostream, ctxt, sk, ea, flags)
//        decryptAndCompare(ctxt, sk, ea, pa);

namespace {

struct Parameters
{
  Parameters(long prm,
             long bitSize,
             long bitSize2,
             long outSize,
             bool bootstrap,
             long seed,
             long nthreads) :
      prm(prm),
      bitSize(bitSize),
      bitSize2(bitSize2),
      outSize(outSize),
      bootstrap(bootstrap),
      seed(seed),
      nthreads(nthreads){};

  long prm;       // parameter size (0-tiny,...,7-huge)
  long bitSize;   // bitSize of input integers (<=32)
  long bitSize2;  // bitSize of 2nd input integer (<=32)
  long outSize;   // bitSize of output integers, as many as needed
  bool bootstrap; // test multiplication with bootstrapping
  long seed;      // PRG seed
  long nthreads;  // number of threads

  friend std::ostream& operator<<(std::ostream& os, const Parameters& params)
  {
    return os << "{"
              << "prm=" << params.prm << ","
              << "bitSize=" << params.bitSize << ","
              << "bitSize2=" << params.bitSize2 << ","
              << "outSize=" << params.outSize << ","
              << "bootstrap=" << params.bootstrap << ","
              << "seed=" << params.seed << ","
              << "nthreads=" << params.nthreads << "}";
  };
};

class GTestBinaryArith :
    public ::testing::TestWithParam<std::tuple<Parameters, int>>
{
protected:
  static std::vector<helib::zzX> unpackSlotEncoding;
  constexpr static long mValues[8][15] = {
      // { p, phi(m),   m,   d, m1, m2, m3,    g1,   g2,   g3, ord1,ord2,ord3,
      // B,c}
      {2, 48, 105, 12, 3, 35, 0, 71, 76, 0, 2, 2, 0, 25, 2},
      {2, 600, 1023, 10, 11, 93, 0, 838, 584, 0, 10, 6, 0, 25, 2},
      {2, 2304, 4641, 24, 7, 3, 221, 3979, 3095, 3760, 6, 2, -8, 25, 3},
      {2, 5460, 8193, 26, 8193, 0, 0, 46, 0, 0, 210, 0, 0, 25, 3},
      {2, 8190, 8191, 13, 8191, 0, 0, 39, 0, 0, 630, 0, 0, 25, 3},
      {2, 10752, 11441, 48, 17, 673, 0, 4712, 2024, 0, 16, -14, 0, 25, 3},
      {2, 15004, 15709, 22, 23, 683, 0, 4099, 13663, 0, 22, 31, 0, 25, 3},
      {2,
       27000,
       32767,
       15,
       31,
       7,
       151,
       11628,
       28087,
       25824,
       30,
       6,
       -10,
       28,
       4}};

  static long correctBitSize(long minimum, long oldBitSize)
  {
    long newBitSize;
    if (oldBitSize <= 0)
      newBitSize = minimum;
    else if (oldBitSize > 32)
      newBitSize = 32;
    else
      newBitSize = oldBitSize;
    return newBitSize;
  };

  // Validates the prm value, throwing if invalid
  static long validatePrm(long prm)
  {
    if (prm < 0 || prm >= 5)
      throw std::invalid_argument("prm must be in the interval [0, 4]");
    return prm;
  };

  static NTL::Vec<long> calculateMvec(const long* vals)
  {
    NTL::Vec<long> mvec;
    NTL::append(mvec, vals[4]);
    if (vals[5] > 1)
      NTL::append(mvec, vals[5]);
    if (vals[6] > 1)
      NTL::append(mvec, vals[6]);
    return mvec;
  };

  static std::vector<long> calculateGens(const long* vals)
  {
    std::vector<long> gens;
    gens.push_back(vals[7]);
    if (vals[8] > 1)
      gens.push_back(vals[8]);
    if (vals[9] > 1)
      gens.push_back(vals[9]);
    return gens;
  };

  static std::vector<long> calculateOrds(const long* vals)
  {
    std::vector<long> ords;
    ords.push_back(vals[10]);
    if (abs(vals[11]) > 1)
      ords.push_back(vals[11]);
    if (abs(vals[12]) > 1)
      ords.push_back(vals[12]);
    return ords;
  };

  static long calculateLevels(bool bootstrap, long outSize, long bitSize)
  {
    long L;
    if (bootstrap)
      L = 900; // that should be enough
    else {
      double nBits =
          (outSize > 0 && outSize < 2 * bitSize) ? outSize : (2 * bitSize);
      double three4twoLvls = log(nBits / 2) / log(1.5);
      double add2NumsLvls = log(nBits) / log(2.0);
      L = (5 + ceil(three4twoLvls + add2NumsLvls)) * 30;
    }
    return L;
  };

  // Returns a reference to the passed-in context once it has been modified to
  // get it ready for test. This is not static because it uses quite a lot of
  // state of the object.
  helib::Context& prepareContext(helib::Context& context)
  {
    if (helib_test::verbose) {
      std::cout << "input bitSizes=" << bitSize << ',' << bitSize2
                << ", output size bound=" << outSize << std::endl;
      if (nthreads > 1)
        std::cout << "  using " << NTL::AvailableThreads() << " threads\n";
      std::cout << "computing key-independent tables..." << std::flush;
    }
    buildModChain(context, L, c, /*willBeBootstrappable=*/bootstrap);
    if (bootstrap) {
      context.makeBootstrappable(mvec, /*t=*/0);
    }
    buildUnpackSlotEncoding(unpackSlotEncoding, *context.ea);
    if (helib_test::verbose) {
      std::cout << " done.\n";
      context.zMStar.printout();
    }
    return context;
  };

  void prepareSecKey(helib::SecKey& secKey)
  {
    if (helib_test::verbose) {
      std::cout << " L=" << L << ", B=" << B << std::endl;
      std::cout << "\ncomputing key-dependent tables..." << std::flush;
    }
    secKey.GenSecKey();
    addSome1DMatrices(secKey); // compute key-switching matrices
    addFrbMatrices(secKey);
    if (bootstrap)
      secKey.genRecryptData();
    if (helib_test::verbose)
      std::cout << " done\n";
  };

  const long prm;
  const long bitSize;
  const long bitSize2;
  const long outSize;
  const bool bootstrap;
  const long seed;
  const long nthreads;

  const long* vals;
  const long p;
  const long m;
  const NTL::Vec<long> mvec;
  const std::vector<long> gens;
  const std::vector<long> ords;
  const long B;
  const long c;
  const long L;
  helib::Context context;
  helib::SecKey secKey;

  GTestBinaryArith() :
      prm(validatePrm(std::get<0>(GetParam()).prm)),
      bitSize(correctBitSize(5, std::get<0>(GetParam()).bitSize)),
      bitSize2(correctBitSize(bitSize, std::get<0>(GetParam()).bitSize2)),
      outSize(std::get<0>(GetParam()).outSize),
      bootstrap(std::get<0>(GetParam()).bootstrap),
      seed(std::get<0>(GetParam()).seed),
      nthreads(std::get<0>(GetParam()).nthreads),
      vals(mValues[prm]),
      p(vals[0]),
      m(vals[2]),
      mvec(calculateMvec(vals)),
      gens(calculateGens(vals)),
      ords(calculateOrds(vals)),
      B(vals[13]),
      c(vals[14]),
      L(calculateLevels(bootstrap, outSize, bitSize)),
      context(m, p, /*r=*/1, gens, ords),
      secKey(prepareContext(context)){};

  void SetUp() override
  {
    if (seed)
      NTL::SetSeed(NTL::ZZ(seed));
    if (nthreads > 1)
      NTL::SetNumThreads(nthreads);

    prepareSecKey(secKey);

    helib::activeContext = &context; // make things a little easier sometimes
#ifdef DEBUG_PRINTOUT
    helib::dbgEa = (helib::EncryptedArray*)context.ea;
    helib::dbgKey = &secKey;
#endif
  }

  virtual void TearDown() override
  {
#ifdef DEBUG_PRINTOUT
    helib::cleanupGlobals();
#endif
  }

  // This gets called once all of the tests have run
public:
  static void TearDownTestCase()
  {
    if (helib_test::verbose)
      helib::printAllTimers(std::cout);
  };
};
constexpr long GTestBinaryArith::mValues[8][15];
std::vector<helib::zzX> GTestBinaryArith::unpackSlotEncoding;

TEST_P(GTestBinaryArith, fifteenForFour)
{
  // Randomly generate up to 15 integers from {0,1} with some entries
  // randomly set to null.  We then encrypt and use the fifteenOrLess4Four
  // function to calculate the binary representation of their sum.  This is
  // then checked against the plaintext calculation.
  // In this case each ciphertext is considered to be the encryption of one
  // bit, however packing more into the slots is possible.
  // LSB is at the start (left) of the vector.

  // vector of ciphertexts corresponding to encrypted input bit vectors.
  std::vector<helib::Ctxt> inBuf(15, helib::Ctxt(secKey));
  std::vector<helib::Ctxt*> inPtrs(15, nullptr);

  // vector of ciphertexts corresponding to the summation of the input bit
  // vectors.
  std::vector<helib::Ctxt> outBuf(5, helib::Ctxt(secKey));

  // Randomly generate and encrypt the input vectors.
  long sum = 0;
  std::string inputBits = "(";
  for (int i = 0; i < 15; i++) {
    if (NTL::RandomBnd(10) > 0) { // Leave empty (null) with small probability.
      inPtrs[i] = &(inBuf[i]);
      long bit = NTL::RandomBnd(2); // Select a randomised bit.
      secKey.Encrypt(inBuf[i], NTL::ZZX(bit));
      inputBits += std::to_string(bit) + ",";
      sum += bit; // Keep track of the plaintext sum.
    } else
      inputBits += "-,"; // This represents a null bit.
  }
  inputBits += ")";

  if (helib_test::verbose) {
    std::cout << std::endl;
    helib::CheckCtxt(inBuf[helib::lsize(inBuf) - 1], "b4 15for4");
  }
  // Add the encrypted bits.
  long numOutputs = fifteenOrLess4Four(helib::CtPtrs_vectorCt(outBuf),
                                       helib::CtPtrs_vectorPt(inPtrs));
  if (helib_test::verbose)
    helib::CheckCtxt(outBuf[helib::lsize(outBuf) - 1], "after 15for4");

  // Check the result.
  long sum2 = 0;
  for (int i = 0; i < numOutputs; i++) {
    NTL::ZZX poly;
    secKey.Decrypt(poly, outBuf[i]);
    sum2 += to_long(ConstTerm(poly)) << i;
  }
  EXPECT_EQ(sum, sum2) << "inputs = " << inputBits << std::endl;
  if (helib_test::verbose) {
    std::cout << "15to4 succeeded, sum" << inputBits << "=" << sum2
              << std::endl;
  }
};

TEST_P(GTestBinaryArith, product)
{
  // Randomly generate a pair of numbers of a specified bit size and then
  // encrypt them in binary representation. Then use multTwoNumbers to
  // multiply the two positive binary numbers and then check against the
  // plaintext calculation. Next, use multTwoNumbers with binary numbers in
  // 2's complement to calculate the product where the multplier is negative
  // and then check against the plaintext calculation.
  //
  // In this case each ciphertext is considered to be the encryption of one
  // bit, however packing more into the slots is possible.
  // LSB is at the start (left) of the vector.

  const helib::EncryptedArray& ea = *context.ea;
  // outSize 1's on the least significant end of mask.
  long mask = (outSize ? ((1L << outSize) - 1) : -1);

  // Choose two random integers with correct bit sizes.
  long multiplicand_data = NTL::RandomBits_long(bitSize);
  long multiplier_data = NTL::RandomBits_long(bitSize2);

  // Encrypt the individual bits.
  NTL::Vec<helib::Ctxt> encrypted_product, encrypted_multiplicand,
      encrypted_multiplier;

  // Resizes the vector of ciphertexts (encrypted_bits) to match the input size.
  helib::resize(encrypted_multiplicand, bitSize, helib::Ctxt(secKey));
  for (long i = 0; i < bitSize; i++) {
    secKey.Encrypt(encrypted_multiplicand[i],
                   NTL::ZZX((multiplicand_data >> i) & 1));
    if (bootstrap) { // If bootstrapping then modulo down to a lower level.
      encrypted_multiplicand[i].bringToSet(context.getCtxtPrimes(5));
    }
  }
  helib::resize(encrypted_multiplier, bitSize2, helib::Ctxt(secKey));
  for (long i = 0; i < bitSize2; i++) {
    secKey.Encrypt(encrypted_multiplier[i],
                   NTL::ZZX((multiplier_data >> i) & 1));
    if (bootstrap) { // If bootstrapping then modulo down to a lower level.
      encrypted_multiplier[i].bringToSet(context.getCtxtPrimes(5));
    }
  }
  if (helib_test::verbose) {
    std::cout << "\n  bits-size " << bitSize << '+' << bitSize2;
    if (outSize > 0)
      std::cout << "->" << outSize;
    helib::CheckCtxt(encrypted_multiplier[0], "b4 multiplication");
  }
  std::vector<long> slots; // Vector that will hold the decrypted result.
  // Test multiplication with two positive numbers.
  // A scope which tests multTwoNumbers using wrappers around the encrypted
  // data.
  {
    helib::CtPtrs_VecCt output_wrapper(
        encrypted_product); // A wrapper around the output vector.
    helib::multTwoNumbers(output_wrapper,
                          helib::CtPtrs_VecCt(encrypted_multiplicand),
                          helib::CtPtrs_VecCt(encrypted_multiplier),
                          /*negative=*/false,
                          outSize,
                          &unpackSlotEncoding);
    helib::decryptBinaryNums(slots, output_wrapper, secKey, ea);
  } // output_wrapper is deleted once out of scope.
  if (helib_test::verbose)
    helib::CheckCtxt(encrypted_product[helib::lsize(encrypted_product) - 1],
                     "after multiplication");

  // Calculate the multiplication in the plain.
  long plaintext_product = multiplicand_data * multiplier_data;
  EXPECT_EQ(slots[0], ((multiplicand_data * multiplier_data) & mask))
      << "Positive product error: multiplicand_data=" << multiplicand_data
      << ", multiplier_data=" << multiplier_data << ", but product=" << slots[0]
      << " (should be " << plaintext_product << '&' << mask << '='
      << (plaintext_product & mask) << ")\n";

  if (helib_test::verbose) {
    std::cout << "positive product succeeded: ";
    if (outSize)
      std::cout << "bottom " << outSize << " bits of ";
    std::cout << multiplicand_data << "*" << multiplier_data << "=" << slots[0]
              << std::endl;
  }

  // Test multiplication of numbers in 2's complement where the multiplier is
  // negative.
  secKey.Encrypt(encrypted_multiplier[bitSize2 - 1], NTL::ZZX(1));
  decryptBinaryNums(slots,
                    helib::CtPtrs_VecCt(encrypted_multiplier),
                    secKey,
                    ea,
                    /*negative=*/true);
  multiplier_data = slots[0];
  encrypted_product.kill(); // Clear the data in encrypted_product.
  // A scope which tests multTwoNumbers (with a negative multiplier) using
  // wrappers around the encrypted data.
  {
    helib::CtPtrs_VecCt output_wrapper(
        encrypted_product); // A wrapper around the output vector.
    multTwoNumbers(output_wrapper,
                   helib::CtPtrs_VecCt(encrypted_multiplicand),
                   helib::CtPtrs_VecCt(encrypted_multiplier),
                   /*negative=*/true,
                   outSize,
                   &unpackSlotEncoding);
    decryptBinaryNums(slots, output_wrapper, secKey, ea, /*negative=*/true);
  } // output_wrapper is deleted once out of scope.

  if (helib_test::verbose)
    helib::CheckCtxt(encrypted_product[helib::lsize(encrypted_product) - 1],
                     "after multiplication");

  // Calculate the multiplication in the plain.
  plaintext_product = multiplicand_data * multiplier_data;
  EXPECT_EQ((slots[0] & mask), (plaintext_product & mask))
      << "Negative product error: multiplicand_data=" << multiplicand_data
      << ", multiplier_data=" << multiplier_data << ", but product=" << slots[0]
      << " (should be " << plaintext_product << '&' << mask << '='
      << (plaintext_product & mask) << ")\n";
  if (helib_test::verbose) {
    std::cout << "negative product succeeded: ";
    if (outSize)
      std::cout << "bottom " << outSize << " bits of ";
    std::cout << multiplicand_data << "*" << multiplier_data << "=" << slots[0]
              << std::endl;
  }

#ifdef DEBUG_PRINTOUT
  // Print out the ciphertext with the lowest level after multiplication
  // if DEBUG_PRINTOUT is defined.
  const helib::Ctxt* minCtxt = nullptr;
  long minLvl = 10000000;
  for (const helib::Ctxt& c : encrypted_product) {
    long lvl = c.logOfPrimeSet();
    if (lvl < minLvl) {
      minCtxt = &c;
      minLvl = lvl;
    }
  }
  decryptAndPrint(
      (std::cout << " after multiplication: "), *minCtxt, secKey, ea, 0);
  std::cout << std::endl;
#endif
};

TEST_P(GTestBinaryArith, add)
{
  // Randomly generate a pair of numbers of a specified bit size and then
  // encrypt them in binary representation. Then use addTwoNumbers to add the
  // two binary numbers and then check against the plaintext calculation.
  //
  // In this case each ciphertext is considered to be the encryption of one
  // bit, however packing more into the slots is possible.
  // LSB is at the start(left) of the vector.

  const helib::EncryptedArray& ea = *context.ea;
  // outSize 1's on the least significant end of mask.
  long mask = (outSize ? ((1L << outSize) - 1) : -1);

  // Choose two random n-bit integers.
  long addend_data = NTL::RandomBits_long(bitSize);
  long augend_data = NTL::RandomBits_long(bitSize2);

  // Encrypt the individual bits.
  NTL::Vec<helib::Ctxt> encrypted_sum, encrypted_addend, encrypted_augend;

  // Resizes the vector of ciphertexts (encrypted_bits) to match the input size.
  helib::resize(encrypted_addend, bitSize, helib::Ctxt(secKey));
  for (long i = 0; i < bitSize; i++) {
    secKey.Encrypt(encrypted_addend[i], NTL::ZZX((addend_data >> i) & 1));
    if (bootstrap) { // If bootstrapping then modulo down to a lower level.
      encrypted_addend[i].bringToSet(context.getCtxtPrimes(5));
    }
  }
  // Resizes the vector of ciphertexts (encrypted_bits) to match the input size.
  helib::resize(encrypted_augend, bitSize2, helib::Ctxt(secKey));
  for (long i = 0; i < bitSize2; i++) {
    secKey.Encrypt(encrypted_augend[i], NTL::ZZX((augend_data >> i) & 1));
    if (bootstrap) { // If bootstrapping then modulo down to a lower level.
      encrypted_augend[i].bringToSet(context.getCtxtPrimes(5));
    }
  }
  if (helib_test::verbose) {
    std::cout << "\n  bits-size " << bitSize << '+' << bitSize2;
    if (outSize > 0)
      std::cout << "->" << outSize;
    std::cout << std::endl;
    helib::CheckCtxt(encrypted_augend[0], "b4 addition");
  }

  std::vector<long>
      decrypted_result; // Vector that will hold the decrypted result.
  // Test addition.
  // A scope which tests addTwoNumbers using wrappers around the encrypted data.
  {
    helib::CtPtrs_VecCt output_wrapper(
        encrypted_sum); // A wrapper around the output vector.
    helib::addTwoNumbers(output_wrapper,
                         helib::CtPtrs_VecCt(encrypted_addend),
                         helib::CtPtrs_VecCt(encrypted_augend),
                         outSize,
                         &unpackSlotEncoding);
    helib::decryptBinaryNums(decrypted_result, output_wrapper, secKey, ea);
  } // output_wrapper is deleted once out of scope.

  if (helib_test::verbose)
    helib::CheckCtxt(encrypted_sum[helib::lsize(encrypted_sum) - 1],
                     "after addition");

  // Calculate the addition in the plain.
  long plaintext_sum = addend_data + augend_data;
  EXPECT_EQ(decrypted_result[0], ((addend_data + augend_data) & mask))
      << "addTwoNums error: addend_data=" << addend_data
      << ", augend_data=" << augend_data
      << ", but plaintext_sum=" << decrypted_result[0]
      << " (should be =" << (plaintext_sum & mask) << ")\n";
  if (helib_test::verbose) {
    std::cout << "addTwoNums succeeded: ";
    if (outSize)
      std::cout << "bottom " << outSize << " bits of ";
    std::cout << addend_data << "+" << augend_data << "=" << decrypted_result[0]
              << std::endl;
  }

#ifdef DEBUG_PRINTOUT
  // Print out the ciphertext with the lowest level after addition if
  // DEBUG_PRINTOUT is defined.
  const helib::Ctxt* minCtxt = nullptr;
  long minLvl = 1000;
  for (const helib::Ctxt& c : encrypted_sum) {
    long lvl = c.logOfPrimeSet();
    if (lvl < minLvl) {
      minCtxt = &c;
      minLvl = lvl;
    }
  }
  decryptAndPrint((std::cout << " after addition: "), *minCtxt, secKey, ea, 0);
  std::cout << std::endl;
#endif
};

TEST_P(GTestBinaryArith, addManyNumbers)
{
  // Randomly generate a vector of numbers of a specified bit size and then
  // encrypt them in binary representation. Then use addManyNumbers to add
  // them all up and then check against the plaintext calculation.
  //
  // In this case each ciphertext is considered to be the encryption of one
  // bit, however packing more into the slots is possible.
  // LSB is at the start(left) of the vector.

  const long num_summands = 5;
  const helib::EncryptedArray& ea = *context.ea;
  // outSize 1's on the least significant end of mask.
  long mask = (outSize ? ((1L << outSize) - 1) : -1);

  // Choose a set of random n-bit integers.
  std::vector<long> summands_data;
  for (long i = 0; i < num_summands; ++i)
    summands_data.push_back(NTL::RandomBits_long(bitSize));

  // Encrypt the individual bits.
  std::vector<helib::Ctxt> encrypted_sum;
  std::vector<std::vector<helib::Ctxt>> encrypted_summands;

  // Utility function for encrypting a number into a binary representation.
  const auto encrypt_binary_number =
      [&](const long num) -> std::vector<helib::Ctxt> {
    std::vector<helib::Ctxt> encrypted_num;
    // Resizes the vector of ciphertexts (encrypted_bits) to match the input
    // size.
    helib::resize(encrypted_num, bitSize, helib::Ctxt(secKey));
    for (long i = 0; i < bitSize; i++) {
      secKey.Encrypt(encrypted_num[i], NTL::ZZX((num >> i) & 1));
      if (bootstrap) { // If bootstrapping then modulo down to a lower level.
        encrypted_num[i].bringToSet(context.getCtxtPrimes(5));
      }
    }
    return encrypted_num;
  };

  // Encrypt the set of numbers into binary representation.
  for (long i = 0; i < num_summands; ++i)
    encrypted_summands.push_back(encrypt_binary_number(summands_data[i]));

  if (helib_test::verbose)
    helib::CheckCtxt(encrypted_summands[0][0], "b4 addition");

  std::vector<long>
      decrypted_result; // Vector that will hold the decrypted result.
  // Test summation.
  // A scope which tests addManyNumbers using wrappers around the encrypted
  // data.
  {
    helib::CtPtrs_vectorCt output_wrapper(
        encrypted_sum); // A wrapper around the output vector.
    helib::CtPtrMat_vectorCt summands_wrapper(
        encrypted_summands); // A wrapper around the output vector.
    helib::addManyNumbers(
        output_wrapper, summands_wrapper, outSize, &unpackSlotEncoding);
    helib::decryptBinaryNums(decrypted_result, output_wrapper, secKey, ea);
  } // output_wrapper is deleted once out of scope.

  if (helib_test::verbose)
    helib::CheckCtxt(encrypted_sum[helib::lsize(encrypted_sum) - 1],
                     "after addition");

  // Calculate the summation in the plain.
  long plaintext_sum = std::accumulate(
      summands_data.begin(), summands_data.end(), 0l, std::plus<long>());
  EXPECT_EQ(decrypted_result[0], plaintext_sum & mask);
  if (helib_test::verbose) {
    std::cout << "addManyNums succeeded: ";
    if (outSize)
      std::cout << "bottom " << outSize << " bits of ";
    std::cout << summands_data[0];
    for (long i = 1; i < num_summands; ++i)
      std::cout << "+" << summands_data[i];
    std::cout << "=" << decrypted_result[0] << std::endl;
  }
}

INSTANTIATE_TEST_SUITE_P(
    smallParameterSizesRepeated,
    GTestBinaryArith,
    ::testing::Combine(
        ::testing::Values(
            // SLOW
            Parameters(1, 5, 0, 0, false, 0, 1)
            // FAST
            // Parameters(0, 5, 0, 0, false, 0, 1)), ::testing::Range(0,2)

            ),
        ::testing::Range(0, 2)) // The range is for repeats
);

} // anonymous namespace
