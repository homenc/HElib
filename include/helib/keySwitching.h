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

#ifndef HELIB_KEY_SWITCHING_H
#define HELIB_KEY_SWITCHING_H
/**
 * @file keySwitching.h
 * @brief - Declaration of key switching matrices
 * @brief - Other key switching-related free functions
 *
 * Copyright IBM Corporation 2019 All rights reserved.
 */

#include <climits>
#include <helib/DoubleCRT.h>
#include <helib/Context.h>
#include <helib/Ctxt.h>

namespace helib {

/**
 * @class KeySwitch
 * @brief Key-switching matrices
 *
 * There are basically two approaches for how to do key-switching: either
 * decompose the mod-q ciphertext into bits (or digits) to make it low-norm,
 * or perform the key-switching operation mod Q>>q. The tradeoff is that when
 * decomposing the (coefficients of the) ciphertext into t digits, we need to
 * increase the size of the key-switching matrix by a factor of t (and the
 * running time similarly grows). On the other hand if we do not decompose at
 * all then we need to work modulo Q>q^2, which means that the bitsize of our
 * largest modulus q0 more than doubles (and hence also the parameter m more
 * than doubles). In general if we decompose into digits of size B then we
 * need to work with Q>q*B.)
 *
 * The part of the spectrum where we expect to find the sweet spot is when
 * we decompose the ciphertext into digits of size B=q0^{1/t} for some small
 * constant t (maybe t=2,3 or so). This means that our largest modulus
 * has to be Q>q0^{1+1/t}, which increases also the parameter m by a factor
 * (1+1/t). It also means that for key-switching in the top levels we would
 * break the ciphertext to t digits, hence the key-switching matrix will have
 * t columns.
 *
 * A key-switch matrix W[s'->s] converts a ciphertext-part with respect to
 * secret-key polynomial s' into a canonical ciphertext (i.e. a two-part
 * ciphertext with respect to (1,s)). The matrix W is a 2-by-t matrix of
 * DoubleCRT objects. The bottom row are just (pseudo)random elements. Then
 * for column j, if the bottom element is aj then the top element is set as
 *     bj = P*Bj*s' + p*ej - s * aj mod P*q0,
 * where p is the plaintext space (i.e. 2 or 2^r, or 1 for CKKS) and Bj
 * is the product of the digits-sizes corresponding to columns 0...i-1.
 * (For example if we have digit sizes 3,5,7 then B0=1, B1=3, B2=15 and
 * B3=105.) Also, q0 is the product of all the "ciphertext primes" and
 * P is roughly the product of all the special primes. (Actually, for BGV,
 * if Q is the product of all the special primes then P=Q*(Q^{-1} mod p).)
 *
 * In this implementation we save some space, by keeping only a PRG seed for
 * generating the pseudo-random elements, rather than the elements themselves.
 *
 * To convert a ciphertext part R, we break R into digits R = sum_j Bj Rj,
 * then set (q0,q1)^T = sum_j Rj * column-j. Note that we have
 * <(1,s),(q0,q1)> = sum_j Rj*(s*aj - s*aj + p*ej +P*Bj*s')
 *       = P * sum_j Bj*Rj * s' + p sum_j Rj*ej
 *       = P * R * s' + p*a-small-element (mod P*q0)
 * where the last element is small since the ej's are small and |Rj|<B.
 * Note that if the ciphertext is encrypted relative to plaintext space p'
 * and then key-switched with matrices W relative to plaintext space p,
 * then we get a mew ciphertext with noise p'*small+p*small, so it is valid
 * relative to plaintext space GCD(p',p).
 *
 * The matrix W is defined modulo Q>t*B*sigma*q0 (with sigma a bound on the
 * size of the ej's), and Q is the product of all the small primes in our
 * moduli chain. However, if p is much smaller than B then is is enough to
 * use W mod Qi with Qi a smaller modulus, Q>p*sigma*q0. Also note that if
 * p<Br then we will be using only first r columns of the matrix W.
 ********************************************************************/
class KeySwitch
{
public:
  /**
   * @brief Class label to be added to JSON serialization as object type
   * information.
   */
  static constexpr std::string_view typeName = "KeySwitch";

  SKHandle fromKey; // A handle for the key s'
  long toKeyID;     // Index of the key s that we are switching into
  long ptxtSpace;   // either p or p^r

  std::vector<DoubleCRT> b; // The top row, consisting of the bi's
  NTL::ZZ prgSeed; // a seed to generate the random ai's in the bottom row
  NTL::xdouble noiseBound; // high probability bound on noise magnitude
  // in each column

  explicit KeySwitch(long sPow = 0,
                     long xPow = 0,
                     long fromID = 0,
                     long toID = 0,
                     long p = 0);
  explicit KeySwitch(const SKHandle& _fromKey,
                     long fromID = 0,
                     long toID = 0,
                     long p = 0);

  bool operator==(const KeySwitch& other) const;
  bool operator!=(const KeySwitch& other) const;

  unsigned long NumCols() const;

  //! @brief returns a dummy static matrix with toKeyId == -1
  static const KeySwitch& dummy();
  bool isDummy() const;

  //! A debugging method
  void verify(SecKey& sk);

  //! @brief Read a key-switching matrix from input
  void readMatrix(std::istream& str, const Context& context);

  //! Raw IO
  /**
   * @brief Write out the `KeySwitch` object in binary format.
   * @param str Output `std::ostream`.
   **/
  void writeTo(std::ostream& str) const;

  /**
   * @brief Read from the stream the serialized `KeySwitch` object in binary
   * format.
   * @param str Input `std::istream`.
   * @return The deserialized `KeySwitch` object.
   **/
  static KeySwitch readFrom(std::istream& str, const Context& context);

  /**
   * @brief Write out the switch key (`KeySwitch`) object to the output
   * stream using JSON format.
   * @param str Output `std::ostream`.
   **/
  void writeToJSON(std::ostream& str) const;

  /**
   * @brief Write out the switch key (`KeySwitch`) object to a `JsonWrapper`.
   * @return The `JsonWrapper`.
   **/
  JsonWrapper writeToJSON() const;

  /**
   * @brief Read from the stream the serialized switch key (`KeySwitch`) object
   * using JSON format.
   * @param str Input `std::istream`.
   * @param context The `Context` to be used.
   * @return The deserialized `KeySwitch` object.
   **/
  static KeySwitch readFromJSON(std::istream& str, const Context& context);

  /**
   * @brief Read from the `JsonWrapper` the serialized switch key (`KeySwitch`)
   * object.
   * @param j The `JsonWrapper` containing the serialized `KeySwitch` object.
   * @param context The `Context` to be used.
   * @return The deserialized `KeySwitch` object.
   **/
  static KeySwitch readFromJSON(const JsonWrapper& j, const Context& context);

  /**
   * @brief In-place read from the stream the serialized switch key
   * (`KeySwitch`) object using JSON format.
   * @param str Input `std::istream`.
   * @param context The `Context` to be used.
   **/
  void readJSON(std::istream& str, const Context& context);

  /**
   * @brief In-place read from the `JsonWrapper` the serialized switch key
   * (`KeySwitch`) object.
   * @param j The `JsonWrapper` containing the serialized `KeySwitch` object.
   * @param context The `Context` to be used.
   **/
  void readJSON(const JsonWrapper& j, const Context& context);
};
std::ostream& operator<<(std::ostream& str, const KeySwitch& matrix);
// We DO NOT have std::istream& operator>>(std::istream& str, KeySwitch&
// matrix); instead must use the readMatrix method above, where you can specify
// context

//! @name Strategies for generating key-switching matrices
//! These functions are implemented in KeySwitching.cpp

//! @brief Constant defining threshold above which a baby-set/giant-step
//! strategy is used
#define HELIB_KEYSWITCH_THRESH (50)

//! @brief Constant defining threshold above which a single
//! giant step matrix is added even in HELIB_KSS_MIN mode.
//! This helps in the matmul routines.
#define HELIB_KEYSWITCH_MIN_THRESH (8)

//! @brief Function that returns number of baby steps.  Used to keep
//! this and matmul routines "in sync".
long KSGiantStepSize(long D);

//! @brief Maximalistic approach:
//! generate matrices s(X^e)->s(X) for all e in Zm*
void addAllMatrices(SecKey& sKey, long keyID = 0);

//! @brief Generate matrices so every s(X^e) can be reLinearized
//! in at most two steps
void addFewMatrices(SecKey& sKey, long keyID = 0);

//! @brief Generate some matrices of the form s(X^{g^i})->s(X), but not all.
//! For a generator g whose order is larger than bound, generate only enough
//! matrices for the giant-step/baby-step procedures (2*sqrt(ord(g))of them).
void addSome1DMatrices(SecKey& sKey,
                       long bound = HELIB_KEYSWITCH_THRESH,
                       long keyID = 0);

//! @brief Generate all matrices s(X^{g^i})->s(X) for generators g of
//! Zm* /(p) and i<ord(g). If g has different orders in Zm* and Zm* /(p)
//! then generate also matrices of the form s(X^{g^{-i}})->s(X)
void add1DMatrices(SecKey& sKey, long keyID = 0);

void addBSGS1DMatrices(SecKey& sKey, long keyID = 0);

//! Generate all/some Frobenius matrices of the form s(X^{p^i})->s(X)
void addSomeFrbMatrices(SecKey& sKey,
                        long bound = HELIB_KEYSWITCH_THRESH,
                        long keyID = 0);

void addFrbMatrices(SecKey& sKey, long keyID = 0);

void addBSGSFrbMatrices(SecKey& sKey, long keyID = 0);

//! These routines just add a single matrix (or two, for bad dimensions)
void addMinimal1DMatrices(SecKey& sKey, long keyID = 0);
void addMinimalFrbMatrices(SecKey& sKey, long keyID = 0);

//! Generate all key-switching matrices for a given permutation network
class PermNetwork;
void addMatrices4Network(SecKey& sKey, const PermNetwork& net, long keyID = 0);

//! Generate specific key-switching matrices, described by the given set
void addTheseMatrices(SecKey& sKey,
                      const std::set<long>& automVals,
                      long keyID = 0);

} // namespace helib

#endif // HELIB_KEY_SWITCHING_H
