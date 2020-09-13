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

#ifndef HELIB_PTXT_H
#define HELIB_PTXT_H

#include <type_traits>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iomanip>

#include <helib/Context.h>
#include <helib/EncryptedArray.h>
#include <helib/assertions.h>
#include <helib/PolyMod.h>

/**
 * @file Ptxt.h
 * @brief Plaintext object parameterised on CKKS and BGV schemes. Also contains
 * definition of `CKKS` and `BGV` structs.
 **/

namespace helib {

/**
 * @brief Type for CKKS scheme, to be used as template parameter.
 **/
struct CKKS
{
  /**
   * @brief Slot type used for CKKS plaintexts: `std::complex<double>`.
   **/
  using SlotType = std::complex<double>;
};

/**
 * @brief Type for BGV scheme, to be used as template parameter.
 **/
struct BGV
{
  /**
   * @brief Slot type used for BGV plaintexts: `helib::PolyMod` i.e. an integer
   * polynomial modulo p^r and G.
   **/
  using SlotType = PolyMod;
};

// Utility functions

/**
 * @brief Converts `std::vector<From>` to `std::vector<Scheme::SlotType>`.
 * @tparam From Type of the element in the input vector.
 * @tparam Scheme The encryption scheme to be used, must be `BGV` or `CKKS`.
 * @param data Vector to be converted.
 * @return Vector of converted values of type `Scheme::SlotType`.
 * @note Only exists for `BGV` and `CKKS`.
 **/
template <typename From, typename Scheme>
inline std::vector<typename Scheme::SlotType> convertDataToSlotVector(
    const std::vector<From>& data,
    const Context& context)
{
  static_assert(std::is_same<Scheme, CKKS>::value ||
                    std::is_same<Scheme, BGV>::value,
                "Can only call convertDataToSlotVector with Scheme equals to "
                "CKKS or BGV");
  using To = typename Scheme::SlotType;
  std::vector<To> res(data.size(), Ptxt<Scheme>::convertToSlot(context, 0l));
  for (std::size_t i = 0; i < data.size(); ++i) {
    res[i] = data[i];
  }
  return res;
}

/**
 * @brief Free function that computes the inner product of two vectors of
 * `Ptxt`.
 * @param result The output `Ptxt` that will hold the result.
 * @param first_vec The first input vector of plaintexts.
 * @param second_vec The second input vector of plaintexts.
 * @note If the two vector sizes differ, the shorter vector will be padded with
 * zeroes.
 **/
template <typename Scheme>
void innerProduct(Ptxt<Scheme>& result,
                  const std::vector<Ptxt<Scheme>>& first_vec,
                  const std::vector<Ptxt<Scheme>>& second_vec)
{
  static_assert(std::is_same<Scheme, CKKS>::value ||
                    std::is_same<Scheme, BGV>::value,
                "Can only call innerProduct with Scheme equals to CKKS or BGV");
  for (std::size_t i = 0; i < first_vec.size(); ++i) {
    assertTrue<RuntimeError>(first_vec[i].isValid(),
                             "Cannot call innerProduct on default-constructed"
                             " Ptxt as first argument at index " +
                                 std::to_string(i));
  }
  for (std::size_t i = 0; i < second_vec.size(); ++i) {
    assertTrue<RuntimeError>(second_vec[i].isValid(),
                             "Cannot call innerProduct on default-constructed"
                             " Ptxt as second argument at index " +
                                 std::to_string(i));
  }
  long n = std::min(first_vec.size(), second_vec.size());
  if (n <= 0) {
    result.clear();
    return;
  }
  result = first_vec[0];
  result *= second_vec[0];
  for (long i = 1; i < n; ++i)
    result += (first_vec[i] * second_vec[i]);
}

/**
 * @brief Deserialize a `std::complex<double>` from the input stream `is`
 * delimited by '[' and ']' (instead of the default '(', ')').
 * @param is The input stream reference.
 * @param num The complex number to deserialize.
 * @throws IOError if the stream contains more than 2 parts.
 *
 * The input stream has to be formatted as '['number']' (this will be
 * deserialized as (number, 0)) or '['real_part', 'imag_part']' (this will be
 * deserialized as (real_part, imag_part)).
 */
void deserialize(std::istream& is, std::complex<double>& num);

/**
 * @brief Serialize a `std::complex<double>` to the output stream `os`
 * delimited by '[' and ']' (instead of the default '(', ')').
 * @param os The output stream reference.
 * @param num The complex number to serialize.
 *
 * The output will be formatted as '['`num.real()`', '`num.imag()`']'.
 */
void serialize(std::ostream& os, const std::complex<double>& num);

// Forward declaration as function is a friend of the templated `Ptxt` class.
/**
 * @brief Function to deserialize a `Ptxt<Scheme>`.
 * @tparam Scheme The `Ptxt` object scheme.  Can be only be `BGV` or `CKKS`.
 * @param is Input `std::istream`.
 * @param ptxt Destination `Ptxt` object.
 * @throws IOError if the stream is badly formatted (i.e. it is not delimited by
 * '[' and ']').
 * @note `ptxt` must be constructed with an appropriate context @b BEFORE
 * calling this function. For example,
 * @code
 * Ptxt my_ptxt(context);
 * deserialize(std::cin, my_ptxt);
 * @endcode
 *
 * The input stream has to be formatted as a comma-separated list surrounded by
 * '[' and ']'.\n
 * Each element of the list will be deserialized as a slot of the
 * type determined by the scheme.\n
 * If the number of tokens in the list is less
 * than the number of slots, the remaining slots will be padded by 0.\n
 * For example '['slot0', 'slot1', 'slot2']' will be deserialized as a
 * plaintext `ptxt` where `ptxt[0]=slot0`, `ptxt[1]=slot1`,
 * `ptxt[2]=slot2` and `ptxt[i]=0` for `i>2`.
 */
template <typename Scheme>
void deserialize(std::istream& is, Ptxt<Scheme>& ptxt);

// Forward declaration as function is a friend of the templated `Ptxt` class.
/**
 * @brief Function to serialize a `Ptxt<Scheme>`.
 * @tparam Scheme The `Ptxt` object scheme.  Can be only be `BGV` or `CKKS`.
 * @param os Output `std::ostream`.
 * @param ptxt `Ptxt` object to be written.
 * @return Input `std::ostream` post serializing.
 * @note `Ptxt` `context` is not serialized, see note of `deserialize`.
 *
 * The output stream will be formatted as a comma-separated list surrounded by
 * '[' and ']'.\n
 * Each slot of `ptxt` will be serialized in an element of such list by the
 * `serialize` function determined by the scheme.\n
 * For example if we have a plaintext `ptxt` such that `ptxt[0]=slot0`,
 * `ptxt[1]=slot1`, `ptxt[2]=slot2`, and `ptxt[i]=0` for `i>2`, it will be
 * serialized as '['slot0', 'slot1', 'slot2', `0`, `0` ...]'.
 */
template <typename Scheme>
void serialize(std::ostream& os, const Ptxt<Scheme>& ptxt);

// Forward declaration as function is a friend of the templated `Ptxt` class.
/**
 * @brief Input shift operator.  Uses the `deserialize` function internally.
 * @param is Input `std::istream`.
 * @param ptxt Destination `Ptxt` object.
 * @return Input `std::istream` post reading.
 * @note `ptxt` must be constructed with an appropriate context @b BEFORE
 * calling this function. For example,
 * @code
 * Ptxt my_ptxt(context);
 * std::cin >> my_ptxt;
 * @endcode
 *
 * The input stream has to be formatted as a comma-separated list surrounded by
 * '[' and ']'.\n
 * Each element of the list will be deserialized as a slot of the
 * type determined by the scheme.\n
 * If the number of tokens in the list is less
 * than the number of slots, the remaining slots will be padded by 0.\n
 * For example '['slot0', 'slot1', 'slot2']' will be deserialized as a
 * plaintext `ptxt` where `ptxt[0]=slot0`, `ptxt[1]=slot1`,
 * `ptxt[2]=slot2` and `ptxt[i]=0` for `i>2`.
 **/
template <typename Scheme>
std::istream& operator>>(std::istream& is, Ptxt<Scheme>& ptxt);

// Forward declaration as function is a friend of the templated `Ptxt` class.
/**
 * @brief Output shift operator.  Uses the `serialize` function internally.
 * @param os Output `std::ostream`.
 * @param ptxt `Ptxt` object to be written.
 * @return Input `std::ostream` post writing.
 * @note `Ptxt` `context` is not serialized, see note of `operator>>`.
 *
 * The output stream will be formatted as a comma-separated list surrounded by
 * '[' and ']'.\n
 * Each slot of `ptxt` will be serialized in an element of such list by the
 * `serialize` function determined by the scheme.\n
 * For example if we have a plaintext `ptxt` such that `ptxt[0]=slot0`,
 * `ptxt[1]=slot1`, `ptxt[2]=slot2`, and `ptxt[i]=0` for `i>2`, it will be
 * serialized as '['slot0', 'slot1', 'slot2', `0`, `0` ...]'.
 **/
template <typename Scheme>
std::ostream& operator<<(std::ostream& is, const Ptxt<Scheme>& ptxt);

/**
 * @class Ptxt
 * @brief An object that mimics the functionality of the `Ctxt` object, and
 * acts as a convenient entry point for inputting/encoding data which is to be
 * encrypted.
 *
 * `Ptxt` is templated on `Scheme`, which may be `CKKS` or `BGV`.
 *
 * In the BGV case, `Ptxt` can be considered to be an element of
 * \f$\mathbb{Z}_p[x]/\Phi(m)\f$ viewed as a vector of slots with values each
 * in \f$\mathbb{Z}_p[x]/G\f$ where G is one of the irreducible factors of
 * \f$\Phi(m)\f$, and all operations are performed entry-wise.
 *
 * General usage:
 * @code
 * helib::Ptxt<BGV> p1(bgv_context, data);
 * helib::Ptxt<BGV> p2(bgv_context, data);
 * p1 += p2;
 * std::cout << p1 << std::endl;
 * @endcode
 *
 * Internally, `Ptxt` objects store their data as a
 * `std::vector<helib::PolyMod>`, where `PolyMod` is a convenience type
 * representing an element of the above ring, \f$\mathbb{Z}_p[x]/G\f$.  The
 * `PolyMod` type can be easily converted via `static_cast` to more convenient
 * types such as `long` and `NTL::ZZX`.
 *
 * In the `CKKS` case, the slot type is `std::complex<double>`, and has
 * sensible operator overloads supporting operations with other `Ptxt<CKKS>`,
 * `Ctxt`, and `std::complex<double>` objects, as well as performing all
 * operations slot-wise.
 *
 * A large number of operator overloads are defined so that `Ptxt` objects
 * should easily inter-operate, as well as providing interoperability with
 * other logically compatible types e.g. `long` and `NTL::ZZX` in the BGV case,
 * `std::complex<double>` in the `CKKS` case, and `helib::Ctxt` in both cases.
 **/
template <typename Scheme>
class Ptxt
{
  static_assert(std::is_same<Scheme, CKKS>::value ||
                    std::is_same<Scheme, BGV>::value,
                "Can only create plaintext object parameterized by the crypto "
                "scheme (CKKS or BGV)");

public:
  /**
   * @brief Alias for type to be stored in the slots.
   *
   * `std::complex<double>` for CKKS, `helib::PolyMod` for BGV.
   **/
  using SlotType = typename Scheme::SlotType;

  /**
   * @brief Default constructor results in invalid `Ptxt` object which throws if
   * used.
   **/
  Ptxt();

  /**
   * @brief Context only constructor, defaults all slots to `0`.
   * @param context `Context` to use.
   **/
  explicit Ptxt(const Context& context);

  /**
   * @brief Single slot constructor, set all slots to `value`.
   * @param context `Context` to use.
   * @param value Value to set all slots to.
   **/
  Ptxt(const Context& context, const SlotType& value);

  /**
   * @brief BGV plaintext polynomial constructor, set all slots to the `value`
   * polynomial.
   * @param context `Context` to use.
   * @param data Polynomial to be converted into slot representation.
   * @note Only exists for `BGV`.
   **/
  template <typename U = Scheme,
            std::enable_if_t<std::is_same<U, BGV>::value>* = nullptr>
  Ptxt(const Context& context, const NTL::ZZX& value);

  /**
   * @brief Slot vector constructor.
   * @param context `Context` to use.
   * @param data Data to populate the slots.
   **/
  Ptxt(const Context& context, const std::vector<SlotType>& data);

  /**
   * @brief Generic slot vector constructor.
   * @param context `Context` to use.
   * @param data Data to populate the slots, must be convertible to `SlotType`.
   **/
  template <typename T>
  Ptxt(const Context& context, const std::vector<T>& data) :
      Ptxt<Scheme>(context, convertDataToSlotVector<T, Scheme>(data, context))
  {}

  /**
   *  @brief Default copy constructor.
   *  @param other `Ptxt` object to copy.
   **/
  Ptxt(const Ptxt<Scheme>& other) = default;

  /**
   * @brief Default move constructor.
   * @param other `Ptxt` to copy.
   **/
  Ptxt(Ptxt<Scheme>&& other) noexcept = default;

  /**
   * @brief Copy assignment operator with other `Ptxt`.
   * @param other `Ptxt` to copy.
   **/
  Ptxt<Scheme>& operator=(const Ptxt<Scheme>& v) = default;

  /**
   * @brief Move assignment operator with other `Ptxt`.
   * @param other `Ptxt` to copy.
   **/
  Ptxt<Scheme>& operator=(Ptxt<Scheme>&& v) noexcept = default;

  /**
   * @brief Default destructor.
   **/
  ~Ptxt() = default;

  /**
   * @brief Check if a `Ptxt` is valid.
   * @return `true` if valid, `false` otherwise.
   **/
  bool isValid() const;

  /**
   * @brief Returns the size (number of slots) of a `Ptxt`.
   * @return Number of slots of the `Ptxt`.
   **/
  size_t size() const;

  /**
   * @brief Returns the size (number of slots) of a `Ptxt` as `long`.
   * @return Number of slots of the `Ptxt`.
   **/
  long lsize() const;

  /**
   * @brief Returns the `context` used to initialize the `Ptxt`
   * @return The `context`.
   **/
  const Context& getContext() const { return *context; }

  /**
   * @brief Set the data.
   * @param data Vector of `SlotType` to populate the slots.
   **/
  void setData(const std::vector<SlotType>& data);

  /**
   * @brief Set the data replicating the input on all slots.
   * @param value `value` to set all slots to.
   **/
  void setData(const SlotType& value);

  /**
   * @brief Set the `Ptxt` data replicating the input polynomial on all slots.
   * @param data Polynomial to be replicate into slots.
   * @note Only works in the `BGV` case.
   **/
  template <typename T = Scheme,
            typename std::enable_if_t<std::is_same<T, BGV>::value>* = nullptr>
  void setData(const NTL::ZZX& value);

  /**
   * @brief Set the `Ptxt` slots using values from decoding `data` to slot
   * representation.
   * @param data Polynomial to be decoded and converted into slot data.
   * @note Only works in the `BGV` case.
   **/
  template <typename T = Scheme,
            typename std::enable_if_t<std::is_same<T, BGV>::value>* = nullptr>
  void decodeSetData(const NTL::ZZX& data);

  /**
   * @brief Sets all slots to `0`.
   **/
  void clear();

  /**
   * @brief Populate slots with random data.
   * @return Reference to `*this` post population.
   **/
  Ptxt<Scheme>& random();

  /**
   * @brief Get the data held in the slots as a `std::vector<SlotType>`.
   * @return Constant reference to the slot vector.
   **/
  const std::vector<SlotType>& getSlotRepr() const;

  /**
   * @brief Converts the slot data in `this` to its single polynomial
   * representation.
   * @return Single encoded polynomial.
   * @note `NTL::ZZX` representation loses some precision in the `CKKS` case.
   **/
  NTL::ZZX getPolyRepr() const;

  /**
   * @brief Square bracket accessor operator.
   * @param i Index of the desired `Ptxt` slot.
   * @return Reference to the data held at slot `i`.
   **/
  SlotType& operator[](long i);

  /**
   * @brief `const` square bracket accessor operator.
   * @param i Index of the desired `Ptxt` slot.
   * @return Copy of the data held at slot `i`.
   **/
  SlotType operator[](long i) const;

  /**
   * @brief `at` accessor operator.
   * @param i Index of the desired `Ptxt` slot.
   * @return Reference to the data held at slot `i`.
   * @note throws if `i` is out of range.
   **/
  SlotType& at(long i);

  /**
   * @brief `const` `at` accessor operator.
   * @param i Index of the desired `Ptxt` slot.
   * @return Copy of the data held at slot `i`.
   * @note throws if `i` is out of range.
   **/
  SlotType at(long i) const;

  /**
   * @brief Equals operator between two `Ptxt` objects.
   * @param other `Ptxt` to compare to.
   * @return `true` if identical, `false` otherwise.
   **/
  bool operator==(const Ptxt<Scheme>& other) const;

  /**
   * @brief Not equals operator between two `Ptxt` objects.
   * @param other `Ptxt` to compare to.
   * @return `true` if differ, `false` otherwise.
   **/
  bool operator!=(const Ptxt<Scheme>& other) const;

  /**
   * @brief Infix multiplication operator.
   * @param rhs Right hand side of multiplication.
   * @return Product of the two `Ptxt` objects.
   **/
  Ptxt<Scheme> operator*(const Ptxt<Scheme>& rhs) const;

  /**
   * @brief Infix addition operator.
   * @param rhs Right hand side of addition.
   * @return Sum of the two `Ptxt` objects.
   **/
  Ptxt<Scheme> operator+(const Ptxt<Scheme>& rhs) const;

  /**
   * @brief Infix subtraction operator.
   * @param rhs Right hand side of subtraction.
   * @return Difference of the two `Ptxt` objects.
   **/
  Ptxt<Scheme> operator-(const Ptxt<Scheme>& rhs) const;

  /**
   * @brief Times equals operator with another `Ptxt`.
   * @param otherPtxt Right hand side of multiplication.
   * @return Reference to `*this` post multiplication.
   **/
  Ptxt<Scheme>& operator*=(const Ptxt<Scheme>& otherPtxt);

  /**
   * @brief Times equals operator with a single `SlotType`.
   * @param scalar Element to be multiplied across all slots.
   * @return Reference to `*this` post multiplication.
   **/
  Ptxt<Scheme>& operator*=(const SlotType& scalar);

  /**
   * @brief Times equals operator with a scalar.
   * @param scalar Element to be added across all slots.
   * @return Reference to `*this` post scalar multiplication.
   **/
  template <typename Scalar>
  Ptxt<Scheme>& operator*=(const Scalar& scalar)
  {
    assertTrue<RuntimeError>(isValid(),
                             "Cannot call operator*= on "
                             "default-constructed Ptxt");
    for (std::size_t i = 0; i < this->slots.size(); i++) {
      this->slots[i] *= scalar;
    }
    return *this;
  }

  /**
   * @brief Plus equals operator with another `Ptxt`.
   * @param otherPtxt Right hand side of addition.
   * @return Reference to `*this` post addition.
   **/
  Ptxt<Scheme>& operator+=(const Ptxt<Scheme>& otherPtxt);

  /**
   * @brief Plus equals operator with a single `SlotType`.
   * @param scalar Element to be added across all slots.
   * @return Reference to `*this` post addition.
   **/
  Ptxt<Scheme>& operator+=(const SlotType& scalar);

  /**
   * @brief Plus equals operator with a scalar.
   * @param scalar Element to be added across all slots.
   * @return Reference to `*this` post scalar addition.
   **/
  template <typename Scalar>
  Ptxt<Scheme>& operator+=(const Scalar& scalar)
  {
    assertTrue<RuntimeError>(isValid(),
                             "Cannot call operator+= on "
                             "default-constructed Ptxt");
    for (std::size_t i = 0; i < this->slots.size(); i++) {
      this->slots[i] += scalar;
    }
    return *this;
  }

  /**
   * @brief Minus equals operator with another `Ptxt`.
   * @param otherPtxt Right hand side of subtraction.
   * @return Reference to `*this` post subtraction.
   **/
  Ptxt<Scheme>& operator-=(const Ptxt<Scheme>& otherPtxt);

  /**
   * @brief Minus equals operator with a single `SlotType`.
   * @param scalar Element to be subtracted across all slots.
   * @return Reference to `*this` post subtraction.
   **/
  Ptxt<Scheme>& operator-=(const SlotType& scalar);

  /**
   * @brief Minus equals operator with a scalar.
   * @param scalar Element to be subtracted across all slots.
   * @return Reference to `*this` post scalar subtraction.
   **/
  template <typename Scalar>
  Ptxt<Scheme>& operator-=(const Scalar& scalar)
  {
    assertTrue<RuntimeError>(isValid(),
                             "Cannot call operator-= on "
                             "default-constructed Ptxt");
    for (std::size_t i = 0; i < this->slots.size(); i++) {
      this->slots[i] -= scalar;
    }
    return *this;
  }

  /**
   * @brief Negate a `Ptxt`.
   * @return Reference to `*this` post negation.
   **/
  Ptxt<Scheme>& negate();

  /**
   * @brief Add a constant to a BGV `Ptxt`.
   * @param scalar Element to be added across all slots.
   * @return Reference to `*this` post scalar addition.
   **/
  template <typename T = Scheme,
            typename Scalar,
            typename std::enable_if_t<std::is_same<T, BGV>::value>* = nullptr>
  Ptxt<Scheme>& addConstant(const Scalar& scalar)
  {
    return *this += scalar;
  }

  /**
   * @brief Add a constant to a CKKS `Ptxt`.
   * @param scalar Element to be added across all slots.
   * @return Reference to `*this` post scalar addition.
   **/
  template <typename T = Scheme,
            typename Scalar,
            typename std::enable_if_t<std::is_same<T, CKKS>::value>* = nullptr>
  Ptxt<Scheme>& addConstantCKKS(const Scalar& scalar)
  {
    return *this += scalar;
  }

  /**
   * @brief Multiplication function between two `Ptxt` objects.
   * @param otherPtxt Right hand side of multiplication.
   * @return Reference to `*this` post multiplication.
   * @note This function is equivalent to operator `*=`.
   **/
  Ptxt<Scheme>& multiplyBy(const Ptxt<Scheme>& otherPtxt);

  /**
   * @brief Multiplication function between three `Ptxt` objects.
   * @param otherPtxt1 First `Ptxt` to multiply with.
   * @param otherPtxt2 Second `Ptxt` to multiply with.
   * @return Reference to `*this` post multiplication.
   **/
  Ptxt<Scheme>& multiplyBy2(const Ptxt& otherPtxt1, const Ptxt& otherPtxt2);

  /**
   * @brief Square operation on a `Ptxt`.
   * @return Reference to `*this` post squaring.
   **/
  Ptxt<Scheme>& square();

  /**
   * @brief Cube operation on a `Ptxt`.
   * @return Reference to `*this` post cube operation.
   **/
  Ptxt<Scheme>& cube();

  /**
   * @brief Power operation to raise a `Ptxt` to an arbitrary non-negative
   *power.
   * @param e Exponent to raise the `Ptxt` by.
   * @return Reference to `*this` post raising to the power e.
   **/
  Ptxt<Scheme>& power(long e);

  /**
   * @brief Rotate slots right by specified amount (slot `i` goes to slot
   * `i+1 mod size`).
   * @param amount Number of slots to rotate by.
   * @return Reference to `*this` post rotation.
   **/
  Ptxt<Scheme>& rotate(long amount);

  /**
   * @brief Rotate slots right by specified amount along a specific dimension.
   * @param dim Dimension in which to rotate.
   * @param amount Number of slots to rotate by.
   * @return Reference to `*this` post rotation.
   **/
  Ptxt<Scheme>& rotate1D(long dim, long amount);

  /**
   * @brief Shifts slots right  by specified amount with `0` fill (slot `i` goes
   * to slot `i+1 mod size`).
   * @param amount Number of slots to shift by.
   * @return Reference to `*this` post shift.
   **/
  Ptxt<Scheme>& shift(long amount);

  /**
   * @brief Shift slots right in one dimension of the hypercube structure with
   * `0` fill.
   * @param dim Dimension in which to shift.
   * @param amount Amount by which to shift.
   * @return Reference to `*this` post shift.
   **/
  Ptxt<Scheme>& shift1D(long dim, long amount);

  /**
   * @brief Apply the automorphism a(X) -> a(X^k) mod Phi_m(X).
   * @param k Exponent of the automorphism to apply.
   * @return Reference to `*this` post automorphism application.
   * @note `k` must be an element of Zm*
   **/
  Ptxt<Scheme>& automorph(long k);

  /**
   * @brief Apply the frobenius automorphism a(X) -> a(X^(p^j)) mod Phi_m(X).
   * @param j Exponent of the automorphism to apply.
   * @return Reference to `*this` post frobenius automorphism application.
   * @note Only valid for the `BGV` scheme.
   **/
  template <typename T = Scheme,
            std::enable_if_t<std::is_same<T, BGV>::value>* = nullptr>
  Ptxt<Scheme>& frobeniusAutomorph(long j);

  /**
   * @brief Replicate single slot across all slots.
   * @param pos Position of the slot to replicate.
   * @return Reference to `*this` post replication.
   **/
  Ptxt<Scheme>& replicate(long pos);

  /**
   * @brief Generate a vector of plaintexts with each slot replicated in each
   * plaintext.
   * @return Vector of replicated plaintext slots.
   * The order of the return vector agrees with the order of the slots. i.e.
   * the i-th plaintext in the return value is a replication of `*this[i]`.
   **/
  std::vector<Ptxt<Scheme>> replicateAll() const;

  /**
   * @brief Apply complex conjugate of complex numbers in slots of a `CKKS`
   * `Ptxt` object.
   * @return Reference to `*this` post complex conjugation.
   * @note Only valid for the `CKKS` scheme.
   **/
  template <typename T = Scheme,
            std::enable_if_t<std::is_same<T, CKKS>::value>* = nullptr>
  Ptxt<Scheme>& complexConj();

  /**
   * @brief Extract the real part of a CKKS plaintext.
   * @return New plaintext containing the real part of each slot.
   * @note Only valid for the `CKKS` scheme.
   **/
  template <typename T = Scheme,
            std::enable_if_t<std::is_same<T, CKKS>::value>* = nullptr>
  Ptxt<Scheme> real() const;

  /**
   * @brief Extract the imaginary part of a CKKS plaintext.
   * @return New plaintext containing the imaginary part of each slot.
   * @note Only valid for the `CKKS` scheme.
   **/
  template <typename T = Scheme,
            std::enable_if_t<std::is_same<T, CKKS>::value>* = nullptr>
  Ptxt<Scheme> imag() const;

  /**
   * @brief Compute the running sum (each slot is the sum of the previous
   * slots).
   * @return Reference to `*this` post summation.
   **/
  Ptxt<Scheme>& runningSums();

  /**
   * @brief Compute the total sum (each slot contains the total sum of every
   * slot).
   * @return Reference to `*this` post summation.
   **/
  Ptxt<Scheme>& totalSums();

  /**
   * @brief Compute the incremental product (each slot is the product of the
   * previous slots).
   * @return Reference to `*this` post multiplication.
   **/
  Ptxt<Scheme>& incrementalProduct();

  /**
   * @brief Compute the total product (each slot contains the total product of
   * every slot).
   * @return Reference to `*this` post multiplication.
   **/
  Ptxt<Scheme>& totalProduct();

  /**
   * @brief Map all non-zero slots to `1`, keeping zero slots as zero.
   * @return Reference to `*this` post mapping.
   **/
  Ptxt<Scheme>& mapTo01();

  // No docs as the four functions have been forward-declared with docs.
  friend void deserialize<Scheme>(std::istream& is, Ptxt& ptxt);

  friend void serialize<Scheme>(std::ostream& os, const Ptxt& ptxt);

  friend std::istream& operator>><Scheme>(std::istream& is, Ptxt& ptxt);

  friend std::ostream& operator<<<Scheme>(std::ostream& os, const Ptxt& ptxt);

  /**
   * @brief Conversion function from `long` to `SlotType`.
   * @param context Context which may be needed to extract algebraic info.
   * @param slot Datum to be converted to a `SlotType`.
   * @return Converted slot.
   **/
  static SlotType convertToSlot(const Context& context, long slot);

  /**
   * @brief To be inline with the `Ctxt` interface. However for `Ptxt` this
   * means do nothing.
   * @return Reference to `*this`.
   **/
  Ptxt<Scheme>& cleanUp() { return *this; }

private:
  const Context* context;

  //! @brief The slot data of the object, where `SlotType` will typically be
  //! `std::complex<double>` (CKKS) or `helib::PolyMod` (BGV).
  std::vector<SlotType> slots;

  /**
   * @brief Helper function to convert between different indexing formats.
   * @param coords Vector of coordinates.
   * @return Index in slot array.
   *
   * Converts a vector of coordinates representing the powers of generators of
   * the Zm* group to the corresponding index in the array of slots.
   * Particularly useful when rotating or shifting.
   **/
  long coordToIndex(const std::vector<long>& coords);

  /**
   * @brief Helper function to convert between different indexing formats.
   * @param index Index in slot array.
   * @return Coordinates as a vector.
   *
   * Converts an index in the array of slots to coordinate form where the
   * coordinates are the powers of the generators of the Zm* group.
   * Particularly useful when rotating or shifting.
   **/
  std::vector<long> indexToCoord(long index);

  /**
   * @brief Convert `BGV` slots into a single polynomial of type `NTL::GF2_X` or
   * `NTL::zz_pX`.
   * @tparam templated Type where the member `RX` is the type to return.
   * @return a single polynomial representing the slots.
   **/
  template <typename type>
  typename type::RX slotsToRX() const;

  /**
   * @brief Verify that all slots have G and p^r which are compatible with
   * those inside `context`.
   * @param slots Slots to be checked against `context`.
   * @throws RuntimeError If `slots` contains anything incompatible with
   * the current `Ptxt` object.
   *
   * @note Does nothing in the CKKS case, since a std::complex<double> is
   * always compatible with `this`.
   **/
  void assertSlotsCompatible(const std::vector<SlotType>& slots) const;

  /**
   * @brief Utility function used to perform the automorphism on a specific
   * type.
   * @tparam type the polynomial type `NTL::GF2_X` or `NTL::zz_pX` to operate
   * on.
   * @param k the exponent of the automorphism, f:(X)->(X^k)
   * @return a single polynomial representing the data to be decoded into the
   * new slots.
   **/
  template <typename type>
  NTL::ZZX automorph_internal(long k);
};
} // namespace helib

#endif // HELIB_PTXT_H
