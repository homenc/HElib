/* Copyright (C) 2019 IBM Corp.
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

#ifndef HELIB_POLYMOD_H
#define HELIB_POLYMOD_H

#include <NTL/ZZX.h>
#include <vector>
#include <memory>
#include <helib/PolyModRing.h>

/**
 * @file PolyMod.h
 * @brief Underlying slot type structure of BGV ptxts
 **/
namespace helib {

/**
 * @class PolyMod
 * @brief An object that contains an `NTL::ZZX` polynomial along with a
 * coefficient modulus `p2r` and a polynomial modulus `G`.
 *
 * A `PolyMod` object can be considered to be an element of
 * \f$\mathbb{Z}_{p^r}[x]/G\f$ where \f$p^r\f$ and \f$ G\f$ are the passed-in
 * parameters `p2r` and `G`.
 *
 * This allows for inter-`PolyMod` operations with the modulo operations
 * performed automatically.
 *
 * General usage:
 * @code
 * helib::PolyMod poly(input_data, p2r, G);
 * @endcode
 * Calling an operation on a default constructed `PolyMod` will throw an
 * `helib::LogicError`.
 **/
class PolyMod
{
public:
  /**
   * @brief Default constructor.
   * @note `PolyMod` objects constructed using this are marked as invalid.
   * If used for any operation whether directly on a `PolyMod` or via another
   * object such as `Ptxt` this will produce an error.
   **/
  PolyMod();

  /**
   * @brief No-data constructor.
   * @param ringDescriptor Descriptor object for the plaintext ring.
   * Contains p^r and G which are to be used for modular reduction.
   * @note This constructor does not take in input data but can be assigned
   * data later via `operator=`.
   **/
  explicit PolyMod(const std::shared_ptr<PolyModRing>& ringDescriptor);

  /**
   * @brief Constant constructor.
   * @param input Input data as a `long`.
   * @param ringDescriptor Descriptor object for the plaintext ring.
   * Contains p^r and G which are to be used for modular reduction.
   * @note This constructor accepts input data as a `long` and converts it
   * into an `NTL::ZZX` polynomial.
   **/
  PolyMod(long input, const std::shared_ptr<PolyModRing>& ringDescriptor);

  /**
   * @brief Coefficient vector constructor.
   * @param input Input data as a vector of `long` (the i'th element of the
   * vector corresponds to the coefficient of x^i).
   * @param ringDescriptor Descriptor object for the plaintext ring.
   * Contains p^r and G which are to be used for modular reduction.
   * @note This constructor accepts input data as a `std::vector<long>` and
   * converts it into an `NTL::ZZX` polynomial.
   **/
  PolyMod(const std::vector<long>& input,
          const std::shared_ptr<PolyModRing>& ringDescriptor);

  /**
   * @brief Polynomial constructor.
   * @param input Input data as an `NTL::ZZX`.
   * @param ringDescriptor Descriptor object for the plaintext ring.
   * Contains p^r and G which are to be used for modular reduction.
   * @note This constructor accepts input data as an `NTL::ZZX`.
   **/
  PolyMod(const NTL::ZZX& input,
          const std::shared_ptr<PolyModRing>& ringDescriptor);

  /**
   * @brief Default copy constructor.
   * @param input `PolyMod` object that is copied.
   **/
  PolyMod(const PolyMod& input) = default;

  /**
   * @brief Default move constructor.
   **/
  PolyMod(PolyMod&& input) noexcept = default;

  /**
   * @brief Default destructor
   */
  ~PolyMod() = default;

  /**
   * @brief Assignment operator.
   * @param input Another `PolyMod` to copy.
   **/
  PolyMod& operator=(const PolyMod& input) = default;

  /**
   * @brief default move assignment operator
   **/
  PolyMod& operator=(PolyMod&& input) =
      default; // noexcept removed as NTL has non-noexcept move constructors
  /**
   * @brief `long` assignment operator, creates a constant polynomial mod G
   * and p2r.
   * @param input `long` datum.
   **/
  PolyMod& operator=(long input);

  /**
   * @brief `std::vector<long>` assignment operator, creates a polynomial mod G
   *and p2r.
   * @param input Coefficient vector of `long` (the i'th element of the vector
   * corresponds to the coefficient of x^i).
   **/
  PolyMod& operator=(const std::vector<long>& input);

  /**
   * @brief `std::initializer_list<long>` assignment operator, creates a
   *polynomial mod G and p2r.
   * @param input coefficient vector as an initializer list of `long` (the
   * i'th element of the vector corresponds to the coefficient of x^i).
   **/
  PolyMod& operator=(const std::initializer_list<long>& input);

  /**
   * @brief `NTL::ZZX` assignment operator, creates a polynomial mod G and p2r.
   * @param input Polynomial.
   **/
  PolyMod& operator=(const NTL::ZZX& input);

  /**
   * @brief Explicit conversion to a `long`.
   * @note This function returns only the constant term even if the polynomial
   * contains higher order terms.
   **/
  explicit operator long() const;

  /**
   * @brief Explicit conversion to `std::vector<long>` (coefficient vector).
   **/
  explicit operator std::vector<long>() const;

  /**
   * @brief Explicit conversion to an `NTL::ZZX`.
   **/
  explicit operator NTL::ZZX() const;

  /**
   * @brief Gets the validity of `this`.  This will be `false` iff `this` was
   * default constructed.
   * @return `true` if `this` is valid, `false` otherwise.
   **/
  bool isValid() const;

  /**
   * @brief Get current p^r value.
   * @return The current p^r value in use.
   **/
  long getp2r() const;

  /**
   * @brief Get current G value.
   * @return The current G value in use.
   **/
  NTL::ZZX getG() const;

  /**
   * @brief Getter function that returns the data of `PolyMod` as an
   * `NTL::ZZX` const reference.
   **/
  const NTL::ZZX& getData() const;

  /**
   * @brief Equals operator between two `PolyMod` objects.
   * @param rhs Other `PolyMod` to compare to.
   * @return `true` if `PolyMod` objects are identical, `false` otherwise.
   **/
  bool operator==(const PolyMod& rhs) const;

  /**
   * @brief Equals operator between a `PolyMod` and a `long`.
   * @param rhs `long` to compare the data against.
   * @return `true` if identical, `false` otherwise.
   * @note Always returns false when called on invalid `PolyMod` objects.
   **/
  bool operator==(long rhs) const;

  /**
   * @brief Equals operator between two `PolyMod` objects.
   * @param rhs Other `PolyMod` to compare to.
   * @return `true` if identical, `false` otherwise.
   * @note Always returns false when called on invalid `PolyMod` objects.
   **/
  bool operator==(const std::vector<long>& rhs) const;

  /**
   * @brief Equals operator between two `PolyMod` objects.
   * @param rhs Other `PolyMod` to compare to.
   * @return `true` if identical, `false` otherwise.
   * @note Always returns false when called on invalid `PolyMod` objects.
   **/
  bool operator==(const NTL::ZZX& rhs) const;

  /**
   * @brief Not equals operator.
   * @param rhs Right-hand side of comparison.
   * @return `true` if not equal, `false` otherwise
   * @note Simply forwards to the correct `operator==` method.
   */
  template <typename T>
  bool operator!=(T&& rhs) const
  {
    return !operator==(std::forward<T>(rhs));
  }

  /**
   * @brief Negate function.
   * @return Reference to `*this` post negation.
   **/
  PolyMod& negate();

  /**
   * @brief Unary minus operator.
   * @return Negation of the `PolyMod`.
   **/
  PolyMod operator-() const;

  /**
   * @brief Infix multiplication operator.
   * @param rhs Right hand side of multiplication.
   * @return Product of the two `PolyMod` objects.
   **/
  PolyMod operator*(const PolyMod& rhs) const;

  /**
   * @brief Infix multiplication operator with `long`.
   * @param rhs Right hand side of multiplication.
   * @return Product of the two values.
   **/
  PolyMod operator*(long rhs) const;

  /**
   * @brief Infix multiplication operator with `NTL::ZZX`.
   * @param rhs Right hand side of multiplication.
   * @return Product of the two objects.
   **/
  PolyMod operator*(const NTL::ZZX& rhs) const;

  /**
   * @brief Infix plus operator.
   * @param rhs Right hand side of addition.
   * @return Sum of the two `PolyMod` objects.
   **/
  PolyMod operator+(const PolyMod& rhs) const;

  /**
   * @brief Infix plus operator with `long`.
   * @param rhs Right hand side of addition.
   * @return Sum of the two values.
   **/
  PolyMod operator+(long rhs) const;

  /**
   * @brief Infix plus operator with `NTL::ZZX`.
   * @param rhs Right hand side of addition.
   * @return Sum of the two objects.
   **/
  PolyMod operator+(const NTL::ZZX& rhs) const;

  /**
   * @brief Infix minus operator.
   * @param rhs Right hand side of subtraction.
   * @return Difference of the two `PolyMod` objects.
   **/
  PolyMod operator-(const PolyMod& rhs) const;

  /**
   * @brief Infix minus operator with `long`.
   * @param rhs Right hand side of subtraction.
   * @return Difference of the two values.
   **/
  PolyMod operator-(long rhs) const;

  /**
   * @brief Infix minus operator with `NTL::ZZX`.
   * @param rhs Right hand side of subtraction.
   * @return Difference of the two objects.
   **/
  PolyMod operator-(const NTL::ZZX& rhs) const;

  /**
   * @brief Times equals operator with `PolyMod` rhs.
   * @param otherPoly Right hand side of multiplication.
   * @return Reference to `*this` post multiplication.
   **/
  PolyMod& operator*=(const PolyMod& otherPoly);

  /**
   * @brief Times equals operator with `long` rhs.
   * @param scalar Right hand side of multiplication.
   * @return Reference to `*this` post multiplication.
   **/
  PolyMod& operator*=(long scalar);

  /**
   * @brief Times equals operator with `NTL::ZZX` rhs.
   * @param otherPoly Right hand side of multiplication.
   * @return Reference to `*this` post multiplication.
   **/
  PolyMod& operator*=(const NTL::ZZX& otherPoly);

  /**
   * @brief Plus equals operator with `PolyMod` rhs.
   * @param otherPoly Right hand side of addition.
   * @return Reference to `*this` post addition.
   **/
  PolyMod& operator+=(const PolyMod& otherPoly);

  /**
   * @brief Plus equals operator with `long` rhs.
   * @param scalar Right hand side of addition.
   * @return Reference to `*this` post addition.
   **/
  PolyMod& operator+=(long scalar);

  /**
   * @brief Plus equals operator with `NTL::ZZX` rhs.
   * @param otherPoly Right hand side of addition.
   * @return Reference to `*this` post addition.
   **/
  PolyMod& operator+=(const NTL::ZZX& otherPoly);

  /**
   * @brief Minus equals operator with `PolyMod` rhs.
   * @param otherPoly Right hand side of subtraction.
   * @return Reference to `*this` post subtraction.
   **/
  PolyMod& operator-=(const PolyMod& otherPoly);

  /**
   * @brief Minus equals operator with `long` rhs.
   * @param scalar Right hand side of subtraction.
   * @return Reference to `*this` post subtraction.
   **/
  PolyMod& operator-=(long scalar);

  /**
   * @brief Minus equals operator with `NTL::ZZX` rhs.
   * @param otherPoly Right hand side of subtraction.
   * @return Reference to `*this` post subtraction.
   **/
  PolyMod& operator-=(const NTL::ZZX& otherPoly);

  /**
   * @brief Deserialize a `PolyMod` object from the input stream `is`.
   * @param is Input `std::istream`.
   * @param poly Destination `PolyMod` object.
   * @throws IOError if the stream is badly formatted (i.e. it is not delimited
   * by '[' and ']').
   * @note `poly` must be constructed with an appropriate p2r and G @b BEFORE
   * calling this function. For example,
   * @code
   * PolyMod my_poly(p2r, G);
   * deserialize(std::cin, my_poly);
   * @endcode
   *
   * The input stream has to be formatted as a comma-separated list surrounded
   * by '[' and ']'.\n
   * Each element of the list will be deserialized as a coefficient of the
   * polynomial.\n
   * For example '['coef0', 'coef1', 'coef2']' will be deserialized as a
   * `PolyMod` object `poly` where `poly[0]=coef0`, `poly[1]=coef1`,
   * `poly[2]=coef2` and `poly[i]=0` for `i>2`.
   **/
  friend void deserialize(std::istream& is, PolyMod& poly);

  /**
   * @brief Serialize a `PolyMod` to the output stream `os`.
   * @param os Output `std::ostream`.
   * @param poly `PolyMod` object to be written.
   * @return Input `std::ostream` post writing.
   * @note p2r and G are not serialized, see note of `deserialize`.
   *
   * The output stream will be formatted as a comma-separated list surrounded by
   * '[' and ']'.\n
   * Each coefficient of `poly` will be serialized in an element of such list by
   * the `>>` operator.\n
   * For example if we have a `PolyMod` object `poly` such that `poly[0]=coef0`,
   * `poly[1]=coef1`, `poly[2]=coef2`, and `poly[i]=0` for `i>2`, it will be
   * serialized as '['coef0', 'coef1', 'coef2']'.
   **/
  friend void serialize(std::ostream& os, const PolyMod& slot);

  /**
   * @brief Input shift operator.
   * @param is Input `std::istream`.
   * @param poly Destination `PolyMod` object.
   * @return Input `std::istream` post reading.
   * @note `poly` must be constructed with an appropriate p2r and G @b BEFORE
   * calling this function. For example,
   * @code
   * PolyMod my_poly(p2r, G);
   * std::cin >> my_poly;
   * @endcode
   *
   * The input stream has to be formatted as a comma-separated list surrounded
   * by '[' and ']'.\n
   * Each element of the list will be deserialized as a coefficient of the
   * polynomial.\n
   * If the number of tokens in the list is less than the number of
   * coefficients, the higher-degree coefficients will be padded by 0.\n
   * For example '['coef0', 'coef1', 'coef2']' will be deserialized as a
   * `PolyMod` object `poly` where `poly[0]=coef0`, `poly[1]=coef1`,
   * `poly[2]=coef2` and `poly[i]=0` for `i>2`.
   **/
  friend std::istream& operator>>(std::istream& is, PolyMod& poly);

  /**
   * @brief Output shift operator.
   * @param os Output `std::ostream`.
   * @param poly `PolyMod` object to be written.
   * @return Input `std::ostream` post writing.
   * @note p2r and G are not serialized, see note of `operator>>`.
   *
   * The output stream will be formatted as a comma-separated list surrounded by
   * '[' and ']'.\n
   * Each coefficient of `poly` will be serialized in an element of such list by
   * the `>>` operator.\n
   * For example if we have a `PolyMod` object `poly` such that `poly[0]=coef0`,
   * `poly[1]=coef1`, `poly[2]=coef2`, and `poly[i]=0` for `i>2`, it will be
   * serialized as '['coef0', 'coef1', 'coef2']'.
   **/
  friend std::ostream& operator<<(std::ostream& os, const PolyMod& poly);

private:
  /**
   * @brief Object containing algebraic information required for modular
   * operations (p^r and G).
   **/
  std::shared_ptr<PolyModRing> ringDescriptor;

  /**
   * @brief The polynomial data this object contains.
   **/
  NTL::ZZX data;

  /**
   * @brief Perform modular reduction on the current data.  To be called
   * after all operations.
   **/
  void modularReduce();

  /**
   * @brief Assert that `poly` is valid, throwing an `helib::LogicError` if
   * not.
   * @param poly Object whose validity is to be verified.
   **/
  static void assertValidity(const PolyMod& poly);

  /**
   * @brief Assert that `lhs` and `rhs` are both valid and interoperable.
   **/
  static void assertInterop(const PolyMod& lhs, const PolyMod& rhs);
};

} // namespace helib
#endif
