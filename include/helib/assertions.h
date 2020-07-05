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

// Header guard
#ifndef HELIB_ASSERTIONS_H
#define HELIB_ASSERTIONS_H

#include <helib/exceptions.h>

/**
 * @file assertions.h
 * @brief Various assertion functions for use within HElib.
 * These are meant as a replacement for C-style asserts, which (undesirably)
 * exit the process without giving the opportunity to clean up and shut down
 * gracefully.  Instead, these functions will throw an exception if their
 * conditions are violated.
 *
 * General usage is of the form:
 *
 * ```
 helib::assertEq(my_variable, 5l, "my_variable must equal 5!");
 * ```
 *
 * Most of these functions will result in an `helib::LogicError` being thrown
 * if their conditions are violated, except for `assertInRange`, which will
 * throw an `helib::OutOfRange`.
 * However, if one wishes to throw some other helib exception, one can specify
 * it as a template argument as follows.  For example:
 *
 * ```
 helib::assertTrue<helib::InvalidArgument>(some_var > 0, "some_var must be
 positive!");
 * ```
 *
 * The full set of helib exceptions can be found in exceptions.h.
 */

namespace helib {

/**
 * Function throwing an exception of type ExceptionTy if the condition is false.
 *
 * @tparam ExceptionTy type of the exception thrown.
 * @tparam T type of the condition being checked (must be a bool).
 * @param value the condition being checked.
 * @param message the message of the exception raised if the condition is false.
 * @throw ExceptionTy exception if condition is false.
 * @note ExceptionTy first and T defaulted to void so that one can specify only
 * ExceptionTy, letting T be inferred from the argument passed.
 */
template <typename ExceptionTy = ::helib::LogicError, typename T = void>
inline void assertTrue(const T& value, const std::string& message)
{
  static_assert(std::is_base_of<::helib::Exception, ExceptionTy>::value,
                "ExceptionTy must inherit from helib::Exception");
  static_assert(std::is_same<bool, T>::value, "Type T is not boolean");
  if (!value) {
    throw ExceptionTy(message);
  }
}

/**
 * Function throwing an exception of type ExceptionTy if the condition is true.
 *
 * @tparam ExceptionTy type of the exception thrown.
 * @tparam T type of the condition being checked (must be a bool).
 * @param value the condition being checked.
 * @param message the message of the exception raised if the condition is true.
 * @throw ExceptionTy exception if condition is true.
 * @note ExceptionTy first and T defaulted to void so that one can specify only
 * ExceptionTy, letting T be inferred from the argument passed.
 */
template <typename ExceptionTy = ::helib::LogicError, typename T = void>
inline void assertFalse(T value, const std::string& message)
{
  static_assert(std::is_same<bool, T>::value, "Type T is not boolean");
  static_assert(std::is_base_of<::helib::Exception, ExceptionTy>::value,
                "ExceptionTy must inherit from helib::Exception");
  if (value) {
    throw ExceptionTy(message);
  }
}

/**
 * Function throwing an exception of type ExceptionTy if the two arguments are
 * not equal.
 *
 * @tparam ExceptionTy type of the exception thrown.
 * @tparam T type of the elements to be compared.
 * @param a the first element to be compared.
 * @param b the second element to be compared.
 * @param message the message of the exception raised if the two values are not
 * equal.
 * @throw ExceptionTy exception if the two values are not equal.
 * @note ExceptionTy first and T defaulted to void so that one can specify only
 * ExceptionTy, letting T be inferred from the argument passed.
 */
template <typename ExceptionTy = ::helib::LogicError, typename T = void>
inline void assertEq(const T& a, const T& b, const std::string& message)
{
  static_assert(std::is_base_of<::helib::Exception, ExceptionTy>::value,
                "ExceptionTy must inherit from helib::Exception");
  bool value = a == b;
  if (!value) {
    throw ExceptionTy(message);
  }
}

/**
 * Function throwing an exception of type ExceptionTy if the two arguments are
 * equal.
 *
 * @tparam ExceptionTy type of the exception thrown.
 * @tparam T type of the elements to be compared.
 * @param a the first element to be compared.
 * @param b the second element to be compared.
 * @param message the message of the exception raised if the two values are
 * equal.
 * @throw ExceptionTy exception if the two values are equal.
 * @note ExceptionTy first and T defaulted to void so that one can specify only
 * ExceptionTy, letting T be inferred from the argument passed.
 */
template <typename ExceptionTy = ::helib::LogicError, typename T = void>
inline void assertNeq(const T& a, const T& b, const std::string& message)
{
  static_assert(std::is_base_of<::helib::Exception, ExceptionTy>::value,
                "ExceptionTy must inherit from helib::Exception");
  bool value = a != b;
  if (!value) {
    throw ExceptionTy(message);
  }
}

/**
 * Function throwing an exception of type ExceptionTy if the argument is
 * nullptr.
 *
 * @tparam ExceptionTy type of the exception thrown.
 * @tparam T type of the element.
 * @param p the element to be tested.
 * @param message the message of the exception raised if the element is nullptr.
 * @throw ExceptionTy exception if p is nullptr.
 * @note ExceptionTy first and T defaulted to void so that one can specify only
 * ExceptionTy, letting T be inferred from the argument passed.
 */
template <typename ExceptionTy = ::helib::LogicError, typename T = void>
inline void assertNotNull(const T& p, const std::string& message)
{
  static_assert(std::is_base_of<::helib::Exception, ExceptionTy>::value,
                "ExceptionTy must inherit from helib::Exception");
  if (p == nullptr) {
    throw ExceptionTy(message);
  }
}

/**
 * Function throwing an exception of type ExceptionTy if the element is in the
 * range [min,max) or [min, max]
 *
 * @tparam ExceptionTy type of the exception thrown.
 * @tparam T type of the element (and of the range).
 * @param elem the element to be tested.
 * @param min the left side of the range (always inclusive).
 * @param max the right side of the range (default exclusive).
 * @param message the message of the exception raised if the element is not in
 * the range.
 * @param right_inclusive flag specifying if the right side is inclusive
 * (default false).
 * @throw ExceptionTy exception if elem is not in the range
 * @note ExceptionTy first and T defaulted to void so that one can specify only
 * ExceptionTy, letting T be inferred from the argument passed.
 */
template <typename ExceptionTy = ::helib::OutOfRangeError, typename T = void>
inline void assertInRange(const T& elem,
                          const T& min,
                          const T& max,
                          const std::string& message,
                          bool right_inclusive = false)
{
  static_assert(std::is_base_of<::helib::Exception, ExceptionTy>::value,
                "ExceptionTy must inherit from helib::Exception");
  bool min_less = min <= elem;
  if (!min_less) {
    throw ExceptionTy(message);
  }
  bool max_great = right_inclusive ? elem <= max : elem < max;
  if (!max_great) {
    throw ExceptionTy(message);
  }
}

} // namespace helib
#endif // End of header guard
