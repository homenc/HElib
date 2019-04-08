//Header guard
#ifndef ASSERTIONS_H
#define ASSERTIONS_H

#include "exceptions.h"

namespace helib {

// Exception type first and T defaulted to void so that one can specify only exception type,
// letting T be inferred from the argument passed.
template <typename ExceptionTy = ::helib::LogicError, typename T = void>
inline void assertTrue(const T& value, const std::string& message) {
  static_assert(std::is_base_of<::helib::Exception, ExceptionTy>::value, "ExceptionTy must inherit from helib::Exception");
  static_assert(std::is_same<bool, T>::value, "Type T is not boolean");
  if (!value) {
    throw ExceptionTy(message);
  }
}

// Exception type first and T defaulted to void so that one can specify only exception type,
// letting T be inferred from the argument passed.
template <typename ExceptionTy = ::helib::LogicError, typename T = void>
inline void assertFalse(T value, const std::string& message) {
  static_assert(std::is_same<bool, T>::value, "Type T is not boolean");
  static_assert(std::is_base_of<::helib::Exception, ExceptionTy>::value, "ExceptionTy must inherit from helib::Exception");
  if (value) {
    throw ExceptionTy(message);
  }
}
  
// Exception type first and T defaulted to void so that one can specify only exception type,
// letting T be inferred from the argument passed.
template <typename ExceptionTy = ::helib::LogicError, typename T = void>
inline void assertEq(const T &a, const T &b, const std::string& message) {
  static_assert(std::is_base_of<::helib::Exception, ExceptionTy>::value, "ExceptionTy must inherit from helib::Exception");
  bool value = a == b;
  if (!value) {
    throw ExceptionTy(message);
  }
}
  
// Exception type first and T defaulted to void so that one can specify only exception type,
// letting T be inferred from the argument passed.
template <typename ExceptionTy = ::helib::LogicError, typename T = void>
inline void assertNeq(const T &a, const T &b, const std::string& message) {
  static_assert(std::is_base_of<::helib::Exception, ExceptionTy>::value, "ExceptionTy must inherit from helib::Exception");
  bool value = a != b;
  if (!value) {
    throw ExceptionTy(message);
  }
}

  
// Exception type first and T defaulted to void so that one can specify only exception type,
// letting T be inferred from the argument passed.
template <typename ExceptionTy = ::helib::LogicError, typename T = void>
inline void assertNotNull(const T &p, const std::string& message) {
  static_assert(std::is_base_of<::helib::Exception, ExceptionTy>::value, "ExceptionTy must inherit from helib::Exception");
  if (p == nullptr) {
    throw ExceptionTy(message);
  }
}
  
/**
 * Function checking if the element is in the range [min,max)
 *
 * @tparam ExceptionTy type of the exception thrown.
 * @tparam T type of the element (and of the range).
 * @param elem the element to be tested.
 * @param min the left side of the range (always inclusive).
 * @param max the right side of the range (default exclusive).
 * @param message the message of the exception raised if the element is not in the range.
 * @param right_inclusive flag specifying if the right side is inclusive (default false).
 * @throw ExceptionTy exception if elem is not in the range
 * @note Exception type first and T defaulted to void so that one can specify only exception type,
 * letting T be inferred from the argument passed.
 */
template <typename ExceptionTy = ::helib::OutOfRangeError, typename T = void>
inline void assertInRange(const T &elem, const T &min, const T &max, const std::string& message, bool right_inclusive = false) {
  static_assert(std::is_base_of<::helib::Exception, ExceptionTy>::value, "ExceptionTy must inherit from helib::Exception");
  bool min_less = min <= elem;
  if (!min_less) {
    throw ExceptionTy(message);
  }
  bool max_great = right_inclusive ? elem <= max : elem < max;
  if (!max_great) {
    throw ExceptionTy(message);
  }
}

} // End of namespace
#endif // End of header guard
