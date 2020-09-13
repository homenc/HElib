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
#ifndef HELIB_EXCEPTIONS_H
#define HELIB_EXCEPTIONS_H

#include <exception>
#include <stdexcept>
#include <sstream>

/**
 * @file exceptions.h
 * @brief Various HElib-specific exception types.
 *
 * This is largely a mirror image of the standard library exceptions, with the
 * added ancestor of `Exception`.  This allows one to distinguish between
 * general exceptions and those specifically thrown by HElib.  For example:
 *
 * ```
 try {
   // Some code including calls to HElib
 }
 catch(const Exception& err) {
   // HElib error handling
 }
 catch(const std::exception& err) {
   // Generic error handling
 }
 * ```
 *
 * To make sure that this is a pattern that can be used, we should only throw
 * exceptions derived from `Exception` wherever possible.
 */

/* @namespace helib*/
namespace helib {

/**
 * @class Exception
 * @brief Base class that other HElib exception classes inherit from.
 */
class Exception
{
public:
  virtual ~Exception() = default;
  /** @fn what returns a pointer to the string of the exception message*/
  virtual const char* what() const noexcept = 0;

protected:
  Exception() = default;
};

/**
 * @class LogicError
 * @brief Inherits from Exception and std::logic_error.
 */
class LogicError : public std::logic_error, public ::helib::Exception
{
public:
  explicit LogicError(const std::string& what_arg) :
      std::logic_error(what_arg){};
  explicit LogicError(const char* what_arg) : std::logic_error(what_arg){};
  virtual ~LogicError(){};
  /** @fn what returns a pointer to the string of the exception message*/
  virtual const char* what() const noexcept override
  {
    return std::logic_error::what();
  };
};

/**
 * @class OutOfRangeError
 * @brief Inherits from Exception and std::out_of_range.
 */
class OutOfRangeError : public std::out_of_range, public ::helib::Exception
{
public:
  explicit OutOfRangeError(const std::string& what_arg) :
      std::out_of_range(what_arg){};
  explicit OutOfRangeError(const char* what_arg) :
      std::out_of_range(what_arg){};
  virtual ~OutOfRangeError(){};
  /** @fn what returns a pointer to the string of the exception message*/
  virtual const char* what() const noexcept override
  {
    return std::out_of_range::what();
  };
};

/**
 * @class RuntimeError
 * @brief Inherits from Exception and std::runtime_error.
 */
class RuntimeError : public std::runtime_error, public ::helib::Exception
{
public:
  explicit RuntimeError(const std::string& what_arg) :
      std::runtime_error(what_arg){};
  explicit RuntimeError(const char* what_arg) : std::runtime_error(what_arg){};
  virtual ~RuntimeError(){};
  /** @fn what returns a pointer to the string of the exception message*/
  virtual const char* what() const noexcept override
  {
    return std::runtime_error::what();
  };
};

/**
 * @class IOError
 * @brief Inherits from Exception and std::runtime_error.
 */
class IOError : public helib::RuntimeError
{
public:
  explicit IOError(const std::string& what_arg) : RuntimeError(what_arg){};
  explicit IOError(const char* what_arg) : RuntimeError(what_arg){};
  virtual ~IOError(){};
  /** @fn what returns a pointer to the string of the exception message */
  virtual const char* what() const noexcept override
  {
    return std::runtime_error::what();
  };
};

/**
 * @class InvalidArgument
 * @brief Inherits from Exception and std::invalid_argument.
 */
class InvalidArgument : public std::invalid_argument, public ::helib::Exception
{
public:
  explicit InvalidArgument(const std::string& what_arg) :
      std::invalid_argument(what_arg){};
  explicit InvalidArgument(const char* what_arg) :
      std::invalid_argument(what_arg){};
  virtual ~InvalidArgument(){};
  /** @fn what returns a pointer to the string of the exception message*/
  virtual const char* what() const noexcept override
  {
    return std::invalid_argument::what();
  };
};

} // namespace helib
#endif // End of header guard
