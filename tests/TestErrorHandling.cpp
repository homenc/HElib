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
#include <helib/assertions.h>
#include <helib/exceptions.h>
#include "test_common.h"
#include "gtest/gtest.h"

namespace {

//  LOGIC_ERROR

TEST(TestErrorHandling, helibLogicErrorCanBeCaughtAsHelibException)
{
  EXPECT_THROW(throw helib::LogicError("Some logic error message"),
               helib::Exception);
}

TEST(TestErrorHandling, helibLogicErrorReturnsWhatStringThroughHelibException)
{
  const std::string what("Some logic error message");
  try {
    throw helib::LogicError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibLogicErrorReturnsCStyleWhatStringThroughHelibException)
{
  const char* what = "Some logic error message";
  try {
    throw helib::LogicError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibLogicErrorCanBeCaughtAsStdException)
{
  EXPECT_THROW(throw helib::LogicError("Some logic error message"),
               std::exception);
}

TEST(TestErrorHandling, helibLogicErrorReturnsWhatStringThroughStdException)
{
  const std::string what("Some logic error message");
  try {
    throw helib::LogicError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibLogicErrorReturnsCStyleWhatStringThroughStdException)
{
  const char* what = "Some logic error message";
  try {
    throw helib::LogicError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibLogicErrorCanBeCaughtAsStdLogicError)
{
  EXPECT_THROW(throw helib::LogicError("Some logic error message"),
               std::logic_error);
}

TEST(TestErrorHandling, helibLogicErrorReturnsWhatStringThroughStdLogicError)
{
  const std::string what("Some logic error message");
  try {
    throw helib::LogicError(what);
  } catch (const std::logic_error& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibLogicErrorReturnsCStyleWhatStringThroughStdLogicError)
{
  const char* what = "Some logic error message";
  try {
    throw helib::LogicError(what);
  } catch (const std::logic_error& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

// RUNTIME_ERROR
TEST(TestErrorHandling, helibRuntimeErrorCanBeCaughtAsHelibException)
{
  EXPECT_THROW(throw helib::RuntimeError("Some runtime error message"),
               helib::Exception);
}

TEST(TestErrorHandling, helibRuntimeErrorReturnsWhatStringThroughHelibException)
{
  const std::string what("Some runtime error message");
  try {
    throw helib::RuntimeError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibRuntimeErrorReturnsCStyleWhatStringThroughHelibException)
{
  const char* what = "Some runtime error message";
  try {
    throw helib::RuntimeError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibRuntimeErrorCanBeCaughtAsStdException)
{
  EXPECT_THROW(throw helib::RuntimeError("Some runtime error message"),
               std::exception);
}

TEST(TestErrorHandling, helibRuntimeErrorReturnsWhatStringThroughStdException)
{
  const std::string what("Some runtime error message");
  try {
    throw helib::RuntimeError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibRuntimeErrorReturnsCStyleWhatStringThroughStdException)
{
  const char* what = "Some runtime error message";
  try {
    throw helib::RuntimeError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibRuntimeErrorCanBeCaughtAsStdRuntimeError)
{
  EXPECT_THROW(throw helib::RuntimeError("Some runtime error message"),
               std::runtime_error);
}

TEST(TestErrorHandling,
     helibRuntimeErrorReturnsWhatStringThroughStdRuntimeError)
{
  const std::string what("Some runtime error message");
  try {
    throw helib::RuntimeError(what);
  } catch (const std::runtime_error& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibRuntimeErrorReturnsCStyleWhatStringThroughStdRuntimeError)
{
  const char* what = "Some runtime error message";
  try {
    throw helib::RuntimeError(what);
  } catch (const std::runtime_error& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

// OutOfRangeError tests

TEST(TestErrorHandling, helibOutOfRangeErrorCanBeCaughtAsHelibOutOfRangeError)
{
  EXPECT_THROW(throw helib::OutOfRangeError("Some out of range error message"),
               helib::OutOfRangeError);
}

TEST(TestErrorHandling,
     helibOutOfRangeErrorReturnsWhatStringThroughHelibOutOfRangeError)
{
  const std::string what("Some out of range error message");
  try {
    throw helib::OutOfRangeError(what);
  } catch (const helib::OutOfRangeError& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibOutOfRangeErrorReturnsCStyleWhatStringThroughHelibOutOfRangeError)
{
  const char* what = "Some out of range error message";
  try {
    throw helib::OutOfRangeError(what);
  } catch (const helib::OutOfRangeError& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibOutOfRangeErrorCanBeCaughtAsHelibException)
{
  EXPECT_THROW(throw helib::OutOfRangeError("Some out of range error message"),
               helib::Exception);
}

TEST(TestErrorHandling,
     helibOutOfRangeErrorReturnsWhatStringThroughHelibException)
{
  const std::string what("Some out of range error message");
  try {
    throw helib::OutOfRangeError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibOutOfRangeErrorReturnsCStyleWhatStringThroughHelibException)
{
  const char* what = "Some out of range error message";
  try {
    throw helib::OutOfRangeError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibOutOfRangeErrorCanBeCaughtAsStdException)
{
  EXPECT_THROW(throw helib::OutOfRangeError("Some out of range error message"),
               std::exception);
}

TEST(TestErrorHandling,
     helibOutOfRangeErrorReturnsWhatStringThroughStdException)
{
  const std::string what("Some out of range error message");
  try {
    throw helib::OutOfRangeError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibOutOfRangeErrorReturnsCStyleWhatStringThroughStdException)
{
  const char* what = "Some out of range error message";
  try {
    throw helib::OutOfRangeError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibOutOfRangeErrorCanBeCaughtAsStdLogicError)
{
  EXPECT_THROW(throw helib::OutOfRangeError("Some out of range error message"),
               std::logic_error);
}

TEST(TestErrorHandling,
     helibOutOfRangeErrorReturnsWhatStringThroughStdLogicError)
{
  const std::string what("Some out of range error message");
  try {
    throw helib::OutOfRangeError(what);
  } catch (const std::logic_error& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibOutOfRangeErrorReturnsCStyleWhatStringThroughStdLogicError)
{
  const char* what = "Some out of range error message";
  try {
    throw helib::OutOfRangeError(what);
  } catch (const std::logic_error& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibOutOfRangeErrorCanBeCaughtAsStdOutOfRange)
{
  EXPECT_THROW(throw helib::OutOfRangeError("Some out of range error message"),
               std::out_of_range);
}

TEST(TestErrorHandling,
     helibOutOfRangeErrorReturnsWhatStringThroughStdOutOfRange)
{
  const std::string what("Some out of range error message");
  try {
    throw helib::OutOfRangeError(what);
  } catch (const std::out_of_range& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibOutOfRangeErrorReturnsCStyleWhatStringThroughStdOutOfRange)
{
  const char* what = "Some out of range error message";
  try {
    throw helib::OutOfRangeError(what);
  } catch (const std::out_of_range& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

// IOError tests

TEST(TestErrorHandling, helibIOErrorCanBeCaughtAsHelibIOError)
{
  EXPECT_THROW(throw helib::IOError("Some IO error message"), helib::IOError);
}

TEST(TestErrorHandling, helibIOErrorReturnsWhatStringThroughHelibIOError)
{
  const std::string what("Some IO error message");
  try {
    throw helib::IOError(what);
  } catch (const helib::IOError& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling, helibIOErrorReturnsCStyleWhatStringThroughHelibIOError)
{
  const char* what = "Some IO error message";
  try {
    throw helib::IOError(what);
  } catch (const helib::IOError& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibIOErrorCanBeCaughtAsHelibException)
{
  EXPECT_THROW(throw helib::IOError("Some IO error message"), helib::Exception);
}

TEST(TestErrorHandling, helibIOErrorReturnsWhatStringThroughHelibException)
{
  const std::string what("Some IO error message");
  try {
    throw helib::IOError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibIOErrorReturnsCStyleWhatStringThroughHelibException)
{
  const char* what = "Some IO error message";
  try {
    throw helib::IOError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibIOErrorCanBeCaughtAsStdException)
{
  EXPECT_THROW(throw helib::IOError("Some IO error message"), std::exception);
}

TEST(TestErrorHandling, helibIOErrorReturnsWhatStringThroughStdException)
{
  const std::string what("Some IO error message");
  try {
    throw helib::IOError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling, helibIOErrorReturnsCStyleWhatStringThroughStdException)
{
  const char* what = "Some IO error message";
  try {
    throw helib::IOError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibIOErrorCanBeCaughtAsStdRuntimeError)
{
  EXPECT_THROW(throw helib::IOError("Some IO error message"),
               std::runtime_error);
}

TEST(TestErrorHandling, helibIOErrorReturnsWhatStringThroughStdRuntimeError)
{
  const std::string what("Some IO error message");
  try {
    throw helib::IOError(what);
  } catch (const std::runtime_error& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibIOErrorReturnsCStyleWhatStringThroughStdRuntimeError)
{
  const char* what = "Some IO error message";
  try {
    throw helib::IOError(what);
  } catch (const std::runtime_error& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibIOErrorCanBeCaughtAsHelibRuntimeError)
{
  EXPECT_THROW(throw helib::IOError("Some IO error message"),
               helib::RuntimeError);
}

TEST(TestErrorHandling, helibIOErrorReturnsWhatStringThroughHelibRuntimeError)
{
  const std::string what("Some IO error message");
  try {
    throw helib::IOError(what);
  } catch (const helib::RuntimeError& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibIOErrorReturnsCStyleWhatStringThroughHelibRuntimeError)
{
  const char* what = "Some IO error message";
  try {
    throw helib::IOError(what);
  } catch (const helib::RuntimeError& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

// InvalidArgument tests

TEST(TestErrorHandling, helibInvalidArgumentCanBeCaughtAsHelibInvalidArgument)
{
  EXPECT_THROW(
      throw helib::InvalidArgument("Some invalid argument error message"),
      helib::InvalidArgument);
}

TEST(TestErrorHandling,
     helibInvalidArgumentReturnsWhatStringThroughHelibInvalidArgument)
{
  const std::string what("Some invalid argument error message");
  try {
    throw helib::InvalidArgument(what);
  } catch (const helib::InvalidArgument& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibInvalidArgumentReturnsCStyleWhatStringThroughHelibInvalidArgument)
{
  const char* what = "Some invalid argument error message";
  try {
    throw helib::InvalidArgument(what);
  } catch (const helib::InvalidArgument& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibInvalidArgumentCanBeCaughtAsHelibException)
{
  EXPECT_THROW(
      throw helib::InvalidArgument("Some invalid argument error message"),
      helib::Exception);
}

TEST(TestErrorHandling,
     helibInvalidArgumentReturnsWhatStringThroughHelibException)
{
  const std::string what("Some invalid argument error message");
  try {
    throw helib::InvalidArgument(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibInvalidArgumentReturnsCStyleWhatStringThroughHelibException)
{
  const char* what = "Some invalid argument error message";
  try {
    throw helib::InvalidArgument(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibInvalidArgumentCanBeCaughtAsStdException)
{
  EXPECT_THROW(
      throw helib::InvalidArgument("Some invalid argument error message"),
      std::exception);
}

TEST(TestErrorHandling,
     helibInvalidArgumentReturnsWhatStringThroughStdException)
{
  const std::string what("Some invalid argument error message");
  try {
    throw helib::InvalidArgument(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibInvalidArgumentReturnsCStyleWhatStringThroughStdException)
{
  const char* what = "Some invalid argument error message";
  try {
    throw helib::InvalidArgument(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibInvalidArgumentCanBeCaughtAsStdLogicError)
{
  EXPECT_THROW(throw helib::InvalidArgument("Some invalid argument message"),
               std::logic_error);
}

TEST(TestErrorHandling,
     helibInvalidArgumentReturnsWhatStringThroughStdLogicError)
{
  const std::string what("Some invalid argument message");
  try {
    throw helib::InvalidArgument(what);
  } catch (const std::logic_error& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibInvalidArgumentReturnsCStyleWhatStringThroughStdLogicError)
{
  const char* what = "Some invalid argument message";
  try {
    throw helib::InvalidArgument(what);
  } catch (const std::logic_error& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(TestErrorHandling, helibInvalidArgumentCanBeCaughtAsStdInvalidArgument)
{
  EXPECT_THROW(throw helib::InvalidArgument("Some invalid argument message"),
               std::invalid_argument);
}

TEST(TestErrorHandling,
     helibInvalidArgumentReturnsWhatStringThroughStdInvalidArgument)
{
  const std::string what("Some invalid argument message");
  try {
    throw helib::InvalidArgument(what);
  } catch (const std::invalid_argument& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(TestErrorHandling,
     helibInvalidArgumentReturnsCStyleWhatStringThroughStdInvalidArgument)
{
  const char* what = "Some invalid argument message";
  try {
    throw helib::InvalidArgument(what);
  } catch (const std::invalid_argument& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

// Cross exception catch disabled
TEST(TestErrorHandling, helibLogicErrorCannotBeCaughtAsRuntimeError)
{
  EXPECT_THROW(
      try {
        throw helib::LogicError("Some logic error message");
      } catch (const helib::RuntimeError& err){},
      helib::LogicError);
}

TEST(TestErrorHandling, helibLogicErrorCannotBeCaughtAsStdRuntimeError)
{
  EXPECT_THROW(
      try {
        throw helib::LogicError("Some logic error message");
      } catch (const std::runtime_error& err){},
      helib::LogicError);
}

TEST(TestErrorHandling, helibRuntimeErrorCannotBeCaughtAsLogicError)
{
  EXPECT_THROW(
      try {
        throw helib::RuntimeError("Some runtime error message");
      } catch (const helib::LogicError& err){},
      helib::RuntimeError);
}

TEST(TestErrorHandling, helibIOErrorCannotBeCaughtAsLogicError)
{
  EXPECT_THROW(
      try {
        throw helib::IOError("Some IO error message");
      } catch (const helib::LogicError& err){},
      helib::RuntimeError);
}

TEST(TestErrorHandling, helibIOErrorCannotBeCaughtAsInvalidArgument)
{
  EXPECT_THROW(
      try {
        throw helib::IOError("Some IO error message");
      } catch (const helib::InvalidArgument& err){},
      helib::RuntimeError);
}

TEST(TestErrorHandling, helibInvalidArgumentCannotBeCaughtAsRuntimeError)
{
  EXPECT_THROW(
      try {
        throw helib::InvalidArgument("Some invalid argument message");
      } catch (const helib::RuntimeError& err){},
      helib::InvalidArgument);
}

TEST(TestErrorHandling, helibInvalidArgumentCannotBeCaughtAsIOError)
{
  EXPECT_THROW(
      try {
        throw helib::InvalidArgument("Some invalid argument message");
      } catch (const helib::IOError& err){},
      helib::InvalidArgument);
}

TEST(TestErrorHandling, helibRuntimeErrorCannotBeCaughtAsStdLogicError)
{
  EXPECT_THROW(
      try {
        throw helib::RuntimeError("Some runtime error message");
      } catch (const std::logic_error& err){},
      helib::RuntimeError);
}

TEST(TestErrorHandling, helibIOErrorCannotBeCaughtAsStdLogicError)
{
  EXPECT_THROW(
      try {
        throw helib::IOError("Some IO error message");
      } catch (const std::logic_error& err){},
      helib::IOError);
}

TEST(TestErrorHandling, helibInvalidArgumentCannotBeCaughtAsStdRuntimeError)
{
  EXPECT_THROW(
      try {
        throw helib::InvalidArgument("Some invalid argument message");
      } catch (const std::runtime_error& err){},
      helib::InvalidArgument);
}

// Testing assertions
// Testing assertTrue
TEST(TestErrorHandling, helibAssertTrueNoThrowsOnTrue)
{
  EXPECT_NO_THROW(helib::assertTrue(true, "Value is false"));
}

TEST(TestErrorHandling, helibAssertTrueThrowsLogicErrorOnFalse)
{
  EXPECT_THROW(helib::assertTrue(false, "Value is false"), helib::LogicError);
}

TEST(TestErrorHandling, helibAssertTrueThrowsNonDefaultExceptions)
{
  EXPECT_THROW(helib::assertTrue<helib::RuntimeError>(false, "Value is false"),
               helib::RuntimeError);
}

// assertFalse
TEST(TestErrorHandling, helibAssertFalseThrowsLogicErrorOnTrue)
{
  EXPECT_THROW(helib::assertFalse(true, "Value is false"), helib::LogicError);
}

TEST(TestErrorHandling, helibAssertFalseNoThrowsOnFalse)
{
  EXPECT_NO_THROW(helib::assertFalse(false, "Value is false"));
}

TEST(TestErrorHandling, helibAssertFalseThrowsNonDefaultExceptions)
{
  EXPECT_THROW(helib::assertFalse<helib::RuntimeError>(true, "Value is false"),
               helib::RuntimeError);
}

// Testing assertEq
TEST(TestErrorHandling, helibAssertEqNoThrowsIfArgsAreEquals)
{
  int a(10), b(10);
  EXPECT_NO_THROW(helib::assertEq(a, b, "Expected a == b"));
}

TEST(TestErrorHandling, helibAssertEqThrowsDefaultLogicError)
{
  int a(10), b(40);
  EXPECT_THROW(helib::assertEq(a, b, "Expected a == b"), helib::LogicError);
}

TEST(TestErrorHandling, helibAssertEqNoThrowsWithCustomEqElemsIdentical)
{
  std::vector<int> a{1, 2};
  std::vector<int> b{1, 2};
  EXPECT_NO_THROW(helib::assertEq(a, b, "Expected a == b"));
}

TEST(TestErrorHandling, helibAssertEqThrowsWithCustomEqElemsDiffer)
{
  std::vector<int> a{1, 2};
  std::vector<int> b{5, 6};
  EXPECT_THROW(helib::assertEq(a, b, "Expected a == b"), helib::LogicError);
}

TEST(TestErrorHandling, helibAssertEqThrowsTemplatedException)
{
  int a(10), b(40);
  EXPECT_THROW(helib::assertEq<helib::RuntimeError>(a, b, "Expected a == b"),
               helib::RuntimeError);
}

// Testing assertNEq
TEST(TestErrorHandling, helibAssertNeqThrowsDefaultLogicErrorIfArgsAreEquals)
{
  int a(10), b(10);
  EXPECT_THROW(helib::assertNeq(a, b, "Expected a != b"), helib::LogicError);
}

TEST(TestErrorHandling, helibAssertNeqNoThrowsDefaultLogicErrorIfArgsDiffer)
{
  int a(10), b(40);
  EXPECT_NO_THROW(helib::assertNeq(a, b, "Expected a != b"));
}

TEST(TestErrorHandling, helibAssertNeqThrowsWithCustomEqElemsIdentical)
{
  std::vector<int> a{1, 2};
  std::vector<int> b{1, 2};
  EXPECT_THROW(helib::assertNeq(a, b, "Expected a != b"), helib::LogicError);
}

TEST(TestErrorHandling, helibAssertNeqNoThrowsWithCustomEqElemsDiffer)
{
  std::vector<int> a{1, 2};
  std::vector<int> b{5, 6};
  EXPECT_NO_THROW(helib::assertNeq(a, b, "Expected a != b"));
}

// assertNotNull
TEST(TestErrorHandling, helibAssertNotNullThrowsDefaultLogicErrorIfArgIsNull)
{
  void* p = nullptr;
  EXPECT_THROW(helib::assertNotNull(p, "Expected not null p"),
               helib::LogicError);
}

TEST(TestErrorHandling, helibAssertNotNullNoThrowsDefaultLogicErrorIfArgsDiffer)
{
  int x;
  int* p = &x;
  EXPECT_NO_THROW(helib::assertNotNull(p, "Expected not null p"));
}

TEST(TestErrorHandling, helibAssertNotNullThrowsWithDefaultConstructedSharedPtr)
{
  std::shared_ptr<int> p;
  EXPECT_THROW(helib::assertNotNull(p, "Expected not null p"),
               helib::LogicError);
}

TEST(TestErrorHandling, helibAssertNotNullNoThrowsWithValidSharedPtr)
{
  std::shared_ptr<int> p = std::make_shared<int>(10);
  EXPECT_NO_THROW(helib::assertNotNull(p, "Expected not null p"));
}

TEST(TestErrorHandling, helibAssertNotNullThrowsWithCNULL)
{
  void* p = NULL;
  EXPECT_THROW(helib::assertNotNull(p, "Expected not null p"),
               helib::LogicError);
}

TEST(TestErrorHandling, helibAssertNotNullThrowsNonStandardHelibException)
{
  int* p = nullptr;
  EXPECT_THROW(
      helib::assertNotNull<helib::RuntimeError>(p, "Expected not null p"),
      helib::RuntimeError);
}

// assertInRange
TEST(TestErrorHandling, helibAssertInRangeNoThrowsOutOfBoundErrorWhenInRangeInt)
{
  int element = 15;
  int min = 10;
  int max = 50;
  EXPECT_NO_THROW(
      helib::assertInRange(element, min, max, "elem not in [min, max)"));
}

TEST(TestErrorHandling,
     helibAssertInRangeNoThrowsOutOfBoundErrorWhenInLeftRange)
{
  int element = 10;
  int min = 10;
  int max = 50;
  EXPECT_NO_THROW(
      helib::assertInRange(element, min, max, "elem not in [min, max)"));
}

TEST(TestErrorHandling, helibAssertInRangeThrowsOutOfBoundErrorWhenLessThanMin)
{
  int element = 5;
  int min = 10;
  int max = 50;
  EXPECT_THROW(
      helib::assertInRange(element, min, max, "elem not in [min, max)"),
      helib::OutOfRangeError);
}

TEST(TestErrorHandling,
     helibAssertInRangeThrowsOutOfBoundErrorWhenGreaterThanMax)
{
  int element = 100;
  int min = 10;
  int max = 50;
  EXPECT_THROW(
      helib::assertInRange(element, min, max, "elem not in [min, max)"),
      helib::OutOfRangeError);
}

TEST(TestErrorHandling, helibAssertInRangeThrowsOutOfBoundErrorWhenEqualsMax)
{
  int element = 50;
  int min = 10;
  int max = 50;
  EXPECT_THROW(
      helib::assertInRange(element, min, max, "elem not in [min, max)"),
      helib::OutOfRangeError);
}

TEST(TestErrorHandling,
     helibAssertInRangeNoThrowsOutOfBoundErrorWhenInRangeIntAndRightInclusive)
{
  int element = 15;
  int min = 10;
  int max = 50;
  EXPECT_NO_THROW(
      helib::assertInRange(element, min, max, "elem not in [min, max)", true));
}

TEST(TestErrorHandling,
     helibAssertInRangeNoThrowsOutOfBoundErrorWhenInLeftRangeAndRightInclusive)
{
  int element = 10;
  int min = 10;
  int max = 50;
  EXPECT_NO_THROW(
      helib::assertInRange(element, min, max, "elem not in [min, max]", true));
}

TEST(TestErrorHandling,
     helibAssertInRangeThrowsOutOfBoundErrorWhenLessThanMinAndRightInclusive)
{
  int element = 5;
  int min = 10;
  int max = 50;
  EXPECT_THROW(
      helib::assertInRange(element, min, max, "elem not in [min, max]", true),
      helib::OutOfRangeError);
}

TEST(TestErrorHandling,
     helibAssertInRangeThrowsOutOfBoundErrorWhenGreaterThanMaxAndRightInclusive)
{
  int element = 100;
  int min = 10;
  int max = 50;
  EXPECT_THROW(
      helib::assertInRange(element, min, max, "elem not in [min, max]", true),
      helib::OutOfRangeError);
}

TEST(TestErrorHandling,
     helibAssertInRangeNoThrowsOutOfBoundErrorWhenEqualsMaxAndRightInclusive)
{
  int element = 50;
  int min = 10;
  int max = 50;
  EXPECT_NO_THROW(
      helib::assertInRange(element, min, max, "elem not in [min, max]", true));
}

TEST(TestErrorHandling, helibAssertInRangeThrowsNonDefaultHelibRuntimeError)
{
  int element = 100;
  int min = 10;
  int max = 50;
  EXPECT_THROW(
      helib::assertInRange<helib::RuntimeError>(element,
                                                min,
                                                max,
                                                "elem not in [min, max]",
                                                true),
      helib::RuntimeError);
}

} // namespace
