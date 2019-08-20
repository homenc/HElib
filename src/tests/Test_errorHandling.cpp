#include "assertions.h"
#include "exceptions.h"
#include "test_common.h"
#include "gtest/gtest.h"

namespace {
  
//  LOGIC_ERROR

TEST(Test_error_handling, helib_LogicError_can_be_caught_as_helib_exception)
{
  EXPECT_THROW(throw helib::LogicError("Some logic error message"), helib::Exception);
}

TEST(Test_error_handling, helib_LogicError_returns_what_string_through_helib_exception)
{
  const std::string what("Some logic error message");
  try {
    throw helib::LogicError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_LogicError_returns_c_style_what_string_through_helib_exception)
{
  const char* what = "Some logic error message";
  try {
    throw helib::LogicError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(Test_error_handling, helib_LogicError_can_be_caught_as_std_exception)
{
  EXPECT_THROW(throw helib::LogicError("Some logic error message"), std::exception);
}

TEST(Test_error_handling, helib_LogicError_returns_what_string_through_std_exception)
{
  const std::string what("Some logic error message");
  try {
    throw helib::LogicError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_LogicError_returns_c_style_what_string_through_std_exception)
{
  const char* what = "Some logic error message";
  try {
    throw helib::LogicError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(Test_error_handling, helib_LogicError_can_be_caught_as_std_logic_error)
{
  EXPECT_THROW(throw helib::LogicError("Some logic error message"), std::logic_error);
}

TEST(Test_error_handling, helib_LogicError_returns_what_string_through_std_logic_error)
{
  const std::string what("Some logic error message");
  try {
    throw helib::LogicError(what);
  } catch (const std::logic_error& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_LogicError_returns_c_style_what_string_through_std_logic_error)
{
  const char* what = "Some logic error message";
  try {
    throw helib::LogicError(what);
  } catch (const std::logic_error& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

  
//RUNTIME_ERROR
TEST(Test_error_handling, helib_RuntimeError_can_be_caught_as_helib_exception)
{
  EXPECT_THROW(throw helib::RuntimeError("Some runtime error message"), helib::Exception);
}

TEST(Test_error_handling, helib_RuntimeError_returns_what_string_through_helib_exception)
{
  const std::string what("Some runtime error message");
  try {
    throw helib::RuntimeError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_RuntimeError_returns_c_style_what_string_through_helib_exception)
{
  const char* what = "Some runtime error message";
  try {
    throw helib::RuntimeError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}
  
TEST(Test_error_handling, helib_RuntimeError_can_be_caught_as_std_exception)
{
  EXPECT_THROW(throw helib::RuntimeError("Some runtime error message"), std::exception);
}

TEST(Test_error_handling, helib_RuntimeError_returns_what_string_through_std_exception)
{
  const std::string what("Some runtime error message");
  try {
    throw helib::RuntimeError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_RuntimeError_returns_c_style_what_string_through_std_exception)
{
  const char* what = "Some runtime error message";
  try {
    throw helib::RuntimeError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(Test_error_handling, helib_RuntimeError_can_be_caught_as_std_runtime_error)
{
  EXPECT_THROW(throw helib::RuntimeError("Some runtime error message"), std::runtime_error);
}

TEST(Test_error_handling, helib_RuntimeError_returns_what_string_through_std_runtime_error)
{
  const std::string what("Some runtime error message");
  try {
    throw helib::RuntimeError(what);
  } catch (const std::runtime_error& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_RuntimeError_returns_c_style_what_string_through_std_runtime_error)
{
  const char* what = "Some runtime error message";
  try {
    throw helib::RuntimeError(what);
  } catch (const std::runtime_error& err) {
    EXPECT_STREQ(err.what(), what);
  }
}
  

//OutOfRangeError tests
  
TEST(Test_error_handling, helib_OutOfRangeError_can_be_caught_as_helib_OutOfRangeError)
{
  EXPECT_THROW(throw helib::OutOfRangeError("Some out of range error message"), helib::OutOfRangeError);
}

TEST(Test_error_handling, helib_OutOfRangeError_returns_what_string_through_helib_OutOfRangeError)
{
  const std::string what("Some out of range error message");
  try {
    throw helib::OutOfRangeError(what);
  } catch (const helib::OutOfRangeError& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_OutOfRangeError_returns_c_style_what_string_through_helib_OutOfRangeError)
{
  const char* what = "Some out of range error message";
  try {
    throw helib::OutOfRangeError(what);
  } catch (const helib::OutOfRangeError& err) {
    EXPECT_STREQ(err.what(), what);
  }
}
  
TEST(Test_error_handling, helib_OutOfRangeError_can_be_caught_as_helib_exception)
{
  EXPECT_THROW(throw helib::OutOfRangeError("Some out of range error message"), helib::Exception);
}

TEST(Test_error_handling, helib_OutOfRangeError_returns_what_string_through_helib_exception)
{
  const std::string what("Some out of range error message");
  try {
    throw helib::OutOfRangeError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_OutOfRangeError_returns_c_style_what_string_through_helib_exception)
{
  const char* what = "Some out of range error message";
  try {
    throw helib::OutOfRangeError(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(Test_error_handling, helib_OutOfRangeError_can_be_caught_as_std_exception)
{
  EXPECT_THROW(throw helib::OutOfRangeError("Some out of range error message"), std::exception);
}

TEST(Test_error_handling, helib_OutOfRangeError_returns_what_string_through_std_exception)
{
  const std::string what("Some out of range error message");
  try {
    throw helib::OutOfRangeError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_OutOfRangeError_returns_c_style_what_string_through_std_exception)
{
  const char* what = "Some out of range error message";
  try {
    throw helib::OutOfRangeError(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}
  
TEST(Test_error_handling, helib_OutOfRangeError_can_be_caught_as_std_logic_error)
{
  EXPECT_THROW(throw helib::OutOfRangeError("Some out of range error message"), std::logic_error);
}
  
TEST(Test_error_handling, helib_OutOfRangeError_returns_what_string_through_std_logic_error)
{
  const std::string what("Some out of range error message");
  try {
    throw helib::OutOfRangeError(what);
  } catch (const std::logic_error& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_OutOfRangeError_returns_c_style_what_string_through_std_logic_error)
{
  const char* what = "Some out of range error message";
  try {
    throw helib::OutOfRangeError(what);
  } catch (const std::logic_error& err) {
    EXPECT_STREQ(err.what(), what);
  }
}
  
  TEST(Test_error_handling, helib_OutOfRangeError_can_be_caught_as_std_out_of_range)
{
  EXPECT_THROW(throw helib::OutOfRangeError("Some out of range error message"), std::out_of_range);
}

  TEST(Test_error_handling, helib_OutOfRangeError_returns_what_string_through_std_out_of_range)
{
  const std::string what("Some out of range error message");
  try {
    throw helib::OutOfRangeError(what);
  } catch (const std::out_of_range& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

  TEST(Test_error_handling, helib_OutOfRangeError_returns_c_style_what_string_through_std_out_of_range)
{
  const char* what = "Some out of range error message";
  try {
    throw helib::OutOfRangeError(what);
  } catch (const std::out_of_range& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

//InvalidArgument tests

TEST(Test_error_handling, helib_InvalidArgument_can_be_caught_as_helib_InvalidArgument)
{
  EXPECT_THROW(throw helib::InvalidArgument("Some invalid argument error message"), helib::InvalidArgument);
}

TEST(Test_error_handling, helib_InvalidArgument_returns_what_string_through_helib_InvalidArgument)
{
  const std::string what("Some invalid argument error message");
  try {
    throw helib::InvalidArgument(what);
  } catch (const helib::InvalidArgument& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_InvalidArgument_returns_c_style_what_string_through_helib_InvalidArgument)
{
  const char* what = "Some invalid argument error message";
  try {
    throw helib::InvalidArgument(what);
  } catch (const helib::InvalidArgument& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(Test_error_handling, helib_InvalidArgument_can_be_caught_as_helib_exception)
{
  EXPECT_THROW(throw helib::InvalidArgument("Some invalid argument error message"), helib::Exception);
}

TEST(Test_error_handling, helib_InvalidArgument_returns_what_string_through_helib_exception)
{
  const std::string what("Some invalid argument error message");
  try {
    throw helib::InvalidArgument(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_InvalidArgument_returns_c_style_what_string_through_helib_exception)
{
  const char* what = "Some invalid argument error message";
  try {
    throw helib::InvalidArgument(what);
  } catch (const helib::Exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(Test_error_handling, helib_InvalidArgument_can_be_caught_as_std_exception)
{
  EXPECT_THROW(throw helib::InvalidArgument("Some invalid argument error message"), std::exception);
}

TEST(Test_error_handling, helib_InvalidArgument_returns_what_string_through_std_exception)
{
  const std::string what("Some invalid argument error message");
  try {
    throw helib::InvalidArgument(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_InvalidArgument_returns_c_style_what_string_through_std_exception)
{
  const char* what = "Some invalid argument error message";
  try {
    throw helib::InvalidArgument(what);
  } catch (const std::exception& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(Test_error_handling, helib_InvalidArgument_can_be_caught_as_std_logic_error)
{
  EXPECT_THROW(throw helib::InvalidArgument("Some invalid argument message"), std::logic_error);
}

TEST(Test_error_handling, helib_InvalidArgument_returns_what_string_through_std_logic_error)
{
  const std::string what("Some invalid argument message");
  try {
    throw helib::InvalidArgument(what);
  } catch (const std::logic_error& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_InvalidArgument_returns_c_style_what_string_through_std_logic_error)
{
  const char* what = "Some invalid argument message";
  try {
    throw helib::InvalidArgument(what);
  } catch (const std::logic_error& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

TEST(Test_error_handling, helib_InvalidArgument_can_be_caught_as_std_invalid_argument)
{
  EXPECT_THROW(throw helib::InvalidArgument("Some invalid argument message"), std::invalid_argument);
}

TEST(Test_error_handling, helib_InvalidArgument_returns_what_string_through_std_invalid_argument)
{
  const std::string what("Some invalid argument message");
  try {
    throw helib::InvalidArgument(what);
  } catch (const std::invalid_argument& err) {
    EXPECT_STREQ(err.what(), what.c_str());
  }
}

TEST(Test_error_handling, helib_InvalidArgument_returns_c_style_what_string_through_std_invalid_argument)
{
  const char* what = "Some invalid argument message";
  try {
    throw helib::InvalidArgument(what);
  } catch (const std::invalid_argument& err) {
    EXPECT_STREQ(err.what(), what);
  }
}

//Cross exception catch disabled
TEST(Test_error_handling, helib_LogicError_cannot_be_caught_as_RuntimeError)
{
  EXPECT_THROW(try { throw helib::LogicError("Some logic error message"); }
               catch (const helib::RuntimeError& err) { },
               helib::LogicError);
}

TEST(Test_error_handling, helib_LogicError_cannot_be_caught_as_std_runtime_error)
{
  EXPECT_THROW(try { throw helib::LogicError("Some logic error message"); }
               catch (const std::runtime_error& err) { },
               helib::LogicError);
}

TEST(Test_error_handling, helib_RuntimeError_cannot_be_caught_as_LogicError)
{
  EXPECT_THROW(try { throw helib::RuntimeError("Some runtime error message"); }
               catch (const helib::LogicError& err) { },
               helib::RuntimeError);
}

TEST(Test_error_handling, helib_RuntimeError_cannot_be_caught_as_std_logic_error)
{
  EXPECT_THROW(try { throw helib::RuntimeError("Some runtime error message"); }
               catch (const std::logic_error& err) { },
               helib::RuntimeError);
}

TEST(Test_error_handling, helib_InvalidArgument_cannot_be_caught_as_std_runtime_error)
{
  EXPECT_THROW(try { throw helib::InvalidArgument("Some invalid argument message"); }
               catch (const std::runtime_error& err) {},
               helib::InvalidArgument);
}

//Testing assertions
//Testing assertTrue
TEST(Test_error_handling, helib_assert_true_no_throws_on_true)
{
  EXPECT_NO_THROW(helib::assertTrue(true, "Value is false"));
}
  
TEST(Test_error_handling, helib_assert_true_throws_LogicError_on_false)
{
  EXPECT_THROW(helib::assertTrue(false, "Value is false"), helib::LogicError);
}
  
TEST(Test_error_handling, helib_assert_true_throws_non_default_exceptions) {
  EXPECT_THROW(helib::assertTrue<helib::RuntimeError>(false, "Value is false"), helib::RuntimeError);
}
  
// assertFalse
TEST(Test_error_handling, helib_assertFalse_throws_LogicError_on_true)
{
  EXPECT_THROW(helib::assertFalse(true, "Value is false"), helib::LogicError);
}
  
TEST(Test_error_handling, helib_assertFalse_no_throws_on_false)
{
  EXPECT_NO_THROW(helib::assertFalse(false, "Value is false"));
}

TEST(Test_error_handling, helib_assertFalse_throws_non_default_exceptions) {
  EXPECT_THROW(helib::assertFalse<helib::RuntimeError>(true, "Value is false"), helib::RuntimeError);
}

//Testing assertEq
TEST(Test_error_handling, helib_assertEq_no_throws_if_args_are_equals)
{
  int a(10), b(10);
  EXPECT_NO_THROW(helib::assertEq(a, b, "Expected a == b"));
}

TEST(Test_error_handling, helib_assertEq_throws_default_LogicError)
{
  int a(10), b(40);
  EXPECT_THROW(helib::assertEq(a, b, "Expected a == b"), helib::LogicError);
}
  
TEST(Test_error_handling, helib_assertEq_no_throws_with_custom_eq_elems_identical) {
  std::vector<int> a{1, 2};
  std::vector<int> b{1, 2};
  EXPECT_NO_THROW(helib::assertEq(a, b, "Expected a == b"));
}

TEST(Test_error_handling, helib_assertEq_throws_with_custom_eq_elems_differ) {
  std::vector<int> a{1, 2};
  std::vector<int> b{5, 6};
  EXPECT_THROW(helib::assertEq(a, b, "Expected a == b"), helib::LogicError);
}

TEST(Test_error_handling, helib_assertEq_throws_templated_exception)
{
  int a(10), b(40);
  EXPECT_THROW(helib::assertEq<helib::RuntimeError>(a, b, "Expected a == b"), helib::RuntimeError);
}
  
//Testing assertNEq
TEST(Test_error_handling, helib_assertNeq_throws_default_LogicError_if_args_are_equals)
{
  int a(10), b(10);
  EXPECT_THROW(helib::assertNeq(a, b, "Expected a != b"), helib::LogicError);
}

TEST(Test_error_handling, helib_assertNeq_no_throws_default_LogicError_if_args_differ)
{
  int a(10), b(40);
  EXPECT_NO_THROW(helib::assertNeq(a, b, "Expected a != b"));
}

TEST(Test_error_handling, helib_assertNeq_throws_with_custom_eq_elems_identical) {
  std::vector<int> a{1, 2};
  std::vector<int> b{1, 2};
  EXPECT_THROW(helib::assertNeq(a, b, "Expected a != b"), helib::LogicError);
}

TEST(Test_error_handling, helib_assertNeq_no_throws_with_custom_eq_elems_differ) {
  std::vector<int> a{1, 2};
  std::vector<int> b{5, 6};
  EXPECT_NO_THROW(helib::assertNeq(a, b, "Expected a != b"));
}
  
//assertNotNull
TEST(Test_error_handling, helib_assertNotNull_throws_default_LogicError_if_arg_is_null)
{
  void *p = nullptr;
  EXPECT_THROW(helib::assertNotNull(p, "Expected not null p"), helib::LogicError);
}

TEST(Test_error_handling, helib_assertNotNull_no_throws_default_LogicError_if_args_differ)
{
  int x;
  int *p = &x;
  EXPECT_NO_THROW(helib::assertNotNull(p, "Expected not null p"));
}

TEST(Test_error_handling, helib_assertNotNull_throws_with_default_constructed_shared_ptr) {
  std::shared_ptr<int> p;
  EXPECT_THROW(helib::assertNotNull(p, "Expected not null p"), helib::LogicError);
}
  
TEST(Test_error_handling, helib_assertNotNull_no_throws_with_valid_shared_ptr) {
  std::shared_ptr<int> p = std::make_shared<int>(10);
  EXPECT_NO_THROW(helib::assertNotNull(p, "Expected not null p"));
}
  
TEST(Test_error_handling, helib_assertNotNull_throws_with_c_NULL) {
  void* p = NULL;
  EXPECT_THROW(helib::assertNotNull(p, "Expected not null p"), helib::LogicError);
}

TEST(Test_error_handling, helib_assertNotNull_throws_non_standard_helib_exception) {
  int *p = nullptr;
  EXPECT_THROW(helib::assertNotNull<helib::RuntimeError>(p, "Expected not null p"), helib::RuntimeError);
}
  
// assertInRange
TEST(Test_error_handling, helib_assertInRange_no_throws_OutOfBoundError_when_in_range_int)
{
  int element = 15;
  int min = 10;
  int max = 50;
  EXPECT_NO_THROW(helib::assertInRange(element, min, max, "elem not in [min, max)"));
}
  
TEST(Test_error_handling, helib_assertInRange_no_throws_OutOfBoundError_when_in_left_range)
{
  int element = 10;
  int min = 10;
  int max = 50;
  EXPECT_NO_THROW(helib::assertInRange(element, min, max, "elem not in [min, max)"));
}
  
TEST(Test_error_handling, helib_assertInRange_throws_OutOfBoundError_when_less_than_min)
{
  int element = 5;
  int min = 10;
  int max = 50;
  EXPECT_THROW(helib::assertInRange(element, min, max, "elem not in [min, max)"), helib::OutOfRangeError);
}
  
TEST(Test_error_handling, helib_assertInRange_throws_OutOfBoundError_when_greater_than_max)
{
  int element = 100;
  int min = 10;
  int max = 50;
  EXPECT_THROW(helib::assertInRange(element, min, max, "elem not in [min, max)"), helib::OutOfRangeError);
}

TEST(Test_error_handling, helib_assertInRange_throws_OutOfBoundError_when_equals_max)
{
  int element = 50;
  int min = 10;
  int max = 50;
  EXPECT_THROW(helib::assertInRange(element, min, max, "elem not in [min, max)"), helib::OutOfRangeError);
}
  
TEST(Test_error_handling, helib_assertInRange_no_throws_OutOfBoundError_when_in_range_int_and_right_inclusive)
{
  int element = 15;
  int min = 10;
  int max = 50;
  EXPECT_NO_THROW(helib::assertInRange(element, min, max, "elem not in [min, max)", true));
}

TEST(Test_error_handling, helib_assertInRange_no_throws_OutOfBoundError_when_in_left_range_and_right_inclusive)
{
  int element = 10;
  int min = 10;
  int max = 50;
  EXPECT_NO_THROW(helib::assertInRange(element, min, max, "elem not in [min, max]", true));
}

TEST(Test_error_handling, helib_assertInRange_throws_OutOfBoundError_when_less_than_min_and_right_inclusive)
{
  int element = 5;
  int min = 10;
  int max = 50;
  EXPECT_THROW(helib::assertInRange(element, min, max, "elem not in [min, max]", true), helib::OutOfRangeError);
}

TEST(Test_error_handling, helib_assertInRange_throws_OutOfBoundError_when_greater_than_max_and_right_inclusive)
{
  int element = 100;
  int min = 10;
  int max = 50;
  EXPECT_THROW(helib::assertInRange(element, min, max, "elem not in [min, max]", true), helib::OutOfRangeError);
}
  
TEST(Test_error_handling, helib_assertInRange_no_throws_OutOfBoundError_when_equals_max_and_right_inclusive)
{
  int element = 50;
  int min = 10;
  int max = 50;
  EXPECT_NO_THROW(helib::assertInRange(element, min, max, "elem not in [min, max]", true));
}
  
TEST(Test_error_handling, helib_assertInRange_throws_non_default_helib_RuntimeError)
{
  int element = 100;
  int min = 10;
  int max = 50;
  EXPECT_THROW(helib::assertInRange<helib::RuntimeError>(element, min, max, "elem not in [min, max]", true), helib::RuntimeError);
}

} // namespace
