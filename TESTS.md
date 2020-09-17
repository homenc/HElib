# Building and running the new tests

1. To build and run the new tests, one must pass `-DENABLE_TEST=on` to cmake
   (see [INSTALL.md](INSTALL.md)).  The google test library will then be
   downloaded and the correct targets will be defined for cmake.
2. To compile them, `make` as normal from a build directory.
3. The entire set of tests can be run with:
    - `./bin/runTests` or `./dependencies/Build/helib_external/bin/runTests` if
      `PACKAGE_BUILD=ON`.  The parameters `noPrint`, `dry` and `verbose` can be
      specified as before, but they will apply to all tests which are run.
    - `ctest`.  In this case the output will be significantly less verbose and
      therefore less helpful in the case of a failure.  The parameter
      `--output-on-failure` can be added to display the output of each test that
      happened to fail.  It is also possible to run tests in parallel by using
      `-j<number of concurrent tests>`.  Detailed output of the tests (as if run
      directly with `runTests`) can be found in
      `Testing/Temporary/LastTest.log`.
    - `make test`, which will run the tests through `ctest`.  In this case no
      arguments can be passed.  Detailed output of the tests (as if run directly
      with `runTests`) can be found in `Testing/Temporary/LastTest.log`.

## Running a subset of the tests

One may often want to run a specific test or set of tests.  This is achieved
using google test's built-in 'filtering' system.  To run only the `PolyEval`
tests, for example, one simply runs
```
./bin/runTests --gtest_filter='*PolyEval*'
```
To list all the tests, use `./bin/runTests --gtest_list_tests`.  More details
on this filtering can be found in Google's documentation [here][1].

It is also possible to run a specific test using `ctest` by adding the argument
`-R "<name of the test>"`.  For example to run only the `PolyEval` test one
simply runs
```
ctest -R "PolyEval"
```
To list all the tests, use `ctest -N`.  More details on this filtering can be
found in the `ctest` documentation `ctest --help`.

## Changing parameters

The majority of the tests currently residing in `tests` are ported from the old
tests which took several parameters.  This maps reasonably well onto google
test's 'Parameterised tests', more of which can be read about [here][2].

Currently, the parameters to run tests with are written in the tests themselves,
and require modification in the code to run with different parameters.  Some
runtime overriding of these parameters would be possible by placing a function
call in the `INSTANTIATE_TEST_CASE_P` macro, but this has not been done yet.

# Writing new tests

## Why google test?

The google test framework gives us a few of main benefits:

- The presence of 'fixtures', which are classes to contain common set-up and
  teardown code.  These are created by inheriting from `::testing::Test` or
  `::testing::TestWithParam<T>`, and typically have a name such as
  `GTest_extractDigits`.
- A rich set of expectation/assertion macros, such as `EXPECT_EQ`, `EXPECT_GT`,
  and many more which can be read about [here][3].  A test will typically
  contain at least one of these, the result of which determines the success or
  failure of the test.
- A convenient runner which provides utilities such as filtering, repeating,
  shuffling, and convenient, consistent printouts.  This adds significant
  convenience when only some tests are failing, since a summary as printed at
  the end specifying which tests failed, with which parameters they failed,
  etc.

## Summary of how to write a google test

1. Create a new file for the tests, let's say called `GTest_foo.cpp`.
2. Add `GTest_foo.cpp` and `GTest_foo` to the appropriate lists in
   `tests/CMakeLists.txt`.
3. In `GTest_foo.cpp`, be sure to include both `"gtest/gtest.h"` and
   `"test_common.h`.  The latter will give access to the parameters
   `helib_test::dry`, `helib_test::noPrint`, `helib_test::special_bits`,
   `helib_test::verbose`, and `helib_test::seed`.
4. Create an anonymous namespace in which to put everything, so that the global
   namespace for the test runner is not polluted.
5. Create a class called `GTest_foo`, inheriting from `::testing::Test` or
   `::testing::TestWithParam<some_type>` for non-parameterized or parameterized
   tests (respectively). Common setup and teardown code/data should be placed in
   this class.
6. Create actual test code inside a `TEST_F` or `TEST_P` for non-parameterized
   or parameterized tests (respectively).
7. If using parameterized tests, create some parameters to run the tests with
   using the `INSTANTIATE_TEST_CASE_P` macro.

## Further reading

More details on how to use these google test utilities can be found by the
helpful google test documentation:

- Basic overview: https://github.com/google/googletest/
- Introduction to using google test:
  https://github.com/google/googletest/blob/master/googletest/docs/primer.md
- More detailed documentation:
  https://github.com/google/googletest/blob/master/googletest/README.md
- Advanced features:
  https://github.com/google/googletest/blob/master/googletest/docs/advanced.md

It may also be helpful to look at the existing examples of the above style of
testing, which can be found in the `src/tests` folder.

[1]: https://github.com/google/googletest/blob/master/googletest/docs/advanced.md#running-a-subset-of-the-tests "Google documentation"
[2]: https://github.com/google/googletest/blob/master/googletest/docs/advanced.md#value-parameterized-tests "Parameterised tests"
[3]: https://github.com/google/googletest/blob/master/googletest/docs/primer.md "expectation/assertion macros"
