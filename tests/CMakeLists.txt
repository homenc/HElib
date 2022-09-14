# Copyright (C) 2019-2020 IBM Corp.
#
# This program is Licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License. See accompanying LICENSE file.

# Gtests configuration file

# If ONLY_ADD_TEST the test targets will not get generated, but only cmake tests
# will be added (see the end of this file).
# Otherwise the target will be generated, tests will be built and added to cmake.
if (NOT ONLY_ADD_TEST)
    # TestVersion will be generated in CMAKE_CURRENT_BINARY_DIR.
    configure_file(TestVersion.in.cpp
                ${CMAKE_CURRENT_BINARY_DIR}/TestVersion.cpp
                @ONLY)

    set(GTEST_SRC
        "test_common.cpp"
        "TestArgMap.cpp"
        "TestBGV.cpp"
        "TestBootstrappingWithMultiplications.cpp"
        "TestCKKS.cpp"
        "TestClonedPtr.cpp"
        "TestContext.cpp"
        "TestCtxt.cpp"
        "TestErrorHandling.cpp"
        "TestHEXL.cpp"
        "TestLogging.cpp"
        "TestMatmulCKKS.cpp"
        "TestMatrix.cpp"
        "TestPartialMatch.cpp"
        "TestPermutations.cpp"
        "TestPolyMod.cpp"
        "TestPolyModRing.cpp"
        "TestPtxt.cpp"
        "TestQuery.cpp"
        "TestSet.cpp"
        "TestBinIO.cpp"
        "TestIO.cpp"
        "${CMAKE_CURRENT_BINARY_DIR}/TestVersion.cpp" # TestVersion.cpp is auto-generated in CMAKE_CURRENT_BINARY_DIR
        )

    set(PORTED_LEGACY_TEST_SRC
        "GTestApproxNums.cpp"
        "GTestBinaryArith.cpp"
        "GTestBinaryCompare.cpp"
        "GTestBootstrapping.cpp"
        "GTestEaCx.cpp"
        "GTestEvalMap.cpp"
        "GTestExtractDigits.cpp"
        "GTestFatboot.cpp"
        "GTestGeneral.cpp"
        "GTestIntraSlot.cpp"
        "GTestMatmul.cpp"
        "GTestPAlgebra.cpp"
        "GTestPGFFT.cpp"
        "GTestPolyEval.cpp"
        "GTestPowerful.cpp"
        "GTestPtrVector.cpp"
        "GTestReplicate.cpp"
        "GTestTableLookup.cpp"
        "GTestThinboot.cpp"
        "GTestThinBootstrapping.cpp"
        "GTestThinEvalMap.cpp")

    # Adding the test binary target
    add_executable(runTests main.cpp ${PORTED_LEGACY_TEST_SRC} ${GTEST_SRC})

    target_link_libraries(runTests gtest)
    target_link_libraries(runTests helib)


    # NOTE: Adding all PRIVATE_HELIB_CXX_FLAGS to the tests will add -fPIC too.
    # At the moment this does not affect compilation of the tests and does not
    # change their running time.
    target_compile_options(runTests PRIVATE ${PRIVATE_HELIB_CXX_FLAGS})

    target_include_directories(
    runTests
    PRIVATE
           # Adding CMAKE_CURRENT_SOURCE_DIR to include files to access
           # test_common.h.
           ${CMAKE_CURRENT_SOURCE_DIR}
           # Adding HELIB_SOURCE_DIR to include local/private headers e.g. io.h
           ${HELIB_SOURCE_DIR}
           "$<BUILD_INTERFACE:${HELIB_DEPENDENCIES_DIR}/json>")

    # Dealing with test resource files (Test_Bin_IO and Test_IO)
    list(APPEND TEST_RESOURCE_FILES "test_resources/iotest_asciiBE.txt"
                                    "test_resources/iotest_asciiLE.txt"
                                    "test_resources/iotest_binBE.bin"
                                    "test_resources/iotest_binLE.bin")

    set(test_resources_dir "${CMAKE_BINARY_DIR}/test_resources/")

    # NOTE: Consider changing add_custom_command with add_custom_target to redo
    # the command if someone deletes the resources without rebuilding the tests
    add_custom_command(TARGET runTests
                    POST_BUILD
                    COMMAND ${CMAKE_COMMAND}
                    ARGS -E make_directory ${test_resources_dir}
                    COMMENT "Copying test resource files")

    foreach (RESOURCE_FILE ${TEST_RESOURCE_FILES})
    add_custom_command(
        TARGET runTests
        POST_BUILD
        # OUTPUT "${CMAKE_BINARY_DIR}/${RESOURCE_FILE}"
        COMMAND ${CMAKE_COMMAND}
        ARGS -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/${RESOURCE_FILE}" "${test_resources_dir}"
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${RESOURCE_FILE}")
    endforeach (RESOURCE_FILE)
endif (NOT ONLY_ADD_TEST)

# Even if not adding test targets we add all tests to cmake tests.
set(TEST_NAMES
    "GTestApproxNums"
    "GTestBinaryArith"
    "GTestBinaryCompare"
    "GTestBootstrapping"
    "GTestEaCx"
    "GTestEvalMap"
    "GTestExtractDigits"
    "GTestFatboot"
    "GTestGeneral"
    "GTestIntraSlot"
    "GTestMatmul"
    "GTestPAlgebra"
    "GTestPGFFT"
    "GTestPolyEval"
    "GTestPowerful"
    "GTestPtrVector"
    "GTestReplicate"
    "GTestTableLookup"
    "GTestThinboot"
    "GTestThinBootstrapping"
    "GTestThinEvalMap"
    "TestArgMap"
    "TestBGV"
    "TestCKKS"
    "TestClonedPtr"
    "TestContext"
    "TestCtxt"
    "TestErrorHandling"
    "TestFatBootstrappingWithMultiplications"
    "TestHEXL"
    "TestLogging"
    "TestMatmulCKKS"
    "TestMatrix"
    "TestPartialMatch"
    "TestPermutations"
    "TestPolyMod"
    "TestPolyModRing"
    "TestPtxt"
    "TestQuery"
    "TestSet"
    "TestThinBootstrappingWithMultiplications"
    "TestBinIO"
    "TestIO"
    "TestVersion"
    )
  message(STATUS "${CMAKE_BINARY_DIR}")
foreach (TEST_NAME ${TEST_NAMES})
    # Tests will be run in ${HELIB_TEST_BIN_DIR} that is the test folder if
    # library build, ${DEPENDENCIES_FOLDER}/Build/helib_external/${HELIB_TEST_BIN_DIR}
    # otherwise.
    # Test output in the XUnit format will be generated in
    # ${CMAKE_BINARY_DIR}/xunit_test_result/TestResult-${TEST_NAME}.xml
    add_test(NAME "${TEST_NAME}"
      COMMAND runTests --gtest_output=xml:${CMAKE_BINARY_DIR}/xunit_test_result/TestResult-${TEST_NAME}.xml --gtest_filter=*${TEST_NAME}*
             WORKING_DIRECTORY ${HELIB_TEST_BIN_DIR})
endforeach (TEST_NAME)
