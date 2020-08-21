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

set(GTEST_SRC
    "test_common.cpp"
    "TestArgMap.cpp"
    "TestBootstrappingWithMultiplications.cpp"
    "TestCKKS.cpp"
    "TestContext.cpp"
    "TestCtxt.cpp"
    "TestErrorHandling.cpp"
    "TestLogging.cpp"
    "TestMatrix.cpp"
    "TestPartialMatch.cpp"
    "TestPolyMod.cpp"
    "TestPolyModRing.cpp"
    "TestPtxt.cpp"
    "TestSet.cpp"
    )

set(PORTED_LEGACY_TEST_SRC
    "GTestApproxNums.cpp"
    "GTestBinaryArith.cpp"
    "GTestBinaryCompare.cpp"
    "GTestBinIO.cpp"
    "GTestBootstrapping.cpp"
    "GTestEaCx.cpp"
    "GTestEvalMap.cpp"
    "GTestExtractDigits.cpp"
    "GTestFatboot.cpp"
    "GTestGeneral.cpp"
    "GTestIntraSlot.cpp"
    "GTestIO.cpp"
    "GTestMatmul.cpp"
    "GTestPAlgebra.cpp"
    "GTestPermutations.cpp"
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

set(TEST_NAMES
    "GTestApproxNums"
    "GTestBinaryArith"
    "GTestBinaryCompare"
    "GTestBinIO"
    "GTestBootstrapping"
    "GTestEaCx"
    "GTestEvalMap"
    "GTestExtractDigits"
    "GTestFatboot"
    "GTestGeneral"
    "GTestIntraSlot"
    "GTestIO"
    "GTestMatmul"
    "GTestPAlgebra"
    "GTestPermutations"
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
    "TestCKKS"
    "TestContext"
    "TestCtxt"
    "TestErrorHandling"
    "TestFatBootstrappingWithMultiplications"
    "TestLogging"
    "TestMatrix"
    "TestPartialMatch"
    "TestPolyMod"
    "TestPolyModRing"
    "TestPtxt"
    "TestSet"
    "TestThinBootstrappingWithMultiplications"
    )

# NOTE: Adding all PRIVATE_HELIB_CXX_FLAGS to the tests will add -fPIC too.
# At the moment this does not affect compilation of the tests and does not
# change their running time.
target_compile_options(runTests PRIVATE ${PRIVATE_HELIB_CXX_FLAGS})

foreach (TEST_NAME ${TEST_NAMES})
  add_test(NAME "${TEST_NAME}" COMMAND runTests --gtest_filter=*${TEST_NAME}*)
endforeach (TEST_NAME)

# Dealing with test resource files (Test_Bin_IO and Test_IO)
list(APPEND TEST_RESOURCE_FILES "test_resources/iotest_asciiBE.txt"
                                "test_resources/iotest_asciiLE.txt"
                                "test_resources/iotest_binBE.bin"
                                "test_resources/iotest_binLE.bin")

set(test_resources_dir "${CMAKE_BINARY_DIR}/test_resources/")

# NOTE: Consider changing add_custom_command with add_custom_target to redo the
# command if someone deletes the resources without rebuilding the tests
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
