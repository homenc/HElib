# Google test configuration file

set(GTEST_SRC
    "test_common.cpp"
    "TestErrorHandling.cpp"
    "TestCKKS.cpp"
    "TestArgMap.cpp"
    "TestBootstrappingWithMultiplications.cpp"
    "TestContext.cpp"
    "TestCtxt.cpp"
    "TestPolyMod.cpp"
    "TestPtxt.cpp"
    "TestPolyModRing.cpp")

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
    "GTestThinEvalMap.cpp"
    )

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
    "TestThinBootstrappingWithMultiplications"
    "TestPolyMod"
    "TestPolyModRing"
    "TestPtxt"
    )

foreach(TEST_NAME ${TEST_NAMES})
  add_test(NAME "${TEST_NAME}" COMMAND runTests --gtest_filter=*${TEST_NAME}*)
endforeach(TEST_NAME)

# Dealing with test resource files (Test_Bin_IO and Test_IO)
list(APPEND TEST_RESOURCE_FILES
    "test_resources/iotest_asciiBE.txt"
    "test_resources/iotest_asciiLE.txt"
    "test_resources/iotest_binBE.bin"
    "test_resources/iotest_binLE.bin")

set(test_resources_dir "${CMAKE_BINARY_DIR}/test_resources/")

# NOTE: Consider changing add_custom_command with add_custom_target to redo the command if someone deletes the resources without rebuilding the tests
add_custom_command(
    TARGET runTests
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} ARGS -E make_directory ${test_resources_dir}
    COMMENT "Copying test resource files"
)

foreach(RESOURCE_FILE ${TEST_RESOURCE_FILES})
    add_custom_command(
      TARGET runTests POST_BUILD
      # OUTPUT "${CMAKE_BINARY_DIR}/${RESOURCE_FILE}"
      COMMAND ${CMAKE_COMMAND} ARGS -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/${RESOURCE_FILE}"
        "${test_resources_dir}"
      DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${RESOURCE_FILE}"
    )
endforeach(RESOURCE_FILE)
