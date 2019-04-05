# Google test configuration file

set(GTEST_SRC
    "test_common.cpp"
    "Test_errorHandling.cpp"
    )

set(PORTED_LEGACY_TEST_SRC
    "GTest_Bin_IO.cpp"
    "GTest_EaCx.cpp"
    "GTest_EvalMap.cpp"
    "GTest_General.cpp"
    "GTest_IO.cpp"
    "GTest_PAlgebra.cpp"
    "GTest_Permutations.cpp"
    "GTest_PolyEval.cpp"
    "GTest_Powerful.cpp"
    "GTest_PtrVector.cpp"
    "GTest_Replicate.cpp"
    "GTest_ThinBootstrapping.cpp"
    "GTest_ThinEvalMap.cpp"
    "GTest_approxNums.cpp"
    "GTest_binaryArith.cpp"
    "GTest_binaryCompare.cpp"
    "GTest_bootstrapping.cpp"
    "GTest_extractDigits.cpp"
    "GTest_intraSlot.cpp"
    "GTest_matmul.cpp"
    "GTest_tableLookup.cpp"
    )

# Adding the test binary target
add_executable(runTests main.cpp ${PORTED_LEGACY_TEST_SRC} ${GTEST_SRC})

target_link_libraries(runTests gtest)
target_link_libraries(runTests helib)

set(TEST_NAMES
    "GTest_Bin_IO"
    "GTest_EaCx"
    "GTest_EvalMap"
    "GTest_General"
    "GTest_IO"
    "GTest_PAlgebra"
    "GTest_Permutations"
    "GTest_PolyEval"
    "GTest_Powerful"
    "GTest_PtrVector"
    "GTest_Replicate"
    "GTest_ThinBootstrapping"
    "GTest_ThinEvalMap"
    "GTest_approxNums"
    "GTest_binaryArith"
    "GTest_binaryCompare"
    "GTest_bootstrapping"
    "GTest_extractDigits"
    "GTest_intraSlot"
    "GTest_matmul"
    "GTest_tableLookup"
    "Test_errorHandling"
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
