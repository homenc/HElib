file(GLOB HELIB_SRC_FILES RELATIVE "${CMAKE_SOURCE_DIR}/src" "*.cpp")
file(GLOB TEST_FILES RELATIVE "${CMAKE_SOURCE_DIR}/src" "Test_*.cpp")

## Remove test files from the source file list.
foreach(prg ${TEST_FILES})
	list(FILTER HELIB_SRC_FILES EXCLUDE REGEX ${prg})
endforeach(prg)

if (BUILD_SHARED)
add_library(fhe SHARED ${HELIB_SRC_FILES})
else (BUILD_SHARED)
add_library(fhe STATIC ${HELIB_SRC_FILES})
endif (BUILD_SHARED)

if (BUILD_AES)
include_directories("${CMAKE_SOURCE_DIR}/src")
add_subdirectory(aes)
endif(BUILD_AES)

if (ENABLE_TEST)
    foreach(TEST_FILE ${TEST_FILES})
        add_executable(${TEST_FILE}_exe ${TEST_FILE})
        target_link_libraries(${TEST_FILE}_exe ${RUN_LIB})
    endforeach(TEST_FILE)
endif(ENABLE_TEST)
