cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(helib
  VERSION 1.0.0
  LANGUAGES CXX)

# We assume enable threads is passed correctly (already checked)
if (HELIB_REQUIRES_PTHREADS)
  find_package(Threads REQUIRED)
endif(HELIB_REQUIRES_PTHREADS)

# Look for the Armadillo package if requested
if (ENABLE_ARMADILLO)
  find_package(Armadillo "${ARMADILLO_MINIMAL_VERSION}" REQUIRED)
else (ENABLE_ARMADILLO)
  # suppress unused variable warning when not needed
  unset(ARMADILLO_MINIMAL_VERSION)
endif(ENABLE_ARMADILLO)

include(GNUInstallDirs)

# Setting up cmake install directories
if (WIN32 AND NOT CYGWIN)
  set(CMAKE_INSTALL_CMAKEDIR "CMake")
else()
  set(CMAKE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_DATAROOTDIR}/cmake/${PROJECT_NAME}")
endif()

if (ENABLE_TEST)
  # Download and unpack googletest at configure time
  message(STATUS "Setting up googletest framework")
  configure_file(${PROJECT_SOURCE_DIR}/../cmake/gtest.cmake ${CMAKE_CURRENT_BINARY_DIR}/googletest-download/CMakeLists.txt)
  execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download
    OUTPUT_QUIET # Remove quiet to log the output (breaking the ncurses ccmake gui)
  )
  if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif()
  execute_process(COMMAND ${CMAKE_COMMAND} --build .
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download
    OUTPUT_QUIET # Remove quiet to log the output (breaking the ncurses ccmake gui)
  )
  if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif()

  # Add googletest directly to our build. This defines
  # the gtest and gtest_main targets.
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                   ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                   EXCLUDE_FROM_ALL)


  enable_testing()
endif(ENABLE_TEST)

set(HELIB_SRCS
        "BenesNetwork.cpp"
        "binaryArith.cpp"
        "binaryCompare.cpp"
        "binio.cpp"
        "bluestein.cpp"
        "CModulus.cpp"
        "Ctxt.cpp"
        "debugging.cpp"
        "DoubleCRT.cpp"
        "EaCx.cpp"
        "EncryptedArray.cpp"
        "eqtesting.cpp"
        "EvalMap.cpp"
        "extractDigits.cpp"
        "fft.cpp"
        "FHEContext.cpp"
        "FHE.cpp"
        "hypercube.cpp"
        "IndexSet.cpp"
        "intraSlot.cpp"
        "KeySwitching.cpp"
        "matching.cpp"
        "matmul.cpp"
        "norms.cpp"
        "NumbTh.cpp"
        "OptimizePermutations.cpp"
        "PAlgebra.cpp"
        "PermNetwork.cpp"
        "permutations.cpp"
        "polyEval.cpp"
        "powerful.cpp"
        "primeChain.cpp"
        "recryption.cpp"
        "replicate.cpp"
        "sample.cpp"
        "tableLookup.cpp"
        "timing.cpp"
        "zzX.cpp"
        )

set(HELIB_HEADERS
        "binaryArith.h"
        "binaryCompare.h"
        "binio.h"
        "bluestein.h"
        "clonedPtr.h"
        "CModulus.h"
        "CtPtrs.h"
        "Ctxt.h"
        "debugging.h"
        "DoubleCRT.h"
        "EncryptedArray.h"
        "EvalMap.h"
        "FHEContext.h"
        "FHE.h"
        "hypercube.h"
        "IndexMap.h"
        "IndexSet.h"
        "intraSlot.h"
        "matching.h"
        "matmul.h"
        "multicore.h"
        "norms.h"
        "NumbTh.h"
        "PAlgebra.h"
        "permutations.h"
        "polyEval.h"
        "powerful.h"
        "primeChain.h"
        "PtrMatrix.h"
        "PtrVector.h"
        "randomMatrices.h"
        "range.h"
        "recryption.h"
        "replicate.h"
        "sample.h"
        "tableLookup.h"
        "timing.h"
        "zzX.h"
        "assertions.h"
        "exceptions.h"
        )

set(LEGACY_TEST_SRCS
        "Test_approxNums.cpp"
        "Test_binaryArith.cpp"
        "Test_binaryCompare.cpp"
        "Test_Bin_IO.cpp"
        "Test_bootstrapping.cpp"
        "Test_EaCx.cpp"
        "Test_EvalMap.cpp"
        "Test_extractDigits.cpp"
        "Test_General.cpp"
        "Test_intraSlot.cpp"
        "Test_IO.cpp"
        "Test_matmul.cpp"
        "Test_PAlgebra.cpp"
        "Test_Permutations.cpp"
        "Test_PolyEval.cpp"
        "Test_Powerful.cpp"
        "Test_PtrVector.cpp"
        "Test_Replicate.cpp"
        "Test_tableLookup.cpp"
        "Test_ThinBootstrapping.cpp"
        "Test_ThinEvalMap.cpp"
        "Test_Timing.cpp"
        )

# Add helib target as a shared/static library
if (BUILD_SHARED)
  add_library(helib SHARED ${HELIB_SRCS})
else (BUILD_SHARED)
  add_library(helib STATIC ${HELIB_SRCS})
endif (BUILD_SHARED)

# Add private/public flags to helib
target_compile_options(helib
                      PRIVATE
                        ${PRIVATE_HELIB_CXX_FLAGS}
                      PUBLIC
                        ${PUBLIC_HELIB_CXX_FLAGS}
                      )

# Define HELIB_THREADS as public helib symbol
target_compile_definitions(helib
  PUBLIC
    $<$<BOOL:${ENABLE_THREADS}>:HELIB_THREADS>
    # NOTE: this should be in a generated configure.h rather than being public
    # NOTE: Generator IF are available form cmake 3.8. Use it in the future.
    # $<IF:$<BOOL:${ENABLE_ARMADILLO}>,FFT_ARMA,FFT_NATIVE>
    $<$<BOOL:${ENABLE_ARMADILLO}>:FFT_ARMA>
    $<$<NOT:$<BOOL:${ENABLE_ARMADILLO}>>:FFT_NATIVE>
)

if (PACKAGE_BUILD)
  # If having a package build export paths as relative to the package root
  file(RELATIVE_PATH NTL_INCLUDE_EXPORTED_PATH "${CMAKE_INSTALL_PREFIX}" "${NTL_INCLUDE_PATHS}")
  file(RELATIVE_PATH NTL_LIBRARIES_EXPORTED_PATH "${CMAKE_INSTALL_PREFIX}" "${NTL_LIBRARIES}")

  # Add escaped import prefix to exported relative paths so it will be copied to the autogenerated target file
  # NOTE: Cmake version 3.13.2 uses variable _IMPORT_PREFIX in the autogenerated target file. If it changes then this will break
  set(NTL_LIBRARIES_EXPORTED_PATH "\${_IMPORT_PREFIX}/${NTL_LIBRARIES_EXPORTED_PATH}")

  # Do the same for GMP iff GMP is not the system one.  If the GMP that we use is on the system,
  # we need to use an absolute path instead.
  if (FETCH_GMP)
    file(RELATIVE_PATH GMP_LIBRARIES_EXPORTED_PATH "${CMAKE_INSTALL_PREFIX}" "${GMP_LIBRARIES}")
    set(GMP_LIBRARIES_EXPORTED_PATH "\${_IMPORT_PREFIX}/${GMP_LIBRARIES_EXPORTED_PATH}")
  else (FETCH_GMP)
    set(GMP_LIBRARIES_EXPORTED_PATH "${GMP_LIBRARIES}")
  endif(FETCH_GMP)
else (PACKAGE_BUILD)
  set(NTL_INCLUDE_EXPORTED_PATH "${NTL_INCLUDE_PATHS}")
  set(NTL_LIBRARIES_EXPORTED_PATH "${NTL_LIBRARIES}")
  set(GMP_LIBRARIES_EXPORTED_PATH "${GMP_LIBRARIES}")
endif(PACKAGE_BUILD)

target_include_directories(helib
  PRIVATE
    "$<$<BOOL:${ENABLE_ARMADILLO}>:${ARMADILLO_INCLUDE_DIRS}>"
  PUBLIC
    "$<BUILD_INTERFACE:${NTL_INCLUDE_PATHS}>"
    "$<INSTALL_INTERFACE:${NTL_INCLUDE_EXPORTED_PATH}>"

    # Headers used from source/build location:
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    # Headers used from the installed location:
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

# Setting the helib properties
set_target_properties(helib
  PROPERTIES
    OUTPUT_NAME "helib"
    SOVERSION "${PROJECT_VERSION}"
    PUBLIC_HEADER "${HELIB_HEADERS}"
)

# Adding gmp and ntl link file as relative/absolute depending on the build/install interfaces
target_link_libraries(helib
  PUBLIC
    "$<BUILD_INTERFACE:${GMP_LIBRARIES}>"
    "$<INSTALL_INTERFACE:${GMP_LIBRARIES_EXPORTED_PATH}>"

    "$<BUILD_INTERFACE:${NTL_LIBRARIES}>"
    "$<INSTALL_INTERFACE:${NTL_LIBRARIES_EXPORTED_PATH}>"

    # Add pthread if required
    $<$<BOOL:${HELIB_REQUIRES_PTHREADS}>:Threads::Threads>
    # Add Armadillo if enabled
    $<$<BOOL:${ENABLE_ARMADILLO}>:${ARMADILLO_LIBRARIES}>
  )


# FIX RPATH
if (BUILD_SHARED)
  # RPATH will be empty since ntl, gmp, and helib are all in PACKAGE_DIR/lib
  set(rel_path "")

  if (APPLE)
    set(helib_rpath "@loader_path/${rel_path}")
  else (APPLE)
    set(helib_rpath "\$ORIGIN/${rel_path}")
  endif(APPLE)

  set_target_properties(helib
    PROPERTIES
      MACOS_RPATH ON
      INSTALL_RPATH "${helib_rpath}"
      INSTALL_RPATH_USE_LINK_PATH ON
      SKIP_BUILD_RPATH OFF
      BUILD_WITH_INSTALL_RPATH OFF
  )
endif(BUILD_SHARED)


if (BUILD_AES)
  add_subdirectory(aes)
endif(BUILD_AES)

if (ENABLE_LEGACY_TEST)
    foreach(TEST_FILE ${LEGACY_TEST_SRCS})
      string(REPLACE "\.cpp" "" TEST_NAME "${TEST_FILE}")
      add_executable(${TEST_NAME} ${TEST_FILE})
      target_compile_options(${TEST_NAME}
                               PRIVATE
                                 ${PRIVATE_HELIB_CXX_FLAGS}
                              )
      target_link_libraries(${TEST_NAME} helib)
    endforeach(TEST_FILE)
endif(ENABLE_LEGACY_TEST)

if (ENABLE_TEST)
  add_subdirectory(tests)
endif(ENABLE_TEST)

# Installing helib
# NOTE: add different lib/dev components
install(TARGETS
          helib
        EXPORT
          helibTargets
        ARCHIVE
          DESTINATION ${CMAKE_INSTALL_LIBDIR}
          COMPONENT lib
        LIBRARY
          DESTINATION ${CMAKE_INSTALL_LIBDIR}
          COMPONENT lib
        PUBLIC_HEADER
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/helib
          COMPONENT lib
      )

# Installing exported targets
install(
  EXPORT
    helibTargets
  # NAMESPACE "" #NOTE: add namespace
  DESTINATION
    ${CMAKE_INSTALL_CMAKEDIR}
  COMPONENT
    lib
  )

# Add auto configuration generating functions
include(CMakePackageConfigHelpers)

# Generating the version files
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/helibConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY ExactVersion #choices AnyNewerVersion SameMajorVersion ExactVersion
)

# Generating the basic helibConfig.cmake file
configure_package_config_file(
  ${PROJECT_SOURCE_DIR}/../cmake/helibConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/helibConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
)

# Installing the remaining files
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/helibConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/helibConfigVersion.cmake
  DESTINATION
    ${CMAKE_INSTALL_CMAKEDIR}
)
