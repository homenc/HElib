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

cmake_minimum_required(VERSION 3.10.2 FATAL_ERROR)

project(helib
        VERSION 1.1.0
        LANGUAGES CXX)

# Location of the header files
set(PROJECT_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/../include")
# Location of the google tests
set(PROJECT_TESTS_DIR "${PROJECT_SOURCE_DIR}/../tests")

# We assume enable threads is passed correctly (already checked)
if (HELIB_REQUIRES_PTHREADS)
  find_package(Threads REQUIRED)
endif (HELIB_REQUIRES_PTHREADS)

include(GNUInstallDirs)

# Setting up cmake install directories
if (WIN32 AND NOT CYGWIN)
  set(CMAKE_INSTALL_CMAKEDIR "CMake")
else ()
  set(CMAKE_INSTALL_CMAKEDIR
      "${CMAKE_INSTALL_DATAROOTDIR}/cmake/${PROJECT_NAME}")
endif ()

if (ENABLE_TEST)
  # Download and unpack googletest at configure time
  message(STATUS "Setting up googletest framework")
  configure_file("${PROJECT_SOURCE_DIR}/../cmake/gtest.cmake"
                 "${CMAKE_CURRENT_BINARY_DIR}/googletest-download/CMakeLists.txt")
  execute_process(
      COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download
      OUTPUT_QUIET # Remove to log the output (breaking the ncurses ccmake gui)
  )
  if (result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
  endif ()
  execute_process(
      COMMAND ${CMAKE_COMMAND} --build .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download
      OUTPUT_QUIET # Remove to log the output (breaking the ncurses ccmake gui)
  )
  if (result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
  endif ()

  # Add googletest directly to our build.  This defines
  # the gtest and gtest_main targets.
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                   ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                   EXCLUDE_FROM_ALL)

  enable_testing()
endif (ENABLE_TEST)

set(HELIB_SRCS
    "BenesNetwork.cpp"
    "binaryArith.cpp"
    "binaryCompare.cpp"
    "binio.cpp"
    "bluestein.cpp"
    "CModulus.cpp"
    "Context.cpp"
    "Ctxt.cpp"
    "debugging.cpp"
    "DoubleCRT.cpp"
    "EaCx.cpp"
    "EncryptedArray.cpp"
    "eqtesting.cpp"
    "EvalMap.cpp"
    "extractDigits.cpp"
    "fhe_stats.cpp"
    "hypercube.cpp"
    "IndexSet.cpp"
    "intraSlot.cpp"
    "keys.cpp"
    "keySwitching.cpp"
    "log.cpp"
    "matching.cpp"
    "matmul.cpp"
    "norms.cpp"
    "NumbTh.cpp"
    "OptimizePermutations.cpp"
    "PAlgebra.cpp"
    "PermNetwork.cpp"
    "permutations.cpp"
    "PGFFT.cpp"
    "polyEval.cpp"
    "PolyMod.cpp"
    "PolyModRing.cpp"
    "powerful.cpp"
    "primeChain.cpp"
    "Ptxt.cpp"
    "randomMatrices.cpp"
    "recryption.cpp"
    "replicate.cpp"
    "sample.cpp"
    "tableLookup.cpp"
    "timing.cpp"
    "zzX.cpp")

set(HELIB_HEADER_DIR "${PROJECT_INCLUDE_DIR}/helib")
set(HELIB_HEADERS
    "${HELIB_HEADER_DIR}/helib.h"
    "${HELIB_HEADER_DIR}/apiAttributes.h"
    "${HELIB_HEADER_DIR}/ArgMap.h"
    "${HELIB_HEADER_DIR}/binaryArith.h"
    "${HELIB_HEADER_DIR}/binaryCompare.h"
    "${HELIB_HEADER_DIR}/binio.h"
    "${HELIB_HEADER_DIR}/bluestein.h"
    "${HELIB_HEADER_DIR}/clonedPtr.h"
    "${HELIB_HEADER_DIR}/CModulus.h"
    "${HELIB_HEADER_DIR}/CtPtrs.h"
    "${HELIB_HEADER_DIR}/Ctxt.h"
    "${HELIB_HEADER_DIR}/debugging.h"
    "${HELIB_HEADER_DIR}/DoubleCRT.h"
    "${HELIB_HEADER_DIR}/EncryptedArray.h"
    "${HELIB_HEADER_DIR}/EvalMap.h"
    "${HELIB_HEADER_DIR}/Context.h"
    "${HELIB_HEADER_DIR}/FHE.h"
    "${HELIB_HEADER_DIR}/keys.h"
    "${HELIB_HEADER_DIR}/keySwitching.h"
    "${HELIB_HEADER_DIR}/log.h"
    "${HELIB_HEADER_DIR}/hypercube.h"
    "${HELIB_HEADER_DIR}/IndexMap.h"
    "${HELIB_HEADER_DIR}/IndexSet.h"
    "${HELIB_HEADER_DIR}/intraSlot.h"
    "${HELIB_HEADER_DIR}/matching.h"
    "${HELIB_HEADER_DIR}/matmul.h"
    "${HELIB_HEADER_DIR}/Matrix.h"
    "${HELIB_HEADER_DIR}/multicore.h"
    "${HELIB_HEADER_DIR}/norms.h"
    "${HELIB_HEADER_DIR}/NumbTh.h"
    "${HELIB_HEADER_DIR}/PAlgebra.h"
    "${HELIB_HEADER_DIR}/partialMatch.h"
    "${HELIB_HEADER_DIR}/permutations.h"
    "${HELIB_HEADER_DIR}/polyEval.h"
    "${HELIB_HEADER_DIR}/PolyMod.h"
    "${HELIB_HEADER_DIR}/PolyModRing.h"
    "${HELIB_HEADER_DIR}/powerful.h"
    "${HELIB_HEADER_DIR}/primeChain.h"
    "${HELIB_HEADER_DIR}/PtrMatrix.h"
    "${HELIB_HEADER_DIR}/PtrVector.h"
    "${HELIB_HEADER_DIR}/Ptxt.h"
    "${HELIB_HEADER_DIR}/randomMatrices.h"
    "${HELIB_HEADER_DIR}/range.h"
    "${HELIB_HEADER_DIR}/recryption.h"
    "${HELIB_HEADER_DIR}/replicate.h"
    "${HELIB_HEADER_DIR}/sample.h"
    "${HELIB_HEADER_DIR}/set.h"
    "${HELIB_HEADER_DIR}/tableLookup.h"
    "${HELIB_HEADER_DIR}/timing.h"
    "${HELIB_HEADER_DIR}/zzX.h"
    "${HELIB_HEADER_DIR}/assertions.h"
    "${HELIB_HEADER_DIR}/exceptions.h"
    "${HELIB_HEADER_DIR}/PGFFT.h"
    "${HELIB_HEADER_DIR}/fhe_stats.h"
    "${HELIB_HEADER_DIR}/zeroValue.h")

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
    "Test_thinboot.cpp"
    "Test_fatboot.cpp"
    "Test_PGFFT.cpp")

# Add helib target as a shared/static library
if (BUILD_SHARED)
  add_library(helib SHARED ${HELIB_SRCS} ${HELIB_HEADERS})
else (BUILD_SHARED)
  add_library(helib STATIC ${HELIB_SRCS} ${HELIB_HEADERS})
endif (BUILD_SHARED)

# Set HElib's properties.  In this case we ask to build HElib with
# POSITION_INDEPENDENT_CODE (-fPIC)
set_target_properties(helib PROPERTIES POSITION_INDEPENDENT_CODE True)

# Add private/public flags to helib
target_compile_options(helib
                       PRIVATE ${PRIVATE_HELIB_CXX_FLAGS}
                       PUBLIC ${PUBLIC_HELIB_CXX_FLAGS})

# Define HELIB_THREADS as public helib symbol
# NOTE: this should be in a generated configure.h rather than being public
target_compile_definitions(helib
                           PUBLIC
                               $<$<BOOL:${ENABLE_THREADS}>:HELIB_THREADS>
                               $<$<BOOL:${ENABLE_THREADS}>:HELIB_BOOT_THREADS>
                               $<$<BOOL:${HELIB_DEBUG}>:HELIB_DEBUG>)

if (PACKAGE_BUILD)
  # If having a package build export paths as relative to the package root
  file(RELATIVE_PATH NTL_INCLUDE_EXPORTED_PATH "${CMAKE_INSTALL_PREFIX}"
                                               "${NTL_INCLUDE_PATHS}")
  file(RELATIVE_PATH NTL_LIBRARIES_EXPORTED_PATH "${CMAKE_INSTALL_PREFIX}"
                                                 "${NTL_LIBRARIES}")

  # Add escaped import prefix to exported relative paths so it will be copied to
  # the autogenerated target file
  # NOTE: Cmake version 3.13.2 uses variable _IMPORT_PREFIX in the autogenerated
  # target file.  If it changes then this will break
  set(NTL_LIBRARIES_EXPORTED_PATH
      "\${_IMPORT_PREFIX}/${NTL_LIBRARIES_EXPORTED_PATH}")

  # Do the same for GMP iff GMP is not the system one.  If the GMP that we use
  # is on the system, we need to use an absolute path instead.
  if (FETCH_GMP)
    file(RELATIVE_PATH
         GMP_LIBRARIES_EXPORTED_PATH "${CMAKE_INSTALL_PREFIX}"
                                     "${GMP_LIBRARIES}")
    set(GMP_LIBRARIES_EXPORTED_PATH
        "\${_IMPORT_PREFIX}/${GMP_LIBRARIES_EXPORTED_PATH}")
  else (FETCH_GMP)
    set(GMP_LIBRARIES_EXPORTED_PATH "${GMP_LIBRARIES}")
  endif (FETCH_GMP)
else (PACKAGE_BUILD)
  set(NTL_INCLUDE_EXPORTED_PATH "${NTL_INCLUDE_PATHS}")
  set(NTL_LIBRARIES_EXPORTED_PATH "${NTL_LIBRARIES}")
  set(GMP_LIBRARIES_EXPORTED_PATH "${GMP_LIBRARIES}")
endif (PACKAGE_BUILD)

target_include_directories(
    helib
    PRIVATE
    PUBLIC # NOTE: The includes must be kept in this order to avoid cmake
           # looking for HElib in /usr/local/include
           # Headers used from source/build location:
           "$<BUILD_INTERFACE:${PROJECT_INCLUDE_DIR}>"
           # Headers used from the installed location:
           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")

target_include_directories(
    helib SYSTEM
    PRIVATE
    PUBLIC
           "$<BUILD_INTERFACE:${NTL_INCLUDE_PATHS}>"
           "$<INSTALL_INTERFACE:${NTL_INCLUDE_EXPORTED_PATH}>")

# Setting the helib properties
set_target_properties(helib
                      PROPERTIES
                          OUTPUT_NAME "helib"
                          SOVERSION "${PROJECT_VERSION}"
                          PUBLIC_HEADER "${HELIB_HEADERS}")

# Adding gmp and ntl link file as relative/absolute depending on the
# build/install interfaces
target_link_libraries(
  helib
  PUBLIC # NTL should be loaded before GMP when NTL is static
         "$<BUILD_INTERFACE:${NTL_LIBRARIES}>"
         "$<INSTALL_INTERFACE:${NTL_LIBRARIES_EXPORTED_PATH}>"
         "$<BUILD_INTERFACE:${GMP_LIBRARIES}>"
         "$<INSTALL_INTERFACE:${GMP_LIBRARIES_EXPORTED_PATH}>"
         # Add pthread if required
         $<$<BOOL:${HELIB_REQUIRES_PTHREADS}>:Threads::Threads>)

# FIX RPATH
if (BUILD_SHARED)
  # RPATH will be empty since ntl, gmp, and helib are all in PACKAGE_DIR/lib
  set(rel_path "")

  if (APPLE)
    set(helib_rpath "@loader_path/${rel_path}")
  else (APPLE)
    set(helib_rpath "\$ORIGIN/${rel_path}")
  endif (APPLE)

  set_target_properties(helib
                        PROPERTIES
                            MACOS_RPATH ON
                            INSTALL_RPATH "${helib_rpath}"
                            INSTALL_RPATH_USE_LINK_PATH ON
                            SKIP_BUILD_RPATH OFF
                            BUILD_WITH_INSTALL_RPATH OFF)
endif (BUILD_SHARED)

if (ENABLE_LEGACY_TEST)
  foreach (TEST_FILE ${LEGACY_TEST_SRCS})
    string(REPLACE "\.cpp" "" TEST_NAME "${TEST_FILE}")
    add_executable(${TEST_NAME} ${TEST_FILE})
    # NOTE: Legacy tests depend only on HElib.
    # HElib's PUBLIC attributes will be added automatically when setting HElib
    # as a dependency.
    target_link_libraries(${TEST_NAME} helib)
  endforeach (TEST_FILE)
endif (ENABLE_LEGACY_TEST)

if (ENABLE_TEST)
  add_subdirectory("${PROJECT_TESTS_DIR}" tests)
endif (ENABLE_TEST)

# Installing helib
# NOTE: add different lib/dev components
install(TARGETS helib
        EXPORT helibTargets
        ARCHIVE       DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT lib
        LIBRARY       DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT lib
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/helib
        COMPONENT lib)

# Installing exported targets
install(EXPORT helibTargets
        # NAMESPACE "" #NOTE: add namespace
        DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
        COMPONENT lib)

# Add auto configuration generating functions
include(CMakePackageConfigHelpers)

# Generating the version files
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/helibConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY
        ExactVersion #choices AnyNewerVersion SameMajorVersion ExactVersion
)

# Generating the basic helibConfig.cmake file
configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/../cmake/helibConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/helibConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_CMAKEDIR})

# Installing the remaining files
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/helibConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/helibConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_CMAKEDIR})
