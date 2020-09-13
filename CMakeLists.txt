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

project(helib_superbuild LANGUAGES C CXX)

# Fail if the target architecture is not 64-bit.
if (NOT (CMAKE_SIZEOF_VOID_P EQUAL 8))
  message(FATAL_ERROR "HElib requires a 64-bit architecture.")
endif ()

# Define standard installation directories (GNU)
include(GNUInstallDirs)

# Use -std=c++14 as default.
set(CMAKE_CXX_STANDARD 14)
# Disable C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
# Require full C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# GMP minimal version to be used if not PACKAGE_BUILD
set(GMP_MINIMAL_VERSION "6.0.0")
# NTL minimal version to be used if NOT PACKAGE_BUILD
set(NTL_MINIMAL_VERSION "11.0.0")

# The -isysroot flag needs set.
if (APPLE OR (CMAKE_CXX_PLATFORM_ID STREQUAL "Darwin"))
  set(macos_isysroot_flag "-isysroot ${CMAKE_OSX_SYSROOT}")
else()
  set(macos_isysroot_flag "")
endif()

# Setting up RelWithDebInfo as default CMAKE_BUILD_TYPE
if (NOT CMAKE_BUILD_TYPE)
  # Setting RelWithDebInfo as it will compile with -O2 -g
  set(CMAKE_BUILD_TYPE RelWithDebInfo
      CACHE
      STRING "Choose the type of build, options are: Debug RelWithDebInfo Release MinSizeRel"
      FORCE)

endif (NOT CMAKE_BUILD_TYPE)

# Adding possible gui values to CMAKE_BUILD_TYPE variable
set_property(CACHE
             CMAKE_BUILD_TYPE
             PROPERTY
             STRINGS "Debug" "RelWithDebInfo" "Release" "MinSizeRel")

# Path containing FindGMP.cmake and FindNTL.cmake
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

option(BUILD_SHARED "Build as shared library" OFF)

option(PACKAGE_BUILD
       "Download dependencies and build as a self-contained package"
       OFF)

option(ENABLE_THREADS
       "Enable threads support.  Requires NTL built with NTL_THREADS=on"
       ON)

option(HELIB_DEBUG
       "Build with HELIB_DEBUG (enables extra debugging info, but needs to be initialized)"
       OFF)
option(ENABLE_TEST "Enable tests" OFF)
option(PEDANTIC_BUILD "Use -Wall -Wpedantic -Wextra -Werror during build" ON)

# Add properties dependent to PACKAGE_BUILD
if (PACKAGE_BUILD)
  # Properties to be enabled when using PACKAGE_BUILD
  set(PACKAGE_DIR
      ""
      CACHE
      STRING "Folder with the compilation output directory (DEFAULT: helib_pack)")

  option(FETCH_GMP "Fetch and compile the GMP library (DEFAULT: ON)" ON)
else (PACKAGE_BUILD)
  # Properties to be enabled when not using PACKAGE_BUILD
  option(ENABLE_LEGACY_TEST
         "Build the legacy test files (does not work with PACKAGE_BUILD)"
         OFF)
endif (PACKAGE_BUILD)

# Allow GMP_DIR and search in it if not using PACKAGE_BUILD or using
# PACKAGE_BUILD and !FETCH_GMP
if ((NOT PACKAGE_BUILD) OR NOT FETCH_GMP)
  set(GMP_DIR
      ""
      CACHE
      STRING "Prefix of the GMP library (ignored when PACKAGE_BUILD is ON)")
endif ((NOT PACKAGE_BUILD) OR NOT FETCH_GMP)

# Allow NTL_DIR only if not PACKAGE_BUILD
if (NOT PACKAGE_BUILD)
  set(NTL_DIR
      ""
      CACHE
      STRING "Prefix of the NTL library (ignored when PACKAGE_BUILD is ON)")
endif (NOT PACKAGE_BUILD)

# Setting flag lists to avoid polluting CMAKE_CXX_FLAGS.
# PUBLIC_HELIB_CXX_FLAGS will be exported to the installed target.
set(PRIVATE_HELIB_CXX_FLAGS "")
set(PUBLIC_HELIB_CXX_FLAGS "")
# Add extra checks during build
if (PEDANTIC_BUILD)
  list(APPEND PRIVATE_HELIB_CXX_FLAGS "-Wall" "-Wpedantic" "-Wextra" "-Werror")
endif (PEDANTIC_BUILD)

if (ENABLE_TEST)
  enable_testing()
endif (ENABLE_TEST)

# Look for pthread using default FindThreads.cmake script
find_package(Threads)
if (ENABLE_THREADS AND Threads_NOTFOUND)
  message(FATAL_ERROR "Cannot find pthreads (ENABLE_THREADS is ON).")
endif ()

# NOTE: Consider reconfiguring everything when PACKAGE_BUILD changes value.
# Options from the previous value will remain otherwise.
# Set up extra properties depending on the value of PACKAGE_BUILD
if (PACKAGE_BUILD)
  # Setting up dependencies versions
  # GMP version to be used (and eventually downloaded) if PACKAGE_BUILD
  set(FETCHED_GMP_VERSION "6.2.0")
  # NTL version to be used (and eventually downloaded) if PACKAGE_BUILD
  set(FETCHED_NTL_VERSION "11.4.3")
  # Setting up default compilation output directory
  if (NOT PACKAGE_DIR)
    set(PACKAGE_DIR "helib_pack")
  endif (NOT PACKAGE_DIR)
  if (NOT IS_ABSOLUTE ${PACKAGE_DIR})
    # Make CMAKE_PACKAGE_DIR absolute
    set(PACKAGE_DIR
        "${CMAKE_BINARY_DIR}/${PACKAGE_DIR}"
        CACHE
        STRING "Folder with the compilation output directory (DEFAULT: helib_pack)"
        FORCE)
  endif ()

  # Setting up download/build path of external dependencies
  set(DEPENDENCIES_FOLDER "${CMAKE_BINARY_DIR}/dependencies")
  set_property(DIRECTORY PROPERTY EP_BASE "${DEPENDENCIES_FOLDER}")

  # Raising warning when PACKAGE_BUILD is ON
  # warn if installing globally
  # NOTE: this is a somewhat fragile check that can be enhanced
  string(FIND "${CMAKE_INSTALL_PREFIX}" "/usr" install_in_usr)
  if ("${install_in_usr}" EQUAL 0)
    message(WARNING
            "CAUTION: Package build should not be installed globally as it will potentially override dependencies.")
  endif ()
  unset(install_in_usr)

  # Warn existing dependencies are ignored and rebuilt
  if (GMP_DIR)
    if (FETCH_GMP)
      message(WARNING "GMP_DIR is ignored when PACKAGE_BUILD is ON.")
    else (FETCH_GMP)
      message(STATUS
              "GMP_DIR is not the system one.  This may prevent relocatability of the package.")
    endif (FETCH_GMP)
  endif (GMP_DIR)

  if (NTL_DIR)
    message(WARNING "NTL_DIR is ignored when PACKAGE_BUILD is ON.")
  endif (NTL_DIR)

  # Add an imported target to propagate the library locations
  add_library(gmp_external SHARED IMPORTED)
  add_library(ntl_external SHARED IMPORTED)

  # RPATH will be empty since ntl, gmp, and helib are all in PACKAGE_DIR/lib
  set(PACKAGE_RPATH "")

  if (NOT FETCH_GMP)
    # find GMP library
    # Try to find the GMP package (using cmake/FindGMP.cmake script)
    # REQUIRED arg make cmake to fail if GMP is not found
    # Checks that at least version GMP_MINIMAL_VERSION is available
    find_package(GMP "${GMP_MINIMAL_VERSION}" REQUIRED)
    set_target_properties(gmp_external
                          PROPERTIES
                          INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_PATHS}"
                          IMPORTED_LOCATION "${GMP_LIBRARIES}")
  endif (NOT FETCH_GMP)

  # Add the external dependencies (done in dependencies/CMakeLists.txt)
  add_subdirectory(dependencies)

  # get gmp and ntl include/link directories
  get_target_property(GMP_INCLUDE_PATHS
                      gmp_external
                      INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(GMP_LIBRARIES gmp_external IMPORTED_LOCATION)
  get_target_property(NTL_INCLUDE_PATHS
                      ntl_external
                      INTERFACE_INCLUDE_DIRECTORIES)
  get_target_property(NTL_LIBRARIES ntl_external IMPORTED_LOCATION)

  # Track if helib requires pthreads as dependency
  set(HELIB_REQUIRES_PTHREADS ${ENABLE_THREADS})
else (PACKAGE_BUILD)
  # find GMP library
  # Try to find the GMP package (using cmake/FindGMP.cmake script)
  # REQUIRED arg make cmake to fail if GMP is not found
  # Checks that at least version GMP_MINIMAL_VERSION is available
  find_package(GMP "${GMP_MINIMAL_VERSION}" REQUIRED)

  # find NTL library
  # Try to find the NTL package (using cmake/FindNTL.cmake script)
  # REQUIRED arg make cmake to fail if NTL is not found
  # Checks that at least version NTL_MINIMAL_VERSION is available
  find_package(NTL "${NTL_MINIMAL_VERSION}" REQUIRED)

  # Setting compiler output directories
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
      ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
      ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
      ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})

  # Thread enabling checks
  # Checking if NTL is built with NTL_THREADS=on
  set(ntl_config_file "${NTL_INCLUDE_PATHS}/NTL/config.h")
  if (EXISTS "${ntl_config_file}")
    include(CheckCXXSymbolExists)
    check_cxx_symbol_exists("NTL_THREADS" "${ntl_config_file}" ntl_has_threads)
  else ()
    message(
      FATAL_ERROR "Cannot locate NTL configuration file (${ntl_config_file}).")
  endif ()

  # Raising errors when threads are misconfigured
  if (ENABLE_THREADS AND NOT ntl_has_threads)
    message(FATAL_ERROR
            "Cannot enable threads since NTL was built without.  Consider re-building NTL with NTL_THREADS=on.")
  endif (ENABLE_THREADS AND NOT ntl_has_threads)

  # This should not really happen
  if (ntl_has_threads AND Threads_NOTFOUND)
    message(FATAL_ERROR "NTL requires pthreads that has not been found.")
  endif (ntl_has_threads AND Threads_NOTFOUND)

  # Track if helib requires pthreads as dependency
  set(HELIB_REQUIRES_PTHREADS ${ntl_has_threads})

  unset(ntl_config_file)
  unset(ntl_has_threads)
endif (PACKAGE_BUILD)

# Building HELIB here
if (PACKAGE_BUILD)
  # Adding HELIB as an external project
  include(ExternalProject)

  # Before building helib_external wait compilation of gmp and ntl.
  list(APPEND helib_external_deps "gmp_external" "ntl_external")

  ExternalProject_Add(
    helib_external
    DEPENDS # await compilation of gmp and ntl
            gmp_external ntl_external
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/src
    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PACKAGE_DIR}
               -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
               -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
               -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD}
               -DCMAKE_CXX_STANDARD_REQUIRED=${CMAKE_CXX_STANDARD_REQUIRED}
               -DCMAKE_CXX_EXTENSIONS=${CMAKE_CXX_EXTENSIONS}
               -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
               -DPRIVATE_HELIB_CXX_FLAGS=${PRIVATE_HELIB_CXX_FLAGS}
               -DPUBLIC_HELIB_CXX_FLAGS=${PUBLIC_HELIB_CXX_FLAGS}
               -DGMP_LIBRARIES=${GMP_LIBRARIES}
               -DNTL_INCLUDE_PATHS=${NTL_INCLUDE_PATHS}
               -DNTL_LIBRARIES=${NTL_LIBRARIES}
               -DENABLE_THREADS=${ENABLE_THREADS}
               -DHELIB_REQUIRES_PTHREADS=${HELIB_REQUIRES_PTHREADS}
               -DBUILD_SHARED=${BUILD_SHARED}
               -DPACKAGE_BUILD=${PACKAGE_BUILD}
               -DFETCH_GMP=${FETCH_GMP}
               -DENABLE_TEST=${ENABLE_TEST}
               -DHELIB_DEBUG=${HELIB_DEBUG}
               -DENABLE_LEGACY_TEST=OFF
    BUILD_ALWAYS ON)

  if (ENABLE_TEST)
    add_test(NAME helib_check
             COMMAND ${CMAKE_MAKE_PROGRAM} test
             WORKING_DIRECTORY "${DEPENDENCIES_FOLDER}/Build/helib_external")
  endif (ENABLE_TEST)

  # To install copy the whole PACKAGE_DIR directory to the defined prefix
  install(
    DIRECTORY ${PACKAGE_DIR}
    DESTINATION
      # this is interpreted relative to the value of CMAKE_INSTALL_PREFIX
      .
    USE_SOURCE_PERMISSIONS)
else (PACKAGE_BUILD)
  # If not compiling as PACKAGE_BUILD then add helib as subfolder and not as an
  # external project
  add_subdirectory(src)
endif (PACKAGE_BUILD)
