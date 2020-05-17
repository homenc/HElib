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

# Use cmake standard find_library package
include(FindPackageHandleStandardArgs)

# Standard `lib` and system specific `lib32/64` directory.
list(APPEND lib_suffixes "lib" "${CMAKE_INSTALL_LIBDIR}")

# Standard `include` system specific directory (usually the same).
list(APPEND header_suffixes "include/NTL" "${CMAKE_INSTALL_INCLUDEDIR}/NTL")

# If CMAKE_LIBRARY_ARCHITECTURE is defined (e.g.: `x86_64-linux-gnu`) add that
# after `lib` and `lib32/64` and after `include` suffixes(required for Ubuntu)
if (CMAKE_LIBRARY_ARCHITECTURE)
  list(APPEND lib_suffixes
       "lib/${CMAKE_LIBRARY_ARCHITECTURE}"
       "${CMAKE_INSTALL_LIBDIR}/${CMAKE_LIBRARY_ARCHITECTURE}")
  list(APPEND header_suffixes
       "include/${CMAKE_LIBRARY_ARCHITECTURE}/NTL"
       "${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_LIBRARY_ARCHITECTURE}/NTL")
endif (CMAKE_LIBRARY_ARCHITECTURE)

if (NTL_DIR)
  # If user-specified folders: look there
  find_library(NTL_LIB
               NAMES ntl libntl
               PATHS ${NTL_DIR}
               PATH_SUFFIXES ${lib_suffixes}
               NO_DEFAULT_PATH
               DOC "NTL library")

  find_path(NTL_HEADERS
            NAMES config.h
            PATHS ${NTL_DIR}
            PATH_SUFFIXES ${header_suffixes}
            NO_DEFAULT_PATH
            DOC "NTL headers")

else (NTL_DIR)
  # Else: look in default paths
  find_library(NTL_LIB
               NAMES ntl libntl
               PATH_SUFFIXES "${lib_suffixes}"
               DOC "NTL library")

  find_path(NTL_HEADERS
            NAMES config.h
            PATH_SUFFIXES ${header_suffixes}
            DOC "NTL headers")
endif (NTL_DIR)

unset(lib_suffixes)
unset(header_suffixes)

if (NTL_HEADERS AND NTL_LIB)
  # Find ntl version
  file(STRINGS "${NTL_HEADERS}/version.h" ntl_version_string
       REGEX "NTL_VERSION[ \t]+\"([0-9.]+)\"")
  string(REGEX REPLACE "[^ \t]*[ \t]+NTL_VERSION[ \t]+\"([0-9.]+)\"" "\\1"
                       ntl_version_string "${ntl_version_string}")

  string(REGEX REPLACE "([0-9]+)\.([0-9]+)\.([0-9]+)" "\\1" ntl_major
                       "${ntl_version_string}")
  string(REGEX REPLACE "([0-9]+)\.([0-9]+)\.([0-9]+)" "\\2" ntl_minor
                       "${ntl_version_string}")
  string(REGEX REPLACE "([0-9]+)\.([0-9]+)\.([0-9]+)" "\\3" ntl_patchlevel
                       "${ntl_version_string}")
  if ((ntl_version_string STREQUAL "")
      OR (ntl_major STREQUAL "")
      OR (ntl_minor STREQUAL "")
      OR (ntl_patchlevel STREQUAL ""))
    # If the version encoding is wrong then it is set to "WRONG VERSION ENCODING" causing find_package_handle_standard_args to fail
    set(NTL_VERSION "WRONG VERSION ENCODING")
    set(NTL_VERSION_MAJOR "WRONG VERSION ENCODING")
    set(NTL_VERSION_MINOR "WRONG VERSION ENCODING")
    set(NTL_VERSION_PATCH "WRONG VERSION ENCODING")
  else ()
    set(NTL_VERSION "${ntl_major}.${ntl_minor}.${ntl_patchlevel}")
    set(NTL_VERSION_MAJOR "${ntl_major}")
    set(NTL_VERSION_MINOR "${ntl_minor}")
    set(NTL_VERSION_PATCH "${ntl_patchlevel}")
  endif ()

  unset(ntl_version_string)
  unset(ntl_major)
  unset(ntl_minor)
  unset(ntl_patchlevel)
endif (NTL_HEADERS AND NTL_LIB)

# Raising an error if ntl with required version has not been found
set(fail_msg
    "NTL required dynamic shared library has not been found.  (Try cmake -DNTL_DIR=<NTL-root-path>).")
if (NTL_DIR)
  set(fail_msg
      "NTL required dynamic shared library has not been found in ${NTL_DIR}.")
endif (NTL_DIR)

# Check the library has been found, handle QUIET/REQUIRED arguments and set
# NTL_FOUND accordingly or raise the error
find_package_handle_standard_args(NTL
                                  REQUIRED_VARS NTL_LIB NTL_HEADERS
                                  VERSION_VAR NTL_VERSION
                                  FAIL_MESSAGE "${fail_msg}")

unset(fail_msg)

# If NTL has been found set the default variables
if (NTL_FOUND)
  set(NTL_LIBRARIES "${NTL_LIB}")
  set(NTL_INCLUDE_PATHS "${NTL_HEADERS}/..")
endif (NTL_FOUND)

# Mark NTL_LIBRARIES NTL_INCLUDE_PATHS as advanced variables
mark_as_advanced(NTL_LIBRARIES NTL_INCLUDE_PATHS)
