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

include(ExternalProject)

set(fetched_gmp_install_prefix "${PACKAGE_DIR}")

# Path of the gmp library file that will be installed (libgmp.so/dylib)
set(fetched_gmp_lib
    "${fetched_gmp_install_prefix}/${CMAKE_INSTALL_LIBDIR}/libgmp${CMAKE_SHARED_LIBRARY_SUFFIX}"
)
# Path of the gmp include directory that will be installed
set(fetched_gmp_includedir
    "${fetched_gmp_install_prefix}/${CMAKE_INSTALL_INCLUDEDIR}")

ExternalProject_Add(
    gmp_fetched
    DEPENDS "setup_install_dir"
    URL "https://gmplib.org/download/gmp/gmp-${FETCHED_GMP_VERSION}.tar.bz2"
    URL_HASH
        # NOTE: hash of version 6.2.0
        SHA256=f51c99cb114deb21a60075ffb494c1a210eb9d7cb729ed042ddb7de9534451ea
    DOWNLOAD_NO_PROGRESS 1
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND
        <SOURCE_DIR>/configure
        "--prefix=${fetched_gmp_install_prefix}"
        "--enable-shared=yes"
        "--disable-static"
        "--with-readline=no"
        # Propagate cmake chosen compiler. 
        # Note must add the macOS isysroot flag here instead of passing
        # CXXFLAGS as that overrides the external project optimize flags.
        "CXX=${CMAKE_CXX_COMPILER} ${macos_isysroot_flag}"
        "CC=${CMAKE_C_COMPILER} ${macos_isysroot_flag}" 
    BUILD_IN_SOURCE ON
    TEST_AFTER_INSTALL OFF
    TEST_EXCLUDE_FROM_MAIN ON
    # TEST_COMMAND
    #     ${CMAKE_MAKE_PROGRAM} check
    # Turn LOG_CONFIGURE ON will log the configure step to a file in the stamp directory
    # LOG_CONFIGURE ON
    # Turn LOG_BUILD ON will log the build step to a file in the stamp directory
    # LOG_BUILD ON
    # No build command specified ==> use the default one
    # BUILD_COMMAND
    # No install command specified ==> build the install target (e.g. make install)
    # INSTALL_COMMAND
    # Turn LOG_INSTALL ON will log the install step to a file in the stamp directory
    # LOG_INSTALL ON
    # Turn LOG_TEST ON will log the test step to a file in the stamp directory
    # LOG_TEST ON
)

if (ENABLE_TEST)
  set(test_path "${DEPENDENCIES_FOLDER}/Source/gmp_fetched")
  add_test(NAME gmp_check
           COMMAND ${CMAKE_MAKE_PROGRAM} check
           WORKING_DIRECTORY ${test_path})
  unset(test_path)
endif (ENABLE_TEST)

# Removed due to rpath delay
# add_dependencies(gmp_external gmp_fetched)

# This is because we are postponing the rpath fix to after both
# the deps are built
add_dependencies(dependencies_built gmp_fetched)

# Add include dirs and library path as properties of gmp_external imported
# target to propagate them above
set_target_properties(
  gmp_external
  PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES
          "${fetched_gmp_install_prefix}/${CMAKE_INSTALL_INCLUDEDIR}"
      IMPORTED_LOCATION "${fetched_gmp_lib}")

# Unset temporary variables
unset(fetched_gmp_install_prefix)
unset(fetched_gmp_byproducts)
unset(fetched_gmp_lib)
unset(fetched_gmp_includedir)
