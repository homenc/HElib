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

set(fetched_ntl_install_prefix "${PACKAGE_DIR}")

# Setting right string for enable_threads
if (ENABLE_THREADS)
  set(build_threads "on")
else (ENABLE_THREADS)
  set(build_threads "off")
endif (ENABLE_THREADS)

# Wait until the directory structure and links are properly created
set(ntl_superbuild_dependencies "setup_install_dir")

# If requiring GMP to be built then in compilation NTL depends to gmp_fetched
# because gmp_external is dependent on this (through dependencies_built) due to
# the rpath delay
if (FETCH_GMP)
  list(APPEND ntl_superbuild_dependencies "gmp_fetched")
endif (FETCH_GMP)

# get the file path of libgmp.(so|dylib) (set beforehand)
get_target_property(gmp_lib_path gmp_external IMPORTED_LOCATION)
# get the dirname of libgmp library (in <prefix>/lib)
get_filename_component(gmp_lib_dir "${gmp_lib_path}" DIRECTORY)
# remove the lib directory to get the prefix
get_filename_component(gmp_prefix "${gmp_lib_dir}" DIRECTORY)

# NOTE: Workaround for the way NTL builds on Linux.
# Without explicitly pointing the linker at the correct GMP, NTL has a bad habit
# of picking up the system's libgmp. This causes erroneous version 'mismatches'.
if ("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
  set(ntl_make_libpath_tweak_cmd
      "export LD_LIBRARY_PATH=${gmp_lib_dir}:$LD_LIBRARY_PATH")
else ()
  set(ntl_make_libpath_tweak_cmd ":")
endif ()

# Path of the ntl library file that will be installed (libntl.so/dylib)
set(fetched_ntl_lib
    "${fetched_ntl_install_prefix}/${CMAKE_INSTALL_LIBDIR}/libntl${CMAKE_SHARED_LIBRARY_SUFFIX}")

# Path of the ntl include directory that will be installed
set(fetched_ntl_includedir
    "${fetched_ntl_install_prefix}/${CMAKE_INSTALL_INCLUDEDIR}")

ExternalProject_Add(
    ntl_fetched
    DEPENDS "${ntl_superbuild_dependencies}"
    URL "https://www.shoup.net/ntl/ntl-${FETCHED_NTL_VERSION}.tar.gz"
    URL_HASH
        # NOTE: hash of version 11.4.3
        SHA256=b7c1ccdc64840e6a24351eb4a1e68887d29974f03073a1941c906562c0b83ad2
    DOWNLOAD_NO_PROGRESS ON
    UPDATE_COMMAND ""
    # Turn LOG_CONFIGURE ON will log the configure step to a file in the stamp directory
    # LOG_CONFIGURE ON
    # SOURCE_SUBDIR <SOURCE_DIR>/src
    CONFIGURE_COMMAND
        cd src && <SOURCE_DIR>/src/configure
            "PREFIX=${PACKAGE_DIR}"
            "SHARED=on"
            "NTL_GMP_LIP=on"
            "NTL_THREADS=${build_threads}"
            "NTL_THREAD_BOOST=${build_threads}"
            "GMP_PREFIX=${gmp_prefix}"
            # Propagate cmake chosen compiler. 
            # Note must add the macOS isysroot flag here instead of passing
            # CXXFLAGS as that overrides the external project optimize flags.
            "CXX=${CMAKE_CXX_COMPILER} ${macos_isysroot_flag}"
    BUILD_IN_SOURCE ON
    BUILD_COMMAND
        cd src &&
        # Adding ${gmp_lib_dir} to LD_LIBRARY_PATH (on Linux) to run program
        # testing GMP during make phase (see note at line 41)
        eval "${ntl_make_libpath_tweak_cmd}" &&
        # Using $(MAKE) set this build as sub-task, so all make
        # arguments (including -j) are forwarded
        $(MAKE)
    INSTALL_COMMAND cd src && $(MAKE) install
    TEST_AFTER_INSTALL OFF
    TEST_EXCLUDE_FROM_MAIN ON
    # Turn LOG_BUILD ON will log the build step to a file in the stamp directory
    # LOG_BUILD ON
    # Turn LOG_INSTALL ON will log the install step to a file in the stamp directory
    # LOG_INSTALL ON
    # Turn LOG_TEST ON will log the test step to a file in the stamp directory
    # LOG_TEST ON
)

if (ENABLE_TEST)
  set(test_path "${DEPENDENCIES_FOLDER}/Source/ntl_fetched/src")
  add_test(NAME ntl_check
           COMMAND ${CMAKE_MAKE_PROGRAM} check
           WORKING_DIRECTORY ${test_path})
  unset(test_path)
endif (ENABLE_TEST)

# Removing static library libntl.a.
set(ntl_static_lib
    "${fetched_ntl_install_prefix}/${CMAKE_INSTALL_LIBDIR}/libntl${CMAKE_STATIC_LIBRARY_SUFFIX}")
add_custom_command(TARGET ntl_fetched
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E remove -f "${ntl_static_lib}")
# Removing also the autotools generated libntl.la to be consistent (it contains references to libntl.a)
set(ntl_static_autotools_conf
    "${fetched_ntl_install_prefix}/${CMAKE_INSTALL_LIBDIR}/libntl.la")
add_custom_command(
    TARGET ntl_fetched
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E remove -f "${ntl_static_autotools_conf}")
unset(ntl_static_lib)
unset(ntl_static_autotools_conf)

# Removed due to rpath delay
# add_dependencies(ntl_external ntl_fetched)

# This is because we are postponing the rpath fix to after both the deps are
# built
add_dependencies(dependencies_built ntl_fetched)

# Add include dirs and library path as properties of ntl_external imported
# target to propagate them above
set_target_properties(
  ntl_external
  PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES
          "${fetched_ntl_install_prefix}/${CMAKE_INSTALL_INCLUDEDIR}"
      IMPORTED_LOCATION "${fetched_ntl_lib}")

# Unset temporary variables
unset(build_threads)
unset(fetched_ntl_install_prefix)
unset(fetched_ntl_byproducts)
unset(fetched_ntl_lib)
unset(fetched_ntl_includedir)
unset(ntl_make_libpath_tweak_cmd)
unset(ntl_superbuild_dependencies)
unset(gmp_prefix)
unset(gmp_lib_dir)
unset(gmp_lib_path)
