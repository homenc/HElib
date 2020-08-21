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
              
## Use -std=c++14 as default.
set(CMAKE_CXX_STANDARD 14)
## Disable C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
## Require full C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(helib_benchmarks
       	LANGUAGES CXX)

# Define standard installation directories (GNU)
include(GNUInstallDirs)

# Setting compiler output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${CMAKE_INSTALL_BINDIR})


# Setting up RelWithDebInfo as default CMAKE_BUILD_TYPE
if (NOT CMAKE_BUILD_TYPE)
  # Setting RelWithDebInfo as it will compile with -02 -g
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: Debug RelWithDebInfo Release MinSizeRel" FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# Adding possible gui values to CMAKE_BUILD_TYPE variable
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "RelWithDebInfo" "Release" "MinSizeRel")

find_package(helib 1.0.0 EXACT REQUIRED)
find_package(benchmark REQUIRED)

# Add each benchmark source file to the executable list
add_executable(helib_benchmark bench_thinboot.cpp bench_fatboot.cpp)
target_link_libraries(helib_benchmark benchmark::benchmark_main helib)
