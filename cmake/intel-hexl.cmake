# Copyright (C) 2021 Intel Corporation
# SPDX-License-Identifier: Apache-2.0
# This program is Licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License. See accompanying LICENSE file.

include(FetchContent)

FetchContent_Declare(
    hexl
    PREFIX hexl
    GIT_REPOSITORY https://github.com/intel/hexl.git
    GIT_TAG v1.0.1
)
FetchContent_GetProperties(hexl) # Get properties used by hexl

if(NOT hexl_POPULATED)
  FetchContent_Populate(hexl)

  set(CMAKE_C_COMPILER ${CMAKE_C_COMPILER} CACHE STRING "" FORCE)
  set(CMAKE_CXX_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING "" FORCE)
  set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/ext_intel_hexl/ CACHE STRING "" FORCE)
  set(HEXL_DEBUG OFF CACHE BOOL "" FORCE)
  set(HEXL_BENCHMARK OFF CACHE BOOL "" FORCE)
  set(HEXL_EXPORT OFF CACHE BOOL "" FORCE)
  set(HEXL_COVERAGE OFF CACHE BOOL "" FORCE)
  set(HEXL_TESTING OFF CACHE BOOL "" FORCE)
  set(HEXL_SHARED_LIB OFF CACHE BOOL "" FORCE)
  set(EXCLUDE_FROM_ALL TRUE)

  add_subdirectory(
    ${hexl_SOURCE_DIR}
    EXCLUDE_FROM_ALL
  )
endif()
