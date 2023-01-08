# Copyright (C) 2022 Intel Corporation
# SPDX-License-Identifier: Apache-2.0

add_executable(TestIO main.cpp TestIO.cpp)

target_include_directories(TestIO PRIVATE "../io"
  "../../../utils/common"
  # "${CMAKE_CURRENT_SOURCE_DIR}/../../../build/src/googletest-src/googletest/include")
  "${CMAKE_CURRENT_SOURCE_DIR}/../../../build/dependencies/Build/helib_external/googletest-src/googletest/include")
target_link_libraries(TestIO PRIVATE
  # "${CMAKE_CURRENT_SOURCE_DIR}/../../../build/lib/libgtest.a")
  "${CMAKE_CURRENT_SOURCE_DIR}/../../../build/dependencies/Build/helib_external/lib/libgtest.a")

target_link_libraries(TestIO PRIVATE helib)
