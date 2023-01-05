# Copyright (C) 2020 IBM Corp.
# Copyright (C) 2022 Intel Corporation
# SPDX-License-Identifier: Apache-2.0
				 
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/test/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/test/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    "${CMAKE_BINARY_DIR}/test/${CMAKE_INSTALL_BINDIR}")

add_executable(test_bootstrap test_bootstrap.cpp)
add_executable(test_bootstrap_keygen test_bootstrap_keygen.cpp)

target_include_directories(test_bootstrap PRIVATE "../common")
target_include_directories(test_bootstrap_keygen PRIVATE "../common")

target_link_libraries(test_bootstrap helib)
target_link_libraries(test_bootstrap_keygen helib)
