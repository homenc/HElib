# Copyright (C) 2020 IBM Corp.
# Copyright (C) 2022 Intel Corporation
# SPDX-License-Identifier: Apache-2.0

add_executable(key-gen key-gen.cpp)
add_executable(gen-data gen-data.cpp)

target_include_directories(key-gen PRIVATE "../common")
target_include_directories(gen-data PRIVATE "../common")

target_link_libraries(key-gen helib)
target_link_libraries(gen-data helib)
