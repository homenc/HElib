# Copyright (C) 2019-2020 IBM Corp.
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

# Use -std=c++14 as default.
set(CMAKE_CXX_STANDARD 14)
# Disable C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
# Require full C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(BGV_binary_arithmetic LANGUAGES CXX)

find_package(helib 1.1.0 EXACT REQUIRED)

add_executable(BGV_binary_arithmetic BGV_binary_arithmetic.cpp)

target_link_libraries(BGV_binary_arithmetic helib)
