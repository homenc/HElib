# Copyright (C) 2020 IBM Corp.
# This program is Licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License. See accompanying LICENSE file.

# Check if the system has fallocate or mkfile (to pre-allocate the ctxt file)
set(alloc_file_name "")
if("${CMAKE_SYSTEM_NAME}" STREQUAL "Linux")
	set(alloc_file_name "fallocate")
elseif (APPLE)
	set(alloc_file_name "mkfile")
else()
	message(FATAL_ERROR "OS not supported.")
endif()
find_program(alloc_file NAMES "${alloc_file_name}")
if(NOT alloc_file)
	message(FATAL_ERROR "Required file allocation command ${alloc_file_name} not found.")
endif()

add_executable(encrypt encrypt.cpp)
add_executable(decrypt decrypt.cpp)

target_include_directories(encrypt PRIVATE "../common")
target_include_directories(decrypt PRIVATE "../common")

target_link_libraries(encrypt helib)
target_link_libraries(decrypt helib)
