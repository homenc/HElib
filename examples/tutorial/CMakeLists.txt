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

add_executable(01_ckks_basics 01_ckks_basics.cpp)
add_executable(02_ckks_depth 02_ckks_depth.cpp)
add_executable(03_ckks_data_movement 03_ckks_data_movement.cpp)
add_executable(04_ckks_matmul 04_ckks_matmul.cpp)
add_executable(05_ckks_multlowlvl 05_ckks_multlowlvl.cpp)
add_executable(06_ckks_complex 06_ckks_complex.cpp)
add_executable(07_ckks_serialization 07_ckks_serialization.cpp)
add_executable(08_ckks_deserialization 08_ckks_deserialization.cpp)

target_link_libraries(01_ckks_basics helib)
target_link_libraries(02_ckks_depth helib)
target_link_libraries(03_ckks_data_movement helib)
target_link_libraries(04_ckks_matmul helib)
target_link_libraries(05_ckks_multlowlvl helib)
target_link_libraries(06_ckks_complex helib)
target_link_libraries(07_ckks_serialization helib)
target_link_libraries(08_ckks_deserialization helib)

