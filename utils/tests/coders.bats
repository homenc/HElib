#!/usr/bin/env bats

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

load "std"

function setup {
  mkdir -p $tmp_folder
  cd $tmp_folder
  print-info-location
  check_python36
  check_locations
  create-bgv-toy-params
  createContext BGV "${prefix_bgv}.params" "${prefix_bgv}"
  genData "${prefix_bgv}.dat" 12 "BGV"
  create-ckks-toy-params
  createContext CKKS "${prefix_ckks}.params" "${prefix_ckks}"
  genData "${prefix_ckks}.dat" 12 "CKKS"
}

function teardown {
  cd -
  remove-test-directory "$tmp_folder"
}

@test "BGV: data == decode(encode(data))" {
  run bash -c "$encode ${prefix_bgv}.dat ${prefix_bgv}.info BGV > ${prefix_bgv}.encoded"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_bgv}.encoded ${prefix_bgv}.info BGV > ${prefix_bgv}.decoded"
  assert [ "$status" -eq 0 ]
  diff "${prefix_bgv}.dat" "${prefix_bgv}.decoded"
}

@test "BGV: data != decode(encode(data)) with padding" {
  # Odd number to create padding
  genData "${prefix_bgv}.pad" 11 "BGV"
  run bash -c "$encode ${prefix_bgv}.pad ${prefix_bgv}.info BGV > ${prefix_bgv}.pad.encoded"
  # techo "${output}"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_bgv}.pad.encoded ${prefix_bgv}.info BGV > ${prefix_bgv}.decoded"
  assert [ "$status" -eq 0 ]
  run diff "${prefix_bgv}.pad" "${prefix_bgv}.decoded"
  assert [ "$status" -ne 0 ]
}

@test "BGV: data == decode(encode(data), nelements) with padding" {
  # Odd number to create padding
  nelements=11
  genData "${prefix_bgv}.pad" ${nelements} "BGV"
  run bash -c "$encode ${prefix_bgv}.pad ${prefix_bgv}.info BGV > ${prefix_bgv}.pad.encoded"
  # techo "${output}"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode --nelements ${nelements} ${prefix_bgv}.pad.encoded ${prefix_bgv}.info BGV > ${prefix_bgv}.decoded"
  assert [ "$status" -eq 0 ]
  run diff "${prefix_bgv}.pad" "${prefix_bgv}.decoded"
  assert [ "$status" -eq 0 ]
}

@test "BGV: matrix data == decode(encode(data))" {
  run bash -c "$encode ${prefix_bgv}.dat ${prefix_bgv}.info BGV --dims 2,3 > ${prefix_bgv}.encoded"
  assert [ "$status" -eq 0 ]
  # Number of lines should be 1 for header and 2 * 3 for data.
  assert [ "$(cat ${prefix_bgv}.encoded | wc -l)" -eq 7 ]
  run bash -c "$decode ${prefix_bgv}.encoded ${prefix_bgv}.info BGV --nelements 12 > ${prefix_bgv}.decoded"
  assert [ "$status" -eq 0 ]
  diff "${prefix_bgv}.dat" "${prefix_bgv}.decoded"
}

@test "BGV: matrix encode(data) does not fit" {
  run bash -c "$encode ${prefix_bgv}.dat ${prefix_bgv}.info BGV --dims 1,3 > ${prefix_bgv}.encoded"
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Number of ptxts 4 > capacity of Matrix: 3 (dims: (1, 3))" ]
}

@test "CKKS: data == decode(encode(data))" {
  run bash -c "$encode ${prefix_ckks}.dat ${prefix_ckks}.info CKKS > ${prefix_ckks}.encoded"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_ckks}.encoded ${prefix_ckks}.info CKKS > ${prefix_ckks}.decoded"
  assert [ "$status" -eq 0 ]
  diff "${prefix_ckks}.dat" "${prefix_ckks}.decoded"
}

@test "CKKS: data != decode(encode(data)) with padding" {
  # Odd number to create padding
  genData "${prefix_ckks}.pad" 11 "CKKS"
  run bash -c "$encode ${prefix_ckks}.pad ${prefix_ckks}.info CKKS > ${prefix_ckks}.pad.encoded"
  # techo "${output}"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_ckks}.pad.encoded ${prefix_ckks}.info CKKS > ${prefix_ckks}.decoded"
  assert [ "$status" -eq 0 ]
  run diff "${prefix_ckks}.pad" "${prefix_ckks}.decoded"
  assert [ "$status" -ne 0 ]
}

@test "CKKS: data == decode(encode(data), nelements) with padding" {
  # Odd number to create padding
  nelements=11
  genData "${prefix_ckks}.pad" ${nelements} "CKKS"
  run bash -c "$encode ${prefix_ckks}.pad ${prefix_ckks}.info CKKS > ${prefix_ckks}.pad.encoded"
  # techo "${output}"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode --nelements ${nelements} ${prefix_ckks}.pad.encoded ${prefix_ckks}.info CKKS > ${prefix_ckks}.decoded"
  assert [ "$status" -eq 0 ]
  run diff "${prefix_ckks}.pad" "${prefix_ckks}.decoded"
  assert [ "$status" -eq 0 ]
}

@test "CKKS: matrix data == decode(encode(data))" {
  run bash -c "$encode ${prefix_ckks}.dat ${prefix_ckks}.info CKKS --dims 2,3 > ${prefix_ckks}.encoded"
  assert [ "$status" -eq 0 ]
  # Number of lines should be 1 for header and 2 * 3 for data.
  assert [ "$(cat ${prefix_ckks}.encoded | wc -l)" -eq 7 ]
  run bash -c "$decode ${prefix_ckks}.encoded ${prefix_ckks}.info CKKS --nelements 12 > ${prefix_ckks}.decoded"
  assert [ "$status" -eq 0 ]
  diff "${prefix_ckks}.dat" "${prefix_ckks}.decoded"
}

@test "CKKS: matrix encode(data) does not fit" {
  run bash -c "$encode ${prefix_ckks}.dat ${prefix_ckks}.info CKKS --dims 1,3 > ${prefix_ckks}.encoded"
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Number of ptxts 6 > capacity of Matrix: 3 (dims: (1, 3))" ]
}

