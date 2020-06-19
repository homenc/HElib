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

# Specialised function for creating BGV data with invalid slot size
function generate-invalid-bgv-data {
  genData "${prefix_bgv}_invalid.dat" 4 "BGV" "--columns 3" # 3 in a slot
  "${encode}" "${prefix_bgv}_invalid.dat" "${prefix_bgv}.info" "BGV" > "${prefix_bgv}_invalid.ptxt"
}

# Specialised function for creating CKKS data with invalid slot size
function generate-invalid-ckks-data {
  genData "${prefix_ckks}_invalid.dat" 4 "CKKS" "--columns 3" # 3 in a slot
  "${encode}" "${prefix_ckks}_invalid.dat" "${prefix_ckks}.info" "CKKS" > "${prefix_ckks}_invalid.ptxt"
}

# Only needed for BGV context with bootstrapping
function generate-bootstrap-data {
  create-bootstrap-toy-params
  createContext BGV "${prefix_bootstrap}.params" "${prefix_bootstrap}" "--bootstrap THIN"
  genData "${prefix_bootstrap}.dat" 12 "BGV" 
  "${encode}" "${prefix_bootstrap}.dat" "${prefix_bootstrap}.info" "BGV" > "${prefix_bootstrap}.ptxt"
}

function setup {
  mkdir -p $tmp_folder
  cd $tmp_folder
  print-info-location
  check_python36
  check_locations
  create-bgv-toy-params
  createContext BGV "${prefix_bgv}.params" "$prefix_bgv"
  genData "${prefix_bgv}.dat" 12 "BGV" 
  "${encode}" "${prefix_bgv}.dat" "${prefix_bgv}.info" "BGV" > "${prefix_bgv}.ptxt"
  create-ckks-toy-params
  createContext CKKS "${prefix_ckks}.params" "$prefix_ckks"
  genData "${prefix_ckks}.dat" 12 "CKKS"
  "${encode}" "${prefix_ckks}.dat" "${prefix_ckks}.info" "CKKS" > "${prefix_ckks}.ptxt"
}

function teardown {
  cd -
  remove-test-directory "$tmp_folder"
}

# Common parameter-related tests
@test "encrypt works with batch size > 1" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt" -b 10
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bgv}.ctxt" ]
}

@test "encrypt works with nthreads > 1" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt" -n 4
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bgv}.ctxt" ]
}

@test "decrypt works with batch size > 1" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "${sk_file_bgv}" "${prefix_bgv}.ctxt" -o result_bgv.decrypt -b 10
  assert [ "$status" -eq 0 ]
  assert [ -f "result_bgv.decrypt" ]
}

@test "decrypt works with nthreads size > 1" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "${sk_file_bgv}" "${prefix_bgv}.ctxt" -o result_bgv.decrypt -n 4
  assert [ "$status" -eq 0 ]
  assert [ -f "result_bgv.decrypt" ]
}

@test "encrypt fails with invalid nthreads < 0" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt" -n -1
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Number of threads must be a positive integer." ]
}

@test "decrypt fails with invalid nthreads < 0" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "${sk_file_bgv}" "${prefix_bgv}.ctxt" -o result_bgv.decrypt -n -1
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Number of threads must a be positive integer." ]
}

@test "encrypt fails with invalid batch size = -1" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt" -b -1
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Batch size must be a positive integer." ]
}

@test "encrypt warns if batch size < nthreads" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt" -b 1 -n 2
  assert [ "$status" -eq 0 ]
  assert [ "${lines[0]}" == "WARNING: Not enough elements in the batch (1) to run with 2 threads." ]
  assert [ "${lines[1]}" == "To achieve better performance set batch size (-b) to be equal to or more than the number of threads (-n). Ideally the batch size should be a multiple of the number of threads." ]
}

@test "decrypt fails with invalid batch size = -1" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "${sk_file_bgv}" "${prefix_bgv}.ctxt" -o result_bgv.decrypt -b -1
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Batch size must be a positive integer." ]
}

@test "decrypt warns if batch size < nthreads" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "${sk_file_bgv}" "${prefix_bgv}.ctxt" -o result_bgv.decrypt -b 1 -n 2
  assert [ "$status" -eq 0 ]
  assert [ "${lines[0]}" == "WARNING: Not enough elements in the batch (1) to run with 2 threads." ]
  assert [ "${lines[1]}" == "To achieve better performance set batch size (-b) to be equal to or more than the number of threads (-n). Ideally the batch size should be a multiple of the number of threads." ]
}

# BGV encryption-related tests TESTS
@test "BGV: encrypt outputs to file with same prefix as input by default" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bgv}.ctxt" ]
}

@test "BGV: encrypt outputs to user specified file" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt" -o result_bgv.ctxt
  assert [ "$status" -eq 0 ]
  assert [ -f "result_bgv.ctxt" ]
}

@test "BGV: data == decrypt(encrypt(data))" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bgv}.ctxt" ]
  run $decrypt "${sk_file_bgv}" "${prefix_bgv}.ctxt" -o result_bgv.decrypt
  assert [ "$status" -eq 0 ]
  assert [ -f result_bgv.decrypt ]
  diff "${prefix_bgv}.ptxt" result_bgv.decrypt
  assert [ "$status" -eq 0 ]
}

@test "BGV: decrypt prints to stdout by default and data == decrypt(encrypt(data))" {
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bgv}.ctxt" ]
  run $decrypt "${sk_file_bgv}" "${prefix_bgv}.ctxt"
  assert [ "$status" -eq 0 ]
  diff <(echo "$output") "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
}

@test "BGV: data == decrypt(encrypt(data)) with bootstrappable context" {
  # This test needs bootstrapping context to work
  generate-bootstrap-data
  run $encrypt "${pk_file_bootstrap}" "${prefix_bootstrap}.ptxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bootstrap}."ctxt ]
  run $decrypt "${sk_file_bootstrap}" "${prefix_bootstrap}.ctxt" -o "${prefix_bootstrap}.decrypt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bootstrap}.decrypt" ]
  diff "${prefix_bootstrap}.ptxt" "${prefix_bootstrap}.decrypt"
  assert [ "$status" -eq 0 ]
}

@test "BGV: encrypt fails if data has more elements than slots" {
  sed 's/nslots = \([0-9][0-9]*\)/nslots = 1\1/' < "${prefix_bgv}.info" > "${prefix_bgv}_more_slots.info"
  "${encode}" "${prefix_bgv}.dat" "${prefix_bgv}_more_slots.info" "BGV" > "${prefix_bgv}_more_slots.ptxt"
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}_more_slots.ptxt"
  assert [ "$status" -ne 0 ]
  assert [ "${lines[0]}" == "Exit due to IOError thrown:" ]
  assert [ "${lines[1]}" == "Cannot deserialize to Ptxt: not enough slots.  Trying to deserialize 13 elements.  Got 3 slots." ]
  assert [ ! -f "${prefix_bgv}_more_slots.ctxt" ]
}

@test "BGV: encrypt fails if data has more elements in slots" {
  generate-invalid-bgv-data
  run $encrypt "${pk_file_bgv}" "${prefix_bgv}_invalid.ptxt"
  assert [ "$status" -ne 0 ]
  assert [ "${lines[0]}" == "Exit due to IOError thrown:" ]
  assert [ "${lines[1]}" == "Cannot deserialize to PolyMod: Degree is too small.  Trying to deserialize 3 coefficients.  Degree is 2." ]
}

# CKKS encryption-related tests TESTS
@test "CKKS: encrypt outputs to file with same prefix as input by default" {
  run $encrypt "${pk_file_ckks}" "${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_ckks}.ctxt" ]
}

@test "CKKS: encrypt outputs to user specified file" {
  run $encrypt "${pk_file_ckks}" "${prefix_ckks}.ptxt" -o "${prefix_ckks}.ctxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_ckks}.ctxt" ]
}

@test "CKKS: data == decrypt(encrypt(data))" {
  run $encrypt "${pk_file_ckks}" "${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_ckks}.ctxt" ]
  run $decrypt "${sk_file_ckks}" "${prefix_ckks}.ctxt" -o "${prefix_ckks}.decrypt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_ckks}.decrypt" ]
  ${diff_threshold} --decrypt "${prefix_ckks}.ptxt" "${prefix_ckks}.decrypt"
  assert [ "$status" -eq 0 ]
}

@test "CKKS: decrypt prints to stdout by default and data == decrypt(encrypt(data))" {
  run $encrypt "${pk_file_ckks}" "${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_ckks}.ctxt" ]
  run $decrypt "${sk_file_ckks}" "${prefix_ckks}.ctxt"
  assert [ "$status" -eq 0 ]
  ${diff_threshold} --decrypt "${prefix_ckks}.ptxt" <(echo "$output")
  assert [ "$status" -eq 0 ]
}

@test "CKKS: encrypt fails if data has more elements than slots" {
  sed 's/nslots = \([0-9][0-9]*\)/nslots = 1\1/' < "${prefix_ckks}.info" > "${prefix_ckks}_more_slots.info"
  "${encode}" "${prefix_ckks}.dat" "${prefix_ckks}_more_slots.info" "CKKS" > "${prefix_ckks}_more_slots.ptxt"
  run $encrypt "${pk_file_ckks}" "${prefix_ckks}_more_slots.ptxt"
  assert [ "$status" -ne 0 ]
  assert [ "${lines[0]}" == "Exit due to IOError thrown:" ]
  assert [ "${lines[1]}" == "Cannot deserialize to Ptxt: not enough slots.  Trying to deserialize 12 elements.  Got 2 slots." ]
  assert [ ! -f "${prefix_ckks}_more_slots.ctxt" ]
}

@test "CKKS: encrypt fails if data has more elements in slots" {
  generate-invalid-ckks-data
  run $encrypt "${pk_file_ckks}" "${prefix_ckks}_invalid.ptxt"
  assert [ "$status" -ne 0 ]
  assert [ "${lines[0]}" == "Exit due to IOError thrown:" ]
  assert [ "${lines[1]}" == "CKKS expects maximum of 2 values per slot (real, imag). Got 3 instead." ]
}
