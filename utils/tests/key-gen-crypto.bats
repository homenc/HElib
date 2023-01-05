#!/usr/bin/env bats

# Copyright (C) 2022 Intel Corporation
# SPDX-License-Identifier: Apache-2.0

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
  keyGeneration BGV "${prefix_bootstrap}.params" "${prefix_bootstrap}" "--bootstrap THIN"
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
  keyGeneration BGV "${prefix_bgv}.params" "$prefix_bgv"
  genData "${prefix_bgv}.dat" 12 "BGV" 
  "${encode}" "${prefix_bgv}.dat" "${prefix_bgv}.info" "BGV" > "${prefix_bgv}.ptxt"
  create-ckks-toy-params
  keyGeneration CKKS "${prefix_ckks}.params" "$prefix_ckks"
  genData "${prefix_ckks}.dat" 12 "CKKS"
  "${encode}" "${prefix_ckks}.dat" "${prefix_ckks}.info" "CKKS" > "${prefix_ckks}.ptxt"
}

function teardown {
  cd $tmp_folder/..
  remove-test-directory "$tmp_folder"
}

# Common parameter-related tests
@test "encrypt with Enc key works with batch size greater than 1" {
  run $encrypt "${enc_pk_file_bgv}" "${prefix_bgv}.ptxt" -b 10
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bgv}.ctxt" ]
}

@test "encrypt with Eval key works with batch size greater than 1" {
  run $encrypt "${eval_pk_file_bgv}" "${prefix_bgv}.ptxt" -b 10
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bgv}.ctxt" ]
}

@test "decrypt (keygen) works with batch size greater than 1" {
  run $encrypt "${enc_pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "${sk_file_bgv}" "${prefix_bgv}.ctxt" -o result_bgv.decrypt -b 10
  assert [ "$status" -eq 0 ]
  assert [ -f "result_bgv.decrypt" ]
}

@test "decrypt (keygen, sk only) works with batch size greater than 1" {
  keyGeneration BGV "${prefix_bgv}.params" "$prefix_bgv" -s
  run $encrypt "${enc_pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "${sk_file_bgv}" "${prefix_bgv}.ctxt" -o result_bgv.decrypt -b 10 -s
  assert [ "$status" -eq 0 ]
  assert [ -f "result_bgv.decrypt" ]
}

@test "decrypt (keygen) works with nthreads size greater than 1" {
  run $encrypt "${enc_pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "${sk_file_bgv}" "${prefix_bgv}.ctxt" -o result_bgv.decrypt -n 4
  assert [ "$status" -eq 0 ]
  assert [ -f "result_bgv.decrypt" ]
}

@test "decrypt (keygen, sk only) works with nthreads size greater than 1" {
  keyGeneration BGV "${prefix_bgv}.params" "$prefix_bgv" -s
  run $encrypt "${enc_pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "${sk_file_bgv}" "${prefix_bgv}.ctxt" -o result_bgv.decrypt -n 4 -s
  assert [ "$status" -eq 0 ]
  assert [ -f "result_bgv.decrypt" ]
}

@test "BGV (keygen): data == decrypt(encrypt(data))" {
  run $encrypt "${enc_pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bgv}.ctxt" ]
  run $decrypt "${sk_file_bgv}" "${prefix_bgv}.ctxt" -o result_bgv.decrypt
  assert [ "$status" -eq 0 ]
  assert [ -f result_bgv.decrypt ]
  diff "${prefix_bgv}.ptxt" result_bgv.decrypt
  assert [ "$status" -eq 0 ]
}

@test "BGV (keygen, sk only): data == decrypt(encrypt(data))" {
  keyGeneration BGV "${prefix_bgv}.params" "$prefix_bgv" -s
  run $encrypt "${enc_pk_file_bgv}" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bgv}.ctxt" ]
  run $decrypt "${sk_file_bgv}" "${prefix_bgv}.ctxt" -o result_bgv.decrypt -s
  assert [ "$status" -eq 0 ]
  assert [ -f result_bgv.decrypt ]
  diff "${prefix_bgv}.ptxt" result_bgv.decrypt
  assert [ "$status" -eq 0 ]
}

@test "BGV (keygen): data == decrypt(encrypt(data)) with bootstrappable context" {
  # This test needs bootstrapping context to work
  generate-bootstrap-data
  run $encrypt "${eval_pk_file_bootstrap}" "${prefix_bootstrap}.ptxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bootstrap}."ctxt ]
  run $decrypt "${sk_file_bootstrap}" "${prefix_bootstrap}.ctxt" -o "${prefix_bootstrap}.decrypt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bootstrap}.decrypt" ]
  diff "${prefix_bootstrap}.ptxt" "${prefix_bootstrap}.decrypt"
  assert [ "$status" -eq 0 ]
}

@test "BGV (keygen, sk only): data == decrypt(encrypt(data)) with bootstrappable context" {
  # This test needs bootstrapping context to work
  create-bootstrap-toy-params
  keyGeneration BGV "${prefix_bootstrap}.params" "${prefix_bootstrap}" "--bootstrap THIN" -s
  genData "${prefix_bootstrap}.dat" 12 "BGV" 
  "${encode}" "${prefix_bootstrap}.dat" "${prefix_bootstrap}.info" "BGV" > "${prefix_bootstrap}.ptxt"
  run $encrypt "${eval_pk_file_bootstrap}" "${prefix_bootstrap}.ptxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bootstrap}."ctxt ]
  run $decrypt "${sk_file_bootstrap}" "${prefix_bootstrap}.ctxt" -o "${prefix_bootstrap}.decrypt" -s
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_bootstrap}.decrypt" ]
  diff "${prefix_bootstrap}.ptxt" "${prefix_bootstrap}.decrypt"
  assert [ "$status" -eq 0 ]
}

# CKKS encryption tests
@test "CKKS: data == decrypt(encrypt(data))" {
  run $encrypt "${enc_pk_file_ckks}" "${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_ckks}.ctxt" ]
  run $decrypt "${sk_file_ckks}" "${prefix_ckks}.ctxt" -o "${prefix_ckks}.decrypt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_ckks}.decrypt" ]
  ${diff_threshold} --decrypt "${prefix_ckks}.ptxt" "${prefix_ckks}.decrypt"
  assert [ "$status" -eq 0 ]
}

@test "CKKS (keygen, sk only): data == decrypt(encrypt(data))" {
  keyGeneration CKKS "${prefix_ckks}.params" "$prefix_ckks" -s
  run $encrypt "${enc_pk_file_ckks}" "${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_ckks}.ctxt" ]
  run $decrypt "${sk_file_ckks}" "${prefix_ckks}.ctxt" -o "${prefix_ckks}.decrypt" -s
  assert [ "$status" -eq 0 ]
  assert [ -f "${prefix_ckks}.decrypt" ]
  ${diff_threshold} --decrypt "${prefix_ckks}.ptxt" "${prefix_ckks}.decrypt"
  assert [ "$status" -eq 0 ]
}