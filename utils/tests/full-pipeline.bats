#!/usr/bin/env bats

# Copyright (C) 2020 IBM Corp.
# Copyright (C) 2022 Intel Corporation
# SPDX-License-Identifier: Apache-2.0
#
# Added tests for separated SK, PK, and Key Switching matrices
#

load "std"

function setup {
  mkdir -p $tmp_folder
  cd $tmp_folder
  print-info-location
  check_python36
  check_locations
  create-bgv-toy-params
  createContext BGV "${prefix_bgv}.params" "$prefix_bgv"
  genData "${prefix_bgv}.dat" 12 "BGV"
  create-ckks-toy-params
  createContext CKKS "${prefix_ckks}.params" "$prefix_ckks"
  genData "${prefix_ckks}.dat" 12 "CKKS"
}

function teardown {
  cd -
  remove-test-directory "$tmp_folder"
}

@test "BGV: full-pipeline encode - encrypt - decrypt - decode" {
  run bash -c "$encode ${prefix_bgv}.dat ${prefix_bgv}.info BGV > ${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $encrypt "$pk_file_bgv" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "$sk_file_bgv" "${prefix_bgv}.ctxt" -o "$prefix_bgv.dec.ptxt"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_bgv}.dec.ptxt ${prefix_bgv}.info BGV > ${prefix_bgv}.decoded"
  assert [ "$status" -eq 0 ]
  diff "${prefix_bgv}.dat" "${prefix_bgv}.decoded"
}

@test "CKKS: full-pipeline encode - encrypt - decrypt - decode" {
  run bash -c "$encode ${prefix_ckks}.dat ${prefix_ckks}.info CKKS > ${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  run $encrypt "$pk_file_ckks" "${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "$sk_file_ckks" "${prefix_ckks}.ctxt" -o "$prefix_ckks.dec.ptxt"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_ckks}.dec.ptxt ${prefix_ckks}.info CKKS > ${prefix_ckks}.decoded"
  assert [ "$status" -eq 0 ]
  ${diff_threshold} ${prefix_ckks}.dat ${prefix_ckks}.decoded
}

@test "BGV: matrix full-pipeline encode - encrypt - decrypt - decode" {
  run bash -c "$encode ${prefix_bgv}.dat ${prefix_bgv}.info BGV --dims 2,3 > ${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $encrypt "$pk_file_bgv" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "$sk_file_bgv" "${prefix_bgv}.ctxt" -o "$prefix_bgv.dec.ptxt"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_bgv}.dec.ptxt ${prefix_bgv}.info BGV --nelements 12 > ${prefix_bgv}.decoded"
  assert [ "$status" -eq 0 ]
  diff "${prefix_bgv}.dat" "${prefix_bgv}.decoded"
}

@test "CKKS: matrix full-pipeline encode - encrypt - decrypt - decode" {
  run bash -c "$encode ${prefix_ckks}.dat ${prefix_ckks}.info CKKS --dims 2,3 > ${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  run $encrypt "$pk_file_ckks" "${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "$sk_file_ckks" "${prefix_ckks}.ctxt" -o "$prefix_ckks.dec.ptxt"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_ckks}.dec.ptxt ${prefix_ckks}.info CKKS --nelements 12 > ${prefix_ckks}.decoded"
  assert [ "$status" -eq 0 ]
  ${diff_threshold} ${prefix_ckks}.dat ${prefix_ckks}.decoded
}

@test "BGV (keygen): full-pipeline encode - encrypt - decrypt - decode" {
  keyGeneration BGV "${prefix_bgv}.params" "$prefix_bgv"
  run bash -c "$encode ${prefix_bgv}.dat ${prefix_bgv}.info BGV > ${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $encrypt "$enc_pk_file_bgv" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "$sk_file_bgv" "${prefix_bgv}.ctxt" -o "$prefix_bgv.dec.ptxt"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_bgv}.dec.ptxt ${prefix_bgv}.info BGV > ${prefix_bgv}.decoded"
  assert [ "$status" -eq 0 ]
  diff "${prefix_bgv}.dat" "${prefix_bgv}.decoded"
}

@test "BGV (keygen, sk only): full-pipeline encode - encrypt - decrypt - decode" {
  keyGeneration BGV "${prefix_bgv}.params" "$prefix_bgv" -s
  run bash -c "$encode ${prefix_bgv}.dat ${prefix_bgv}.info BGV > ${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $encrypt "$enc_pk_file_bgv" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "$sk_file_bgv" "${prefix_bgv}.ctxt" -o "$prefix_bgv.dec.ptxt" -s
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_bgv}.dec.ptxt ${prefix_bgv}.info BGV > ${prefix_bgv}.decoded"
  assert [ "$status" -eq 0 ]
  diff "${prefix_bgv}.dat" "${prefix_bgv}.decoded"
}

@test "CKKS (keygen): full-pipeline encode - encrypt - decrypt - decode" {
  keyGeneration CKKS "${prefix_ckks}.params" "$prefix_ckks"
  run bash -c "$encode ${prefix_ckks}.dat ${prefix_ckks}.info CKKS > ${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  run $encrypt "$enc_pk_file_ckks" "${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "$sk_file_ckks" "${prefix_ckks}.ctxt" -o "$prefix_ckks.dec.ptxt"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_ckks}.dec.ptxt ${prefix_ckks}.info CKKS > ${prefix_ckks}.decoded"
  assert [ "$status" -eq 0 ]
  ${diff_threshold} ${prefix_ckks}.dat ${prefix_ckks}.decoded
}

@test "CKKS (keygen, sk only): full-pipeline encode - encrypt - decrypt - decode" {
  keyGeneration CKKS "${prefix_ckks}.params" "$prefix_ckks" -s
  run bash -c "$encode ${prefix_ckks}.dat ${prefix_ckks}.info CKKS > ${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  run $encrypt "$enc_pk_file_ckks" "${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "$sk_file_ckks" "${prefix_ckks}.ctxt" -o "$prefix_ckks.dec.ptxt" -s
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_ckks}.dec.ptxt ${prefix_ckks}.info CKKS > ${prefix_ckks}.decoded"
  assert [ "$status" -eq 0 ]
  ${diff_threshold} ${prefix_ckks}.dat ${prefix_ckks}.decoded
}

@test "BGV (keygen): matrix full-pipeline encode - encrypt - decrypt - decode" {
  keyGeneration BGV "${prefix_bgv}.params" "$prefix_bgv"
  run bash -c "$encode ${prefix_bgv}.dat ${prefix_bgv}.info BGV --dims 2,3 > ${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $encrypt "$enc_pk_file_bgv" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "$sk_file_bgv" "${prefix_bgv}.ctxt" -o "$prefix_bgv.dec.ptxt"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_bgv}.dec.ptxt ${prefix_bgv}.info BGV --nelements 12 > ${prefix_bgv}.decoded"
  assert [ "$status" -eq 0 ]
  diff "${prefix_bgv}.dat" "${prefix_bgv}.decoded"
}

@test "BGV (keygen, sk only): matrix full-pipeline encode - encrypt - decrypt - decode" {
  keyGeneration BGV "${prefix_bgv}.params" "$prefix_bgv" -s
  run bash -c "$encode ${prefix_bgv}.dat ${prefix_bgv}.info BGV --dims 2,3 > ${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $encrypt "$enc_pk_file_bgv" "${prefix_bgv}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "$sk_file_bgv" "${prefix_bgv}.ctxt" -o "$prefix_bgv.dec.ptxt" -s
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_bgv}.dec.ptxt ${prefix_bgv}.info BGV --nelements 12 > ${prefix_bgv}.decoded"
  assert [ "$status" -eq 0 ]
  diff "${prefix_bgv}.dat" "${prefix_bgv}.decoded"
}

@test "CKKS (keygen): matrix full-pipeline encode - encrypt - decrypt - decode" {
  keyGeneration CKKS "${prefix_ckks}.params" "$prefix_ckks"
  run bash -c "$encode ${prefix_ckks}.dat ${prefix_ckks}.info CKKS --dims 2,3 > ${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  run $encrypt "$enc_pk_file_ckks" "${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "$sk_file_ckks" "${prefix_ckks}.ctxt" -o "$prefix_ckks.dec.ptxt"
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_ckks}.dec.ptxt ${prefix_ckks}.info CKKS --nelements 12 > ${prefix_ckks}.decoded"
  assert [ "$status" -eq 0 ]
  ${diff_threshold} ${prefix_ckks}.dat ${prefix_ckks}.decoded
}

@test "CKKS (keygen, sk only): matrix full-pipeline encode - encrypt - decrypt - decode" {
  keyGeneration CKKS "${prefix_ckks}.params" "$prefix_ckks" -s
  run bash -c "$encode ${prefix_ckks}.dat ${prefix_ckks}.info CKKS --dims 2,3 > ${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  run $encrypt "$enc_pk_file_ckks" "${prefix_ckks}.ptxt"
  assert [ "$status" -eq 0 ]
  run $decrypt "$sk_file_ckks" "${prefix_ckks}.ctxt" -o "$prefix_ckks.dec.ptxt" -s
  assert [ "$status" -eq 0 ]
  run bash -c "$decode ${prefix_ckks}.dec.ptxt ${prefix_ckks}.info CKKS --nelements 12 > ${prefix_ckks}.decoded"
  assert [ "$status" -eq 0 ]
  ${diff_threshold} ${prefix_ckks}.dat ${prefix_ckks}.decoded
}
