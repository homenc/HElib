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
  create-ckks-toy-params
}

function teardown {
  cd -
  remove-test-directory "$tmp_folder"
}

@test "params file does not exist" {
  run "${create_context}" abcd.params 
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Could not open file 'abcd.params'." ]  
}

@test "create BGV toy params" {
  run "${create_context}" "${prefix_bgv}.params"
  assert [ "$status" -eq 0 ]
  filePrefix=$(sed -n 's,File prefix:.*\(test_bgv[0-9]*\).*,\1,p' <<< "${lines[0]}")
  # techo "Detected file prefix: '$filePrefix'"
  assert [ -f "${filePrefix}.sk" ]
  assert [ -f "${filePrefix}.pk" ]
}

@test "creates the correct files given a prefix" {
  run "${create_context}" "${prefix_bgv}.params" -o abcd
  assert [ "$status" -eq 0 ]
  assert [ -f "abcd.pk" ]
  assert [ -f "abcd.sk" ]
}

@test "create BGV toy params with info file" {
  run "${create_context}" "${prefix_bgv}.params" --info-file -o bcde
  assert [ "$status" -eq 0 ]
  assert [ -f "bcde.info" ]
}

@test "key switching matrices are created by default and indicated to cmd line" {
  run "${create_context}" "${prefix_bgv}.params" -o a
  assert [ "$status" -eq 0 ]
  assert [ "${lines[0]}" == "Key switching matrices created." ]
}

@test "key switching matrices are created by default and indicated in the info-file" {
  run "${create_context}" "${prefix_bgv}.params" --info-file -o defg
  assert [ "$status" -eq 0 ]
  assert [ -f "defg.info" ]
  line=$(head -n 1 "defg.info")
  assert [ "$line" == "Key switching matrices created." ]
}

@test "no-skm flag disables generation of key switching matrices" {
  run "${create_context}" "${prefix_bgv}.params" --info-file --no-skm -o cdef
  assert [ "$status" -eq 0 ]
  assert [ -f "cdef.info" ]
  line=$(awk '/Key switching matrices created/ { print $1 }' cdef.info)
  assert [ "$line" == '' ]
}

# Frobenius matrices created.\n

@test "Frobenius switching matrices are not created by default and indicated to cmd line" {
  run "${create_context}" "${prefix_bgv}.params" -o a
  assert [ "$status" -eq 0 ]
  line=$(echo $output | awk '/Frobenius matrices created./ { print $1 }')
  assert [ "${line}" == "" ]
}

@test "Frobenius switching matrices are not created by default and not indicated in the info-file" {
  run "${create_context}" "${prefix_bgv}.params" --info-file -o defg
  assert [ "$status" -eq 0 ]
  assert [ -f "defg.info" ]
  line=$(awk '/Frobenius matrices created./ { print $1 }' defg.info)
  assert [ "${line}" == "" ]
}

@test "frob-skm creates Frobenius matrices and indicates it to cmd line" {
  run "${create_context}" "${prefix_bgv}.params" -o a --frob-skm
  assert [ "$status" -eq 0 ]
  assert [ "${lines[1]}" == "Frobenius matrices created." ]
}

@test "frob-skm creates Frobenius matrices and indicates it to the info-file" {
  run "${create_context}" "${prefix_bgv}.params" --info-file --frob-skm -o defg
  assert [ "$status" -eq 0 ]
  assert [ -f "defg.info" ]
  line="$(sed '2q;d' defg.info)"
  assert [ "$line" == "Frobenius matrices created." ]
}

@test "no-skm and frob-skm flags are not a valid combination" {
  run "${create_context}" "${prefix_bgv}.params" --no-skm --frob-skm -o cdef
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Frobenius matrices reqires switch-key matrices to be generated." ]
}

@test "default scheme fails if p < 2" {
  run "${create_context}" "${prefix_ckks}.params"
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "BGV invalid plaintext modulus. In BGV it must be a prime number greater than 1." ]
}

@test "BGV: fails if p < 2" {
  run "${create_context}" "${prefix_ckks}.params" --scheme BGV
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "BGV invalid plaintext modulus. In BGV it must be a prime number greater than 1." ]
}

@test "CKKS: fails if p != -1" {
  run "${create_context}" "${prefix_bgv}.params" --scheme CKKS
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "CKKS invalid plaintext modulus. In CKKS it must be set to -1." ]
}

@test "BGV: fails on non-prime p" {
  sed -e 's/p=[0-9][0-9]*/p=4/' < "${prefix_bgv}.params" > ${prefix_bgv}_non_prime_p.params
  run "${create_context}" "${prefix_bgv}_non_prime_p.params" --scheme BGV
  assert [ "$status" -ne 0 ]
  assert [ "${lines[1]}" == "Modulus pp is not prime (nor -1)" ]
}

@test "scheme flag works correctly for BGV and CKKS" {
  run "${create_context}" "${prefix_bgv}.params" --scheme BGV -o bgv
  assert [ "$status" -eq 0 ]
  p=$(sed -n 's,.*p = \([0-9]*\).*,\1,p' <<< "${lines[1]}")
  # techo "${lines[1]}"
  # techo "p = '$p'"
  assert [ "$p" != "-1" ]
  run "${create_context}" "${prefix_ckks}.params" --scheme CKKS -o ckks
  assert [ "$status" -eq 0 ]
  p=$(sed -n 's,.*p = \(-[0-9]*\).*,\1,p' <<< "${lines[1]}")
  # techo "${lines[1]}"
  # techo "p = '$p'"
  assert [ "$p" == "-1" ]
}

@test "CKKS: fails on odd m" {
  sed -e 's/m=[0-9][0-9]*/m=33/' < "${prefix_ckks}.params" > ${prefix_ckks}_odd_m.params
  run "${create_context}" "${prefix_ckks}_odd_m.params" --scheme CKKS
  assert [ "$status" -ne 0 ]
  assert [ "${lines[1]}" == "CKKS scheme only supports m as a power of two." ]
}

@test "CKKS: fails on m non-power of 2" {
  sed -e 's/m=[0-9][0-9]*/m=34/' < "${prefix_ckks}.params" > ${prefix_ckks}_nonp2_m.params
  run "${create_context}" "${prefix_ckks}_nonp2_m.params" --scheme CKKS
  assert [ "$status" -ne 0 ]
  assert [ "${lines[1]}" == "CKKS scheme only supports m as a power of two." ]
}

@test "fails on invalid scheme option" {
  run "${create_context}" "${prefix_bgv}.params" --scheme BFV -o bfv
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Unrecognized scheme 'BFV'." ]
}

@test "bootstrap = NONE option works as normal" {
  run "${create_context}" "${prefix_bgv}.params" --bootstrap NONE -o test
  assert [ "$status" -eq 0 ]
  assert [ -f "test.pk" ]
  assert [ -f "test.sk" ]
}

@test "bootstrap = THIN option works with bootstrap params" {
  create-bootstrap-toy-params
  run "${create_context}" "${prefix_bootstrap}.params" --bootstrap THIN -o thin
  assert [ "$status" -eq 0 ]
  assert [ -f "thin.pk" ]
  assert [ -f "thin.sk" ]
}

@test "bootstrap = THICK option works with bootstrap params" {
  create-bootstrap-toy-params
  run "${create_context}" "${prefix_bootstrap}.params" --bootstrap THICK -o thick
  assert [ "$status" -eq 0 ]
  assert [ -f "thick.pk" ]
  assert [ -f "thick.sk" ]
}

@test "fails when passing invalid bootstrap option" {
  create-bootstrap-toy-params
  run "${create_context}" "${prefix_bootstrap}.params" --bootstrap INVALID
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Bad boostrap option: INVALID.  Allowed options are NONE, THIN, THICK." ]
}

@test "bootstrap fails when --no-skm is specified" {
  create-bootstrap-toy-params
  run "${create_context}" "${prefix_bootstrap}.params" --bootstrap THIN --no-skm
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Cannot generate bootstrappable context without switch-key and frobenius matrices." ]
  run "${create_context}" "${prefix_bootstrap}.params" --bootstrap THICK --no-skm
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Cannot generate bootstrappable context without switch-key and frobenius matrices." ]
}

@test "bootstrap fails when missing parameters (mvec, gens, ords)" {
  create-bootstrap-toy-params
  awk '!/mvec/' "${prefix_bootstrap}.params" > "${prefix_bootstrap}_mvec.params"
  run "${create_context}" "${prefix_bootstrap}_mvec.params" --bootstrap THICK
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Missing mvec parameter for bootstrapping in ${prefix_bootstrap}_mvec.params." ]
  awk '!/gens/' "${prefix_bootstrap}.params" > "${prefix_bootstrap}_gens.params"
  run "${create_context}" "${prefix_bootstrap}_gens.params" --bootstrap THICK
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Missing gens parameter for bootstrapping in ${prefix_bootstrap}_gens.params." ]
  awk '!/ords/' "${prefix_bootstrap}.params" > "${prefix_bootstrap}_ords.params"
  run "${create_context}" "${prefix_bootstrap}_ords.params" --bootstrap THICK
  assert [ "$status" -ne 0 ]
  assert [ "$output" == "Missing ords parameter for bootstrapping in ${prefix_bootstrap}_ords.params." ]
}
# bootstrapping creates SKM and frobenius
@test "bootstrap creates SKM, Frobenius matrices and recrypt data and indicates it to cmd line" {
  create-bootstrap-toy-params
  run "${create_context}" "${prefix_bootstrap}.params" --bootstrap THICK -o defg
  assert [ "$status" -eq 0 ]
  assert [ "${lines[0]}" == "Key switching matrices created." ]
  assert [ "${lines[1]}" == "Frobenius matrices created." ]
  assert [ "${lines[2]}" == "Recrypt data created." ]
}

@test "bootstrap creates SKM, Frobenius matrices and recrypt data and indicates it to the info-file" {
  create-bootstrap-toy-params
  run "${create_context}" "${prefix_bootstrap}.params" --info-file --bootstrap THICK -o defg
  assert [ "$status" -eq 0 ]
  assert [ -f "defg.info" ]
  line="$(sed '1q;d' defg.info)"
  assert [ "$line" == "Key switching matrices created." ]
  line="$(sed '2q;d' defg.info)"
  assert [ "$line" == "Frobenius matrices created." ]
  line="$(sed '3q;d' defg.info)"
  assert [ "$line" == "Recrypt data created." ]
}

@test "bootstrap thick context can perform bootstrapping" {
  create-bootstrap-toy-params
  cat "${prefix_bootstrap}.params"
  run "${create_context}" "${prefix_bootstrap}.params" --bootstrap THICK -o boots
  assert [ "$status" -eq 0 ]
  run "${test_bootstrap}" --thick boots.sk --quiet
  assert [ "$status" -eq 0 ]
}

@test "bootstrap thin context can perform bootstrapping" {
 create-bootstrap-toy-params
 run "${create_context}" "${prefix_bootstrap}.params" --bootstrap THIN -o boots
 assert [ "$status" -eq 0 ]
 run "${test_bootstrap}" boots.sk --quiet
 assert [ "$status" -eq 0 ]
}
