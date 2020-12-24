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

BGV_packed_arithmetic="$examples_bin/BGV_packed_arithmetic"

function setup {
  mkdir -p $tmp_folder
  cd $tmp_folder
  print-info-location
}

function teardown {
  cd -
  remove-test-directory "$tmp_folder"
}

results=("{\"HElibVersion\":\"$(helib_version)\",\"content\":{\"scheme\":\"BGV\",\"slots\":[[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0],[0]]},\"serializationVersion\":\"0.0.1\",\"type\":\"Ptxt\"}"
         "{\"HElibVersion\":\"$(helib_version)\",\"content\":{\"scheme\":\"BGV\",\"slots\":[[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2],[2]]},\"serializationVersion\":\"0.0.1\",\"type\":\"Ptxt\"}")

@test "BGV_packed_arithmetic works" {
  run $BGV_packed_arithmetic 
  result1=$(echo -e "$output" | awk '/Decrypted Result:/{ print substr($0, index($0, $3)) ; exit 0 }')
  result2=$(echo -e "$output" | awk '/Decrypted Result:/{ i++; if(i==2) { print substr($0, index($0, $3)) } }')
  assert [ "$result1" == "${results[0]}" ]
  assert [ "$result2" == "${results[1]}" ]
}
