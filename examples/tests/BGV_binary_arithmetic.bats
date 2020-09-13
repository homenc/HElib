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

BGV_binary_arithmetic="$examples_bin/BGV_binary_arithmetic"

function setup {
  mkdir -p $tmp_folder
  cd $tmp_folder
  print-info-location
}

function teardown {
  cd -
  remove-test-directory "$tmp_folder"
}

# Calculate the pop count
function popcount {
  declare -i n="$1"
  declare -i i="$2"
  local count=0
  while [ $n -ne 0 ] && [ $i -gt 0 ]; do
    if [ $((n & 1)) -eq 1 ]; then
      ((count++))
    fi
    n=$((n >> 1))
    ((i--))
  done
  echo "$count"
}

@test "BGV_binary_arithmetic works" {
  run $BGV_binary_arithmetic 
  a=$(echo -e "$output" | awk '/^a =/{ printf $3 }')
  b=$(echo -e "$output" | awk '/^b =/{ printf $3 }')
  c=$(echo -e "$output" | awk '/^c =/{ printf $3 }')
  output1=$(echo -e "$output" | awk '/a\*b\+c/{ printf $3 }')
  output2=$(echo -e "$output" | awk '/a\+b\+c/{ printf $3 }')
  output3=$(echo -e "$output" | awk '/popcnt/{ printf $3 }')

  result1=$(($a*$b+$c))
  result2=$(($a+$b+$c))

  assert [ "$output1" == "$result1" ]
  assert [ "$output2" == "$result2" ]
  # Note: The HE function for this has a maximum of 15 bits
  assert [ "$output3" == "$(popcount $a 15)" ]
}
