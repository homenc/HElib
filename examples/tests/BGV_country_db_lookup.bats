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

BGV_country_db_lookup="$examples_bin/BGV_country_db_lookup"
db_file="$examples_root/BGV_country_db_lookup/countries_dataset.csv"

function setup {
  mkdir -p $tmp_folder
  cd $tmp_folder
  print-info-location
}

function teardown {
  cd -
  remove-test-directory "$tmp_folder"
}

# Grabs the result of the query and assigns it to result
function getResult {
  result=$(echo -e "$output" | awk '/Query result:/ { printf substr($0, index($0,$3)) }')
}

@test "BGV_country_db_lookup returns London for England" {
  run $BGV_country_db_lookup db_filename=$db_file <<< "England"
  getResult
  assert [ "$result" == "London" ]
}

@test "BGV_country_db_lookup returns message if country not found" {
  run $BGV_country_db_lookup db_filename=$db_file <<< "Antarctica"
  getResult
  assert [ "$result" == "Country name not in the database." ]
}

@test "BGV_country_db_lookup returns Sarajevo for Bosnia and Herzegovina" {
  run $BGV_country_db_lookup db_filename=$db_file <<< "Bosnia and Herzegovina"
  getResult
  assert [ "$result" == "Sarajevo" ]
}

@test "BGV_country_db_lookup returns Andorra la Vella for Andorra" {
  run $BGV_country_db_lookup db_filename=$db_file <<< "Andorra"
  getResult
  assert [ "$result" == "Andorra la Vella" ]
}

@test "BGV_country_db_lookup does not find match if country has trailing space" {
  run $BGV_country_db_lookup db_filename=$db_file <<< "Scotland "
  getResult
  assert [ "$result" == "Country name not in the database." ]
}
