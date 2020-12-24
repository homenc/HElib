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

function random-char-string {
  local N=${1:-8}
  echo $(head /dev/urandom | LC_ALL=C tr -dc A-Za-z0-9 | head -c $N)
}

tmp_folder="tmp_$(random-char-string)"
helib_root="../.."
examples_root="../.."
examples_bin="$examples_root/build/bin"

function helib_version {
  cat "${helib_root}/VERSION"
}

function assert {
  if "$@"; then 
    return 0
  else
    techo "Output: $output"
    techo "Status: $status"
    return 1
  fi
}

function print-info-location {
  if [ "$DEBUG" == "true" ] || [ "$DEBUG" == "1" ]; then
    # Whitespace after because in `-p` flag mode printout (BATS bug)
    # seems to overwrite test name which then appear afterwards anyway.
    techo "DEBUG info in ${tmp_folder}.                "
  fi
}

function remove-test-directory {
  if [ "$DEBUG" == "true" ] || [ "$DEBUG" == "1" ]; then
    : # Don't delete.
  elif [ -z "$DEBUG" ] || [ "$DEBUG" == "false" ] || [ "$DEBUG" == "0" ]; then
    rm -rf ./$1
  else            
    techo "Teardown: unrecognized value for DEBUG ${DEBUG}, assume false."
    rm -rf ./$1
  fi
}

function techo {
  local line=
  while IFS= read -r line; do
    echo "# $line" >&3
  done <<< "$1"
}
