#!/bin/bash

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

# You can pass your own clang-format as first arg.
FORMAT_PRG="${1:-"clang-format"}"
FORMAT_CMD="${FORMAT_PRG} -style=file -i"

# Check we can find clang-format
if [ -z $(which $FORMAT_PRG) ]; then
  echo "$FORMAT_PRG not found."
  exit 1
fi

# Check format program is correct major version
if [ ! $($FORMAT_PRG --version | cut -d' ' -f3 | cut -d'.' -f1) -ge 9 ]; then
  >&2 echo "Clang-format version below 9. Require 9+."
  exit 2
fi

function echo_header {
  # Truncate if required.
  local str=${1:0:72}

  # If odd length append a space.
  if [ $(( ${#str} % 2 )) -eq 1 ]; then str+=" "; fi

  # Build header.
  local n=$(( (80 - ${#str}) / 2 ))
  half_padding=$(printf "%0.s#" $(seq 1 $n) )

  echo "$half_padding $str $half_padding"
}

echo "Using '${FORMAT_CMD}' to perform formatting."

previous_dir=""
for file in $(find -E . -type f -regex '.*\.(c|h|cpp|hpp)' \
              ! -path '*/misc/*' \
              ! -path '*/build/*' \
              ! -name "PGFFT.*"); do
  current_dir=$(dirname $file)
  if [ "$current_dir" != "$previous_dir" ]; then
    previous_dir="$current_dir"
    echo_header "Formatting $current_dir"
  fi
  $FORMAT_CMD $file && echo $file
done

