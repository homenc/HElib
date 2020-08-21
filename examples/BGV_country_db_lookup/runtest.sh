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

queries=('Italy:Rome'
         'Germany:Berlin'
         'Spain:Madrid'
         'France:Paris'
         'Andorra:Andorra la Vella'
         'Antarctica:Country name not in the database.')

# Make sure we are running in the correct place
dir=$(dirname $0)

# Number of queries
rt=${#queries[@]}

for query in "${queries[@]}"; do
  # Capture the result value for comparison
  value=$( $dir/build/BGV_country_db_lookup "db_filename=$dir/countries_dataset.csv" <<< "${query%:*}" \
    | awk '/Query result:/{ $1=$2=""; print $0 }' \
    | awk '{ $1=$1; print $0 }') # Trim leading whitespace
  echo "${query%:*} gives '$value'"
  # If the query result matches then decrement the return value rt
  if [ "$value" = "${query#*:}" ]; then
    rt=$((rt-1))
  fi
done

# This will return the overall pass or fail
exit "$rt"
