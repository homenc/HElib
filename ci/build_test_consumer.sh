#!/usr/bin/env bash

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

set -xe

if [ "$#" -ne 3 ]; then
  echo "Wrong parameter number. Usage ./${0} <CONSUMER_FOLDER> <PACKAGE_BUILD> <CMAKE_INSTALL_PREFIX>"
  exit 1
fi

CONSUMER_FOLDER="${1}"
PACKAGE_BUILD="${2}"
CMAKE_INSTALL_PREFIX="${3}"
ROOT_DIR="$(pwd)"

cd "${ROOT_DIR}/${CONSUMER_FOLDER}"
mkdir build
cd build

# This test is quite brittle, but we can assume PACKAGE_BUILD is ON or OFF
if [ "${PACKAGE_BUILD}" == "ON" ]; then
  CMAKE_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}/helib_pack"
fi

cmake -Dhelib_DIR="${CMAKE_INSTALL_PREFIX}/share/cmake/helib" ..
make -j4 VERBOSE=1
cd "${ROOT_DIR}/${CONSUMER_FOLDER}/tests"
bats -t .

cd "${ROOT_DIR}"
