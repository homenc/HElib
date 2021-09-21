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

# Copyright (C) 2021 Intel Corporation
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

set -xe

if [ "$#" -ne 6 ]; then
  echo "Wrong parameter number. Usage ./${0} <PACKAGE_BUILD> <CMAKE_INSTALL_PREFIX> <C_COMPILER> <CXX_COMPILER> <USE_INTEL_HEXL> <INTEL_HEXL_DIR>"
  exit 1
fi

PACKAGE_BUILD="${1}"
CMAKE_INSTALL_PREFIX="${2}"
C_COMPILER="${3}"
CXX_COMPILER="${4}"
USE_INTEL_HEXL="${5}"
INTEL_HEXL_DIR="${6}"
ROOT_DIR="$(pwd)"

mkdir build
cd build
# We assume PACKAGE_BUILD argument to be a valid cmake option
cmake -DPACKAGE_BUILD="${PACKAGE_BUILD}" \
      -DCMAKE_INSTALL_PREFIX="${CMAKE_INSTALL_PREFIX}" \
      -DUSE_INTEL_HEXL="${USE_INTEL_HEXL}" \
      -DCMAKE_C_COMPILER="${C_COMPILER}" \
      -DCMAKE_CXX_COMPILER="${CXX_COMPILER}" \
      -DHEXL_DIR="${INTEL_HEXL_DIR}" \
      -DBUILD_SHARED=OFF \
      -DENABLE_TEST=ON \
      -DTARGET_ARCHITECTURE=x86-64 \
      ..
make -j4 VERBOSE=1
make install
cd "${ROOT_DIR}"
