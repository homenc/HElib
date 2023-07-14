#!/usr/bin/env bash

# Copyright (C) 2020-2021 IBM Corp.
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

if [ "$#" -ne 5 ]; then
  echo "Wrong parameter number. Usage ./${0} <PACKAGE_BUILD> <RUNNER_OS> <C_COMPILER> <CXX_COMPILER> <USE_INTEL_HEXL>"
  exit 1
fi

PACKAGE_BUILD="${1}"
RUNNER_OS="${2}"
C_COMPILER="${3}"
CXX_COMPILER="${4}"
USE_INTEL_HEXL="${5}"

if [ "${USE_INTEL_HEXL}" == "ON" ]; then
  git clone https://github.com/intel/hexl.git -b v1.2.4
  cd hexl
  cmake -B build \
    -DCMAKE_INSTALL_PREFIX=./ \
    -DHEXL_SHARED_LIB=OFF \
    -DCMAKE_C_COMPILER="${C_COMPILER}" \
    -DCMAKE_CXX_COMPILER="${CXX_COMPILER}"
  cmake --build build -j4 --target install
  cd ../
fi

cd $HOME

# Remove substr after dash char e.g. ubuntu20.04 -> ubuntu
if [ "${RUNNER_OS%-*}" == "ubuntu" ]; then
    sudo apt-get -yq --no-install-suggests --no-install-recommends install libgmp-dev libntl-dev bats
fi

if [ "${RUNNER_OS}" == "macos-latest" ]  || [ "${RUNNER_OS}" == "macos-12" ] || [ "${RUNNER_OS}" == "macos-13" ] ; then
  brew install ntl bats-core
fi

cd "$HOME"
