#!/bin/bash

set -ev

# if possible, ask for the precise number of processors,
# otherwise take 2 processors as reasonable default.
if [ -x /usr/bin/getconf ]; then
    NPROCESSORS=$(/usr/bin/getconf _NPROCESSORS_ONLN)
else
    NPROCESSORS=2
fi

MAKEFLAGS="j${NPROCESSORS}"
export MAKEFLAGS

mkdir build
cd build

cmake -DENABLE_TEST=ON ..
make 
make test
