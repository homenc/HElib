#!/bin/bash

# Script for building and testing HElib in various configurations.
# Default behaviour is build homenc/HElib master using the package build
# without running any tests.

# Default arguments
repo=https://github.com/homenc/HElib.git
branch=master
package=true
tests=false
examples=false
utils=false
gbench=false
helib_install="$HOME/helib_build_install"
helib_DIR="${helib_install}/helib_pack/share/cmake/helib/"

# Function for printing usage info.
function printUsage {
  echo "Usage: CMD [-h] [-r <repo>] [-b <branch>] [-p] [-l] [-t] [-e] [-u] [-g] [-a]" 
  echo "    -h             Displays this help message."
  echo "    -r <repo>      HElib repo to clone (Default = https://github.com/IBM-HElib/HElib.git)."
  echo "    -b <branch>    Branch of HElib to checkout (Default = master)."
  echo "    -p             Flag to indicate a package build (This is the default build type)."
  echo "    -l             Flag to indicate a library build."
  echo "    -t             Run the HElib Google tests."
  echo "    -e             Build and test the HElib examples directory."
  echo "    -u             Build and test the HElib utils directory."
  echo "    -g             Build the HElib Google benchmark directory."
  echo "    -a             Build and test all of the above."
}

while getopts ":hr:b:plteuga" opt; do
  case ${opt} in
    h ) printUsage
        exit 0
        ;;
    r ) repo=$OPTARG
        ;;
    b ) branch=$OPTARG
        ;;
    p ) package=true
        ;;
    l ) package=false
        helib_DIR="${helib_install}/share/cmake/helib/"
        ;;
    t ) tests=true 
        ;;
    e ) examples=true
        ;;
    u ) utils=true
        ;;
    g ) gbench=true
        ;;
    a ) tests=true
        examples=true
        utils=true
        gbench=true
        ;;
    : ) echo "Invalid option: $OPTARG requires an argument." 1>&2
        printUsage
        exit 1
        ;;
    \?) echo "Invalid option: $OPTARG" 1>&2
        printUsage
        exit 1
        ;;
  esac
done


# Clone and build HElib in the home directory
cd
git clone ${repo} # Default = https://github.com/IBM-HElib/HElib.git
cd HElib
git checkout ${branch} # Default = master
mkdir build
cd build
cmake -DPACKAGE_BUILD=${package} -DBUILD_SHARED=ON -DCMAKE_INSTALL_PREFIX=${helib_install} -DENABLE_TEST=${tests} ..
make -j4 VERBOSE=1
make install

# Run the google tests
if [[ ${tests} = "true" ]]; then
  ctest -j4 --output-on-failure --no-compress-output --test-output-size-passed 32768 --test-output-size-failed 262144 -T Test
fi

# Build and test the examples
if [[ ${examples} = "true" ]]; then
  cd ../examples
  mkdir build
  cd build
  cmake -Dhelib_DIR=${helib_DIR} ..
  make -j4 VERBOSE=1
  cd ../tests
  bats -j 4 .
  cd ..
fi

# Build and test the utilities
if [[ ${utils} = "true" ]]; then
  cd ../utils
  mkdir build
  cd build
  cmake -Dhelib_DIR=${helib_DIR} ..
  make -j4 VERBOSE=1
  cd ../tests
  bats -j 4 .
  cd ..
fi

# Build benchmarks
if [[ ${gbench} = "true" ]]; then
  cd ../benchmarks
  mkdir build
  cd build
  cmake -Dhelib_DIR=${helib_DIR} ..
  make -j4 VERBOSE=1
fi
