# Build and Tests

## Introduction
This directory provides dockerfiles and a script for building and testing
HElib. This can be used to test various versions and configurations of HElib in
self-contained docker containers.

## What is provided
Currently provided:
- Dockerfiles:
    - Ubuntu 20.04
    - Fedora 33
    - CentOS 8
- Test script for cloning, build and testing HElib and all subprojects

The Dockerfiles provided contain recipes for installing all prerequisites for
HElib as well as the test script found in `build_scripts`.

## Running the build script

The script `build_and_test_helib.sh` has a help method which prints the usage
to the console and can be run with
```
$ ./build_and_test_helib.sh -h
Usage: CMD [-h] [-r <repo>] [-b <branch>] [-p] [-l] [-t] [-e] [-u] [-g] [-a]
    -h             Displays this help message.
    -r <repo>      HElib repo to clone (Default = https://github.com/IBM-HElib/HElib.git).
    -b <branch>    Branch of HElib to checkout (Default = master).
    -p             Flag to indicate a package build (This is the default build type).
    -l             Flag to indicate a library build.
    -t             Run the HElib Google tests.
    -e             Build and test the HElib examples directory.
    -u             Build and test the HElib utils directory.
    -g             Build the HElib Google benchmark directory.
    -a             Build and test all of the above.
```

By default, the script clones
[IBM-HElib/HElib](https://github.com/IBM-HElib/HElib), performs a package build
and only compiles and installs the library. To run the tests or use a lib
build, extra flags must be passed to the script.

**NOTE:** HElib is always cloned in the user's home directory.

It is possible to tell the script what repository to clone using the `-r
<repo>` option, where `<repo>` is the HTTPS link to perform the clone from.
Additionally, one can select a specific branch to build and test using the `-b
<branch>` option.

## Running the docker containers

1. Using the dockerfiles provided, first build the images by running
   ```
   docker build -t <name:tag of image> -f <path to Dockerfile> .
   ```
   For example from the `HElib/docker/build_and_test` directory, if you want to
   build the Ubuntu image
   ```
   docker build -t he-ready-ubuntu:20.04 -f Dockerfile.UBUNTU .
   ```
   will build the Ubuntu 20.04 image and tag it with the name
   `he-ready-ubuntu:20.04`.
   
2. Once the images are built you can run the containers by running
   ```
   docker run --name <name of container> <name:tag of image>
   ```
   The containers will by default run the bash script described above. For
   example using the example image above
   ```
   docker run --name default_ubuntu_test he-ready-ubuntu:20.04
   ```
   will spin a container called `default_ubuntu_test` and run internally
   `./root/build_and_test_helib.sh -a`. This will clone `IBM-HElib/HElib`,
   perform a package build and install, run the Google tests, and then build
   and run the tests for all subprojects.

   **OPTIONAL:** It is possible to override the default build of each docker
   container by specifying your own inputs during the run command by running
   ```
   docker run --name <name of container> <name:tag of image> ./root/build_and_test_helib.sh <optional arguments or flags>
   ```
   For example running this command
   ```
   docker run --name my_ubuntu_test he-ready-ubuntu:20.04 ./root/build_and_test_helib.sh -r https://github.com/homenc/HElib.git -lt
   ```
   will run an Ubuntu 20.04 container called `my_ubuntu_test` with the
   following custom options:
   - the `-r` flag will instead clone [homenc/HElib](https://github.com/homenc/HElib)
   - because the `-b` was not specified the master branch will be tested by default
   - the `-l` flag tells the script to build HElib using the library build
   - the `-t` flag tells the script to only run the HElib Google tests and nothing else

   To see the additional options available look at the usage info for the script.
