# Building and installing HElib

The current HElib build, install, and regression tests suite have been built
and tested on Ubuntu 22.04, macOS Montrey >=12.6.7, and macOS Ventura >=13.4.1.
Previous HElib versions included build, install, and regression tests on 
Ubuntu 18.04, ubuntu 20.04, Fedora 33, CentOS 8.2, macOS Mojave >=10.14.6, 
macOS Catalina >=10.15.7, and macOS Big Sur >=11.7.8

There are two different ways to build and install HElib. The first one will
automatically download and build the GMP and NTL dependencies and pack the
libraries in a relocatable folder. The second way, instead, requires the
dependencies to be installed by you and available in the system.

This release of HElib has experimental support for the Intel® [HEXL](https://github.com/intel/hexl) 
acceleration library for homomorphic encryption that exploits Intel® Advanced Vector Extensions 512. 
Instructions to enable and link to HEXL are given 
[below](#enabling-and-linking-to-intel-hexl).

```diff
- Please read these instructions in full to better choose the type of build that is best for you.
```

## General prerequisites

- pthreads
- git >= 2.36 (required to build and run the HElib test suite)

**Default Linux environment:**

- Ubuntu 22.04 LTS
- GNU make >= 4.3
- g++ >= 11.3.0
- cmake >= 3.22

**macOS environment:**

- Apple clang >= 14.0.0 (available with the latest Xcode for the tested versions of macOS)
- Xcode Command Line Tools (can be installed with the command `xcode-select
  --install` in a terminal)
- cmake >= 3.22 (available from [CMake](https://cmake.org/) or [MacPorts
  Project](https://www.macports.org/) and [Homebrew](https://brew.sh/) as
  packages)
- GNU make >= 4.3

**For HElib development:**

- clang-format >= 14.0.0 (available with your linux distribution and for macOS
  from [MacPorts Project](https://www.macports.org/) and
  [Homebrew](https://brew.sh/) as packages)

## Option 1: package build (recommended for most users)

This option bundles HElib and its dependencies (NTL and GMP) in one directory
which can then be moved around freely on the system.  NTL and GMP will be
automatically fetched and compiled.  It can be installed globally (i.e. under
`/usr/local`), which is the default option if no `CMAKE_INSTALL_PREFIX` is
specified, but this should only be done with caution as existing versions of
NTL, GMP, or HElib will be overwritten.  These additional two prerequisites are
required in this case:

- m4 >= 1.4.18
- patchelf >= 0.14.3 (if building on Linux)

Please note that if changing from library build to package build, it is safer to
use a clean build directory.

### Instructions

1. Create a build directory, typically as a sibling of `src`:

```bash
cd HElib
mkdir build
cd build
```

2. Run the cmake configuration step, specifying that you want a package build
   (via -DPACKAGE_BUILD=ON) and saying where you would like the installation to
   be. To install in `/home/alice/helib_install`, for example:

```bash
cmake -DPACKAGE_BUILD=ON -DCMAKE_INSTALL_PREFIX=/home/alice/helib_install ..
```

Extra options can be specified here, such as enabling HElib tests with
`-DENABLE_TEST=ON`.  See later section entitled "HElib build options" for
details.

3. Compile, with an optional number of threads specified (16 in this example).
   The output of this will be in the relocatable folder `helib_pack`:

```bash
make -j16
```

4. (optional) If step 2 was performed with `-DENABLE_TEST=ON`, HElib tests can
   be run as follows:

```bash
ctest
```

Detailed HElib-specific test logs can be found in
`Testing/Temporary/LastTest.log`.

5. (optional) Run the install step, to copy the folder `helib_pack` to
   `${CMAKE_INSTALL_PREFIX}` (in this example `/home/alice/helib_install`):

```bash
make install
```

of course, if the `CMAKE_INSTALL_PREFIX` was kept as the default `/usr/local` or
some other system-wide path, step 5 may require `sudo` privileges.

## Option 2: library build (advanced)

This option involves building HElib on its own, linking against pre-existing
dependencies (NTL and GMP) on the system.  In this way, the HElib library can be
moved around, but its dependencies (NTL and GMP) cannot, as they are absolute
paths.  For this option, you must build GMP >=6.2.1 and NTL >=11.5.1 yourself.
For details on how to do this, please see the section on building dependencies
later.  It is assumed throughout this installation option that the environment
variables `$GMPDIR` and `$NTLDIR` are set to point to the installation
directories of GMP and NTL respectively.

Please note that if changing from package build to library build, it is safer to
use a clean build directory.

1. Create a build directory, typically as a sibling of `src`:

```bash
cd HElib
mkdir build
cd build
```

2. Run the cmake configuration step, specifying where to find NTL and GMP.  If
   not specified, system-wide locations such as `/usr/local/lib` will be
   searched. To install in `/home/alice/helib_install`, for example:

```bash
cmake -DGMP_DIR="${GMPDIR}" -DNTL_DIR="${NTLDIR}" -DCMAKE_INSTALL_PREFIX=/home/alice/helib_install ..
```

Extra options can be specified here, such as enabling HElib tests with
`-DENABLE_TEST=ON`.  See later section entitled "HElib build options" for
details.

3. Compile, with an optional number of threads specified (16 in this example):

```bash
make -j16
```

4. (optional) If step 2 was performed with `-DENABLE_TEST=ON`, tests can be run
   as follows:

```bash
ctest
```

Detailed HElib test logs can be found in `Testing/Temporary/LastTest.log`.

5. Run the install step, to copy the files to `${CMAKE_INSTALL_PREFIX}` (in this
   example `/home/alice/helib_install`):

```bash
make install
```

of course, if the `CMAKE_INSTALL_PREFIX` was kept as the default `/usr/local` or
some other system-wide path, step 5 may require `sudo` privileges.

## Building dependencies (for option 2)

### GMP

Many distributions come with GMP pre-installed. If not, you can install GMP as
follows.

1. Download GMP from [http://www.gmplib.org](http://www.gmplib.org) -- make sure
   that you get GMP >=6.2.0 (current version is 6.2.1).
2. Decompress and cd into the gmp directory (e.g., `gmp-6.2.1`).
3. GMP is compiled in the standard unix way:

```bash
      ./configure
      make
      sudo make install
```

This will install GMP into `/usr/local` by default.

**NOTE:** For further options when building GMP, run `./configure --help` in
step 3.

### NTL

You can install NTL as follows:

1. Download NTL >=11.5.1 from
   [https://libntl.org/download.html](https://libntl.org/download.html)
2. Decompress and cd into the directory, e.g., `ntl-11.5.1/src`
3. NTL is configured, built and installed in the standard Unix way (but remember
   to specify the following flags to `configure`):

```bash
      ./configure NTL_GMP_LIP=on SHARED=on  NTL_THREADS=on NTL_THREAD_BOOST=on
      make
      sudo make install
```

This should install NTL into `/usr/local`.

**NOTE:** For further options when building NTL, run `./configure --help` in
step 3.

**NOTE**: if linking against a non-system GMP, pass `GMP_PREFIX=<path/to/gmp>`
to the `./configure` step.

## Enabling and linking to Intel® HEXL
**NOTE:** HElib with HEXL acceleration is only supported on the processors with AVX512DQ and 
AVX512-IFMA such as the 3rd generation Intel® Xeon® or the 11th generation Intel® Core®

**NOTE:** It is currently only possible to use HEXL with HElib when using the
library build and when building HElib as a static library. i.e.
`-DPACKAGE_BUILD=OFF` and `-DBUILD_SHARED=OFF`.

First you must download and build HEXL from source.  Currently, HElib only
works with HEXL version >= 1.2.1 Using git this would be

```bash
git clone https://github.com/intel/hexl --branch 1.2.1
```
Follow the instructions for HEXL installation in the README.md for all
available options.  Note previous versions of HEXL requires the deprecated
`-DENABLE_EXPORT=ON` otherwise the cmake metadata for linking a cmake project
is not created. Modern versions do not have this flag and the metadata is
created by default.  For a quick start most people will want,

```bash
cd hexl
cmake -S . -B build/ [-DCMAKE_INSTALL_PREFIX=<install-location-for-HEXL>/lib/cmake/hexl-1.2.1] 
cmake --build build -j [<parallel-jobs>]
cmake --install build
```
If you do not provide an optional install location for HEXL the default is
`/usr/local`.

To enable and link HEXL in HElib, you must configure cmake for HElib with
`-DUSE_INTEL_HEXL=ON`. If HEXL is not in the default system location or you
wish to use another installation tell HElib where to find it using
`-DHEXL_DIR=<install-location-for-hexl>`.
There is no requirement to provide any HElib subprojects with the location of
HEXL.

## HElib build options

### Generic options

- `BUILD_SHARED=ON/OFF` (default is `OFF`): Build as a shared library. Note that
  building HElib (regardless of `BUILD_SHARED`) will fail if NTL is not built as
  a shared library. The default for NTL is static library, to build NTL as a
  shared library use `./configure SHARED=on` in step 1. 
- `CMAKE_BUILD_TYPE`: (default is `RelWithDebInfo`): Choose the type of build,
  options are: `Debug`, `RelWithDebInfo`, `Release`, `MinSizeRel`.
- `CMAKE_INSTALL_PREFIX`: Desired installation directory for HElib.
- `ENABLE_TEST=ON/OFF` (default is `OFF`): Enable building of tests. This will
  include an automatic download step for the google test framework stable
  release (googletest v1.10.0)
- `ENABLE_THREADS=ON/OFF` (default is `ON`): Enable threading support. This must
  be on if and only if NTL was built with `NTL_THREADS=ON`.
- `PEDANTIC_BUILD=ON/OFF` (default is `ON`): Use `-Wall -Wpedantic -Wextra
  -Werror` during build.
- `HELIB_DEBUG=ON/OFF` (default is `OFF`): Activate the debug module when
  building HElib (by defining the `HELIB_DEBUG` macro). When the debug module is
  active, this generates extra information used for debugging purposes.
  `HELIB_DEBUG` will propagate to programs using HElib, when using cmake. When
  this is enabled, programs using HElib will generate a warning during
  configuration.  This is to remind the user that use of the debug module can
  cause issues, such as `sigsegv`, if initialized incorrectly.

### Parameters specific to option 1 (package build)

- `PACKAGE_DIR`: Location that a package build will be installed to.  Defaults
  to `${CMAKE_INSTALL_PREFIX}/helib_pack`.
- `FETCH_GMP`: Whether or not to fetch and build GMP.  Defaults to `ON`.  If set
  to `OFF`, there should either exist a system-installed GMP library, or
  `GMP_DIR` should point to a valid GMP prefix.
- `GMP_DIR`: Prefix of the GMP library.  Ignored if `FETCH_GMP=ON`.

### Parameters specific to option 2 (library build)

- `GMP_DIR`: Prefix of the GMP library.
- `NTL_DIR`: Prefix of the NTL library.
- `USE_INTEL_HEXL`: Enable the Intel HEXL library.
- `HEXL_DIR`: Prefix of the Intel HEXL library.

# Using HElib in a project

## Standard method

After `make install` has been run in either option 1 or option 2, one can find
the required shared library files to link against in `lib` and the header files
in `include`.  These can be used in the preferred way with your build system of
choice.

## Package build with cmake

Another, easier way is possible if you are using HElib in a cmake project.

1. Include the following line in your `CMakeLists.txt`:

```cmake
find_package(helib)
```

2. Run your `cmake` step with 
  `-Dhelib_DIR=<helib install prefix>/share/cmake/helib`.

## Example

Full working examples of cmake-based projects which uses HElib can be found in
the `examples` directory.
