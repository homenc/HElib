# Building and running HElib's benchmarks

## Instructions

### Build and install HElib

Build and install HElib as in [INSTALL.md](../INSTALL.md).

### Build and install google benchmark

Download the library https://github.com/google/benchmark from github, build and
then install it.
```
git clone https://github.com/google/benchmark
cd benchmark
mkdir build
cd build
cmake -DBENCHMARK_ENABLE_GTEST_TESTS=OFF -DBENCHMARK_ENABLE_INSTALL=ON -DCMAKE_INSTALL_PREFIX=<installation dir> ..
make
make install
```

### Write a benchmark test and link to both installed libraries.

Write the benchmark in an existing source file or a new one in
`HElib/benchmarks`, adding the new source file to `CMakeLists.txt` if
appropriate, then build the project.
```
cd HElib/benchmarks
mkdir build
cd build
cmake -Dhelib_DIR=<helib installation dir>/share/cmake/helib -Dbenchmark_DIR=<benchmark installation dir>/lib/cmake/benchmark/ ..
make
```

When writing benchmarks use `BENCHMARK_CAPTURE(custom params)` for custom
params and the extra option `BENCHMARK_CAPTURE()->Iterations(iterations)` for a
specific number of iterations or `BENCHMARK_CAPTURE()->MinTime(time in
seconds)` to set a minimum time for the benchmark to run for.

NOTE: Both options cannot be used at the same time.
