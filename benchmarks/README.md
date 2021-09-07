# Some HElib benchmarks

## Prerequisites 

### Build and install HElib

Build and install HElib as described in [INSTALL.md](../INSTALL.md).

### Build and install google benchmark

Download the library https://github.com/google/benchmark from github, build and
then install it. We tested with Google Benchmark v1.5.6

```
git clone https://github.com/google/benchmark.git -b v1.5.6
cd benchmark
mkdir build
cd build
cmake -DBENCHMARK_ENABLE_GTEST_TESTS=OFF -DBENCHMARK_ENABLE_INSTALL=ON -DCMAKE_INSTALL_PREFIX=<installation dir> ..
make
make install
```

## Build benchmark

Benchmarks are found in `HElib/benchmarks` and are grouped by themes; the file
names reflect the themes they belong to. By default, the benchmarks are built as 
separate executables. If you prefer to build a single executable then set the
`-DSINGLE_EXEC=ON`. Compiled tests appear in the `bin` directory.  Now build
the project.

```
cd HElib/benchmarks
mkdir build
cd build
cmake -Dhelib_DIR=<helib installation dir>/share/cmake/helib -Dbenchmark_DIR=<benchmark installation dir>/lib/cmake/benchmark/ ..
make
```

When writing benchmarks please use `BENCHMARK_CAPTURE(<custom params>)` for
custom parameters and do not forget to add the new target to `CMakeLists.txt`.
Use the extra option `BENCHMARK_CAPTURE()->Iterations(<iterations>)` for a
specific number of iterations or `BENCHMARK_CAPTURE()->MinTime(<time in
seconds>)` to set a minimum time for the benchmark to run for.

NOTE: Both `Iterations` and `MinTime` cannot be used together.

## Run benchmark

To execute individual tests run the following

```
./bin/<benchmark to run>
```

If you have compiled the tests as a single executable using the
`-DSINGLE_EXEC=ON` flag, you can execute the tests by running the following

```
./bin/helib_benchmark
```
