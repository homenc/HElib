# Building and Running the lookup/scoring tests

## Prequisites

- You must have bats-core installed.
- Build and install HElib.
- Build the utilities directory.

## Build the tests

1. In the directory `<path-to-HElib>/HElib/misc/psi` run
```
mkdir build
cd build
cmake -Dhelib_DIR=<absolute-path-to-HElib-install-dir> ..
make -j3
```

## Run the tests

Make sure you are in the directory `<path-to-HElib>/HElib/misc/psi/test`
directory.

Make sure you set which tests you want to run by commenting out the skip in
`lookup.bats` and `scoring.bats`.

Pay particular attention to setting the correct `nslots`, `modulus`, `p`, and
`m` in the files `lookup.bats`, `scoring.bats`, and `gen-params.batx`. 

1. Create the parameters
```
./gen-params.batx -f "generate params"
```
This places the parameters, context, keys and encrypted database/query in the
directory `data_and_params`.

2. Run the tests
```
DEBUG=1 bats .
```
You should see a lot of tmp directories, one for each test.

3. View the results
```
./report.sh <number-of-threads-you-want-to-view>
```
Note: It scans the directories for the results you have asked for but you will
be prompted to press enter to see the next result.
