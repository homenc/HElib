# Building and Running the lookup/scoring tests

## Prerequisites

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

Note that the test data is not provided. Data can be generated and encoded
using the utilities [here](../../utils/).

Due to how the tests have been written, please ensure both the database and
query have at least four columns. For example, a database of 2 rows and 4
columns, where the data can be seen in bold

[comment]: <> (Use these formatting and HTML tags instead of the conventional)
[comment]: <> (code block to allow formatted text within what looks like a)
[comment]: <> (code block)

<pre>
2 4
{"HElibVersion":"2.2.0","content":{"scheme":"BGV","slots":<b>[[0,1],[2,3],[4,5]]</b>},"serializationVersion":"0.0.1","type":"Ptxt"}
{"HElibVersion":"2.2.0","content":{"scheme":"BGV","slots":<b>[[0],[0],[0]]</b>},"serializationVersion":"0.0.1","type":"Ptxt"}
{"HElibVersion":"2.2.0","content":{"scheme":"BGV","slots":<b>[[6,7],[8,9],[10,11]]</b>},"serializationVersion":"0.0.1","type":"Ptxt"}
{"HElibVersion":"2.2.0","content":{"scheme":"BGV","slots":<b>[[1],[1],[1]]</b>},"serializationVersion":"0.0.1","type":"Ptxt"}
{"HElibVersion":"2.2.0","content":{"scheme":"BGV","slots":<b>[[12,13],[14,15],[16,17]]</b>},"serializationVersion":"0.0.1","type":"Ptxt"}
{"HElibVersion":"2.2.0","content":{"scheme":"BGV","slots":<b>[[2],[2],[2]]</b>},"serializationVersion":"0.0.1","type":"Ptxt"}
{"HElibVersion":"2.2.0","content":{"scheme":"BGV","slots":<b>[[18,19],[0],[0]]</b>},"serializationVersion":"0.0.1","type":"Ptxt"}
{"HElibVersion":"2.2.0","content":{"scheme":"BGV","slots":<b>[[3],[3],[3]]</b>},"serializationVersion":"0.0.1","type":"Ptxt"}
</pre>
and an example query
<pre>
1 4
{"HElibVersion":"2.2.0","content":{"scheme":"BGV","slots":<b>[[0,1],[2,3],[4,5]]</b>},"serializationVersion":"0.0.1","type":"Ptxt"}
{"HElibVersion":"2.2.0","content":{"scheme":"BGV","slots":<b>[[6,7],[8,9],[10,11]]</b>},"serializationVersion":"0.0.1","type":"Ptxt"}
{"HElibVersion":"2.2.0","content":{"scheme":"BGV","slots":<b>[[12,13],[14,15],[16,17]]</b>},"serializationVersion":"0.0.1","type":"Ptxt"}
{"HElibVersion":"2.2.0","content":{"scheme":"BGV","slots":<b>[[18,19],[0],[0]]</b>},"serializationVersion":"0.0.1","type":"Ptxt"}
</pre>

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
You should see a lot of `tmp` directories, one for each test.

3. View the results
```
./report.sh <number-of-threads-you-want-to-view>
```
Note: It scans the directories for the results you have asked for but you will
be prompted to press enter to see the next result.
