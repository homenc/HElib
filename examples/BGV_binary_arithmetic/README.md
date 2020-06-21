# To run BGV_binary_arithmetic example
There are two ways to build and run this example.
1. build and run the example on a local machine
2. build and run the example in a docker container

## build and run the example on a local machine
1. assume the HElib is installed
2. run the following command

        $ mkdir build && cd build
        $ cmake -DCMAKE_PREFIX_PATH="helib installed root" ..

## build and run the example in a docker container
1. go to the root folder of HElib
2. run `docker build -f examples/BGV_binary_arithmetic/Dockerfile .`
