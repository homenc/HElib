#!/bin/bash
set -e # Fail on error
run='[run.sh]:'

echo $run make HElib
make

echo $run compile test prog
g++ -v  -g -O2 -Wfatal-errors -Wshadow -Wall -I/Users/hamish/HElibBin/local/include -std=c++11  -o Test_Bin_IO_x Test_Bin_IO.cpp fhe.a -L/Users/hamish/HElibBin/local/lib -lntl -lgmp  -lm

echo $run Exec test prog
./Test_Bin_IO_x

echo $run Compare files
diff iotest_ascii1.txt iotest_ascii2.txt 
