cmake_minimum_required(VERSION 3.2)
project(fhe)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(ENABLE_TEST "Build the test files" OFF)
option(ENABLE_THREADS "To use multi-thread. Make Sure the NTL is built with NTL_THREADS=on" OFF)
option(BUILD_SHARED "To build shared objects" OFF)
option(BUILD_AES "To build homomorphic AES" OFF)

## Might modify the following path to NTL according to your environment
set(NTL_HEADER /usr/local/include) 
set(NTL_LIB /usr/local/lib) 
include_directories(${NTL_HEADER})
link_directories(${NTL_LIB})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
## Use -std=c++11 as dedeault. 
set(CMAKE_CXX_STANDARD 11)
## On old platform, might change to -std=c++0x
#set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

## To use multi-threads.
## Make sure the NTL is built with NTL_THREADS=on
## For NTL's verison older than 7.0, should turn FHE_THREADS off.
if (ENABLE_THREADS)
set(FHE_THREADS ON) 
add_definitions(-DFHE_THREADS)
set(RUN_LIB fhe ntl gmp pthread)
else (ENABLE_THREADS)
set(FHE_THREADS OFF) 
set(RUN_LIB fhe ntl gmp)
endif (ENABLE_THREADS)

add_subdirectory(src)
