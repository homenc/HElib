cmake_minimum_required(VERSION 3.2)
project(fhe)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(NTL_HEADER /usr/local/include) # might modify this path according to your environment
include_directories(${NTL_HEADER})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
set(FHE_THREADS ON) # make sure the NTL is built with NTL_THREADS=on
add_subdirectory(src)
