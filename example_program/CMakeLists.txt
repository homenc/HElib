cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

## Use -std=c++11 as default.
set(CMAKE_CXX_STANDARD 11)
## Disable C++ extensions
set(CMAKE_CXX_EXTENSIONS OFF)
## Require full C++ standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(helib_example
       	LANGUAGES CXX)

find_package(helib 1.0.0 EXACT REQUIRED)

add_executable(helib_example helib_example.cpp)

target_link_libraries(helib_example helib)

