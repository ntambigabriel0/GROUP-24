cmake_minimum_required(VERSION 3.14)
project(numerical-computing LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# The library: add each new .cpp file from src/ here
add_library(numcomp src/gram_schmidt.cpp)
target_include_directories(numcomp PUBLIC include)

# Tests (run with: ctest --test-dir build)
enable_testing()
add_executable(testgramschmidt tests/test_gram_schmidt.cpp)
target_link_libraries(testgramschmidt numcomp)
add_test(NAME gramschmidt COMMAND testgramschmidt)

# Examples
add_executable(gramschmidtexample examples/gram_schmidt_example.cpp)
target_link_libraries(gramschmidtexample numcomp)
