cmake_minimum_required(VERSION 3.14)
project(numerical-computing LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# 
# Library: every algorithm in src/ (one .cpp per module).
# When you add a module, add its .cpp here AND one add_module(...) line below.
#
add_library(numcomp
  src/gaussian.cpp
  src/lu.cpp
  src/jacobi.cpp
  src/gauss_seidel.cpp
  src/determinant.cpp
  src/matrix_inversion.cpp
  src/cholesky.cpp
  src/gram_schmidt.cpp
  src/qr_decomposition.cpp
  src/eigenvalues.cpp
  src/power_iteration.cpp
  src/least_squares.cpp
  src/polynomial_fitting.cpp
  src/lagrange.cpp
  src/newton_interpolation.cpp
  src/polynomial_interpolation.cpp
  src/romberg_integration.cpp
  src/gaussian_quadrature.cpp
  src/runge_kutta2.cpp
  src/Runge_Kutta_algorithm.cpp
)

# Public headers are included as <numerical-computing/...>.
# include/numerical-computing is also exposed so older bare includes such as
#   #include "romberg_integration.hpp"
# still compile. Once every file uses the <numerical-computing/...> form,
# that second directory can be removed.
target_include_directories(numcomp PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
  ${CMAKE_CURRENT_SOURCE_DIR}/include/numerical-computing
)

if(MSVC)
  target_compile_options(numcomp PRIVATE /W4)
  target_compile_definitions(numcomp PUBLIC _USE_MATH_DEFINES)  # M_PI on MSVC
else()
  target_compile_options(numcomp PRIVATE -Wall -Wextra)
endif()

# ---------------------------------------------------------------------------
# Tests and examples
# ---------------------------------------------------------------------------
enable_testing()

# add_module(<name> <test source> <example source>)
#   builds test_<name> (registered with CTest) and example_<name>.
# Paths are relative to the project root; pass an empty string to skip one.
function(add_module name test_src example_src)
  if(NOT "${test_src}" STREQUAL "")
    add_executable(test_${name} "${test_src}")
    target_link_libraries(test_${name} PRIVATE numcomp)
    # CMake's Release mode defines NDEBUG, which would silently disable any
    # assert() checks in a test. Undefine it so tests can still fail.
    if(MSVC)
      target_compile_options(test_${name} PRIVATE /UNDEBUG)
    else()
      target_compile_options(test_${name} PRIVATE -UNDEBUG)
    endif()
    add_test(NAME ${name} COMMAND test_${name})
  endif()
  if(NOT "${example_src}" STREQUAL "")
    add_executable(example_${name} "${example_src}")
    target_link_libraries(example_${name} PRIVATE numcomp)
  endif()
endfunction()

# Gaussian elimination: prefer the renamed test, fall back to the old file name
# so the build works before and after the rename.
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_gaussian.cpp")
  set(GAUSSIAN_TEST tests/test_gaussian.cpp)
else()
  set(GAUSSIAN_TEST "tests/example testing_guassian.cpp")
endif()

#          name                      test                                   example
add_module(gaussian                  "${GAUSSIAN_TEST}"                     examples/example_gaussian.cpp)
add_module(lu                        tests/test_lu.cpp                      examples/lu_example.cpp)
add_module(jacobi                    tests/test_jacobi.cpp                  examples/jacobi_example.cpp)
add_module(gauss_seidel              tests/gauss_seidel_test.cpp            examples/example_gauss_seidel.cpp)
add_module(determinant               tests/test_determinant.cpp             examples/example_determinant.cpp)
add_module(matrix_inversion          tests/test_matrix_inversion.cpp        examples/matrix_inversion_example.cpp)
add_module(cholesky                  tests/test_cholesky.cpp                examples/example_cholesky.cpp)
add_module(gram_schmidt              tests/test_gram_schmidt.cpp            examples/gram_schmidt_example.cpp)
add_module(qr_decomposition          tests/qrDecomposition_test.cpp         examples/qrDecomposition_example-1.cpp)
add_module(eigenvalues               tests/test_eigenvalues.cpp             examples/eigenvalues_example.cpp)
add_module(power_iteration           tests/test_power_iteration.cpp         examples/example_power_iteration.cpp)
add_module(least_squares             tests/test_least_squares.cpp           examples/example_least_squares.cpp)
add_module(polynomial_fitting        tests/polynomial_fitting_test.cpp      examples/polynomial_fitting_example.cpp)
add_module(lagrange                  tests/test_lagrange.cpp                examples/lagrange_example.cpp)
add_module(newton_interpolation      tests/newton_interpolation_test.cpp    examples/newton_example.cpp)
add_module(polynomial_interpolation  tests/test_polynomial_interpolation.cpp examples/polynomial_interpolation_example.cpp)
add_module(romberg_integration       tests/romberg_integration_test.cpp     examples/romberg_integration_example.cpp)
add_module(gaussian_quadrature       tests/test_gaussian_quadrature.cpp     examples/gaussian_quadrature_example.cpp)
add_module(runge_kutta2              tests/runge_kutta2_test.cpp            examples/runge_kutta2_example.cpp)
add_module(runge_kutta4              tests/Runge_Kutta_test.cpp             examples/Runge_Kutta_example.cpp)
