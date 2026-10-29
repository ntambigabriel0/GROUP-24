
# Numerical Computing Library in C++ (GROUP-24)

## Group Members

| No. | Name | Registration Number |
|---|---|---|
| 1 | NTAMBI GABRIEL | 25/U/0855 |
| 2 | MUWONGE EMMANUEL | 25/U/08755/PS |
| 3 | KIRYA GABRIEL SSESE | 25/U/08729/PS |
| 4 | MUSIIME NUWAGABA | 25/U/08600/PS |
| 5 | CARGAN HAGAN | 25/U/08707/PS |
| 6 | NABIRYE FORTUNATE | 25/U/08759/PS |
| 7 | NAKALULE LEILA | 25/U/03499/PSA |

A reusable C++17 numerical computing library written from scratch, with no external libraries. It provides tools for solving linear systems, matrix factorisation, eigenvalues, interpolation, curve fitting, integration and ODE solving. Each algorithm lives in its own header/source pair with its own test and example. The library covers 21 algorithms, the last one added being Gaussian quadrature (Gauss-Legendre).

Everything in the library is in the `numcomp` namespace, and public headers are included as `<numerical-computing/...>`. A few modules are still being aligned to this layout; the [status table](#4-module-status) shows which.

## Table of Contents
- [Group members](#group-members)

1. [Requirements](#1-requirements)
2. [What the library does](#2-what-the-library-does)
3. [Project structure](#3-project-structure)
4. [Module status](#4-module-status)
5. [Build](#5-build)
6. [Use the library in your own program](#6-use-the-library-in-your-own-program)
7. [Examples: numerical use cases](#7-examples-numerical-use-cases)
8. [Run the tests](#8-run-the-tests)
9. [Template candidates](#9-template-candidates)
10. [Design notes](#10-design-notes)
11. [Weekly reports](#11-weekly-reports)
12. [Team workflow](#12-team-workflow)

## 1. Requirements

| Tool | Version | Purpose |
|---|---|---|
| C++ compiler (GCC, Clang, MSVC) | C++17 | Compiling the library |
| CMake | 3.14 or newer | Build system (`CMakeLists.txt`) |
| Git | any recent version | Cloning and version control |

No third-party libraries are needed.

## 2. What the library does

| Area | Algorithm | Header | Main function(s) |
|---|---|---|---|
| Linear systems | Gaussian elimination | `gaussian.hpp` | `gaussianElimination` |
| | LU decomposition (partial pivoting) | `lu.h` | `luDecompose`, `luSolve` |
| | Jacobi iteration | `jacobi.hpp` | `jacobi` |
| | Gauss-Seidel iteration | `gauss_seidel.hpp` | `GaussSeidel::solve` |
| Matrix operations | Determinant (via LU) | `determinant.hpp` | `determinant` |
| | Matrix inversion (via LU) | `matrix_inversion.hpp` | `invertMatrix` |
| Factorisation | Cholesky | `cholesky.hpp` | `choleskyDecomposition` |
| | Gram-Schmidt orthogonalisation | `gram_schmidt.hpp` | `gramSchmidt`, `dot`, `norm` |
| | QR decomposition (built on Gram-Schmidt) | `qr_decomposition.hpp` | `qrDecomposition` |
| Eigenvalues | QR algorithm | `eigenvalues.hpp` | `eigenvalues` |
| | Power iteration | `power_iteration.hpp` | `powerIteration` |
| Fitting | Least-squares solution (via QR) | `least_squares.hpp` | `leastSquares` |
| | Least-squares polynomial fitting | `polynomial_fitting.hpp` | `polynomialFit` |
| Interpolation | Lagrange | `lagrange.hpp` | `lagrangeInterpolate` |
| | Newton (divided differences) | `newton_interpolation.hpp` | `newtonCoefficients`, `newtonEvaluate` |
| | Polynomial interpolation | `polynomial_interpolation.hpp` | `polynomialInterpolation`, `evaluatePolynomial` |
| Integration | Romberg | `romberg_integration.hpp` | `rombergIntegration` |
| | Gaussian quadrature (Gauss-Legendre) | `gaussian_quadrature.hpp` | `gaussLegendreRule`, `gaussLegendre`, `compositeGaussLegendre` |
| ODE solvers | Runge-Kutta 2nd order | `runge_kutta2.hpp` | `rungeKutta2` |
| | Runge-Kutta 4th order | `Runge_Kutta.hpp` | `RK4Solver` |

Several modules build on each other, so they reuse results instead of repeating code:

```
gram_schmidt  ->  qr_decomposition  ->  least_squares  ->  polynomial_fitting
                         \
                          ->  eigenvalues
lu  ->  determinant
    ->  matrix_inversion

gaussian_quadrature   (standalone: uses only the C++ standard library)
```

## 3. Project structure

```
GROUP-24/
├── include/numerical-computing/   # Public headers (one per algorithm)
├── src/                           # Implementations (.cpp)
├── examples/                      # One demo program per module
├── tests/                         # One test program per module
├── reports/
│   ├── week-01/                   # Gaussian, Gram-Schmidt, LU
│   ├── week-02/                   # Determinant, matrix inversion, polynomial interpolation, QR
│   ├── week-03/                   # Cholesky, eigenvalues, Jacobi, Lagrange, least squares,
│   │                              # polynomial fitting, power iteration, Romberg, Runge-Kutta 2
│   ├── week-04/                   # Gauss-Seidel, Newton interpolation
│   └── week-05/                   # Gaussian quadrature, inheritance/polymorphism/encapsulation report
├── .github/CODEOWNERS
├── CMakeLists.txt
├── LICENSE
└── README.md
```

## 4. Module status

This table reflects what was checked by compiling each module's example and test with `g++ -std=c++17`.

| Module | Builds | Tests | Notes |
|---|---|---|---|
| Gram-Schmidt | Yes | Pass | Also built by `CMakeLists.txt` |
| LU decomposition | Yes | Pass | |
| Determinant | Yes | Pass | Needs `src/lu.cpp` |
| Matrix inversion | Yes | Pass | Needs `src/lu.cpp` |
| Cholesky | Yes | Pass | |
| QR decomposition | Yes | Runs | Needs `src/gram_schmidt.cpp`; test prints the Q and R matrices |
| Least squares | Yes | Pass (6/6) | Needs QR and Gram-Schmidt sources |
| Polynomial fitting | Yes | Pass | Needs least squares, QR and Gram-Schmidt sources |
| Lagrange interpolation | Yes | Pass | |
| Newton interpolation | Yes | Pass | |
| Gaussian quadrature | Yes | Pass (35/35) | Standalone: needs only `src/gaussian_quadrature.cpp`. Also builds with no warnings under `-Wall -Wextra -Wpedantic` |
| Jacobi iteration | Yes | Pass | Test also needs `src/lu.cpp` |
| Eigenvalues | Example builds | Test does not compile | `tests/test_eigenvalues.cpp` line 46 has a stray `1` after a statement. Needs QR and Gram-Schmidt sources |
| Gaussian elimination | Not yet | Not yet | `src/gaussian.cpp` has its own `main()`, which clashes with the example and test |
| Gauss-Seidel | Not yet | Not yet | `src/gauss_seidel.cpp` is a standalone program with `main()`, not an implementation of `GaussSeidel::solve` |
| Polynomial interpolation | Not yet | Not yet | Includes `gaussian.hpp`, which contains two copies of the header (stray `#endif`) |
| Power iteration | Example builds | Test does not compile | The test includes `power_iteration_example.h`, which is not in the repo |
| Romberg integration | Not yet | Not yet | Source, example and test include `romberg.hpp`, but the header is `romberg_integration.hpp` |
| Runge-Kutta 2 | Not yet | Not yet | Source, example and test include `rk2.hpp`, but the header is `runge_kutta2.hpp` |
| Runge-Kutta 4 | Not yet | Not yet | Source and tests include `algorithm.hpp` with curly quotes, and the header uses `std::functional` instead of `std::function` |

## 5. Build

### With CMake

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

`CMakeLists.txt` currently registers only Gram-Schmidt. It builds the `numcomp` library from `src/gram_schmidt.cpp`, the `testgramschmidt` test and the `gramschmidtexample` program. Each new module is added by listing its `.cpp` file in `add_library(numcomp ...)` and adding its test and example the same way.

### Building another module directly

Compile a module's example together with its source and the sources it depends on. From the repository root:

```bash
# LU decomposition
g++ -std=c++17 -Iinclude examples/lu_example.cpp src/lu.cpp -o lu_example

# Least squares (needs QR and Gram-Schmidt)
g++ -std=c++17 -Iinclude examples/example_least_squares.cpp \
    src/least_squares.cpp src/qr_decomposition.cpp src/gram_schmidt.cpp -o least_squares_example

# Polynomial fitting (needs least squares, QR and Gram-Schmidt)
g++ -std=c++17 -Iinclude examples/polynomial_fitting_example.cpp src/polynomial_fitting.cpp \
    src/least_squares.cpp src/qr_decomposition.cpp src/gram_schmidt.cpp -o polyfit_example

# Gaussian quadrature (standalone, no other sources needed)
g++ -std=c++17 -Iinclude examples/gaussian_quadrature_example.cpp src/gaussian_quadrature.cpp -o gq_example
```

Source dependencies per module:

| Module | Sources to compile together |
|---|---|
| `determinant`, `matrix_inversion` | module + `lu.cpp` |
| `qr_decomposition` | module + `gram_schmidt.cpp` |
| `least_squares` | module + `qr_decomposition.cpp` + `gram_schmidt.cpp` |
| `polynomial_fitting` | module + `least_squares.cpp` + `qr_decomposition.cpp` + `gram_schmidt.cpp` |
| `eigenvalues` | module + `qr_decomposition.cpp` + `gram_schmidt.cpp` |
| `gaussian_quadrature` | module only |

## 6. Use the library in your own program

Include the header of the module you need and call its function inside the `numcomp` namespace:

```cpp
#include <numerical-computing/lu.h>
#include <iostream>

int main() {
    numcomp::Matrix A = {{2, 1}, {1, 3}};
    numcomp::LUResult f = numcomp::luDecompose(A);   // A = L * U with row swaps recorded
    numcomp::Vector x = numcomp::luSolve(f, {3, 5}); // solve A x = b
    std::cout << "x = " << x[0] << ", y = " << x[1] << "\n";
}
```

Most modules share the same type names:

```cpp
using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;
```

Invalid input is reported with exceptions from the `std::invalid_argument` family, for example `SingularMatrixError` in `lu.h`.

## 7. Examples: numerical use cases

The outputs below were produced by running the programs in `examples/`.

### Solving a linear system with LU decomposition (`lu_example.cpp`)
Factor once, then reuse the factors for several right-hand sides.
```
L:
1 0
0.5 1
U:
2 1
0 2.5
swaps: 0
x = 0.8, y = 1.4
second solve: x = 0.6, y = -0.2
Error caught: luDecompose: zero pivot in column 1, matrix is singular
```

### Determinant and inverse (`example_determinant.cpp`, `matrix_inversion_example.cpp`)
For A = [[2, 1], [1, 3]]:
```
Determinant = 5

Inverse Matrix A^-1:
0.6   -0.2
-0.2   0.4
```

### Cholesky factorisation (`example_cholesky.cpp`)
For the symmetric positive-definite matrix [[4, 12, -16], [12, 37, -43], [-16, -43, 98]]:
```
Cholesky factor L:
 2  0  0
 6  1  0
-8  5  3
```

### Orthogonalising vectors with Gram-Schmidt (`gram_schmidt_example.cpp`)
```
Orthonormal vectors:
0.707107 0.707107 0
0.408248 -0.408248 0.816497
-0.57735 0.57735 0.57735
q0 . q1 = 1.11022e-16
Error caught: gramSchmidt: vectors are linearly dependent
```

### Fitting a line by least squares (`example_least_squares.cpp`)
```
a = 1.167
b = 0.500
Best-fit line: y = 1.167 + 0.500x
```

### Polynomial fitting (`polynomial_fitting_example.cpp`)
Data points (0, 1), (1, 3), (2, 5), (3, 8):
```
[Degree 1 - Line of Best Fit]
y = 0.8000 + 2.3000x

[Degree 2 - Best Fit Parabola]
y = 1.0500 + 1.5500x + 0.2500x^2
```

### Interpolation (`lagrange_example.cpp`, `newton_example.cpp`)
Newton interpolation through (0, 1), (1, 3), (2, 7):
```
Divided Differences (Coefficients):
  c0 = 1
  c1 = 2
  c2 = 1

P(3.0000) = 13.0000
```
Lagrange interpolation evaluates the same kind of data at one point or many:
```
y at x = 3 is 13  (expected 13)
x = 0.5 -> y = 1.75
x = 1.5 -> y = 4.75
x = 2.5 -> y = 9.75
```

### Numerical integration with Gaussian quadrature (`gaussian_quadrature_example.cpp`)
An n-point Gauss-Legendre rule integrates every polynomial of degree 2n - 1 or lower exactly. Selected lines from the output:
```
3-point rule on [-1, 1]:
  node -0.774596669241   weight 0.555555555556
  node 0   weight 0.888888888889
  node 0.774596669241   weight 0.555555555556

Integral of sin(x) from 0 to pi (exact value 2):
  n = 2  result = 1.93581957465  error = 0.0641804253489
  n = 5  result = 2.00000011028  error = 1.10284472221e-07
  n = 8  result = 2  error = 4.21884749358e-15

3 points on x^k over [0, 1]:
  x^5: 0.166666666667  exact 0.166666666667  (degree <= 5, exact)
  x^6: 0.1425  exact 0.142857142857  (degree 6 > 5, not exact)

Normal distribution P(-1 <= Z <= 1), 10 points: 0.682689492137
Work to stretch a spring (F = 200x, 0 to 0.3 m), 2 points: 9 J

Integral of sin(x) from 0 to 20:
  one 5-point rule:         error = 4.52018893146
  10 panels, 5 points each: error = 2.78347456195e-10
Error caught: gaussLegendre: number of points must be between 1 and 100
```

### Iterative solving with Jacobi (`jacobi_example.cpp`)
```
First three rounds from (0,0):
  round 1: (1.5, 2.4)
  round 2: (0.9, 1.8)
  round 3: (1.05, 2.04)
Solution: x = 1, y = 2 after 22 rounds (converged: yes)
Not diagonally dominant: converged: no
```

### Eigenvalues with the QR algorithm (`eigenvalues_example.cpp`)
```
[[2,1],[1,2]] (expect 3 and 1):
  eigenvalues: 3 1
  rounds: 22, converged: yes
rotation matrix (complex eigenvalues, expect converged: no):
  eigenvalues: 0 0
  rounds: 50, converged: no
```

## 8. Run the tests

Each module has a test program in `tests/` that prints one `PASS` or `FAIL` line per check and ends with a summary. Build a test the same way as an example, for instance:

```bash
g++ -std=c++17 -Iinclude tests/test_lu.cpp src/lu.cpp -o test_lu && ./test_lu
```

The tests cover normal use, edge cases and invalid inputs. For example, the LU tests check that P·A equals L·U, that solutions have a near-zero residual, and that singular matrices and wrong-sized right-hand sides throw.

| Test file | Module |
|---|---|
| `test_gram_schmidt.cpp` | Gram-Schmidt |
| `test_lu.cpp` | LU decomposition |
| `test_determinant.cpp` | Determinant |
| `test_matrix_inversion.cpp` | Matrix inversion |
| `test_cholesky.cpp` | Cholesky |
| `qrDecomposition_test.cpp` | QR decomposition |
| `test_least_squares.cpp` | Least squares |
| `polynomial_fitting_test.cpp` | Polynomial fitting |
| `test_lagrange.cpp` | Lagrange interpolation |
| `newton_interpolation_test.cpp` | Newton interpolation |
| `test_gaussian_quadrature.cpp` | Gaussian quadrature (35 checks) |
| `test_jacobi.cpp` | Jacobi iteration |
| `test_eigenvalues.cpp` | Eigenvalues |
| `test_power_iteration.cpp` | Power iteration |
| `test_polynomial_interpolation.cpp` | Polynomial interpolation |
| `gauss_seidel_test.cpp` | Gauss-Seidel |
| `romberg_integration_test.cpp` | Romberg integration |
| `runge_kutta2_test.cpp` | Runge-Kutta 2 |
| `Runge_Kutta_test.cpp` | Runge-Kutta 4 |
| `example testing_guassian.cpp` | Gaussian elimination |

## 9. Template candidates

The library headers currently use no C++ templates. The only templates in the repository are the `template <typename F>` check helpers in the test files (`test_eigenvalues.cpp`, `test_least_squares.cpp`, `test_lagrange.cpp`, `newton_interpolation_test.cpp`, `test_matrix_inversion.cpp`, `test_gram_schmidt.cpp`, `test_jacobi.cpp`, `polynomial_fitting_test.cpp`, `test_lu.cpp`, `test_gaussian_quadrature.cpp`), which accept any callable.

Good candidates for templating in the library:

| Where | Current code | Candidate |
|---|---|---|
| Type aliases in `lu.h`, `jacobi.hpp`, `determinant.hpp`, `eigenvalues.hpp`, `gram_schmidt.hpp` | `using Vector = std::vector<double>;` | An alias template `template <typename T> using Vector = std::vector<T>;` so the algorithms work with `float` or `long double` |
| `Matrix` class in `gaussian.hpp` | `std::vector<std::vector<double>> data;` | `template <typename T> class Matrix` |
| Result structs `LUResult`, `QRResult`, `JacobiResult`, `EigenResult`, `QuadratureRule` | Fields of type `Matrix`, `Vector` and `double` | `template <typename T>` structs |
| Free functions such as `luDecompose`, `choleskyDecomposition`, `determinant`, `gramSchmidt`, `dot`, `norm`, `lagrangeInterpolate`, `newtonCoefficients`, `polynomialFit`, `gaussLegendreRule` | Parameters and return values fixed to `double` | Function templates over the number type |
| Default tolerances | `1e-12` in `lu.h`, `cholesky.hpp`, `determinant.hpp`; `1e-10` in `gram_schmidt.hpp`, `jacobi.hpp`, `eigenvalues.hpp`; `1e-9` in `qr_decomposition.hpp`, `least_squares.hpp` | Variable templates such as `template <typename T> constexpr T defaultTol = T(1e-12);` |
| `rombergIntegration`, `gaussLegendre`, `compositeGaussLegendre` and `rungeKutta2` | `std::function<double(double)>` and `std::function<double(double, double)>` | A template parameter for the callable, as the test helpers already do |
| `powerIteration` | Fixed `double A[20][20]` and `double x[20]` arrays | A size and type template parameter, or `std::vector` |

## 10. Design notes

- **Separation of interface and logic:** declarations are in `include/numerical-computing/*.hpp`, implementations in `src/*.cpp`.
- **Layered modules:** `qrDecomposition` calls `gramSchmidt` and `dot`; `leastSquares` calls `qrDecomposition`; `polynomialFit` calls `leastSquares`; `determinant` and `invertMatrix` use `luDecompose`.
- **Rules as data:** `gaussLegendreRule` returns a `QuadratureRule` (nodes and weights); `compositeGaussLegendre` builds it once and reuses it on every panel.
- **Result structs:** `LUResult`, `QRResult`, `JacobiResult` and `EigenResult` bundle several outputs (for example L, U, the permutation and the swap count) into one return value.
- **Floating-point care:** pivots and rank checks are compared against a tolerance parameter with a default, not against exact zero.
- **Error handling:** invalid input throws `std::invalid_argument` or a derived class such as `SingularMatrixError`.
- **Object-oriented features:** the library is mostly free functions in the `numcomp` namespace. The classes are `Matrix` (`gaussian.hpp`) and `RK4Solver` (`Runge_Kutta.hpp`), with private data behind public member functions, and `SingularMatrixError` (`lu.h`), the one derived class, which inherits from `std::invalid_argument`. Overloading, `std::function` parameters and helper hiding are covered in `reports/week-05/inheritance_polymorphism_encapsulation_report.md`.
- **Known inconsistencies:** some modules still use global functions or other include names (see the [status table](#4-module-status)); the aim is for all of them to use the `numcomp` namespace and the `<numerical-computing/...>` include path.

## 11. Weekly reports

Progress reports, one per module, are in `reports/week-01` to `reports/week-05`. They describe how each module was built, the issues faced and the design decisions.

| Week | Reports |
|---|---|
| `week-01` | Gaussian elimination, Gram-Schmidt, LU decomposition |
| `week-02` | Determinant, matrix inversion, polynomial interpolation, QR decomposition |
| `week-03` | Cholesky, eigenvalues, Jacobi, Lagrange, least squares, polynomial fitting, power iteration, Romberg, Runge-Kutta 2 |
| `week-04` | Gauss-Seidel, Newton interpolation, Gaussian quadrature |


## 12. Team workflow

- `main` is protected by a branch ruleset (`protect-main`) that requires a pull request with one approval before merging.
- Work is done on a feature branch per module and merged by pull request.
- `.github/CODEOWNERS` assigns review ownership.
- The project is released under the licence in `LICENSE`.
