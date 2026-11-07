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

### Algorithm allocation

| Member | Algorithms |
|---|---|
| Ntambi | Gram-Schmidt Orthogonalization (Tier 1), Lagrange Interpolation, |
| Muwonge | Matrix Inversion (Tier 2), Least-Squares Polynomial Fitting (Tier 4), Newton Interpolation |
| Kirya | Polynomial Interpolation (Tier 2), Romberg Integration, Runge-Kutta 2nd Order (RK2) |
| Leila | Gaussian Elimination (Tier 1), Gauss-Seidel Iteration, Power Iteration |
| Cargan | LU Decomposition (Tier 1), Eigenvalue Computation (Tier 3), Jacobi Iteration |
| Fortunate | Determinant Computation (Tier 2), Least-Squares Solution (Tier 3), Cholesky Decomposition |
| Musiime | QR Decomposition (Tier 2), Runge-Kutta 4th Order (RK4) |

Total: 6 members x 3 + 1 member x 2 = 20 algorithms.

A reusable C++17 numerical computing library written from scratch, with no external libraries. It provides tools for solving linear systems, matrix factorisation, eigenvalues, interpolation, curve fitting, integration and ODE solving. Each algorithm lives in its own header/source pair with its own test and example. The library covers 20 algorithms, the last one added being Gaussian quadrature (Gauss-Legendre).

Everything in the library is in the `numcomp` namespace, and public headers are included as `<numerical-computing/...>`. A few modules still use the older layout (global functions or a different include name); the [status table](#4-module-status) lists them.

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
13. [Limitations](#13-limitations)
14. [Commit history: proof of contribution](#14-commit-history-proof-of-contribution)

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

### Module dependencies

20 algorithms: 8 depend on other modules (in 3 chains) and 12 are standalone. In the diagram, an arrow means "is used by" (it points from a module to the module built on top of it).

```
CHAIN 1: orthogonalisation and least squares

                 ..............................
                 : Gram-Schmidt               :
                 ..............................
                               |
                               v
                 ..............................
                 : QR decomposition           :
                 ..............................
                               |
               +---------------+---------------+
               |                               |
               v                               v
..............................    ..............................
: Least squares              :    : Eigenvalues                :
..............................    ..............................
               |
               v
..............................
: Polynomial fitting         :
..............................


CHAIN 2: LU

                 ..............................
                 : LU decomposition           :
                 ..............................
                               |
               +---------------+---------------+
               |                               |
               v                               v
..............................    ..............................
: Determinant                :    : Matrix inversion           :
..............................    ..............................


CHAIN 3: Gaussian elimination

..............................
: Gaussian elimination       :
..............................
               |
               v
..............................
: Polynomial interpolation   :
..............................


STANDALONE (standard library only, no arrows)

..............................    ..............................
: jacobi                     :    : gauss_seidel               :
..............................    ..............................

..............................    ..............................
: cholesky                   :    : power_iteration            :
..............................    ..............................

..............................    ..............................
: lagrange                   :    : newton_interpolation       :
..............................    ..............................

..............................    ..............................
: romberg_integration        :    : gaussian_quadrature        :
..............................    ..............................

..............................    ..............................
: runge_kutta2               :    : runge_kutta4               :
..............................    ..............................

Legend:  [ box ] = module     --->  = "is used by"
```

| Module | Uses | Used by |
|---|---|---|
| `gram_schmidt` | none | `qr_decomposition` |
| `qr_decomposition` | `gram_schmidt` | `least_squares`, `eigenvalues` |
| `least_squares` | `qr_decomposition` | `polynomial_fitting` |
| `polynomial_fitting` | `least_squares` | none |
| `eigenvalues` | `qr_decomposition` | none |
| `lu` | none | `determinant`, `matrix_inversion` |
| `determinant` | `lu` | none |
| `matrix_inversion` | `lu` | none |
| `gaussian` | none | `polynomial_interpolation` |
| `polynomial_interpolation` | `gaussian` | none |

Standalone modules (use only the C++ standard library, so each is compiled with just its own `.cpp`):

| Area | Standalone modules |
|---|---|
| Linear systems | `jacobi`, `gauss_seidel` |
| Factorisation | `cholesky` |
| Eigenvalues | `power_iteration` |
| Interpolation | `lagrange`, `newton_interpolation` |
| Integration | `romberg_integration`, `gaussian_quadrature` |
| ODE solvers | `runge_kutta2`, `Runge_Kutta` (RK4) |

Notes:
- `gram_schmidt`, `lu` and `gaussian` are roots: other modules build on them, but they build on nothing.
- `jacobi`'s test includes `lu.h` only to cross-check answers, so the library module itself is standalone.
- `gaussian` is standalone as a module, but it is a dependency of `polynomial_interpolation`, so it is not a leaf.
- `polynomial_interpolation` needs `gaussian.cpp` (it builds its Vandermonde system with `numcomp::gaussianElimination`). `CMakeLists.txt` compiles every source into one `numcomp` library, so this is handled automatically when you build with CMake.

## 3. Project structure

```
GROUP-24/
├── include/numerical-computing/   # Public headers (one per algorithm)
├── src/                           # Implementations (.cpp)
├── examples/                      # One demo program per module
├── tests/                         # One test program per module
├── reports/
│   ├── week-01/                   # Gaussian, Gram-Schmidt, LU
│   ├── week-02/                   # Determinant, matrix inversion, polynomial interpolation, QR (report(2).md)
│   ├── week-03/                   # Cholesky, eigenvalues, Jacobi, Lagrange, least squares,
│   │                              # polynomial fitting, power iteration, Romberg, Runge-Kutta 2
│   └── week-04/                   # Gauss-Seidel, Gaussian quadrature, Newton interpolation
├── .github/CODEOWNERS
├── CMakeLists.txt
├── LICENSE
└── README.md
```

## 4. Module status

All 20 modules build with no warnings under `-Wall -Wextra` and all 20 tests pass (`ctest`, 20/20).

| Module | Builds | Tests | Notes |
|---|---|---|---|
| Gaussian elimination | Yes | Pass | |
| LU decomposition | Yes | Pass | |
| Jacobi iteration | Yes | Pass | Test also uses `lu.h` to cross-check answers |
| Gauss-Seidel | Yes | Pass | Header still uses the `GaussSeidel` namespace, not `numcomp` |
| Determinant | Yes | Pass | Needs `src/lu.cpp` |
| Matrix inversion | Yes | Pass | Needs `src/lu.cpp` |
| Cholesky | Yes | Pass | |
| Gram-Schmidt | Yes | Pass | |
| QR decomposition | Yes | Pass | Needs `src/gram_schmidt.cpp` |
| Eigenvalues | Yes | Pass | Needs QR and Gram-Schmidt sources |
| Power iteration | Yes | Pass | Header is not in the `numcomp` namespace |
| Least squares | Yes | Pass | Needs QR and Gram-Schmidt sources |
| Polynomial fitting | Yes | Pass | Needs least squares, QR and Gram-Schmidt sources |
| Lagrange interpolation | Yes | Pass | |
| Newton interpolation | Yes | Pass | |
| Polynomial interpolation | Yes | Pass | Needs `src/gaussian.cpp`; header is not in the `numcomp` namespace |
| Romberg integration | Yes | Pass | Header is not in the `numcomp` namespace |
| Gaussian quadrature | Yes | Pass | Standalone; builds with no warnings under `-Wall -Wextra -Wpedantic` |
| Runge-Kutta 2 | Yes | Pass | Header is not in the `numcomp` namespace |
| Runge-Kutta 4 | Yes | Pass | |

## 5. Build

### With CMake

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

`CMakeLists.txt` builds one `numcomp` library from every `.cpp` file in `src/`. For each module it also builds a test (`test_<name>`, registered with CTest) and an example (`example_<name>`). To add a module, list its `.cpp` file in `add_library(numcomp ...)` and add one `add_module(<name> <test> <example>)` line.

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
| `polynomial_interpolation` | module + `gaussian.cpp` |

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

Each module has a test program in `tests/`. The simplest way to run them all is through CMake:

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

To build a single test by hand, compile it the same way as an example, for instance:

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
| `test_gaussian.cpp` | Gaussian elimination |

## 9. Template candidates

The library headers currently use no C++ templates. The only templates in the repository are the `template <typename F>` check helpers in the test files (`test_eigenvalues.cpp`, `test_least_squares.cpp`, `test_lagrange.cpp`, `newton_interpolation_test.cpp`, `test_matrix_inversion.cpp`, `test_gram_schmidt.cpp`, `test_jacobi.cpp`, `polynomial_fitting_test.cpp`, `test_lu.cpp`, `test_gaussian_quadrature.cpp`), which accept any callable.

Good candidates for templating in the library:

| Where | Current code | Candidate |
|---|---|---|
| Type aliases in `lu.h`, `jacobi.hpp`, `determinant.hpp`, `eigenvalues.hpp`, `gram_schmidt.hpp` | `using Vector = std::vector<double>;` | An alias template `template <typename T> using Vector = std::vector<T>;` so the algorithms work with `float` or `long double` |
| `Matrix` alias in `gaussian.hpp` | `using Matrix = std::vector<std::vector<double>>;` | `template <typename T> using Matrix = std::vector<std::vector<T>>;` |
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
- **Object-oriented features:** the library is mostly free functions in the `numcomp` namespace. The main class is `RK4Solver` (`Runge_Kutta.hpp`), with private data behind public member functions. `SingularMatrixError` (`lu.h`) is the one derived class and inherits from `std::invalid_argument`.
- **Known inconsistencies:** the headers for `gauss_seidel` (namespace `GaussSeidel`), `polynomial_interpolation`, `power_iteration`, `romberg_integration` and `runge_kutta2` are not yet in the `numcomp` namespace. They work as they are; moving them would also mean updating their tests and examples. The CMake include path also still exposes `include/numerical-computing` so that older bare includes keep compiling.

## 11. Weekly reports

Progress reports, one per module, are in `reports/week-01` to `reports/week-04`. They describe how each module was built, the issues faced and the design decisions.

| Week | Reports |
|---|---|
| `week-01` | Gaussian elimination, Gram-Schmidt, LU decomposition |
| `week-02` | Determinant, matrix inversion, polynomial interpolation, QR decomposition (`report(2).md`) |
| `week-03` | Cholesky, eigenvalues, Jacobi, Lagrange, least squares, polynomial fitting, power iteration, Romberg, Runge-Kutta 2 |
| `week-04` | Gauss-Seidel, Newton interpolation, Gaussian quadrature |

## 12. Team workflow

- `main` is protected by a branch ruleset (`protect-main`) that requires a pull request with one approval before merging.
- Work is done on a feature branch per module and merged by pull request.
- `.github/CODEOWNERS` assigns review ownership.
- The project is released under the licence in `LICENSE`.

### 12.1 Review flow

1. A member builds their module on a feature branch (header, source, test, example, weekly report).
2. The member opens a pull request into `main`.
3. A reviewer approves it. The ruleset needs one approval before the **Merge** button works.
4. The PR is merged. The merge commit records who pressed merge, and the PR page on GitHub shows the reviewers and their approvals.

### 12.2 Who did what: per member

| Member | Modules | Commits on `main` | PRs they authored | PRs they merged |
|---|---|---|---|---|
| Ntambi Gabriel | Gram-Schmidt, Lagrange interpolation, Gaussian quadrature | 47 | 23 | 55 |
| Muwonge Emmanuel | Matrix inversion, Polynomial fitting, Newton interpolation | 16 | 3 | 2 |
| Kirya Gabriel | Polynomial interpolation, Romberg integration, Runge-Kutta 2 | 58 | 13 | 3 |
| Musiime Nuwagaba | QR decomposition, Runge-Kutta 4 | 25 | 8 | 0 |
| Cargan Hagan | LU decomposition, Eigenvalues, Jacobi iteration | 18 | 5 | 2 |
| Nabirye Fortunate | Determinant, Cholesky, Least squares | 19 | 3 | 1 |
| Nakalule Leila | Gaussian elimination, Power iteration, Gauss-Seidel | 46 | 10 | 2 |

The team has 65 merged pull requests in total. Commit history per member is in [section 14](#14-commit-history-proof-of-contribution).


## 13. Limitations

Every method has conditions under which it is slow, inaccurate or unsafe to use. This section lists them per algorithm, so that a user can pick the right tool. Most points come from the headers and the weekly reports in `reports/`.

**General limits (apply to the whole library)**

- Everything uses `double`. There is no support for `float`, `long double` or complex numbers (see [section 9](#9-template-candidates)).
- Matrices are dense `std::vector<std::vector<double>>`. There is no sparse storage, so large sparse systems waste memory and time.
- Tolerances are mostly **absolute**, not relative. A matrix scaled by a very small factor (for example `1e-13 * I`) can be reported as singular.
- Some modules take a matrix as **rows** (`luDecompose`, `jacobi`) and others as **columns** (`gramSchmidt`, `qrDecomposition`, `leastSquares`). Check the header before calling.
- Five headers are not yet in the `numcomp` namespace (see [section 10](#10-design-notes)).

### Linear systems

| Algorithm | Limitations |
|---|---|
| Gaussian elimination | Dense only, cost grows as n³. A pivot below `1e-12` is treated as singular and throws `std::runtime_error`, so badly scaled but valid systems may be rejected. Takes `A` and `b` by value, so large systems are copied. |
| LU decomposition | Square matrices only. Pivots below the absolute tolerance (default `1e-12`) throw `SingularMatrixError`. Cannot factor singular or nearly singular matrices. |
| Jacobi iteration | Guaranteed to converge only for strictly diagonally dominant matrices. A zero on the diagonal throws, and the rows must be reordered by the caller. If it does not converge, `converged` is `false` and `x` must not be trusted. Usually needs more rounds than Gauss-Seidel. |
| Gauss-Seidel iteration | Same convergence condition as Jacobi. The defaults are loose (`tolerance = 0.0001`, `maxIterations = 100`). It returns `false` when it fails, and `x` is modified in place. Uses the `GaussSeidel` namespace, not `numcomp`. |

### Matrix operations and factorisation

| Algorithm | Limitations |
|---|---|
| Determinant | Built on LU, so it inherits the absolute pivot tolerance. A matrix LU calls singular returns exactly `0`, even if the true determinant is tiny but non-zero. Large matrices can overflow or underflow, because the result is a product of the pivots. |
| Matrix inversion | Square, non-singular matrices only. Inverting is slower and less accurate than solving `Ax = b` directly with `luSolve`, so avoid it when only a solution is needed. |
| Cholesky | The matrix must be square, symmetric and positive definite. Anything else is rejected. |
| Gram-Schmidt | The modified variant is used, but it still loses orthogonality for nearly dependent vectors. Linearly dependent vectors throw. It needs at least as many dimensions as vectors. |
| QR decomposition | Built on Gram-Schmidt, so it has the same loss of orthogonality on ill-conditioned input. It is not based on Householder reflections. Rank-deficient matrices throw. |

### Eigenvalues

| Algorithm | Limitations |
|---|---|
| Eigenvalues (QR algorithm) | Reliable for symmetric matrices and matrices with all-real eigenvalues. With complex eigenvalues it never settles, returns `converged = false`, and the values must not be trusted. Returns eigenvalues only, not eigenvectors. Uses an unshifted iteration (except for singular matrices), so convergence can be slow when eigenvalues are close in size. |
| Power iteration | Finds only the dominant eigenvalue. The matrix is limited to **20 × 20** (fixed `double A[20][20]` arrays). It stops after 1000 iterations. It can fail when two eigenvalues have the same magnitude, and throws `std::runtime_error` if it does not converge. Not in the `numcomp` namespace. |

### Fitting and interpolation

| Algorithm | Limitations |
|---|---|
| Least squares | Needs at least as many rows as columns (`m >= n`) and full column rank. It solves via QR rather than the normal equations, which is stable, but rank-deficient systems throw. |
| Polynomial fitting | Needs at least `degree + 1` points. A high degree gives an ill-conditioned fit, and the coefficients are in the monomial basis. |
| Lagrange interpolation | Recomputes everything on each call, at O(n²) per point. Evenly spaced points with a high degree cause **Runge's phenomenon** (wild swings between points). Duplicate x values throw. No extrapolation warning is given. |
| Newton interpolation | Same Runge's phenomenon as Lagrange. Duplicate x values throw. Divided differences must be recomputed if a point is added in the middle. |
| Polynomial interpolation | Solves a Vandermonde system, which is badly conditioned as the number of points grows, so use it only for small data sets. The duplicate-x check uses exact equality (`==`), so near-duplicates are not caught. Not in the `numcomp` namespace. |

### Integration

| Algorithm | Limitations |
|---|---|
| Romberg integration | There is **no input validation**: `n <= 0` silently returns `0.0`, and an empty function, `NaN` or infinite limits are not rejected. It uses a fixed number of levels (default 5) and has no tolerance or error estimate. Needs 2^(n-1) + 1 function evaluations, and the shift `1 << (i - 1)` overflows for very large `n`. Rounding error appears at deep levels (`n > 10`). Works best on smooth functions. Not in the `numcomp` namespace. |
| Gaussian quadrature | Works best for smooth functions. A corner or jump inside a panel lowers the accuracy (put a panel edge on the corner). Functions that blow up at an endpoint are not handled. The number of points is capped at 100. `gaussLegendre` rebuilds the rule on every call, whereas `compositeGaussLegendre` builds it once. |

### ODE solvers

| Algorithm | Limitations |
|---|---|
| Runge-Kutta 2 | Fixed step size, with no error control. Only second-order accurate, so it needs small steps for good accuracy. Solves one scalar equation `dy/dt = f(t, y)`, not systems. Not suited to stiff problems. Not in the `numcomp` namespace. |
| Runge-Kutta 4 | Fixed step size, with no error control. Scalar equations only. Stores every step in `t_history` and `y_history`, so very small steps use a lot of memory. Not suited to stiff problems. |

---

## 14. Commit history: proof of contribution

Each member's work can be checked on GitHub through the commit history filtered by their account. The pull request history shows who opened and who approved each module (see [section 12](#12-team-workflow)).

| Member | Reg. No. | Modules | Commit history |
|---|---|---|---|
| Ntambi Gabriel | 25/U/0855 | Gram-Schmidt, Lagrange, Gaussian quadrature | [View commits](https://github.com/ntambigabriel0/GROUP-24/commits/main/?author=ntambigabriel0) |
| Muwonge Emmanuel | 25/U/08755/PS | Matrix inversion, Least-squares polynomial fitting, Newton interpolation | [View commits](https://github.com/ntambigabriel0/GROUP-24/commits/main/?author=Muwonge-Emmanuel) |
| Kirya Gabriel Ssese | 25/U/08729/PS | Polynomial interpolation, Romberg, Runge-Kutta 2 | [View commits](https://github.com/ntambigabriel0/GROUP-24/commits/main/?author=kirya-gabriel) |
| Musiime Nuwagaba | 25/U/08600/PS | QR decomposition, Runge-Kutta 4 | [View commits](https://github.com/ntambigabriel0/GROUP-24/commits/main/?author=Mus1ime) |
| Cargan Hagan | 25/U/08707/PS | LU decomposition, Eigenvalues, Jacobi | [View commits](https://github.com/ntambigabriel0/GROUP-24/commits/main/?author=hcargan) |
| Nabirye Fortunate | 25/U/08759/PS | Determinant, Least squares, Cholesky | [View commits](https://github.com/ntambigabriel0/GROUP-24/commits/main/?author=fortunate-design) |
| Nakalule Leila | 25/U/03499/PSA | Gaussian elimination, Gauss-Seidel, Power iteration | [View commits](https://github.com/ntambigabriel0/GROUP-24/commits/main/?author=leilanakalule-debug) |

The repository-wide history is at `https://github.com/ntambigabriel0/GROUP-24/commits/main`.
