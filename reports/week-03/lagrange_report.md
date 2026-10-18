# Week 2 Progress Report: Lagrange Interpolation

**Module:** Lagrange Interpolation
**Files:** `include/numerical-computing/lagrange.hpp`, `src/lagrange.cpp`, `tests/test_lagrange.cpp`, `examples/lagrange_example.cpp`

## Completed
- Implemented the Lagrange Interpolation module in the `numcomp` namespace: `lagrangeInterpolate` (one x value) and an overload that takes many x values at once
- Added input validation: empty data, `xs` and `ys` of different sizes, duplicate x values, and NaN or infinite values
- Wrote 17 tests covering normal use, edge cases (single point, unsorted points, empty query list), and invalid inputs. All pass with no compiler warnings
- Wrote a usage example
- Updated `CMakeLists.txt` so the library, test, and example build with Lagrange alongside the existing modules; tests run with `ctest`

## How the Lagrange Code Was Built
1. **Defined the problem.** Given n known points (x, y) with different x values, find the one polynomial of degree n-1 that passes through all of them, and return its value at a chosen x.
2. **Designed the interface first.** Put the declarations in the header and kept the implementation in `src/`, so the public API is separate from the logic. The function takes `xs`, `ys`, and the x to evaluate at.
3. **Validated input before any maths.** Reject empty data, mismatched sizes, NaN or infinite values, and duplicate x values, using `std::invalid_argument`. Duplicates must be rejected because the formula divides by `xs[i] - xs[j]`, which would be zero.
4. **Built each basis piece.** For each point i, the piece `L_i(x)` is the product over every other point j of `(x - xs[j]) / (xs[i] - xs[j])`. It equals 1 at `xs[i]` and 0 at every other given x.
5. **Added the pieces together.** The answer is the sum of `ys[i] * L_i(x)`. This is why the polynomial passes exactly through every given point.
6. **Added a many-values version.** It checks the data once and then evaluates each query, so the checks are not repeated for every x.
7. **Tested as we went.** Used answers worked out by hand (points on `y = x² + x + 1`), confirmed the result returns the original y at every given x, confirmed a cubic is reproduced exactly, compared against `sin(x)`, then tested the invalid inputs.

## Program Organization (Classes)

### 1. What classes did you create?
None. This module has no classes. It is organized as free functions inside the `numcomp` namespace.

### 2. Why did you organize it this way?
Lagrange interpolation is a calculation: points go in and a number comes out. Nothing has to be remembered between calls, so there is no state for a class to hold. Plain functions keep the code simple, and the namespace stops our names from clashing with other code.

### 3. What data does it hold?
Nothing. There are no data members and no new types. The data points are passed in by `const` reference as `std::vector<double>`, and the working values (the basis piece and the running total) live only inside the function while it runs.

### 4. What functions operate on that data?
| Function | Public / Private | What it does |
|----------|------------------|--------------|
| `lagrangeInterpolate(xs, ys, x)` | Public | Value of the interpolating polynomial at one x |
| `lagrangeInterpolate(xs, ys, queries)` | Public | Same, for a list of x values; returns a list of results |
| `validate(xs, ys)` | Private (`static` in `lagrange.cpp`) | Checks the data and throws `std::invalid_argument` if it is bad |
| `evaluate(xs, ys, x)` | Private (`static` in `lagrange.cpp`) | Does the actual sum-of-products calculation |

### 5. Why are some members public and others hidden?
The two `lagrangeInterpolate` functions are declared in the header, so other modules can call them. `validate` and `evaluate` are marked `static` in `lagrange.cpp`, so they are visible only inside that file. Users only see the function signatures and their documented behavior, so the maths can change later without breaking their code.

```

## How This Module Connects to Others
- **Uses:** nothing from the library; only the C++ standard library.
- **Used by:** nothing yet. Newton Interpolation (Muwonge) and Polynomial Interpolation (Kirya) fit the same polynomial through the same points, so their results can be compared against this module as a cross-check.

## Issues Faced
- **Duplicate x values.** The formula divides by `xs[i] - xs[j]`, so two equal x values would divide by zero. They are rejected up front with a clear error message.
- **Floating-point comparison.** Decimal results are rarely exactly equal to the expected value, so the tests compare with a small tolerance instead of `==`.
- **NaN and infinite input.** These would silently spread through the sum and give a meaningless answer, so they are rejected at the start.
- **File names with spaces.** The Gram-Schmidt files were first saved with spaces instead of underscores, which stopped other modules (QR) from including them. The Lagrange files were named with underscores from the start.
- **Known limitation.** With many points, especially evenly spaced ones, a single high-degree polynomial can swing wildly between the points (Runge's phenomenon). The module is accurate for a small number of points; this is a property of the method, not a bug.

## In Progress
- Adding the Lagrange files to the repository and opening a pull request
- Reviewing each other's code through pull requests

## Challenges / Blockers
- Making sure every member has visible commits
- Keeping `CMakeLists.txt` in step as each member adds their own files

## Next Week
- Implement Gaussian Quadrature (the next module assigned to me) with matching tests and an example
- Compare Lagrange results with Newton Interpolation (Muwonge) on the same points
- Verify a clean build from scratch following the README steps
- Merge pull requests so every member has commits in the history


