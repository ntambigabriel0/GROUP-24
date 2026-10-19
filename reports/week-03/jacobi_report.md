# Week 3 Progress Report: Jacobi Iteration

**Module:** Jacobi Iteration (standalone iterative solver)
**Files:** `include/numerical-computing/jacobi.hpp`, `src/jacobi.cpp`, `tests/test_jacobi.cpp`, `examples/jacobi_example.cpp`

## Completed
- Implemented the Jacobi Iteration module in the `numcomp` namespace: `jacobi`, which solves A·x = b approximately by repeated rounds
- Returns a `JacobiResult` holding the solution, the number of rounds, a `converged` flag, and the size of the last change
- Added an optional starting guess and a round limit, because the method only converges for some matrices
- Added input validation: empty matrix, non-square matrix, wrong-sized `b` or starting guess, a zero on the diagonal, non-positive tolerance, round limit below 1
- Detects values that blow up and stops early instead of looping on garbage
- Wrote 25 tests covering normal use, edge cases, and invalid inputs. All pass with no compiler warnings
- Wrote a usage example that prints the first three rounds
- Documented in the header when the method is guaranteed to work (strictly diagonally dominant matrices)

## How the Jacobi Code Was Built
1. **Defined the problem.** Solve A·x = b without factoring A: guess, plug in, and repeat until the values stop changing.
2. **Designed the interface first.** Put the declarations in the header (`JacobiResult`, `jacobi`) and kept the logic in `src/`. The function takes `A`, `b`, an optional starting guess `x0`, a tolerance `tol`, and a round limit `maxIterations`.
3. **Validated input before any maths.** Reject bad sizes, a zero on the diagonal (the method divides by it), a tolerance that is not positive, and a round limit below 1, using `std::invalid_argument`.
4. **Rearranged each equation.** Each equation is solved for its own unknown: the new value is `(b[i] minus the other terms) / A[i][i]`.
5. **Kept an old and a new copy.** Every new value uses only the previous round's values, never a value from the same round. This is what makes it Jacobi and not Gauss-Seidel.
6. **Added the stopping rule.** After each round, find the biggest change in any value. When it is at most `tol`, stop. The round limit stops matrices that never settle.
7. **Handled divergence.** If any value becomes infinite or not a number, stop and report `converged = false`.
8. **Tested as we went.** Reproduced the plan's rounds by hand (4x + y = 6 and 2x + 5y = 12 from (0,0) give (1.5, 2.4), then (0.9, 1.8), then (1.05, 2.04)), cross-checked a 3×3 answer against the LU solver, then tested the edge cases and the invalid inputs.

## Program Organization (Classes)

### 1. What classes did you create?
None. This module has no classes. It has one result struct, `JacobiResult`, and one free function, `jacobi`, inside the `numcomp` namespace.

### 2. Why did you organize it this way?
Solving a system is a calculation: a matrix and a vector go in and a solution comes out. Nothing has to be remembered between calls, so there is no state for a class to hold. A struct bundles the solution with the information needed to judge it (rounds, whether it converged, and the last change).

### 3. What data does it hold?
| Item | Type | Public / Private | Purpose |
|---|---|---|---|
| `x` | `Vector` | Public (in `JacobiResult`) | The approximate solution |
| `iterations` | `int` | Public (in `JacobiResult`) | How many rounds were done |
| `converged` | `bool` | Public (in `JacobiResult`) | Whether the values stopped changing |
| `change` | `double` | Public (in `JacobiResult`) | Biggest change in any value during the last round |

### 4. What functions operate on that data?
| Function | Public / Private | What it does |
|---|---|---|
| `jacobi` | Public | Checks the input, runs the rounds, and returns the result |

All the work is inside `jacobi`, which keeps this standalone module easy to read.

### 5. Why are some members public and others hidden?
The header declares only `jacobi` and its result, because that is all other code needs. There are no helper functions to hide.

### 6. How does the rest of the program use it?
Other code calls `numcomp::jacobi(A, b)` and reads `.x`. If `.converged` is false, the solution should not be trusted. The matrix is a list of rows, the same as the LU module.

## How This Module Connects to Others
- **Uses:** only the C++ standard library (the tests also use the LU module to cross-check answers)
- **Used by:** nothing; it is a standalone solver
- **Related modules:** Gauss-Seidel is the same idea but uses each new value straight away, so it usually needs fewer rounds. Gaussian elimination and LU are the exact (direct) alternatives

## Issues Faced
- **Old and new copies.** It is easy to update the values in place by mistake, which turns the method into Gauss-Seidel. A test checks that round 2 uses only round 1's values.
- **Zero on the diagonal.** The method divides by the diagonal entry, so a matrix such as `[[0,1],[1,0]]` is rejected even though it has a solution. The error tells the user to reorder the rows.
- **Divergence.** A matrix that is not diagonally dominant can make the values grow without limit, and they can overflow. The code stops and reports `converged = false`.
- **Floating-point comparison.** Results are never exactly the expected numbers, so the tests compare with a small tolerance instead of `==`.
- **Comparing with Gauss-Seidel.** The plan asks to compare convergence speed, but Gauss-Seidel is not yet in the repository, so this has not been done.
- **Git/GitHub workflow.** Getting the files into the shared repository through a pull request.

## In Progress
- Adding `src/jacobi.cpp` to the library list in `CMakeLists.txt`
- Opening the pull request for review

## Challenges / Blockers
- `CMakeLists.txt` currently builds only the Gram-Schmidt module
- The Gaussian elimination header does not compile at the moment, so the answers were cross-checked against LU instead

## Next Week
- Verify a clean build from scratch following the README steps
- Review teammates' code through pull requests
- Compare the number of rounds with Gauss-Seidel on the same system once it is available

## AI Use
- Tool: Claude
- Purpose: Help with writing and testing the Jacobi code, comparing it with the repository, and drafting this report
- Reason: To save time and to check that the code compiles and the tests pass
