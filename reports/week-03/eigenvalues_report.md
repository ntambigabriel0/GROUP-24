# Week 3 Progress Report: Eigenvalue Computation

**Module:** Eigenvalue Computation (Tier 3), the QR algorithm
**Files:** `include/numerical-computing/eigenvalues.hpp`, `src/eigenvalues.cpp`, `tests/test_eigenvalues.cpp`, `examples/eigenvalues_example.cpp`

## Completed
- Implemented the Eigenvalue Computation module in the `numcomp` namespace: `eigenvalues`, which finds **all** the eigenvalues of a square matrix
- It repeatedly uses `qrDecomposition` from the QR module: factor a shifted copy of the matrix as Q·R, then multiply back in reverse order (R·Q) until the last row is tiny, record that eigenvalue, and repeat on the smaller matrix
- Returns an `EigenResult` holding the eigenvalues (largest first), the number of rounds, and a `converged` flag
- Uses a Wilkinson shift with deflation, so it also converges when two eigenvalues are nearly equal in size or exactly +x and -x (the plain unshifted method slows down or never finishes there)
- Handles singular matrices (an eigenvalue of 0): QR refuses a nearly singular matrix, so the shift is nudged slightly and the round is retried
- The tolerance is relative to the size of the matrix, so it works at any scale
- Added input validation: empty matrix, non-square matrix, non-positive tolerance, negative round limit
- Wrote 26 tests covering normal use, edge cases, and invalid inputs. All pass with no compiler warnings
- Verified against numpy on 1,140 random matrices (symmetric, singular and non-symmetric with real eigenvalues): every one converged and matched to about 1e-10
- Wrote a usage example
- Documented in the header that it is reliable for symmetric matrices or matrices whose eigenvalues are all real

## How the Eigenvalue Code Was Built
1. **Defined the problem.** An eigenvalue tells how much a matrix stretches along a special direction. The QR algorithm finds all of them at once, not only the largest one.
2. **Designed the interface first.** Put the declarations in the header (`EigenResult`, `eigenvalues`) and kept the logic in `src/`. The function takes the matrix, a stopping tolerance `tol`, and a round limit `maxIterations`.
3. **Validated input before any maths.** Reject an empty or non-square matrix, a tolerance that is not positive, and a negative round limit, using `std::invalid_argument`.
4. **Built one round.** Subtract a shift from the diagonal, call `qrDecomposition` to get Q and R, then multiply them in reverse order to get R·Q and add the shift back. This is a similarity transformation, so the eigenvalues stay the same. Both Q and R are stored as lists of columns, so entry (row i, column j) of R·Q is the sum over k of `R[k][i] * Q[j][k]`.
5. **Chose the shift.** The first version had no shift. Testing against numpy showed that it did not converge for 8% of random symmetric matrices (eigenvalues nearly equal in size but opposite in sign) and never for exact ties such as `[[0,1],[1,0]]`. The Wilkinson shift, taken from the bottom-right 2x2 block, fixes this and makes each eigenvalue appear in a few rounds.
6. **Added deflation and the stopping rule.** When every entry left of the diagonal in the last row is at most `tol` times the biggest entry of A, the last diagonal entry is an eigenvalue. It is recorded and the matrix shrinks by one row and column. A round limit stops matrices that never settle.
7. **Handled nearly singular shifts.** A good shift is close to an eigenvalue, which makes the shifted matrix nearly singular, and the QR module (through Gram-Schmidt) refuses those. When that happens the shift is nudged slightly and the round is retried. QR is asked to refuse at 1e-6, because Gram-Schmidt loses accuracy on nearly singular matrices.
8. **Tested as we went.** Used matrices with known eigenvalues (the plan example gives 3 and 1), checked that the eigenvalues add up to the trace, compared the largest value with a power iteration, then tested the edge cases and the invalid inputs.

## Program Organization (Classes)

### 1. What classes did you create?
None. This module has no classes. It has one result struct, `EigenResult`, and one free function, `eigenvalues`, inside the `numcomp` namespace.

### 2. Why did you organize it this way?
Finding eigenvalues is a calculation: a matrix goes in and numbers come out. Nothing has to be remembered between calls, so there is no state for a class to hold. A struct bundles the three things the caller needs (the values, the round count, and whether it converged), because a C++ function can return only one value.

### 3. What data does it hold?
| Item | Type | Public / Private | Purpose |
|---|---|---|---|
| `values` | `Vector` | Public (in `EigenResult`) | The eigenvalues, largest first |
| `iterations` | `int` | Public (in `EigenResult`) | How many QR rounds were done |
| `converged` | `bool` | Public (in `EigenResult`) | Whether the entries below the diagonal became tiny |

### 4. What functions operate on that data?
| Function | Public / Private | What it does |
|---|---|---|
| `eigenvalues` | Public | Runs the QR algorithm and returns the eigenvalues |
| `chooseShift` | Private (in `src/eigenvalues.cpp`) | Picks the Wilkinson shift from the bottom-right 2x2 block |
| `shiftedRound` | Private | One round: factor the shifted matrix as Q·R, return R·Q plus the shift |
| `roundWithNudging` | Private | Retries a round with a slightly changed shift if QR refuses the matrix |
| `lastRowLeftOfDiagonal` | Private | Checks whether the last eigenvalue has appeared |
| `maxAbsEntry` | Private | Sizes the relative tolerance and the nudge |

### 5. Why are some members public and others hidden?
Only `eigenvalues` and its result are in the header, because they are all other modules need. The helpers are in an unnamed namespace inside the `.cpp` file, so they can change without affecting anyone else's code.

### 6. How does the rest of the program use it?
Other code calls `numcomp::eigenvalues(A)` and reads `.values`. If `.converged` is false, the values should not be trusted.

## How This Module Connects to Others
- **Uses:** `qrDecomposition` (QR module), which uses `gramSchmidt` and `dot` (Gram-Schmidt module)
- **Used by:** nothing yet; it is the "all eigenvalues" tool
- **Related module:** Power Iteration finds only the largest eigenvalue, so the largest value from this module should agree with it
- **Rows or columns:** a matrix and its transpose have the same eigenvalues, so it does not matter whether the inner vectors are rows (like LU) or columns (like Gram-Schmidt and QR). A test confirms both give the same answer

## Issues Faced
- **Depending on another module.** The QR module had to exist first. Its function takes and returns lists of columns, so the R·Q product had to be written for column storage.
- **Two different tolerances.** QR has its own tolerance for deciding that columns are dependent, and the eigenvalue loop has its own for deciding when to stop. They mean different things, so the loop's `tol` is not passed to QR.
- **Singular matrices.** A matrix with an eigenvalue of 0 made the QR step throw. The first fix added a large shift to the diagonal, but that squeezed the eigenvalues together and a third of random singular matrices then ran out of rounds. Nudging the Wilkinson shift instead fixed this.
- **Slow or no convergence.** Checking the first version against numpy on random matrices showed that nearly-equal-size eigenvalues with opposite signs converge very slowly, and exact ties never do. The shifted version converged on all 1,140 random matrices, in at most 24 rounds.
- **Accuracy and nearly singular matrices.** Gram-Schmidt loses accuracy on nearly singular matrices, and a good shift makes the matrix nearly singular. Asking QR to refuse at 1e-6 and nudging the shift by 1e-4 of the matrix size gave about 1e-10 accuracy.
- **Scale.** An absolute tolerance treated a matrix with entries around 1e-11 as already finished, so the tolerance is now relative to the biggest entry.
- **Known limit.** Very badly scaled matrices (eigenvalues spanning a factor of a million or more) can still run out of rounds. The function then reports `converged = false`.
- **Complex eigenvalues.** A matrix such as a 90-degree rotation never settles into triangular form. The function reports `converged = false` instead of returning wrong numbers silently, and the header says so.
- **Checking against Power Iteration.** In the repository the Power Iteration source contains only a `main` program and no callable `powerIteration` function, so the tests use a small power iteration of their own. This should be swapped for the real one once it is available.
- **Git/GitHub workflow.** Getting the files into the shared repository through a pull request.

## In Progress
- Adding `src/eigenvalues.cpp` (and `src/qr_decomposition.cpp`) to the library list in `CMakeLists.txt`
- Opening the pull request for review
- Replacing the test's own power iteration with the repository's Power Iteration

## Challenges / Blockers
- `CMakeLists.txt` currently builds only the Gram-Schmidt module
- Power Iteration is not yet available as a function that other code can call

## Next Week
- Verify a clean build from scratch following the README steps
- Review teammates' code through pull requests
- Check the largest eigenvalue against the real Power Iteration once it can be called

## AI Use
- Tool: Claude
- Purpose: Help with writing and testing the eigenvalue code (including the shifted QR method with deflation), verifying it against numpy, comparing it with the repository, and drafting this report
- Reason: To save time and to check that the code compiles and the tests pass
