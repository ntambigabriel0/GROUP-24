# Progress Report: LU Decomposition (Tier 1)

**Author:** Cargan

## What LU decomposition is
LU decomposition splits a square matrix **A** into two triangular matrices so that **A = L·U** (with row swaps recorded when needed, so strictly P·A = L·U).

- **L** is lower triangular with 1s on its diagonal. Below the diagonal it stores the multipliers used during elimination.
- **U** is upper triangular. It is the result of Gaussian elimination.

It is the reusable factorization of the library. Once A is split, solving for any new right-hand side `b` is cheap: solve forward through L, then backward through U. The Determinant and Inversion modules both build on it, so the function also returns the row permutation and the swap count.

**Worked example.** For `[[2,1],[1,3]]` the only multiplier is `1/2 = 0.5`. L holds it, U is `[[2,1],[0,2.5]]`, and multiplying L·U returns the original matrix.

## Completed
- Implemented LU decomposition with partial pivoting in the `numcomp` namespace: `luDecompose` and `luSolve`
- `luDecompose` returns an `LUResult` holding **L, U, the permutation and the swap count**
- `luSolve` reuses a stored factorization to solve A·x = b for any `b`
- Added input validation: empty matrix, non-square matrix, singular matrix (`SingularMatrixError`), and a right-hand side of the wrong length
- Wrote tests (`tests/test_lu.cpp`) covering normal use, edge cases and invalid inputs (24 checks, all passing)
- Wrote a usage example (`examples/lu_example.cpp`)
- Matched the project format used for Gram-Schmidt: header in `include/numerical-computing/`, implementation in `src/`, `numcomp` namespace, `std::invalid_argument` family of errors

## How the LU code was built
1. **Defined the problem.** Take a square matrix and produce L and U with A = L·U, keeping the multipliers so that the work of elimination is not thrown away.
2. **Designed the interface first.** The header declares `Matrix`, `LUResult`, `SingularMatrixError`, `luDecompose` and `luSolve`. The logic lives in `src/lu.cpp`.
3. **Validated input before any maths.** Reject an empty or non-square matrix with `std::invalid_argument`.
4. **Processed one column at a time.** For each column `k`:
   - find the row at or below `k` with the largest absolute value in that column (partial pivoting) and swap it into place;
   - if even the best pivot is smaller than `tol`, the matrix is singular, so throw `SingularMatrixError`;
   - for each row below the pivot, compute the multiplier `m = U[i][k] / U[k][k]`, store it in `L[i][k]`, and subtract `m` times the pivot row from that row.
5. **Recorded the swaps.** Each swap is applied to U, to the permutation list, and to the multipliers already stored in L, and the swap counter goes up by one.
6. **Finished L.** Put 1s on L's diagonal and return `L`, `U`, `perm` and `swaps` together in one `LUResult`.
7. **Wrote the solver.** `luSolve` reorders `b` using `perm`, then does forward substitution through L (no division needed, the diagonal is 1) and backward substitution through U.
8. **Tested as we went.** Checked that P·A equals L·U, that L and U have the right shape, that solutions have a near-zero residual, and then tested the invalid inputs.

## Issues faced
- **Missing closing braces.** A missing `}` after the row-swap block put the elimination loop inside the `if`, so elimination only ran when a swap happened and the code did not compile. A missing `};` after the error class caused a similar failure. Recompiling after every edit and auto-formatting the file made these easy to spot.
- **Swapping the multipliers too.** When two rows are swapped, the multipliers already saved in L for those rows must be swapped as well, otherwise they stay attached to the wrong rows.
- **Floating-point comparison.** Results are never exactly 0 or 1, so a pivot is compared against a tolerance (`tol`) rather than `== 0`, the cleared entry below each pivot is set to exactly 0, and the tests compare with a small tolerance instead of `==`.
- **Counting down with an unsigned type.** The backward substitution loop must count down to 0. With `std::size_t`, `i >= 0` is always true, so the loop is written as `for (i = n; i-- > 0;)`.
- **Matching the team format.** The first version used plain functions with no namespace. It was reworked to match the Gram-Schmidt module: `numcomp` namespace, camelCase names, `Matrix` type, doc comments on every function, and errors that derive from `std::invalid_argument`.

## Design decisions
- **Partial pivoting.** Choosing the largest pivot avoids dividing by zero or by tiny numbers. A test with a pivot of `1e-20` confirms the answer stays accurate.
- **`SingularMatrixError` inherits from `std::invalid_argument`.** Singular matrices can be caught specifically, and code that catches `std::invalid_argument` still catches them.
- **The permutation and swap count are returned.** The determinant needs the swap count to set its sign, and inversion needs the permutation to reorder each right-hand side.

## In Progress
- Agreeing the output format (`LUResult`: `L`, `U`, `perm`, `swaps`) with the owners of Inversion and Determinant
- Confirming where the shared `Matrix` type should be defined, so modules do not define it differently
- Opening a pull request for review

## Challenges/Blockers
- Making sure `Matrix` is defined once, identically, across modules
- Learning the pull request workflow under the branch protection rules

## Next Week
- Review teammates' code through pull requests
- Check that Inversion and Determinant work correctly with the `LUResult` returned here
- Add `LU` to the build (`CMakeLists.txt`) so its tests run with `ctest`

## AI Use
- Tool: Claude
- Purpose: Help with writing and debugging the LU code (including finding the missing braces), restructuring it to match the team's format, and drafting the tests, the example and this report
- Reason: To save time and to check that the code compiles and the tests pass
