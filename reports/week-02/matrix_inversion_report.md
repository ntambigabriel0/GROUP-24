# Progress Report: Matrix Inversion (Tier 2)

**Author:** Muwonge

## What Matrix Inversion is
The inverse of a square matrix **A** is a matrix **A⁻¹** such that multiplying them together produces the Identity matrix (**I**). This module computes the inverse by generating the columns of the identity matrix and solving the system $A \cdot x = e_i$ for each column. By utilizing LU decomposition, the matrix is factored exactly once, allowing for fast, repeated solving.

## Completed
- Implemented the matrix inversion logic in the `numcomp` namespace: `invertMatrix`.
- Integrated directly with the Tier 1 `luDecompose` and `luSolve` functions via `<numerical-computing/lu.h>`.
- Validated inputs to reject empty or non-square matrices and catch singular matrices.
- Documented the absolute tolerance caveat where tiny matrices (e.g., `1e-13 * I`) are evaluated as singular.
- Wrote tests verifying $A \cdot A^{-1} \approx I$ for 2x2, 3x3, and pivoting cases.

## How the Matrix Inversion code was built
1. **Designed the interface:** Set up `matrix_inversion.hpp` to rely exclusively on the shared `Matrix` and `Vector` aliases defined in `lu.h`.
2. **Validated input:** Rejected empty matrices and ragged dimensions using `std::invalid_argument`.
3. **Factored once:** Called `numcomp::luDecompose(A)` at the start of the function to obtain the **L**, **U**, and **P** factors.
4. **Solved for columns:** Created a loop that builds identity columns and passes them to `numcomp::luSolve(f, e)`.
5. **Tested the logic:** Built a multiplication helper in the test suite to mathematically prove that multiplying the original matrix by the computed inverse yields the Identity matrix.

## Issues faced
- **Header conflicts:** Early versions duplicated the `Matrix` and `Vector` definitions. This was resolved by including `lu.h` directly in the header to ensure type consistency across the library.
- **Redundant checks:** Removed an impossible size check that tested for a dimension of zero after already verifying the matrix was not empty.

## Design decisions
- **Free function over class:** Elected not to create a custom class for this module to preserve compatibility with the shared `numcomp::Matrix` alias and the existing direct solvers.