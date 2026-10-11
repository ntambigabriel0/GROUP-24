# Project Report: Class Structure of the QR Decomposition Program

## Overview

This report is about QR decomposition, which is a Tier 2 component of the project. It factors a matrix A into Q and R. The matrix is stored as a list of column vectors. All code is in the `numcomp` namespace.

Because this is Tier 2, my work is built on the Tier 1 Gram-Schmidt implementation (`gram_schmidt.hpp`), which was already made by my groupmate Ntambi. I did not write Gram-Schmidt. My `qrDecomposition` function calls Ntambi's `gramSchmidt` to get Q, and uses `dot` to compute R.

The program does not use traditional classes with private data and member functions. It uses one struct (`QRResult`) to hold results and free functions to do the work. This report explains that design and why it was chosen.

## What class did you create?

I created one struct: `QRResult`. In C++, a struct is a class whose members are public by default.

The program also uses `numcomp::Vector`, which comes from Ntambi's Tier 1 file `gram_schmidt.hpp`. It appears to be a vector of doubles. I did not define it.

## Why did you create it?

`qrDecomposition` has to return two results, Q and R. A function can only return one value, so `QRResult` packages both into a single object. The caller reads `result.Q` and `result.R`.

## What data does it hold?

- `Q` (`std::vector<Vector>`): the orthonormal columns of the matrix.
- `R` (`std::vector<Vector>`): the upper-triangular scaling matrix, stored as columns.

## What functions operate on that data?

`QRResult` has no member functions. Its data is filled in by the free function:

`QRResult qrDecomposition(const std::vector<Vector>& A_columns, double tol = 1e-9)`

This function:
1. Throws `std::invalid_argument` if the input is empty.
2. Calls `gramSchmidt` to get the orthonormal columns Q.
3. Creates R filled with zeros.
4. Fills the entries on and above the diagonal of R, using `dot(Q[i], A[j])` for each i <= j.
5. Returns a `QRResult` containing Q and R.

The function uses `gramSchmidt` and `dot` from Tier 1 (Ntambi's work).

## Which members are private and which are public?

All members of `QRResult` are public. There are no private members.

## Why are all members public?

`QRResult` is only a data container. It has no rules to protect: any pair of Q and R matrices is a valid value for it. Public members let later code (the tests, the example program, and Tier 3) read Q and R directly without getter functions.

If the struct had to keep its data consistent (for example, guarantee that R is always upper triangular), private members with controlled access would be the better design.

## How is the class used by the rest of the program?

1. The caller builds the matrix A as a `std::vector<Vector>` of columns.
2. The caller calls `numcomp::qrDecomposition(A)`.
3. The function returns a `QRResult`.
4. The caller reads `result.Q` and `result.R` to print them or pass them to Tier 3.

The example file (3x3 matrix) and the test file (2x2 matrix) both follow this pattern.

## How the code is organized

- `numerical-computing/qr_decomposition.hpp`: declares `QRResult` and `qrDecomposition`.
- `numerical-computing/qr_decomposition.cpp`: implements `qrDecomposition`.
- `qrDecomposition_example.cpp`: example program using a 3x3 matrix.
- `qrDecomposition_test.cpp`: test program using a 2x2 matrix.
