#include <numerical-computing/lu.h>

#include <cmath>      // std::fabs
#include <stdexcept>  // std::invalid_argument
#include <string>     // std::to_string
#include <utility>    // std::swap, std::move

namespace numcomp {

LUResult luDecompose(const Matrix& A, double tol) {
    // Step 1: check the input before doing any maths \
    const std::size_t n = A.size();  // n = number of rows
    if (n == 0) {
        throw std::invalid_argument("luDecompose: matrix must not be empty");
    }
    for (const Vector& row : A) {
        if (row.size() != n) {  // every row must be as long as the matrix is tall
            throw std::invalid_argument("luDecompose: matrix must be square");
        }
    }

    // Step 2: set up the working space 
    Matrix U = A;                  // copy of A: it is reduced into U, and A is untouched
    Matrix L(n, Vector(n, 0.0));   // n x n of zeros: multipliers are stored here
    std::vector<int> perm(n);      // row swap record
    for (std::size_t i = 0; i < n; ++i) {
        perm[i] = static_cast<int>(i);  // starts as [0, 1, 2, ...]
    }
    int swaps = 0;                 // counts the row swaps

    // Step 3: one column at a time ----
    for (std::size_t k = 0; k < n; ++k) {  // k = current column and pivot position

        // Find the row (from k downward) with the biggest value in column k
        std::size_t pivotRow = k;
        for (std::size_t i = k + 1; i < n; ++i) {
            if (std::fabs(U[i][k]) > std::fabs(U[pivotRow][k])) {
                pivotRow = i;
            }
        }

        // If even the best pivot is about zero, the matrix is singular
        if (std::fabs(U[pivotRow][k]) < tol) {
            throw SingularMatrixError("luDecompose: zero pivot in column " +
                                      std::to_string(k) + ", matrix is singular");
        }

        // Swap the best row into place (in U, in perm, and in the multipliers so far)
        if (pivotRow != k) {
            std::swap(U[k], U[pivotRow]);
            std::swap(perm[k], perm[pivotRow]);
            for (std::size_t j = 0; j < k; ++j) {
                std::swap(L[k][j], L[pivotRow][j]);  // keep multipliers with their rows
            }
            ++swaps;
        }

        // Elimination: clear the entries below the pivot and keep the multipliers
        for (std::size_t i = k + 1; i < n; ++i) {
            const double m = U[i][k] / U[k][k];  // the multiplier
            L[i][k] = m;                         // keep it in L
            for (std::size_t j = k; j < n; ++j) {
                U[i][j] -= m * U[k][j];          // subtract m times the pivot row
            }
            U[i][k] = 0.0;                       // force the cleared entry to exactly 0
        }
    }

    // Step 4: L has 1s on its diagonal
    for (std::size_t i = 0; i < n; ++i) {
        L[i][i] = 1.0;
    }

    return LUResult{std::move(L), std::move(U), std::move(perm), swaps};
}

Vector luSolve(const LUResult& f, const Vector& b) {
    const std::size_t n = f.L.size();
    if (b.size() != n) {
        throw std::invalid_argument(
            "luSolve: b must have one entry per row of the matrix");
    }

    // Reorder b to match the row swaps
    Vector pb(n);
    for (std::size_t i = 0; i < n; ++i) {
        pb[i] = b[static_cast<std::size_t>(f.perm[i])];
    }

    // Forward substitution: solve L*y = pb, top to bottom (no division: L's diagonal is 1)
    Vector y(n);
    for (std::size_t i = 0; i < n; ++i) {
        double s = pb[i];
        for (std::size_t j = 0; j < i; ++j) {
            s -= f.L[i][j] * y[j];
        }
        y[i] = s;
    }

    // Backward substitution: solve U*x = y, bottom to top.
    // (i-- > 0 is used because an unsigned size_t can never go below 0)
    Vector x(n);
    for (std::size_t i = n; i-- > 0;) {
        double s = y[i];
        for (std::size_t j = i + 1; j < n; ++j) {
            s -= f.U[i][j] * x[j];
        }
        x[i] = s / f.U[i][i];
    }
    return x;
}

}  // namespace numcomp
