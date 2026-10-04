#ifndef NUMERICAL_COMPUTING_LU_H  // include guard: stops the file being included twice
#define NUMERICAL_COMPUTING_LU_H

#include <stdexcept>
#include <string>
#include <vector>

namespace numcomp {  // everything in the library lives in this namespace

using Vector = std::vector<double>;  // nickname: one vector = a list of numbers
using Matrix = std::vector<Vector>;  // nickname: a matrix = a list of rows

// Error thrown when a matrix is singular (it has no usable pivot, so it
// cannot be factored). It inherits from std::invalid_argument, so code that
// catches std::invalid_argument will catch this one too.
class SingularMatrixError : public std::invalid_argument {
public:
    explicit SingularMatrixError(const std::string& msg)
        : std::invalid_argument(msg) {}
};

// Everything LU decomposition returns, bundled together.
// The permutation and swap count are included because Determinant and
// Inversion both need them.
struct LUResult {
    Matrix L;               // lower triangle (1s on the diagonal, multipliers below)
    Matrix U;               // upper triangle
    std::vector<int> perm;  // row swap record: row i of L*U is row perm[i] of A
    int swaps;              // how many row swaps happened (sign of the determinant)
};

// LU decomposition with partial pivoting: P*A = L*U.
//
// Takes a square matrix A and splits it into a lower triangular matrix L
// and an upper triangular matrix U. The largest available pivot is swapped
// into place in each column for numerical stability. A itself is not changed.
//
// Throws std::invalid_argument if:
//   - the matrix is empty,
//   - the matrix is not square.
// Throws SingularMatrixError if the matrix is singular (within the tolerance).
//
// `tol` is an absolute tolerance: a pivot smaller than this is treated as zero.
LUResult luDecompose(const Matrix& A,
                     double tol = 1e-12);  // default tolerance if the user gives none

// Solve A*x = b using a factorization already computed by luDecompose.
// Reuse the same factorization for many different b: forward substitution
// through L, then backward substitution through U.
//
// Throws std::invalid_argument if b does not have one entry per row of A.
Vector luSolve(const LUResult& f, const Vector& b);

}  // namespace numcomp

#endif  // NUMERICAL_COMPUTING_LU_H
