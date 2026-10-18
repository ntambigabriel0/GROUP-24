#ifndef NUMERICAL_COMPUTING_EIGENVALUES_HPP  // include guard: stops the file being included twice
#define NUMERICAL_COMPUTING_EIGENVALUES_HPP

#include <vector>

namespace numcomp {  // everything in the library lives in this namespace

using Vector = std::vector<double>;  // nickname: one vector = a list of numbers
using Matrix = std::vector<Vector>;  // nickname: a matrix = a list of vectors

// Rows or columns: a matrix and its transpose have the same eigenvalues, so it
// does not matter whether the inner vectors are rows (like luDecompose) or
// columns (like gramSchmidt and qrDecomposition). Either gives the same answer.

/// Everything the eigenvalue computation returns, bundled together.
struct EigenResult {
    Vector values;   // the eigenvalues, largest first
    int iterations;  // how many QR rounds were done
    bool converged;  // true if the entries below the diagonal became tiny
};

/// All eigenvalues of a square matrix, using the QR algorithm.
///
/// Repeatedly factors A = Q*R (using qrDecomposition) and replaces A with R*Q.
/// The eigenvalues stay the same, but the matrix approaches triangular form,
/// and the diagonal then holds the eigenvalues. A itself is not changed.
///
/// Reliable for symmetric matrices, or matrices whose eigenvalues are all real.
/// Matrices with complex eigenvalues never settle into triangular form: the
/// result then has converged = false and its values should not be trusted.
///
/// If the matrix is singular (it has an eigenvalue of 0), the QR step cannot
/// factor it. The algorithm then adds a shift to the diagonal, runs again, and
/// subtracts the shift from the answers, so singular matrices still work.
///
/// Throws std::invalid_argument if:
///   - the matrix is empty,
///   - the matrix is not square,
///   - tol is not positive,
///   - maxIterations is negative.
///
/// `tol` is an absolute tolerance: the algorithm stops when every entry below
/// the diagonal is at most this. `maxIterations` stops it looping forever.
EigenResult eigenvalues(const Matrix& A,
                        double tol = 1e-10,
                        int maxIterations = 1000);  // defaults if the user gives none

}  // namespace numcomp

#endif  // NUMERICAL_COMPUTING_EIGENVALUES_HPP
