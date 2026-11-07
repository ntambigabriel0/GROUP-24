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

/// All eigenvalues of a square matrix, using the shifted QR algorithm.
///
/// Repeatedly factors a shifted copy of the matrix as Q*R (using qrDecomposition)
/// and multiplies the factors back in reverse order, R*Q, which keeps the
/// eigenvalues the same while pushing the matrix toward triangular form. The
/// shift (a Wilkinson shift, taken from the bottom-right 2x2 block) makes the
/// last eigenvalue appear in a handful of rounds, even when two eigenvalues are
/// nearly equal in size or exactly +x and -x (where the plain unshifted method
/// slows down or never finishes). Each time the last row has become tiny, that
/// eigenvalue is recorded and the matrix shrinks by one. A itself is not changed.
///
/// Reliable for symmetric matrices, or matrices whose eigenvalues are all real.
/// Matrices with complex eigenvalues never settle into triangular form: the
/// result then has converged = false and its values should not be trusted.
/// Very badly scaled matrices (eigenvalues spanning a factor of a million or
/// more) may also run out of rounds; converged = false then says so.
///
/// Singular matrices (an eigenvalue of 0) are handled: the QR step refuses a
/// nearly singular matrix, and the shift is then nudged slightly and tried again.
///
/// Throws std::invalid_argument if:
///   - the matrix is empty,
///   - the matrix is not square,
///   - tol is not positive,
///   - maxIterations is negative.
///
/// `tol` is a RELATIVE tolerance: an entry counts as tiny when it is at most
/// tol times the biggest entry of A, so the matrix may be scaled by any factor.
/// `maxIterations` is the total number of QR rounds allowed across all the
/// eigenvalues; it stops the loop running forever.
EigenResult eigenvalues(const Matrix& A,
                        double tol = 1e-10,
                        int maxIterations = 1000);  // defaults if the user gives none

}  // namespace numcomp

#endif  // NUMERICAL_COMPUTING_EIGENVALUES_HPP
