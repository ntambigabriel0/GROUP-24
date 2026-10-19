#ifndef NUMERICAL_COMPUTING_JACOBI_HPP  // include guard: stops the file being included twice
#define NUMERICAL_COMPUTING_JACOBI_HPP

#include <vector>

namespace numcomp {  // everything in the library lives in this namespace

using Vector = std::vector<double>;  // nickname: one vector = a list of numbers
using Matrix = std::vector<Vector>;  // nickname: a matrix = a list of rows

// Matrix convention for this function: a Matrix is a list of ROWS, so A[i][j] is
// row i, column j. That matches how a system of equations is written, and it is
// the same as luDecompose.

/// Everything the Jacobi iteration returns, bundled together.
struct JacobiResult {
    Vector x;         // the approximate solution
    int iterations;   // how many rounds were done
    bool converged;   // true if the last round changed every value by at most tol
    double change;    // the biggest change in any value during the last round
};

/// Solves A*x = b with the Jacobi iteration.
///
/// Each equation is rearranged to solve for its own unknown. Starting from a
/// guess, every new value is worked out using only the PREVIOUS round's values,
/// so an old copy and a new copy are kept. Rounds repeat until the values stop
/// changing. This is the simpler sibling of Gauss-Seidel, which uses each new
/// value straight away and usually needs fewer rounds.
///
/// Example: 4x + y = 6 and 2x + 5y = 12 from (0,0) give (1.5, 2.4), then
/// (0.9, 1.8), then (1.05, 2.04), closing in on (1, 2).
///
/// It is only guaranteed to work when A is strictly diagonally dominant (on every
/// row the diagonal entry is larger than the sum of the other entries in absolute
/// value). For other matrices it may converge slowly or not at all: it then
/// returns converged = false, and x should not be trusted. A and b are not changed.
///
/// Throws std::invalid_argument if:
///   - A is empty,
///   - A is not square,
///   - b does not have one entry per row of A,
///   - x0 is given but does not have one entry per row of A,
///   - a diagonal entry of A is zero (reorder the rows so it is not),
///   - tol is not positive,
///   - maxIterations is less than 1.
///
/// `x0` is the starting guess. If it is empty, every value starts at 0.
/// `tol` is an absolute tolerance: it stops when no value changes by more than this.
/// `maxIterations` stops it looping forever on matrices that do not converge.
JacobiResult jacobi(const Matrix& A,
                    const Vector& b,
                    const Vector& x0 = Vector(),
                    double tol = 1e-10,
                    int maxIterations = 1000);  // defaults if the user gives none

}  // namespace numcomp

#endif  // NUMERICAL_COMPUTING_JACOBI_HPP
