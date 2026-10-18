#ifndef NUMERICAL_COMPUTING_CHOLESKY_HPP
#define NUMERICAL_COMPUTING_CHOLESKY_HPP

#include <vector>

namespace numcomp {

using Matrix = std::vector<std::vector<double>>;

// Returns the lower triangular matrix L such that A = L * L^T.
// The input matrix must be square, symmetric, and positive definite.
Matrix choleskyDecomposition(
    const Matrix& A,
    double tol = 1e-12
);

} // namespace numcomp

#endif
