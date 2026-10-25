#ifndef NUMERICAL_COMPUTING_GAUSS_SEIDEL_HPP
#define NUMERICAL_COMPUTING_GAUSS_SEIDEL_HPP

#include <vector>
#include <stdexcept>

// Solves a system of linear equations using the
// Gauss-Seidel iteration method.
//
// Parameters:
// A              - Coefficient matrix
// b              - Right-hand-side vector
// x              - Initial guess; updated with the solution
// maxIterations  - Maximum number of iterations
// tolerance      - Convergence tolerance
//
// Returns:
// true if the method converges; false otherwise.
//
// Throws:
// std::invalid_argument if the matrix or vectors are invalid.

namespace GaussSeidel {

bool solve(
    const std::vector<std::vector<double>>& A,
    const std::vector<double>& b,
    std::vector<double>& x,
    int maxIterations = 100,
    double tolerance = 0.0001
);

} // namespace GaussSeidel

#endif // NUMERICAL_COMPUTING_GAUSS_SEIDEL_HPP
