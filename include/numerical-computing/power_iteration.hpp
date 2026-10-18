#ifndef NUMERICAL_COMPUTING_POWER_ITERATION_HPP
#define NUMERICAL_COMPUTING_POWER_ITERATION_HPP

#include <stdexcept>

// Finds the dominant eigenvalue of a square matrix
// using the Power Iteration method.
//
// Parameters:
// n       - Order of the square matrix (maximum 20)
// A       - Input matrix
// x       - Initial vector; updated to the eigenvector approximation
// error   - Convergence tolerance
//
// Returns:
// Approximation of the dominant eigenvalue.
//
// Throws:
// std::invalid_argument if the matrix size or initial vector is invalid.
// std::runtime_error if the iteration fails or does not converge.

double powerIteration(
    int n,
    double A[20][20],
    double x[20],
    double error = 0.0001
);

#endif // NUMERICAL_COMPUTING_POWER_ITERATION_HPP
