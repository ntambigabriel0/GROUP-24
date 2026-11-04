#ifndef NUMERICAL_COMPUTING_GAUSSIAN_HPP
#define NUMERICAL_COMPUTING_GAUSSIAN_HPP

#include <vector>

namespace numcomp {

using Vector = std::vector<double>;  // nickname: one vector = a list of numbers
using Matrix = std::vector<Vector>;  // nickname: a matrix = a list of rows

// Solves Ax = b using Gaussian elimination with partial pivoting.
// Throws std::invalid_argument if A is empty, not square, or b does not match A.
// Throws std::runtime_error if A is singular (a pivot is below 1e-12).
Vector gaussianElimination(Matrix A, Vector b);

}  // namespace numcomp

#endif // NUMERICAL_COMPUTING_GAUSSIAN_HPP
