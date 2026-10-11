#ifndef NUMERICAL_COMPUTING_MATRIX_INVERSION_HPP
#define NUMERICAL_COMPUTING_MATRIX_INVERSION_HPP

#include <numerical-computing/lu.h>

namespace numcomp {

// Computes the inverse of a square matrix using LU decomposition.
//
// Throws std::invalid_argument if the matrix is empty, contains empty rows, 
// or is not square.
// Throws SingularMatrixError if the matrix is singular.
//
// Note: Due to the absolute tolerance caveat in the underlying LU decomposition,
// matrices scaled by very small factors (e.g., 1e-13 * I) may be reported as singular.
Matrix invertMatrix(const Matrix& A);

}  // namespace numcomp

#endif  // NUMERICAL_COMPUTING_MATRIX_INVERSION_HPP