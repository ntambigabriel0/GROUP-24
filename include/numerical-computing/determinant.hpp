#ifndef NUMERICAL_COMPUTING_DETERMINANT_HPP
#define NUMERICAL_COMPUTING_DETERMINANT_HPP

#include <vector>

namespace numcomp {

using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;

// Calculate the determinant of a square matrix.
double determinant(const Matrix& A, double tol = 1e-12);

}  // namespace numcomp

#endif
