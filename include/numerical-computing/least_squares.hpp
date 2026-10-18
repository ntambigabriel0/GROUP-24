#ifndef NUMCOMP_LEAST_SQUARES_HPP
#define NUMCOMP_LEAST_SQUARES_HPP

#include <numerical-computing/qr_decomposition.hpp>
#include <vector>

namespace numcomp {

// A_columns contains the columns of A.
// b contains the right-hand-side vector.
// Returns x that minimizes ||Ax - b|| squared.
Vector leastSquares(
    const std::vector<Vector>& A_columns,
    const Vector& b,
    double tol = 1e-9
);

} // namespace numcomp

#endif
