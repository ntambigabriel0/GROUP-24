#ifndef NUMCOMP_QR_DECOMPOSITION_HPP
#define NUMCOMP_QR_DECOMPOSITION_HPP

#include <vector>

// Provides numcomp::Vector, numcomp::gramSchmidt and numcomp::dot
#include <numerical-computing/gram_schmidt.hpp>

namespace numcomp {

// Packages the matrices so callers can read .Q and .R
struct QRResult {
    std::vector<Vector> Q; // Orthonormal columns matrix
    std::vector<Vector> R; // Upper-triangular scaling matrix
};

// QR decomposition of a matrix given as a list of column vectors.
// Throws std::invalid_argument if A_columns is empty.
QRResult qrDecomposition(const std::vector<Vector>& A_columns, double tol = 1e-9);

} // namespace numcomp

#endif // NUMCOMP_QR_DECOMPOSITION_HPP
