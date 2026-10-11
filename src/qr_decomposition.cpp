#include <numerical-computing/qr_decomposition.hpp>

#include <cstddef>
#include <stdexcept>
#include <utility>

namespace numcomp {

QRResult qrDecomposition(const std::vector<Vector>& A_columns, double tol) {
    if (A_columns.empty()) {
        throw std::invalid_argument("qrDecomposition: Input matrix A cannot be empty");
    }

    const std::size_t num_cols = A_columns.size();

    // Extract the orthonormal vector basis
    std::vector<Vector> Q_columns = numcomp::gramSchmidt(A_columns, tol);

    // Pre-allocate R with zeros
    std::vector<Vector> R_columns(num_cols, Vector(num_cols, 0.0));

    // Compute the upper-triangular entries
    for (std::size_t j = 0; j < num_cols; ++j) {
        for (std::size_t i = 0; i <= j; ++i) { // Entries below the diagonal remain 0.0
            R_columns[j][i] = numcomp::dot(Q_columns[i], A_columns[j]);
        }
    }

    return QRResult{std::move(Q_columns), std::move(R_columns)};
}

} // namespace numcomp
