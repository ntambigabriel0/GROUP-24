#include "numerical-computing/cholesky.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace numcomp {

Matrix choleskyDecomposition(const Matrix& A, double tol) {
    if (A.empty()) {
        throw std::invalid_argument(
            "Cholesky: matrix must not be empty."
        );
    }

    if (!std::isfinite(tol) || tol <= 0.0) {
        throw std::invalid_argument(
            "Cholesky: tolerance must be positive and finite."
        );
    }

    const std::size_t n = A.size();

    // Check that the matrix is square and contains finite values.
    double scale = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        if (A[i].size() != n) {
            throw std::invalid_argument(
                "Cholesky: matrix must be square."
            );
        }

        for (std::size_t j = 0; j < n; ++j) {
            if (!std::isfinite(A[i][j])) {
                throw std::invalid_argument(
                    "Cholesky: matrix entries must be finite."
                );
            }

            scale = std::max(scale, std::abs(A[i][j]));
        }
    }

    if (scale == 0.0) {
        throw std::invalid_argument(
            "Cholesky: zero matrix is not positive definite."
        );
    }

    // Check symmetry.
    const double symmetryTolerance = tol * scale;

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            if (std::abs(A[i][j] - A[j][i])
                > symmetryTolerance) {
                throw std::invalid_argument(
                    "Cholesky: matrix must be symmetric."
                );
            }
        }
    }

    // Initialize L as a zero matrix.
    Matrix L(n, std::vector<double>(n, 0.0));

    const double pivotTolerance = tol * scale;

    // Compute the lower triangular matrix L.
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            double sum = 0.0;

            for (std::size_t k = 0; k < j; ++k) {
                sum += L[i][k] * L[j][k];
            }

            if (i == j) {
                const double diagonal = A[i][i] - sum;

                if (!std::isfinite(diagonal) ||
                    diagonal <= pivotTolerance) {
                    throw std::invalid_argument(
                        "Cholesky: matrix is not sufficiently "
                        "positive definite for the chosen tolerance."
                    );
                }

                L[i][j] = std::sqrt(diagonal);
            } else {
                L[i][j] = (A[i][j] - sum) / L[j][j];

                if (!std::isfinite(L[i][j])) {
                    throw std::runtime_error(
                        "Cholesky: numerical calculation failed."
                    );
                }
            }
        }
    }

    return L;
}

} // namespace numcomp
