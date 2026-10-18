#include <numerical-computing/least_squares.hpp>
#include <numerical-computing/qr_decomposition.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace numcomp {

Vector leastSquares(
    const std::vector<Vector>& A_columns,
    const Vector& b,
    double tol
) {
    // Check the tolerance.
    if (!std::isfinite(tol) || tol < 0.0) {
        throw std::invalid_argument(
            "leastSquares: tolerance must be finite and non-negative"
        );
    }

    // Check that the matrix and right-hand side are not empty.
    if (A_columns.empty() || b.empty()) {
        throw std::invalid_argument(
            "leastSquares: matrix and vector must not be empty"
        );
    }

    const std::size_t m = b.size();
    const std::size_t n = A_columns.size();

    // Least squares requires at least as many equations as unknowns.
    if (m < n) {
        throw std::invalid_argument(
            "leastSquares: fewer equations than unknowns"
        );
    }

    // Check that every column has the correct number of entries.
    for (const Vector& column : A_columns) {
        if (column.size() != m) {
            throw std::invalid_argument(
                "leastSquares: matrix columns have inconsistent sizes"
            );
        }

        for (double value : column) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument(
                    "leastSquares: matrix contains a non-finite value"
                );
            }
        }
    }

    for (double value : b) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument(
                "leastSquares: vector contains a non-finite value"
            );
        }
    }

    // Factor A into Q and R using the group's QR implementation.
    QRResult qr = qrDecomposition(A_columns, tol);

    const std::vector<Vector>& Q = qr.Q;
    const std::vector<Vector>& R = qr.R;

    // Compute y = Q^T b.
    Vector y(n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t k = 0; k < m; ++k) {
            y[i] += Q[i][k] * b[k];
        }
    }

    // Solve Rx = y by back substitution.
    // R is stored as a list of columns.
    Vector x(n, 0.0);

    for (std::size_t i = n; i-- > 0;) {
        double sum = y[i];

        for (std::size_t j = i + 1; j < n; ++j) {
            sum -= R[j][i] * x[j];
        }

        const double diagonal = R[i][i];

        if (std::abs(diagonal) <= tol) {
            throw std::invalid_argument(
                "leastSquares: matrix is rank deficient"
            );
        }

        x[i] = sum / diagonal;
    }

    return x;
}

} // namespace numcomp
