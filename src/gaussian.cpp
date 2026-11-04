#include <numerical-computing/gaussian.hpp>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace numcomp {

Vector gaussianElimination(Matrix A, Vector b) {

    const std::size_t n = A.size();

    // Handle an empty matrix first, so nothing below reads out of bounds
    if (n == 0) {
        throw std::invalid_argument("Matrix cannot be empty.");
    }

    // Check matrix dimensions: A must be square and b must match
    if (b.size() != n) {
        throw std::invalid_argument(
            "Matrix must be square and match vector size."
        );
    }
    for (std::size_t i = 0; i < n; i++) {
        if (A[i].size() != n) {
            throw std::invalid_argument(
                "Matrix must be square and match vector size."
            );
        }
    }

    // Forward elimination with partial pivoting
    for (std::size_t k = 0; k + 1 < n; k++) {

        std::size_t pivotRow = k;

        // Find the row with the largest pivot
        for (std::size_t i = k + 1; i < n; i++) {
            if (std::fabs(A[i][k]) > std::fabs(A[pivotRow][k])) {
                pivotRow = i;
            }
        }

        // Check for a zero pivot
        if (std::fabs(A[pivotRow][k]) < 1e-12) {
            throw std::runtime_error(
                "Matrix is singular or has no unique solution."
            );
        }

        // Swap rows when necessary
        if (pivotRow != k) {
            std::swap(A[k], A[pivotRow]);
            std::swap(b[k], b[pivotRow]);
        }

        // Eliminate entries below the pivot
        for (std::size_t i = k + 1; i < n; i++) {

            double factor = A[i][k] / A[k][k];

            for (std::size_t j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }

            b[i] -= factor * b[k];
        }
    }

    // Check the last pivot
    if (std::fabs(A[n - 1][n - 1]) < 1e-12) {
        throw std::runtime_error(
            "Matrix is singular or has no unique solution."
        );
    }

    // Back substitution (i counts down from n to 1, row index is i - 1)
    Vector x(n);

    for (std::size_t i = n; i-- > 0;) {

        double sum = b[i];

        for (std::size_t j = i + 1; j < n; j++) {
            sum -= A[i][j] * x[j];
        }

        x[i] = sum / A[i][i];
    }

    return x;
}

}  // namespace numcomp
