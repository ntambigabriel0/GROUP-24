#include "numerical-computing/gaussian.hpp"
#include <cmath>
#include <stdexcept>
#include <vector>
#include <utility>

std::vector<double> gaussianElimination(Matrix A, std::vector<double> b) {

    int n = A.rowCount();

    // Check matrix dimensions
    if (A.colCount() != n || b.size() != static_cast<std::size_t>(n)) {
        throw std::invalid_argument(
            "Matrix must be square and match vector size."
        );
    }

    // Handle an empty matrix
    if (n == 0) {
        throw std::invalid_argument("Matrix cannot be empty.");
    }

    // Forward elimination with partial pivoting
    for (int k = 0; k < n - 1; k++) {

        int pivotRow = k;

        // Find the row with the largest pivot
        for (int i = k + 1; i < n; i++) {
            if (std::fabs(A.at(i, k)) > std::fabs(A.at(pivotRow, k))) {
                pivotRow = i;
            }
        }

        // Check for a zero pivot
        if (std::fabs(A.at(pivotRow, k)) < 1e-12) {
            throw std::runtime_error(
                "Matrix is singular or has no unique solution."
            );
        }

        // Swap rows when necessary
        if (pivotRow != k) {
            A.swapRows(k, pivotRow);
            std::swap(b[k], b[pivotRow]);
        }

        // Eliminate entries below the pivot
        for (int i = k + 1; i < n; i++) {

            double factor = A.at(i, k) / A.at(k, k);

            for (int j = k; j < n; j++) {
                A.at(i, j) -= factor * A.at(k, j);
            }

            b[i] -= factor * b[k];
        }
    }

    // Check the last pivot
    if (std::fabs(A.at(n - 1, n - 1)) < 1e-12) {
        throw std::runtime_error(
            "Matrix is singular or has no unique solution."
        );
    }

    // Back substitution
    std::vector<double> x(n);

    for (int i = n - 1; i >= 0; i--) {

        double sum = b[i];

        for (int j = i + 1; j < n; j++) {
            sum -= A.at(i, j) * x[j];
        }

        x[i] = sum / A.at(i, i);
    }

    return x;
}
