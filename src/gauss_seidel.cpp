#include "numerical-computing/gauss_seidel.hpp"
#include <cmath>
#include <stdexcept>
#include <vector>

namespace GaussSeidel {

bool solve(
    const std::vector<std::vector<double>>& A,
    const std::vector<double>& b,
    std::vector<double>& x,
    int maxIterations,
    double tolerance
) {

    int n = static_cast<int>(A.size());

    // Check matrix and vector dimensions
    if (n == 0 || b.size() != A.size() || x.size() != A.size()) {
        throw std::invalid_argument(
            "Matrix and vector dimensions are invalid."
        );
    }

    // Check iteration count and tolerance
    if (maxIterations <= 0 || tolerance <= 0 ||
        !std::isfinite(tolerance)) {
        throw std::invalid_argument(
            "Invalid maximum iterations or tolerance."
        );
    }

    // Check that the matrix is square and has valid values
    for (int i = 0; i < n; i++) {

        if (A[i].size() != A.size()) {
            throw std::invalid_argument(
                "Coefficient matrix must be square."
            );
        }

        if (std::fabs(A[i][i]) < 1e-12) {
            throw std::invalid_argument(
                "A diagonal element is zero."
            );
        }

        if (!std::isfinite(b[i]) || !std::isfinite(x[i])) {
            throw std::invalid_argument(
                "Input vectors must contain finite values."
            );
        }

        for (int j = 0; j < n; j++) {
            if (!std::isfinite(A[i][j])) {
                throw std::invalid_argument(
                    "Matrix must contain finite values."
                );
            }
        }
    }

    // Gauss-Seidel iteration
    for (int k = 0; k < maxIterations; k++) {

        double maxError = 0.0;

        for (int i = 0; i < n; i++) {

            double sum = b[i];

            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum -= A[i][j] * x[j];
                }
            }

            double newX = sum / A[i][i];

            if (!std::isfinite(newX)) {
                return false;
            }

            double error = std::fabs(newX - x[i]);

            if (error > maxError) {
                maxError = error;
            }

            x[i] = newX;
        }

        // Check convergence
        if (maxError < tolerance) {
            return true;
        }
    }

    // Maximum iterations reached without convergence
    return false;
}

} // namespace GaussSeidel
