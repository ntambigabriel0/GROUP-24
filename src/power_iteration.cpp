#include "numerical-computing/power_iteration.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace {

const int MAX_ITERATIONS = 1000;
const double ZERO_LIMIT = 1e-12;

// Runs power iteration using one starting vector.
bool runIteration(
    int n,
    double A[20][20],
    const double start[20],
    double tolerance,
    double& eigenvalue,
    double result[20]
) {
    double current[20] = {0.0};
    double next[20] = {0.0};

    for (int i = 0; i < n; i++) {
        current[i] = start[i];
    }

    double largest = 0.0;

    for (int i = 0; i < n; i++) {
        largest = std::max(largest, std::fabs(current[i]));
    }

    if (largest < ZERO_LIMIT) {
        return false;
    }

    for (int i = 0; i < n; i++) {
        current[i] /= largest;
    }

    for (int iteration = 0;
         iteration < MAX_ITERATIONS;
         iteration++) {

        double largestNext = 0.0;

        for (int i = 0; i < n; i++) {
            next[i] = 0.0;

            for (int j = 0; j < n; j++) {
                next[i] += A[i][j] * current[j];
            }

            largestNext = std::max(
                largestNext, std::fabs(next[i])
            );
        }

        if (largestNext < ZERO_LIMIT) {
            return false;
        }

        // Calculate the Rayleigh quotient.
        double numerator = 0.0;
        double denominator = 0.0;

        for (int i = 0; i < n; i++) {
            numerator += current[i] * next[i];
            denominator += current[i] * current[i];
        }

        if (denominator < ZERO_LIMIT) {
            return false;
        }

        double lambda = numerator / denominator;

        // Check the eigenvalue residual.
        double residualSquared = 0.0;

        for (int i = 0; i < n; i++) {
            double difference =
                next[i] - lambda * current[i];

            residualSquared += difference * difference;
        }

        double residual = std::sqrt(residualSquared);

        if (residual <= tolerance *
                        (1.0 + std::fabs(lambda))) {

            eigenvalue = lambda;

            for (int i = 0; i < n; i++) {
                result[i] = current[i];
            }

            return true;
        }

        // Normalize the next vector.
        for (int i = 0; i < n; i++) {
            current[i] = next[i] / largestNext;
        }
    }

    return false;
}

}

double powerIteration(
    int n,
    double A[20][20],
    double x[20],
    double error
) {
    if (n < 1 || n > 20) {
        throw std::invalid_argument(
            "Matrix size must be between 1 and 20."
        );
    }

    if (error <= 0.0 || !std::isfinite(error)) {
        throw std::invalid_argument(
            "Tolerance must be positive and finite."
        );
    }

    bool hasNonZeroValue = false;

    for (int i = 0; i < n; i++) {

        if (!std::isfinite(x[i])) {
            throw std::invalid_argument(
                "Initial vector must contain finite values."
            );
        }

        if (std::fabs(x[i]) > ZERO_LIMIT) {
            hasNonZeroValue = true;
        }

        for (int j = 0; j < n; j++) {
            if (!std::isfinite(A[i][j])) {
                throw std::invalid_argument(
                    "Matrix must contain finite values."
                );
            }
        }
    }

    if (!hasNonZeroValue) {
        throw std::invalid_argument(
            "Initial vector cannot be zero."
        );
    }

    std::vector<double> eigenvalues;
    std::vector<std::vector<double>> eigenvectors;

    double value;
    double result[20] = {0.0};

    // First try the initial vector supplied by the user.
    if (runIteration(n, A, x, error, value, result)) {
        eigenvalues.push_back(value);
        eigenvectors.push_back(
            std::vector<double>(result, result + n)
        );
    }

    // Try alternative starting vectors to avoid missing
    // the dominant eigenvalue because of an orthogonal start.
    for (int k = 0; k < n; k++) {

        double start[20] = {0.0};
        start[k] = 1.0;

        if (runIteration(n, A, start, error, value, result)) {
            eigenvalues.push_back(value);
            eigenvectors.push_back(
                std::vector<double>(result, result + n)
            );
        }
    }

    if (eigenvalues.empty()) {
        throw std::runtime_error(
            "Power iteration failed to converge."
        );
    }

    // Select the eigenvalue with the greatest magnitude.
    int best = 0;

    for (int i = 1;
         i < static_cast<int>(eigenvalues.size());
         i++) {

        if (std::fabs(eigenvalues[i]) >
            std::fabs(eigenvalues[best])) {
            best = i;
        }
    }

    // Detect distinct eigenvalues with equal maximum magnitude.
    for (int i = 0;
         i < static_cast<int>(eigenvalues.size());
         i++) {

        if (i == best) {
            continue;
        }

        double scale = std::max(
            1.0,
            std::max(
                std::fabs(eigenvalues[i]),
                std::fabs(eigenvalues[best])
            )
        );

        if (std::fabs(
                std::fabs(eigenvalues[i]) -
                std::fabs(eigenvalues[best])
            ) <= error * scale
            &&
            std::fabs(
                eigenvalues[i] - eigenvalues[best]
            ) > error * scale) {

            throw std::runtime_error(
                "There is no unique dominant eigenvalue."
            );
        }
    }

    // Return the corresponding eigenvector through x.
    for (int i = 0; i < n; i++) {
        x[i] = eigenvectors[best][i];
    }

    return eigenvalues[best];
}
