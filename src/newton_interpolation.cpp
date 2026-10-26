#include <numerical-computing/newton_interpolation.hpp>
#include <stdexcept>
#include <cmath>

namespace numcomp {

std::vector<double> newtonCoefficients(const std::vector<double>& xs, const std::vector<double>& ys) {
    if (xs.empty() || ys.empty()) {
        throw std::invalid_argument("newtonCoefficients: data vectors cannot be empty");
    }
    if (xs.size() != ys.size()) {
        throw std::invalid_argument("newtonCoefficients: x and y vectors must have the same size");
    }

    // We use a 1D array to save memory, overwriting it from the back to compute the divided differences
    std::vector<double> coeffs = ys;
    const std::size_t n = xs.size();

    for (std::size_t j = 1; j < n; ++j) {
        for (std::size_t i = n - 1; i >= j; --i) {
            double denominator = xs[i] - xs[i - j];
            if (std::fabs(denominator) < 1e-14) {
                throw std::invalid_argument("newtonCoefficients: x values must be distinct");
            }
            coeffs[i] = (coeffs[i] - coeffs[i - 1]) / denominator;
        }
    }

    return coeffs;
}

double newtonEvaluate(const std::vector<double>& xs, const std::vector<double>& coeffs, double target_x) {
    if (xs.empty() || coeffs.empty()) {
        throw std::invalid_argument("newtonEvaluate: data vectors cannot be empty");
    }
    if (xs.size() != coeffs.size()) {
        throw std::invalid_argument("newtonEvaluate: x and coefficients vectors must have the same size");
    }

    // Evaluate using a nested multiplication loop (Horner's method)
    double result = coeffs.back();
    for (int i = static_cast<int>(coeffs.size()) - 2; i >= 0; --i) {
        result = result * (target_x - xs[i]) + coeffs[i];
    }

    return result;
}

} 