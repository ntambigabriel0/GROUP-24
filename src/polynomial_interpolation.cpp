#include <numerical-computing/polynomial_interpolation.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace numcomp {

Vector polynomialInterpolation(const Vector& xs, const Vector& ys) {
    const std::size_t n = xs.size();

    // ---- Check the input ----
    if (n == 0) {
        throw std::invalid_argument("polynomialInterpolation: no points given");
    }
    if (ys.size() != n) {
        throw std::invalid_argument("polynomialInterpolation: xs and ys must be the same size");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(xs[i]) || !std::isfinite(ys[i])) {
            throw std::invalid_argument("polynomialInterpolation: inputs must be finite numbers");
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {
            const double scale = std::max({1.0, std::fabs(xs[i]), std::fabs(xs[j])});
            if (std::fabs(xs[i] - xs[j]) <= 1e-12 * scale) {
                throw std::invalid_argument(
                    "polynomialInterpolation: x values must be distinct");
            }
        }
    }

    // ---- Build the Vandermonde matrix: row i is [1, x_i, x_i^2, ...] ----
    Matrix V(n, Vector(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        double power = 1.0;
        for (std::size_t j = 0; j < n; ++j) {
            V[i][j] = power;
            power *= xs[i];
        }
    }

    // ---- Solve V * coeffs = ys with the library's Gaussian elimination ----
    try {
        return gaussianElimination(V, ys);
    } catch (const std::runtime_error&) {
        throw std::runtime_error(
            "polynomialInterpolation: system is numerically singular "
            "(too many points, or x values too close together)");
    }
}

double evaluatePolynomial(const Vector& coeffs, double x) {
    double result = 0.0;
    for (std::size_t i = coeffs.size(); i-- > 0;) {
        result = result * x + coeffs[i];  // Horner's method
    }
    return result;
}

}  // namespace numcomp
