#include <numerical-computing/polynomial_fitting.hpp>
#include <cmath>
#include <stdexcept>

namespace numcomp {

Vector polynomialFit(const std::vector<double>& xs, const std::vector<double>& ys, int degree) {
    if (xs.empty() || ys.empty()) {
        throw std::invalid_argument("polynomialFit: data vectors cannot be empty");
    }
    if (xs.size() != ys.size()) {
        throw std::invalid_argument("polynomialFit: x and y vectors must have the same size");
    }
    if (degree < 0) {
        throw std::invalid_argument("polynomialFit: degree cannot be negative");
    }
    if (xs.size() < static_cast<std::size_t>(degree + 1)) {
        throw std::invalid_argument("polynomialFit: not enough data points for the requested degree");
    }

    const std::size_t num_points = xs.size();
    
    // We need (degree + 1) columns: [x^0, x^1, x^2, ..., x^degree].
    std::vector<Vector> A_columns;
    A_columns.reserve(degree + 1);

    for (int j = 0; j <= degree; ++j) {
        Vector current_column(num_points, 0.0);
        for (std::size_t i = 0; i < num_points; ++i) {
            current_column[i] = std::pow(xs[i], j);
        }
        A_columns.push_back(current_column);
    }

    // Pass the constructed matrix to the Tier 3 solver
    return leastSquares(A_columns, ys);
}

} 