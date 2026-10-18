#include <numerical-computing/lagrange.hpp>

#include <cmath>     
#include <stdexcept>  

namespace numcomp {

// Checks the data once, before any maths is done
static void validate(const std::vector<double>& xs, const std::vector<double>& ys) {
    if (xs.empty()) {
        throw std::invalid_argument("lagrangeInterpolate: no data points");
    }
    if (xs.size() != ys.size()) {  // every x needs exactly one y
        throw std::invalid_argument(
            "lagrangeInterpolate: xs and ys must have the same size");
    }
    for (std::size_t i = 0; i < xs.size(); ++i) {
        if (!std::isfinite(xs[i]) || !std::isfinite(ys[i])) {
            throw std::invalid_argument(
                "lagrangeInterpolate: values must be finite numbers");
        }
    }
    // Two equal x values would make us divide by zero below
    for (std::size_t i = 0; i < xs.size(); ++i) {
        for (std::size_t j = i + 1; j < xs.size(); ++j) {
            // "Equal" is judged relative to the size of the numbers, so
            // very small or very large x values are treated fairly
            const double scale = std::fmax(std::fabs(xs[i]), std::fabs(xs[j]));
            if (std::fabs(xs[i] - xs[j]) <= 1e-14 * scale) {
                throw std::invalid_argument(
                    "lagrangeInterpolate: x values must be distinct");
            }
        }
    }
}

// Value at x, assuming the data has already been validated
static double evaluate(const std::vector<double>& xs,
                       const std::vector<double>& ys,
                       double x) {
    const std::size_t n = xs.size();
    double result = 0.0;

    // The polynomial is a sum of n pieces: ys[i] * L_i(x)
    for (std::size_t i = 0; i < n; ++i) {
        // L_i(x) equals 1 at xs[i] and 0 at every other xs[j]:
        //   L_i(x) = product over j != i of (x - xs[j]) / (xs[i] - xs[j])
        double basis = 1.0;
        for (std::size_t j = 0; j < n; ++j) {
            if (j != i) {
                basis *= (x - xs[j]) / (xs[i] - xs[j]);
            }
        }
        result += ys[i] * basis;  // add this piece
    }
    return result;
}

double lagrangeInterpolate(const std::vector<double>& xs,
                           const std::vector<double>& ys,
                           double x) {
    validate(xs, ys);
    return evaluate(xs, ys, x);
}

std::vector<double> lagrangeInterpolate(const std::vector<double>& xs,
                                        const std::vector<double>& ys,
                                        const std::vector<double>& queries) {
    validate(xs, ys);  // check once, then reuse for every query
    std::vector<double> out;
    out.reserve(queries.size());
    for (double q : queries) {
        out.push_back(evaluate(xs, ys, q));
    }
    return out;
}

}  // namespace numcomp
