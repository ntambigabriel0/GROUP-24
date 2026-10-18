#ifndef NUMERICAL_COMPUTING_LAGRANGE_HPP  // include guard: stops the file being included twice
#define NUMERICAL_COMPUTING_LAGRANGE_HPP

#include <vector>

namespace numcomp {  // everything in the library lives in this namespace

/// Lagrange interpolation: finds the one polynomial of degree n-1 that passes
/// through n known points (xs[i], ys[i]), and returns its value at `x`.
///
/// Example: points (0,1), (1,3), (2,7) lie on y = x^2 + x + 1,
///          so lagrangeInterpolate({0,1,2}, {1,3,7}, 3) returns 13.
///
/// Throws std::invalid_argument if:
///   - xs is empty,
///   - xs and ys have different sizes,
///   - any two x values are equal (or equal to within rounding error), because one x
///     cannot have two different y values on a polynomial,
///   - any input value is NaN or infinite.
double lagrangeInterpolate(const std::vector<double>& xs,
                           const std::vector<double>& ys,
                           double x);

/// Same as above, but evaluates the polynomial at many x values at once.
std::vector<double> lagrangeInterpolate(const std::vector<double>& xs,
                                        const std::vector<double>& ys,
                                        const std::vector<double>& queries);

}  // namespace numcomp

#endif  // NUMERICAL_COMPUTING_LAGRANGE_HPP
