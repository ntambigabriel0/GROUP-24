#ifndef NUMERICAL_COMPUTING_POLYNOMIAL_INTERPOLATION_HPP  // include guard
#define NUMERICAL_COMPUTING_POLYNOMIAL_INTERPOLATION_HPP

#include <vector>

#include <numerical-computing/gaussian.hpp>  

namespace numcomp {

// Polynomial interpolation via a Vandermonde system.
//
// Finds the one polynomial of degree n-1 that passes through the n points
// (xs[i], ys[i]). Each point gives one equation about the unknown
// coefficients; the equations are stacked into a Vandermonde system and
// solved with gaussianElimination.
//
// Returns the coefficients in increasing order of power:
//   result[0] + result[1]*x + result[2]*x^2 + ...
//
// Example: points (0,1), (1,3), (2,7) give {1, 1, 1}, i.e. x^2 + x + 1.
//
// Throws std::invalid_argument if:
//   - xs is empty,
//   - xs and ys have different sizes,
//   - any input value is NaN or infinite,
//   - any two x values are equal (or equal to within rounding error).
// Throws std::runtime_error if the Vandermonde system is numerically singular.
//
// WARNING: Vandermonde matrices become badly conditioned as the number of
// points grows (roughly beyond 10 to 15 points) or when x values are close
// together, so the coefficients lose accuracy. For many points, use
// lagrangeInterpolate (same polynomial, no linear system) instead.
Vector polynomialInterpolation(const Vector& xs, const Vector& ys);

// Evaluates coeffs[0] + coeffs[1]*x + coeffs[2]*x^2 + ... using Horner's
// method. An empty coefficient list is the zero polynomial and returns 0.
double evaluatePolynomial(const Vector& coeffs, double x);

}  // namespace numcomp

#endif  // NUMERICAL_COMPUTING_POLYNOMIAL_INTERPOLATION_HPP
