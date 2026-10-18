#ifndef NUMERICAL_COMPUTING_POLYNOMIAL_FITTING_HPP
#define NUMERICAL_COMPUTING_POLYNOMIAL_FITTING_HPP

#include <numerical-computing/least_squares.hpp>
#include <vector>

namespace numcomp {

// Fits a polynomial of the specified degree to a set of x and y data points
// using the least-squares method.
//
// Returns a vector of coefficients [a0, a1, a2, ...] such that:
// y = a0 + a1*x + a2*x^2 + ...
//
// Throws std::invalid_argument if:
// - The data vectors are empty or mismatched in size.
// - The requested degree is negative.
// - There are fewer data points than (degree + 1).
Vector polynomialFit(const std::vector<double>& xs, const std::vector<double>& ys, int degree);

} 

#endif // NUMERICAL_COMPUTING_POLYNOMIAL_FITTING_HPP