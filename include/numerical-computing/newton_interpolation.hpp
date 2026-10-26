#ifndef NUMERICAL_COMPUTING_NEWTON_INTERPOLATION_HPP
#define NUMERICAL_COMPUTING_NEWTON_INTERPOLATION_HPP

#include <vector>

namespace numcomp {

// Computes the divided differences (coefficients) for Newton's interpolation polynomial.
// Throws std::invalid_argument if vectors are empty, mismatched, or contain duplicate x values.
std::vector<double> newtonCoefficients(const std::vector<double>& xs, const std::vector<double>& ys);

// Evaluates the Newton polynomial at a specific x value using Horner's method.
// Throws std::invalid_argument if the input vectors are empty or mismatched.
double newtonEvaluate(const std::vector<double>& xs, const std::vector<double>& coeffs, double target_x);

} 

#endif 