#ifndef NUMERICAL_COMPUTING_POLYNOMIAL_INTERPOLATION_HPP
#define NUMERICAL_COMPUTING_POLYNOMIAL_INTERPOLATION_HPP

#include <vector>

std::vector<double> polynomialInterpolation(const std::vector<double>& xs,
                                            const std::vector<double>& ys);

double evaluatePolynomial(const std::vector<double>& coeffs, double x);

#endif  // NUMERICAL_COMPUTING_POLYNOMIAL_INTERPOLATION_HPP
