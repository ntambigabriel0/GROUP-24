#ifndef NUMERICAL_COMPUTING_GAUSSIAN_QUADRATURE_HPP  // include guard: stops the file being included twice
#define NUMERICAL_COMPUTING_GAUSSIAN_QUADRATURE_HPP

#include <functional>
#include <vector>

namespace numcomp {  // everything in the library lives in this namespace

// Largest number of points the rule generator accepts.
constexpr int kMaxQuadraturePoints = 100;

// A Gauss-Legendre rule on [-1, 1]: integral of f ~ sum of weights[i] * f(nodes[i]).
struct QuadratureRule {
    std::vector<double> nodes;    // the sample points, smallest first
    std::vector<double> weights;  // the matching weights (they add up to 2)
};

// Builds the n-point Gauss-Legendre rule on [-1, 1].
// The nodes are the roots of the Legendre polynomial P_n, found with Newton's method.
// An n-point rule integrates every polynomial of degree 2n - 1 or lower exactly.
// Throws std::invalid_argument if n < 1 or n > kMaxQuadraturePoints.
QuadratureRule gaussLegendreRule(int n);

// Integral of f over [a, b] using one n-point Gauss-Legendre rule.
// Works when a > b (the result changes sign) and returns 0 when a == b.
// Throws std::invalid_argument if f is empty, a or b is not finite, or n is out of range.
double gaussLegendre(const std::function<double(double)>& f,
                     double a, double b,
                     int n = 5);  // default number of points if the user gives none

// Integral of f over [a, b]: splits [a, b] into equal panels and applies an
// n-point Gauss-Legendre rule on each one. Use it when f is not smooth enough
// for a single rule over the whole interval.
// Throws std::invalid_argument if panels < 1, plus the same checks as gaussLegendre.
double compositeGaussLegendre(const std::function<double(double)>& f,
                              double a, double b,
                              int panels,
                              int n = 3);  // default points per panel if the user gives none

}  // namespace numcomp

#endif  // NUMERICAL_COMPUTING_GAUSSIAN_QUADRATURE_HPP
