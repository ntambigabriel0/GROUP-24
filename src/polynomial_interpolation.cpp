#include "polynomial_interpolation.hpp"
#include <cmath>
#include <stdexcept>

#include "gaussian.hpp"

using namespace std;

vector<double> polynomialInterpolation(const vector<double>& xs,
                                       const vector<double>& ys) {
    if (xs.size() != ys.size()) {
        throw invalid_argument("xs and ys must have the same size.");
    }
    if (xs.empty()) {
        throw invalid_argument("At least one data point is required.");
    }

    int n = static_cast<int>(xs.size());

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (xs[i] == xs[j]) {
                throw invalid_argument("x-values must be distinct.");
            }
        }
    }

    Matrix V(n, n);
    for (int i = 0; i < n; i++) {
        double power = 1.0;
        for (int j = 0; j < n; j++) {
            V.at(i, j) = power;
            power *= xs[i];
        }
    }

    return gaussianElimination(V, ys);
}

double evaluatePolynomial(const vector<double>& coeffs, double x) {
    double result = 0.0;
    for (int i = static_cast<int>(coeffs.size()) - 1; i >= 0; i--) {
        result = result * x + coeffs[i];
    }
    return result;
}
