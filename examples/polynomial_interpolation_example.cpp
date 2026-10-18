#include <numerical-computing/polynomial_interpolation.hpp>

#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace numcomp;

int main() {
    Vector xs = {0, 1, 2};
    Vector ys = {1, 3, 7};

    try {
        Vector coeffs = polynomialInterpolation(xs, ys);

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Coefficients [c0 + c1*x + c2*x^2 + ...]:\n";
        for (std::size_t i = 0; i < coeffs.size(); ++i) {
            std::cout << "  c" << i << " = " << coeffs[i] << "\n";
        }

        double x = 3.0;
        std::cout << "P(" << x << ") = " << evaluatePolynomial(coeffs, x) << "\n";
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
    }

    
    try {
        polynomialInterpolation({1, 1}, {2, 3});
    } catch (const std::invalid_argument& e) {
        std::cout << "Error caught: " << e.what() << "\n";
    }
    return 0;
}
