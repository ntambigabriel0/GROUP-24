#include "polynomial_interpolation.hpp"
#include <iostream>
#include <vector>
#include <iomanip>

using namespace std;

int main() {
    vector<double> x_vals = {0, 1, 2};
    vector<double> y_vals = {1, 3, 7};

    try {
        vector<double> coeffs = polynomialInterpolation(x_vals, y_vals);

        cout << fixed << setprecision(4);
        cout << "Calculated Polynomial Coefficients [c0 + c1*x + c2*x^2 + ...]:\n";
        for (size_t i = 0; i < coeffs.size(); i++) {
            cout << "c" << i << " = " << coeffs[i] << "\n";
        }

        double test_x = 3.0;
        double test_y = evaluatePolynomial(coeffs, test_x);
        
        cout << "\nEvaluating polynomial at x = " << test_x << ":\n";
        cout << "P(" << test_x << ") = " << test_y << "\n";

    } catch (const exception& e) {
        cout << "Error: " << e.what() << endl;
    }

    return 0;
}