#include "polynomial_interpolation.hpp"
#include <iostream>
#include <vector>
#include <iomanip>
#include <cassert>

using namespace std;

void runTests() {
    cout << "Running Polynomial Interpolation Tests...\n";

    // Test case: points (0, 1), (1, 3), (2, 7)
    // Expected polynomial: P(x) = 1 + 2x + x^2 (since P(0)=1, P(1)=3, P(2)=7)
    vector<double> x_vals = {0.0, 1.0, 2.0};
    vector<double> y_vals = {1.0, 3.0, 7.0};

    vector<double> coeffs = polynomialInterpolation(x_vals, y_vals);

    // Check coefficient values with tolerance
    assert(coeffs.size() == 3);
    assert(fabs(coeffs[0] - 1.0) < 1e-6); // c0
    assert(fabs(coeffs[1] - 2.0) < 1e-6); // c1
    assert(fabs(coeffs[2] - 1.0) < 1e-6); // c2

    // Test evaluation
    double test_x = 3.0;
    double expected_y = 1.0 + 2.0 * (3.0) + pow(3.0, 2); // 1 + 6 + 9 = 16
    double calculated_y = evaluatePolynomial(coeffs, test_x);

    assert(fabs(calculated_y - expected_y) < 1e-6);

    cout << "All tests passed successfully!\n\n";
}

int main() {
    try {
        runTests();

        // Standard execution example
        vector<double> x_vals = {0, 1, 2};
        vector<double> y_vals = {1, 3, 7};
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
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}