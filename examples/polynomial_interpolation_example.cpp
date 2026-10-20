#include <iomanip>
#include <iostream>
#include <vector>

#include "polynomial_interpolation.hpp"

using namespace std;

int main() {
    vector<double> xs = {0.0, 1.0, 2.0};
    vector<double> ys = {1.0, 3.0, 7.0};

    try {
        vector<double> c = polynomialInterpolation(xs, ys);

        cout << fixed << setprecision(4);
        cout << "Coefficients (lowest power first): ";
        for (double v : c) cout << v << " ";
        cout << endl;

        cout << "p(3.0) = " << evaluatePolynomial(c, 3.0)
             << "  (expected 13)" << endl;
    } catch (const exception& e) {
        cout << "Error: " << e.what() << endl;
    }

    return 0;
}
