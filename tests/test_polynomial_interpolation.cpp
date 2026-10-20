#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "polynomial_interpolation.hpp"

using namespace std;

static int failures = 0;

static bool close(double a, double b, double tol = 1e-8) {
    return fabs(a - b) < tol;
}

static void check(bool ok, const char* name) {
    cout << (ok ? "  PASS: " : "  FAIL: ") << name << endl;
    if (!ok) failures++;
}

static bool coeffsMatch(const vector<double>& got, const vector<double>& expected) {
    if (got.size() != expected.size()) return false;
    for (size_t i = 0; i < got.size(); i++) {
        if (!close(got[i], expected[i])) return false;
    }
    return true;
}

int main() {
    cout << "--- POLYNOMIAL INTERPOLATION TESTS ---\n" << endl;

    cout << "[Test 0] Documented example" << endl;
    {
        vector<double> c = polynomialInterpolation({0.0, 1.0, 2.0}, {1.0, 3.0, 7.0});
        check(coeffsMatch(c, {1.0, 1.0, 1.0}), "coefficients are {1, 1, 1}");
        check(close(evaluatePolynomial(c, 3.0), 13.0), "p(3) = 13");
    }

    cout << "[Test 1] Single point" << endl;
    {
        vector<double> c = polynomialInterpolation({5.0}, {7.0});
        check(coeffsMatch(c, {7.0}), "constant polynomial p(x) = 7");
    }

    cout << "[Test 2] Linear (2 points)" << endl;
    {
        vector<double> c = polynomialInterpolation({0.0, 1.0}, {1.0, 3.0});
        check(coeffsMatch(c, {1.0, 2.0}), "coefficients are {1, 2}");
    }

    cout << "[Test 3] Quadratic (3 points)" << endl;
    {
        vector<double> c = polynomialInterpolation({0.0, 1.0, 3.0}, {2.0, 0.0, 2.0});
        check(coeffsMatch(c, {2.0, -3.0, 1.0}), "coefficients are {2, -3, 1}");
        check(close(evaluatePolynomial(c, 2.0), 0.0), "p(2) = 0");
    }

    cout << "[Test 4] Cubic, unsorted x-values" << endl;
    {
        vector<double> xs = {2.0, -1.0, 0.0, 1.0};
        vector<double> ys;
        for (double x : xs) ys.push_back(x * x * x - 2 * x + 1);
        vector<double> c = polynomialInterpolation(xs, ys);
        check(coeffsMatch(c, {1.0, -2.0, 0.0, 1.0}), "coefficients are {1, -2, 0, 1}");
        check(close(evaluatePolynomial(c, 3.0), 22.0), "predicts p(3) = 22 (outside data)");
    }

    cout << "[Test 5] Passes through all data points" << endl;
    {
        vector<double> xs = {-2.0, -0.5, 1.0, 2.5, 4.0};
        vector<double> ys = {3.0, -1.0, 4.0, 1.5, -2.0};
        vector<double> c = polynomialInterpolation(xs, ys);
        bool ok = true;
        for (size_t i = 0; i < xs.size(); i++) {
            if (!close(evaluatePolynomial(c, xs[i]), ys[i], 1e-6)) ok = false;
        }
        check(ok, "p(xs[i]) == ys[i] for all i");
    }

    cout << "[Test 6] Duplicate x-values" << endl;
    {
        bool caught = false;
        try {
            polynomialInterpolation({1.0, 2.0, 1.0}, {1.0, 2.0, 3.0});
        } catch (const invalid_argument&) {
            caught = true;
        }
        check(caught, "throws invalid_argument");
    }

    cout << "[Test 7] Mismatched sizes" << endl;
    {
        bool caught = false;
        try {
            polynomialInterpolation({1.0, 2.0}, {1.0});
        } catch (const invalid_argument&) {
            caught = true;
        }
        check(caught, "throws invalid_argument");
    }

    cout << "[Test 8] Empty input" << endl;
    {
        bool caught = false;
        try {
            polynomialInterpolation({}, {});
        } catch (const invalid_argument&) {
            caught = true;
        }
        check(caught, "throws invalid_argument");
    }

    cout << "[Test 9] Evaluate empty polynomial" << endl;
    check(close(evaluatePolynomial({}, 3.0), 0.0), "returns 0");

    cout << "\nAll tests finished running. Failures: " << failures << endl;
    return failures == 0 ? 0 : 1;
}
