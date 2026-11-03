#include <iostream>
#include <cmath>
#include <stdexcept>

#include "numerical-computing/power_iteration.hpp"

using namespace std;

bool approximatelyEqual(double a, double b, double tolerance = 0.001) {
    return fabs(a - b) < tolerance;
}

// Test 1: Check that an orthogonal initial vector does not
// prevent the algorithm from finding the dominant eigenvalue.
bool testOrthogonalInitialVector() {
    double A[20][20] = {};
    A[0][0] = 2;
    A[0][1] = 1;
    A[1][0] = 1;
    A[1][1] = 2;

    double x[20] = {-1, 1};

    try {
        double result = powerIteration(2, A, x);

        if (!approximatelyEqual(result, 3.0)) {
            cout << "FAIL: Expected eigenvalue 3, got "
                 << result << endl;
            return false;
        }

        cout << "PASS: Orthogonal initial vector test" << endl;
        return true;
    }
    catch (const exception& e) {
        cout << "FAIL: Unexpected exception: "
             << e.what() << endl;
        return false;
    }
}

// Test 2: Check that equal-magnitude eigenvalues
// are not incorrectly reported as zero.
bool testEqualMagnitudeEigenvalues() {
    double A[20][20] = {};
    A[0][0] = 2;
    A[1][1] = -2;

    double x[20] = {1, 1};

    try {
        powerIteration(2, A, x);

        cout << "FAIL: Expected runtime_error for "
             << "non-unique dominant eigenvalue." << endl;
        return false;
    }
    catch (const runtime_error&) {
        cout << "PASS: Equal-magnitude eigenvalues test" << endl;
        return true;
    }
    catch (const exception& e) {
        cout << "FAIL: Wrong exception: "
             << e.what() << endl;
        return false;
    }
}

// Test 3: Check that non-convergence produces runtime_error.
bool testNonConvergence() {
    double A[20][20] = {};
    A[0][1] = -1;
    A[1][0] = 1;

    double x[20] = {1, 0};

    try {
        powerIteration(2, A, x);

        cout << "FAIL: Expected runtime_error "
             << "when iteration does not converge." << endl;
        return false;
    }
    catch (const runtime_error&) {
        cout << "PASS: Non-convergence test" << endl;
        return true;
    }
    catch (const exception& e) {
        cout << "FAIL: Wrong exception: "
             << e.what() << endl;
        return false;
    }
}

// Test 4: Check a normal dominant eigenvalue calculation.
bool testNormalCalculation() {
    double A[20][20] = {};
    A[0][0] = 4;
    A[1][1] = 2;

    double x[20] = {1, 1};

    try {
        double result = powerIteration(2, A, x);

        if (!approximatelyEqual(result, 4.0)) {
            cout << "FAIL: Expected eigenvalue 4, got "
                 << result << endl;
            return false;
        }

        cout << "PASS: Normal calculation test" << endl;
        return true;
    }
    catch (const exception& e) {
        cout << "FAIL: Unexpected exception: "
             << e.what() << endl;
        return false;
    }
}

int main() {
    cout << "Running Power Iteration Tests..." << endl << endl;

    int passed = 0;
    int total = 4;

    if (testOrthogonalInitialVector()) passed++;
    if (testEqualMagnitudeEigenvalues()) passed++;
    if (testNonConvergence()) passed++;
    if (testNormalCalculation()) passed++;

    cout << endl;
    cout << "Tests passed: " << passed
         << " out of " << total << endl;

    if (passed == total) {
        cout << "All tests passed!" << endl;
        return 0;
    }

    cout << "Some tests failed." << endl;
    return 1;
}
