#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <numerical-computing/gauss_seidel.hpp>

using namespace std;

// Test 1: Check whether the solver gives a correct solution
void testSuccessfulSolution() {
    vector<vector<double>> A = {
        {4, 1, 1},
        {2, 5, 1},
        {1, 1, 3}
    };

    vector<double> b = {6, 12, 10};
    vector<double> x(3, 0.0);

    bool success = GaussSeidel::solve(A, b, x);

    double tolerance = 0.001;

    if (success &&
        fabs(4*x[0] + x[1] + x[2] - 6) < tolerance &&
        fabs(2*x[0] + 5*x[1] + x[2] - 12) < tolerance &&
        fabs(x[0] + x[1] + 3*x[2] - 10) < tolerance) {

        cout << "Successful solution test: PASSED" << endl;
    } else {
        cout << "Successful solution test: FAILED" << endl;
    }
}

// Test 2: Check whether invalid input throws an exception
void testInvalidInput() {
    vector<vector<double>> A = {
        {4, 1},
        {2, 5}
    };

    vector<double> b = {6};
    vector<double> x(2, 0.0);

    try {
        GaussSeidel::solve(A, b, x);
        cout << "Invalid input test: FAILED" << endl;
    }
    catch (const invalid_argument&) {
        cout << "Invalid input test: PASSED" << endl;
    }
}

// Test 3: Check whether the solver returns false when it cannot converge
void testNonConvergence() {
    vector<vector<double>> A = {
        {1, 2},
        {2, 1}
    };

    vector<double> b = {3, 3};
    vector<double> x(2, 0.0);

    bool success = GaussSeidel::solve(A, b, x, 1, 0.000001);

    if (!success) {
        cout << "Non-convergence test: PASSED" << endl;
    } else {
        cout << "Non-convergence test: FAILED" << endl;
    }
}

int main() {
    cout << "--- GAUSS-SEIDEL TESTS ---" << endl << endl;

    testSuccessfulSolution();
    testInvalidInput();
    testNonConvergence();

    return 0;
}
