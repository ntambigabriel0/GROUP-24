#include <iostream>
#include <vector>
#include <cmath>
#include "gauss_seidel.hpp"

using namespace std;

int main() {
    vector<vector<double>> A = {
        {4, 1, 1},
        {2, 5, 1},
        {1, 1, 3}
    };

    vector<double> b = {6, 12, 10};
    vector<double> x(3, 0.0);

    bool success = GaussSeidel::solve(A, b, x);

    if (!success) {
        cout << "TEST FAILED: Solver did not converge."
             << endl;
        return 1;
    }

    double tolerance = 0.001;

    // Check the first equation
    double r1 = 4*x[0] + x[1] + x[2];

    // Check the second equation
    double r2 = 2*x[0] + 5*x[1] + x[2];

    // Check the third equation
    double r3 = x[0] + x[1] + 3*x[2];

    if (fabs(r1 - 6) < tolerance &&
        fabs(r2 - 12) < tolerance &&
        fabs(r3 - 10) < tolerance) {

        cout << "TEST PASSED!" << endl;
        cout << "All equations are satisfied." << endl;

    } else {
        cout << "TEST FAILED!" << endl;
        cout << "The calculated solution is inaccurate."
             << endl;
        return 1;
    }

    cout << "x1 = " << x[0] << endl;
    cout << "x2 = " << x[1] << endl;
    cout << "x3 = " << x[2] << endl;

    return 0;
}
