#include <iostream>
#include <iomanip>
#include <exception>

#include "numerical-computing/power_iteration.hpp"

using namespace std;

int main() {
    const int n = 2;

    double A[20][20] = {
        {2, 1},
        {1, 2}
    };

    double x[20] = {1, 0};

    try {
        double eigenvalue = powerIteration(n, A, x);

        cout << fixed << setprecision(4);

        cout << "Power Iteration Example" << endl;
        cout << "-----------------------" << endl;

        cout << "Dominant Eigenvalue: "
             << eigenvalue << endl;

        cout << "Corresponding Eigenvector: ";

        for (int i = 0; i < n; i++) {
            cout << x[i] << " ";
        }

        cout << endl;
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }

    return 0;
}
