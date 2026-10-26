#include <iostream>
#include <vector>
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

    if (success) {
        cout << "Gauss-Seidel Solution:" << endl;

        for (int i = 0; i < 3; i++) {
            cout << "x" << i + 1 << " = "
                 << x[i] << endl;
        }
    } else {
        cout << "The method did not converge "
             << "or the input is invalid." << endl;
    }

    return 0;
}
