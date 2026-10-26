#include <iostream>
#include <cmath>
using namespace std;

int main() {
    int n;
    const int MAX = 20;
    const int MAX_ITER = 100;
    const double TOL = 0.0001;

    double A[MAX][MAX];
    double b[MAX];
    double x[MAX] = {0};

    cout << "Enter the number of equations: ";
    cin >> n;

    if (n < 1 || n > MAX) {
        cout << "Invalid number of equations." << endl;
        return 1;
    }

    cout << "Enter the coefficients of the matrix A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    cout << "Enter the constants vector b:" << endl;
    for (int i = 0; i < n; i++) {
        cin >> b[i];
    }

    for (int i = 0; i < n; i++) {
        if (fabs(A[i][i]) < 1e-12) {
            cout << "Error: A diagonal element is zero."
                 << endl;
            return 1;
        }
    }

    bool converged = false;

    for (int k = 0; k < MAX_ITER; k++) {
        double maxError = 0.0;

        for (int i = 0; i < n; i++) {
            double sum = b[i];

            for (int j = 0; j < n; j++) {
                if (j != i) {
                    sum -= A[i][j] * x[j];
                }
            }

            double newX = sum / A[i][i];

            double error = fabs(newX - x[i]);

            if (error > maxError) {
                maxError = error;
            }

            x[i] = newX;
        }

        if (maxError < TOL) {
            converged = true;
            break;
        }
    }

    cout << "\nSolution:" << endl;

    for (int i = 0; i < n; i++) {
        cout << "x" << i + 1 << " = "
             << x[i] << endl;
    }

    if (converged) {
        cout << "Convergence achieved." << endl;
    } else {
        cout << "Maximum iterations reached." << endl;
    }

    return 0;
}
