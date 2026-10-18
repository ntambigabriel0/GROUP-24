#include <iostream>
#include <cmath>
using namespace std;

int main() {
    int n, i, j, k;
    double A[20][20], x[20], y[20];
    double max, lambda, oldLambda;
    double error = 0.0001;

    cout << "Enter the order of the matrix: ";
    cin >> n;

    if (n < 1 || n > 20) {
        cout << "Invalid matrix size." << endl;
        return 1;
    }

    cout << "Enter the matrix elements:" << endl;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    cout << "Enter the initial vector:" << endl;
    for (i = 0; i < n; i++) {
        cin >> x[i];
    }

    max = 0;
    for (i = 0; i < n; i++) {
        if (fabs(x[i]) > max)
            max = fabs(x[i]);
    }

    if (max == 0) {
        cout << "Initial vector cannot be zero." << endl;
        return 1;
    }

    for (i = 0; i < n; i++)
        x[i] = x[i] / max;

    oldLambda = 0;

    for (k = 1; k <= 100; k++) {
        for (i = 0; i < n; i++) {
            y[i] = 0;
            for (j = 0; j < n; j++)
                y[i] += A[i][j] * x[j];
        }

        max = 0;
        for (i = 0; i < n; i++) {
            if (fabs(y[i]) > max)
                max = fabs(y[i]);
        }

        if (max == 0) {
            cout << "Iteration failed." << endl;
            return 1;
        }

        for (i = 0; i < n; i++)
            x[i] = y[i] / max;

        double top = 0, bottom = 0;

        for (i = 0; i < n; i++) {
            double sum = 0;
            for (j = 0; j < n; j++)
                sum += A[i][j] * x[j];

            top += x[i] * sum;
            bottom += x[i] * x[i];
        }

        lambda = top / bottom;

        cout << "Iteration " << k
             << " Eigenvalue = " << lambda << endl;

        if (k > 1 && fabs(lambda - oldLambda) < error)
            break;

        oldLambda = lambda;
    }

    cout << "Dominant eigenvalue = " << lambda << endl;
    cout << "Eigenvector = ";

    for (i = 0; i < n; i++)
        cout << x[i] << " ";

    cout << endl;

    return 0;
}

