#include <iostream>
#include <cmath>
#include "power_iteration_example.h"
using namespace std;

double powerIteration(const PowerIterationTest &test) {
    int i, j;
    double x[2], y[2];
    double lambda = 0, oldLambda = 0;

    for (i = 0; i < test.n; i++)
        x[i] = test.initial[i];

    double largest = 0;

    for (i = 0; i < test.n; i++) {
        if (fabs(x[i]) > largest)
            largest = fabs(x[i]);
    }

    if (largest == 0)
        return NAN;

    for (i = 0; i < test.n; i++)
        x[i] /= largest;

    for (int k = 1; k <= 100; k++) {
        for (i = 0; i < test.n; i++) {
            y[i] = 0;

            for (j = 0; j < test.n; j++)
                y[i] += test.A[i][j] * x[j];
        }

        largest = 0;

        for (i = 0; i < test.n; i++) {
            if (fabs(y[i]) > largest)
                largest = fabs(y[i]);
        }

        if (largest == 0)
            return NAN;

        for (i = 0; i < test.n; i++)
            x[i] = y[i] / largest;

        double top = 0, bottom = 0;

        for (i = 0; i < test.n; i++) {
            double sum = 0;

            for (j = 0; j < test.n; j++)
                sum += test.A[i][j] * x[j];

            top += x[i] * sum;
            bottom += x[i] * x[i];
        }

        lambda = top / bottom;

        cout << "Iteration " << k
             << ": Eigenvalue = " << lambda << endl;

        if (k > 1 && fabs(lambda - oldLambda) < 0.0001)
            return lambda;

        oldLambda = lambda;
    }

    return lambda;
}

int main() {
    cout << "Running: " << example.name << endl;

    double result = powerIteration(example);

    cout << "\nExpected eigenvalue: "
         << example.expectedEigenvalue << endl;

    cout << "Calculated eigenvalue: "
         << result << endl;

    if (isfinite(result) &&
        fabs(result - example.expectedEigenvalue)
        < example.tolerance) {
        cout << "TEST PASSED" << endl;
    } else {
        cout << "TEST FAILED" << endl;
    }

    return 0;
}
