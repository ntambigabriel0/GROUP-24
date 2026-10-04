#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iomanip>

using namespace std;


// ===============================
// Matrix Class
// ===============================

class Matrix {
private:
    int rows;
    int cols;
    vector<vector<double>> data;

public:
    Matrix(int r, int c) : rows(r), cols(c) {
        data.resize(rows, vector<double>(cols, 0.0));
    }

    double& at(int i, int j) {
        return data[i][j];
    }

    double at(int i, int j) const {
        return data[i][j];
    }

    int rowCount() const {
        return rows;
    }

    int colCount() const {
        return cols;
    }

    void swapRows(int i, int j) {
        swap(data[i], data[j]);
    }
};


// ===============================
// Gaussian Elimination Algorithm
// ===============================

vector<double> gaussianElimination(Matrix A, vector<double> b) {

    int n = A.rowCount();

    if (A.colCount() != n || b.size() != n) {
        throw invalid_argument(
            "Matrix must be square and match vector size."
        );
    }

    // Forward elimination
    for (int k = 0; k < n - 1; k++) {

        // Find largest pivot
        int pivotRow = k;

        for (int i = k + 1; i < n; i++) {

            if (fabs(A.at(i, k)) > fabs(A.at(pivotRow, k))) {
                pivotRow = i;
            }
        }

        // Check for zero pivot
        if (fabs(A.at(pivotRow, k)) < 1e-12) {
            throw runtime_error(
                "Matrix is singular or has no unique solution."
            );
        }

        // Swap rows if necessary
        if (pivotRow != k) {
            A.swapRows(k, pivotRow);
            swap(b[k], b[pivotRow]);
        }

        // Eliminate values below pivot
        for (int i = k + 1; i < n; i++) {

            double factor = A.at(i, k) / A.at(k, k);

            for (int j = k; j < n; j++) {
                A.at(i, j) -= factor * A.at(k, j);
            }

            b[i] -= factor * b[k];
        }
    }

    // Check final pivot
    if (fabs(A.at(n - 1, n - 1)) < 1e-12) {
        throw runtime_error(
            "Matrix is singular or has no unique solution."
        );
    }

    // Back substitution
    vector<double> x(n);

    for (int i = n - 1; i >= 0; i--) {

        double sum = b[i];

        for (int j = i + 1; j < n; j++) {
            sum -= A.at(i, j) * x[j];
        }

        x[i] = sum / A.at(i, i);
    }

    return x;
}


// ===============================
// Compare Results
// ===============================

bool approximatelyEqual(double a, double b) {

    const double tolerance = 1e-9;

    return fabs(a - b) < tolerance;
}


// ===============================
// Test Function
// ===============================

void runTest(
    string testName,
    Matrix A,
    vector<double> b,
    vector<double> expected
) {

    cout << "TEST: " << testName << endl;

    try {

        vector<double> result =
            gaussianElimination(A, b);

        bool passed = true;

        for (int i = 0; i < result.size(); i++) {

            if (!approximatelyEqual(result[i], expected[i])) {
                passed = false;
            }
        }

        cout << fixed << setprecision(6);

        cout << "Expected: ";

        for (double value : expected) {
            cout << value << " ";
        }

        cout << endl;

        cout << "Actual:   ";

        for (double value : result) {
            cout << value << " ";
        }

        cout << endl;

        if (passed) {
            cout << "RESULT: PASS" << endl;
        }
        else {
            cout << "RESULT: FAIL" << endl;
        }

    }
    catch (const exception& e) {

        cout << "ERROR: " << e.what() << endl;
        cout << "RESULT: FAIL" << endl;
    }

    cout << "--------------------------------" << endl;
}


// ===============================
// Main Testing Program
// ===============================

int main() {

    // --------------------------------
    // TEST 1: 2 x 2 system
    // --------------------------------

    Matrix A1(2, 2);

    A1.at(0, 0) = 2;
    A1.at(0, 1) = 1;

    A1.at(1, 0) = 1;
    A1.at(1, 1) = 3;

    vector<double> b1 = {3, 5};

    vector<double> expected1 = {0.8, 1.4};

    runTest(
        "Test 1 - 2x2 System",
        A1,
        b1,
        expected1
    );


    // --------------------------------
    // TEST 2: 3 x 3 system
    // --------------------------------

    Matrix A2(3, 3);

    A2.at(0, 0) = 3;
    A2.at(0, 1) = 2;
    A2.at(0, 2) = -1;

    A2.at(1, 0) = 2;
    A2.at(1, 1) = -2;
    A2.at(1, 2) = 4;

    A2.at(2, 0) = -1;
    A2.at(2, 1) = 0.5;
    A2.at(2, 2) = -1;

    vector<double> b2 = {1, -2, 0};

    vector<double> expected2 = {1, -2, -2};

    runTest(
        "Test 2 - 3x3 System",
        A2,
        b2,
        expected2
    );


    // --------------------------------
    // TEST 3: Partial Pivoting
    // --------------------------------

    Matrix A3(2, 2);

    A3.at(0, 0) = 0;
    A3.at(0, 1) = 2;

    A3.at(1, 0) = 1;
    A3.at(1, 1) = 3;

    vector<double> b3 = {4, 7};

    vector<double> expected3 = {1, 2};

    runTest(
        "Test 3 - Partial Pivoting",
        A3,
        b3,
        expected3
    );


    // --------------------------------
    // TEST 4: Singular Matrix
    // --------------------------------

    cout << "TEST: Test 4 - Singular Matrix" << endl;

    Matrix A4(2, 2);

    A4.at(0, 0) = 1;
    A4.at(0, 1) = 2;

    A4.at(1, 0) = 2;
    A4.at(1, 1) = 4;

    vector<double> b4 = {3, 6};

    try {

        vector<double> result =
            gaussianElimination(A4, b4);

        cout << "RESULT: FAIL" << endl;
        cout << "The singular matrix was not detected."
             << endl;

    }
    catch (const exception& e) {

        cout << "Expected error detected:" << endl;
        cout << e.what() << endl;
        cout << "RESULT: PASS" << endl;
    }

    cout << "--------------------------------" << endl;

    cout << endl;
    cout << "All Gaussian Elimination tests completed."
         << endl;

    return 0;
}