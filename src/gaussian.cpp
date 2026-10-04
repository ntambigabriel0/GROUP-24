#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iomanip>

using namespace std;

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


vector<double> gaussianElimination(Matrix A, vector<double> b) {

    int n = A.rowCount();

    if (A.colCount() != n || b.size() != n) {
        throw invalid_argument("Matrix must be square and match vector size.");
    }
    for (int k = 0; k < n - 1; k++) {

        
        int pivotRow = k;

        for (int i = k + 1; i < n; i++) {
            if (fabs(A.at(i, k)) > fabs(A.at(pivotRow, k))) {
                pivotRow = i;
            }
        }

        
        if (fabs(A.at(pivotRow, k)) < 1e-12) {
            throw runtime_error("Matrix is singular or has no unique solution.");
        }

        
        if (pivotRow != k) {
            A.swapRows(k, pivotRow);
            swap(b[k], b[pivotRow]);
        }

        
        for (int i = k + 1; i < n; i++) {

            double factor = A.at(i, k) / A.at(k, k);

            for (int j = k; j < n; j++) {
                A.at(i, j) -= factor * A.at(k, j);
            }

            b[i] -= factor * b[k];
        }
    }
    if (fabs(A.at(n - 1, n - 1)) < 1e-12) {
        throw runtime_error("Matrix is singular or has no unique solution.");
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


int main() {

    // Example:
    //
    // 2x + y = 3

    Matrix A(2, 2);

    A.at(0, 0) = 2;
    A.at(0, 1) = 1;

    A.at(1, 0) = 1;
    A.at(1, 1) = 3;

    vector<double> b = {3, 5};

    try {

        vector<double> x = gaussianElimination(A, b);

        cout << fixed << setprecision(4);

        cout << "Solution:" << endl;
        cout << "x = " << x[0] << endl;
        cout << "y = " << x[1] << endl;

    }
    catch (const exception& e) {
        cout << "Error: " << e.what() << endl;
    }

    return 0;
}
