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


// Gaussian Elimination with Partial Pivoting
vector<double> gaussianElimination(Matrix A, vector<double> b) {

    int n = A.rowCount();

    if (A.colCount() != n || b.size() != n) {
        throw invalid_argument(
            "Matrix must be square and match vector size."
        );
    }

    // Forward elimination
    for (int k = 0; k < n - 1; k++) {

        // Find the largest pivot
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

        // Eliminate values below the pivot
        for (int i = k + 1; i < n; i++) {

            double factor = A.at(i, k) / A.at(k, k);

            for (int j = k; j < n; j++) {
                A.at(i, j) -= factor * A.at(k, j);
            }
            b[i] -= factor * b[k];
        }
    }

    // Back substitution
    vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < n; j++) {
            sum += A.at(i, j) * x[j];
        }
        x[i] = (b[i] - sum) / A.at(i, i);
    }

    return x;
}

int main() {
    Matrix A(3, 3);
    A.at(0, 0) = 2.0; A.at(0, 1) = 1.0; A.at(0, 2) = -1.0;
    A.at(1, 0) = -3.0; A.at(1, 1) = -1.0; A.at(1, 2) = 2.0;
    A.at(2, 0) = -2.0; A.at(2, 1) = 1.0; A.at(2, 2) = 2.0;

    vector<double> b = {8.0, -11.0, -3.0};

    vector<double> x = gaussianElimination(A, b);

    cout << fixed << setprecision(4);
    for (double value : x) {
        cout << value << " ";
    }
    cout << endl;

    return 0;
}