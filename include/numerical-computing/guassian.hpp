#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

class Matrix {
private:
    int rows, cols;
    vector<vector<double>> data;

public:
    Matrix(int r, int c) {
        rows = r;
        cols = c;
        data.resize(rows, vector<double>(cols, 0));
    }

    double& at(int i, int j) {
        return data[i][j];
    }

    int rowCount() {
        return rows;
    }

    int colCount() {
        return cols;
    }

    void swapRows(int i, int j) {
        swap(data[i], data[j]);
    }
};

vector<double> gaussianElimination(Matrix A, vector<double> b) {
    int n = A.rowCount();

    for (int k = 0; k < n - 1; k++) {

        int pivot = k;

        for (int i = k + 1; i < n; i++) {
            if (fabs(A.at(i, k)) > fabs(A.at(pivot, k))) {
                pivot = i;
            }
        }

        if (fabs(A.at(pivot, k)) < 0.000001) {
            cout << "No unique solution." << endl;
            return vector<double>();
        }

        if (pivot != k) {
            A.swapRows(k, pivot);
            swap(b[k], b[pivot]);
        }

        for (int i = k + 1; i < n; i++) {
            double factor = A.at(i, k) / A.at(k, k);

            for (int j = k; j < n; j++) {
                A.at(i, j) = A.at(i, j) - factor * A.at(k, j);
            }

            b[i] = b[i] - factor * b[k];
        }
    }

    vector<double> x(n);

    for (int i = n - 1; i >= 0; i--) {
        double sum = b[i];

        for (int j = i + 1; j < n; j++) {
            sum = sum - A.at(i, j) * x[j];
        }

        x[i] = sum / A.at(i, i);
    }

    return x;
}

int main() {

    Matrix A(2, 2);

    A.at(0, 0) = 2;
    A.at(0, 1) = 1;

    A.at(1, 0) = 1;
    A.at(1, 1) = 3;

    vector<double> b = {3, 5};

    vector<double> x = gaussianElimination(A, b);

    if (!x.empty()) {
        cout << "Solution:" << endl;
        cout << "x = " << x[0] << endl;
        cout << "y = " << x[1] << endl;
    }

    return 0;
}