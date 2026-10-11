#include "polynomial_interpolation.hpp"

// Matrix Member Functions
Matrix::Matrix(int r, int c) : rows(r), cols(c) {
    data.resize(rows, vector<double>(cols, 0.0));
}

double& Matrix::at(int i, int j) {
    return data[i][j];
}

double Matrix::at(int i, int j) const {
    return data[i][j];
}

int Matrix::rowCount() const {
    return rows;
}

int Matrix::colCount() const {
    return cols;
}

void Matrix::swapRows(int i, int j) {
    swap(data[i], data[j]);
}

// Numerical Methods Functions
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

vector<double> polynomialInterpolation(const vector<double>& x_vals, const vector<double>& y_vals) {
    if (x_vals.size() != y_vals.size() || x_vals.empty()) {
        throw invalid_argument("Vectors must be non-empty and of the same size.");
    }

    int n = x_vals.size();
    Matrix A(n, n);
    vector<double> b = y_vals;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A.at(i, j) = pow(x_vals[i], j);
        }
    }

    return gaussianElimination(A, b);
}

double evaluatePolynomial(const vector<double>& coeffs, double x) {
    double result = 0.0;
    for (size_t i = 0; i < coeffs.size(); i++) {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}