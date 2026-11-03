#ifndef NUMERICAL_COMPUTING_GAUSSIAN_HPP
#define NUMERICAL_COMPUTING_GAUSSIAN_HPP

#include <utility>
#include <vector>

// A dense matrix of doubles, stored as a vector of rows.
class Matrix {
private:
    int rows;
    int cols;
    std::vector<std::vector<double>> data;

public:
    // Creates an r x c matrix filled with zeros.
    Matrix(int r, int c) : rows(r), cols(c) {
        data.resize(rows, std::vector<double>(cols, 0.0));
    }

    // Read/write access to matrix elements.
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

    // Swaps two rows.
    void swapRows(int i, int j) {
        std::swap(data[i], data[j]);
    }
};

// Solves Ax = b using Gaussian elimination with partial pivoting.
std::vector<double> gaussianElimination(Matrix A, std::vector<double> b);

#endif // NUMERICAL_COMPUTING_GAUSSIAN_HPP
