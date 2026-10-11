
#define GAUSSIAN_HPP

#include <vector>
#include <stdexcept>
#include <cmath>

class Matrix {
private:
    int rows;
    int cols;
    std::vector<std::vector<double>> data;

public:
    Matrix(int r, int c);

    double& at(int i, int j);
    double at(int i, int j) const;

    int rowCount() const;
    int colCount() const;

    void swapRows(int i, int j);
};

std::vector<double> gaussianElimination(Matrix A, std::vector<double> b);

#endif // GAUSSIAN_HPP#ifndef NUMERICAL_COMPUTING_GAUSSIAN_HPP
#define NUMERICAL_COMPUTING_GAUSSIAN_HPP

#include <utility>
#include <vector>

// A dense matrix of doubles, stored as a vector of rows.
class Matrix {
private:
    int rows;                                  // number of rows
    int cols;                                  // number of columns
    std::vector<std::vector<double>> data;    // the entries, row by row

public:
    // Creates an r x c matrix filled with zeros.
    Matrix(int r, int c) : rows(r), cols(c) {
        data.resize(rows, std::vector<double>(cols, 0.0));
    }

    // Read/write access to the entry in row i, column j (no bounds checking).
    double& at(int i, int j) { return data[i][j]; }
    double at(int i, int j) const { return data[i][j]; }

    int rowCount() const { return rows; }
    int colCount() const { return cols; }

    // Swaps row i and row j.
    void swapRows(int i, int j) { std::swap(data[i], data[j]); }
};

// Solves A*x = b by Gaussian elimination with partial pivoting.
// A and b are taken by value, so the caller's copies are not changed.
//
// Throws std::invalid_argument if A is not square or b's size does not
// match A.
// Throws std::runtime_error if the matrix is singular (a pivot is smaller
// than 1e-12), meaning there is no unique solution.
std::vector<double> gaussianElimination(Matrix A, std::vector<double> b);

#endif  // NUMERICAL_COMPUTING_GAUSSIAN_HPP

