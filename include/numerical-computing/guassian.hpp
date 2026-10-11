#ifndef GAUSSIAN_HPP
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

#endif // GAUSSIAN_HPP
