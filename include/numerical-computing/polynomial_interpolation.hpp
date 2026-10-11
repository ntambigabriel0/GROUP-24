#ifndef POLYNOMIAL_INTERPOLATION_HPP
#define POLYNOMIAL_INTERPOLATION_HPP

#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

class Matrix {
private:
    int rows;
    int cols;
    vector<vector<double>> data;

public:
    Matrix(int r, int c);
    double& at(int i, int j);
    double at(int i, int j) const;
    int rowCount() const;
    int colCount() const;
    void swapRows(int i, int j);
};

vector<double> gaussianElimination(Matrix A, vector<double> b);
vector<double> polynomialInterpolation(const vector<double>& x_vals, const vector<double>& y_vals);
double evaluatePolynomial(const vector<double>& coeffs, double x);

#endif 