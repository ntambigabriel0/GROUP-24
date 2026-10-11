# Polynomial Interpolation

### 1 How I built it?

I used the gaussian elimination direct linear solver. It uses row operations to zero out everything below the diagonal and swaps the largest pivot correctly first to avoid dividing by zero or tiny numbers.

Then from the polynomial interpolation where data points are stacked into Vandermonde system and solves for the coefficients using Gaussian elimination. 

The logic was validated using specific points such as (0,1), (1,3) and (2,7). The output is correct according to the polynomial and points to 13 for x = 3.

There is added error detection. It's guarded against singular matrices and validated that input vectors are non empty and of equal size.

### 2 What classes were created?
One core class: Matrix. Manages the 2D grid structure required for linear operations. It was organized this way to encapsulate 2D vectors inside a Matrix class and keeps the row swapping and indexing clean, preventing  raw index errors and abstractiing away the underlying std::vector<std::vector<double>> storage from the solver logic.

### 3 What data does it hold?
The rows and cols are int type and are private, they hold the number of rows in the matrix and number of columns in the matrix respectively.
Data is vector<vector<double>>, it is private and it is a 2D container for matrix elements.

### 4 What functions operate on that data?
Matrix::at(i, j), Matrix::rowCount(), colCount(), Matrix::swapRows(i, j), gaussianElimination(A, b), polynomialInterpolation(x_vals, y_vals), polynomialInterpolation(x_vals, y_vals) are all public and they perform operations on this data.

### 5 Why some members are public and others hidden?
Matrix dimensions and raw data vectors are private to prevent unintended corruption. Access is done through public functions such as at() and swapRows() making sure there is safe encapsulation while exposing the solver functions.

### 6 How the rest of the program uses it?
vector<double> x_vals = {0, 1, 2};
vector<double> y_vals = {1, 3, 7};
vector<double> coeffs = polynomialInterpolation(x_vals, y_vals);
double test_y = evaluatePolynomial(coeffs, 3.0);
It connects to others by using the standard C++ library headers ( <vector>, <cmath>, <stdexcept>).
It can be integrated into larger numerical computing libraries for curve fitting, data approximation or numerical analysis.

### 7 what issues were faced?
Instability in vandermonde Matrices.
Singular Matrix detection.

Other alternative formulations such as Newton interpolation and Lagrange interpolation are to be added to avoid Vandermonde ill-conditioned issues for larger datasets.


