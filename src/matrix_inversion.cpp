#include <numerical-computing/matrix_inversion.hpp>
#include <stdexcept>

namespace numcomp {

Matrix invertMatrix(const Matrix& A) {
    if (A.empty()) {
        throw std::invalid_argument("invertMatrix: matrix is empty");
    }
    
    const std::size_t n = A.size();
    if (A[0].empty()) {
        throw std::invalid_argument("invertMatrix: matrix contains empty rows");
    }

    for (const Vector& row : A) {
        if (row.size() != n) {
            throw std::invalid_argument("invertMatrix: matrix must be square");
        }
    }

    // Factor the matrix once
    LUResult f = luDecompose(A);

    Matrix inverseA(n, Vector(n, 0.0));

    // Solve for each column of the identity matrix
    for (std::size_t col = 0; col < n; ++col) {
        Vector e(n, 0.0);
        e[col] = 1.0;

        Vector x = luSolve(f, e);

        for (std::size_t row = 0; row < n; ++row) {
            inverseA[row][col] = x[row];
        }
    }

    return inverseA;
}

}  // namespace numcomp