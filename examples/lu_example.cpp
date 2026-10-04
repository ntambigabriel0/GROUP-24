#include <numerical-computing/lu.h>

#include <iostream>
#include <stdexcept>

// Prints a matrix one row per line
static void printMatrix(const numcomp::Matrix& m) {
    for (const auto& row : m) {
        for (double x : row) std::cout << x << " ";
        std::cout << "\n";
    }
}

int main() {
    // The example from the project plan
    numcomp::Matrix A = {{2, 1}, {1, 3}};

    // Split A into a lower triangle L and an upper triangle U
    numcomp::LUResult f = numcomp::luDecompose(A);

    std::cout << "L:\n";
    printMatrix(f.L);
    std::cout << "U:\n";
    printMatrix(f.U);
    std::cout << "swaps: " << f.swaps << "\n";

    // Solve A*x = b using the stored factors: forward through L, backward through U
    numcomp::Vector x = numcomp::luSolve(f, {3, 5});
    std::cout << "x = " << x[0] << ", y = " << x[1] << "\n";

    // The point of LU: a new right-hand side reuses the same factors
    numcomp::Vector x2 = numcomp::luSolve(f, {1, 0});
    std::cout << "second solve: x = " << x2[0] << ", y = " << x2[1] << "\n";

    // Bad input: the second row is twice the first, so the matrix is singular
    try {
        numcomp::luDecompose({{1, 2}, {2, 4}});
    } catch (const numcomp::SingularMatrixError& e) {
        std::cout << "Error caught: " << e.what() << "\n";
    }
    return 0;
}
