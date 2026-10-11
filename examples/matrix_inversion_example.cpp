#include <numerical-computing/matrix_inversion.hpp>
#include <iostream>

static void printMatrix(const numcomp::Matrix& m) {
    for (const auto& row : m) {
        for (double x : row) std::cout << x << "\t";
        std::cout << "\n";
    }
}

int main() {
    numcomp::Matrix A = {{2, 1}, {1, 3}};

    std::cout << "Original Matrix A:\n";
    printMatrix(A);
    
    try {
        numcomp::Matrix inv = numcomp::invertMatrix(A);
        std::cout << "\nInverse Matrix A^-1:\n";
        printMatrix(inv);
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << "\n";
    }

    try {
        std::cout << "\nAttempting to invert a singular matrix...\n";
        numcomp::invertMatrix({{1, 2}, {2, 4}});
    } catch (const numcomp::SingularMatrixError& e) {
        std::cout << "Caught expected error: " << e.what() << "\n";
    }

    return 0;
}