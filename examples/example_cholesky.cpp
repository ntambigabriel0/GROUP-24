#include "numerical-computing/cholesky.hpp"

#include <iomanip>
#include <iostream>

int main() {
    numcomp::Matrix A = {
        {  4,  12, -16 },
        { 12,  37, -43 },
        { -16, -43,  98 }
    };

    try {
        numcomp::Matrix L =
            numcomp::choleskyDecomposition(A);

        std::cout << "Input matrix A:\n";

        for (const auto& row : A) {
            for (double value : row) {
                std::cout << std::setw(10)
                          << value << ' ';
            }
            std::cout << '\n';
        }

        std::cout << "\nCholesky factor L:\n";
        std::cout << std::fixed << std::setprecision(4);

        for (const auto& row : L) {
            for (double value : row) {
                std::cout << std::setw(10)
                          << value << ' ';
            }
            std::cout << '\n';
        }

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
