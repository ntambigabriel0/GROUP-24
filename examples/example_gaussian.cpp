#include <iostream>
#include <vector>
#include <iomanip>
#include <stdexcept>
#include "numerical-computing/gaussian.hpp"

int main() {

    Matrix A(3, 3);

    A.at(0, 0) = 2.0;
    A.at(0, 1) = 1.0;
    A.at(0, 2) = -1.0;

    A.at(1, 0) = -3.0;
    A.at(1, 1) = -1.0;
    A.at(1, 2) = 2.0;

    A.at(2, 0) = -2.0;
    A.at(2, 1) = 1.0;
    A.at(2, 2) = 2.0;

    std::vector<double> b = {8.0, -11.0, -3.0};

    try {

        std::vector<double> x = gaussianElimination(A, b);

        std::cout << std::fixed << std::setprecision(4);

        std::cout << "Solution:" << std::endl;

        for (int i = 0; i < static_cast<int>(x.size()); i++) {
            std::cout << "x" << i + 1 << " = " << x[i] << std::endl;
        }

    }
    catch (const std::exception& e) {

        std::cout << "Error: " << e.what() << std::endl;

    }

    return 0;
}
