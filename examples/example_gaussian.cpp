#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <numerical-computing/gaussian.hpp>

int main() {
    using numcomp::Matrix;
    using numcomp::Vector;
    using numcomp::gaussianElimination;

    Matrix A = {{ 2.0,  1.0, -1.0},
                {-3.0, -1.0,  2.0},
                {-2.0,  1.0,  2.0}};

    Vector b = {8.0, -11.0, -3.0};

    try {
        Vector x = gaussianElimination(A, b);

        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Solution:" << std::endl;
        for (std::size_t i = 0; i < x.size(); i++) {
            std::cout << "x" << i + 1 << " = " << x[i] << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }

    // A singular system is reported with an exception, not a wrong answer
    Matrix S = {{1, 2},
                {2, 4}};
    try {
        gaussianElimination(S, {3.0, 6.0});
    } catch (const std::exception& e) {
        std::cout << "Singular system: " << e.what() << std::endl;
    }

    return 0;
}
