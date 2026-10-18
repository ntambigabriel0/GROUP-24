#include <numerical-computing/least_squares.hpp>

#include <iomanip>
#include <iostream>

int main() {
    using namespace numcomp;

    // Fit y = a + bx to the three points.
    // The matrix is represented as column vectors.
    std::vector<Vector> A_columns = {
        {1.0, 1.0, 1.0},  // Column for a
        {0.0, 1.0, 2.0}   // Column for b
    };

    Vector b = {1.0, 2.0, 2.0};

    try {
        Vector x = leastSquares(A_columns, b);

        std::cout << std::fixed << std::setprecision(3);
        std::cout << "Least-squares solution:\n";
        std::cout << "a = " << x[0] << '\n';
        std::cout << "b = " << x[1] << '\n';
        std::cout << "Best-fit line: y = "
                  << x[0] << " + " << x[1] << "x\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
        return 1;
    }

    return 0;
}
