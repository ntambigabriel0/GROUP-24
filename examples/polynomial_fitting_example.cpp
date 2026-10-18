#include <numerical-computing/polynomial_fitting.hpp>
#include <iostream>
#include <vector>
#include <iomanip>

int main() {
    std::cout << "--- Least-Squares Polynomial Fitting Example ---\n\n";

    // Dataset from the numerical methods project plan
    std::vector<double> xs = {0.0, 1.0, 2.0, 3.0};
    std::vector<double> ys = {1.0, 3.0, 5.0, 8.0};

    std::cout << "Data Points (x, y):\n";
    for (std::size_t i = 0; i < xs.size(); ++i) {
        std::cout << "  (" << xs[i] << ", " << ys[i] << ")\n";
    }

    try {
        // Fit a line (degree 1)
        numcomp::Vector line_coeffs = numcomp::polynomialFit(xs, ys, 1);
        
        std::cout << "\n[Degree 1 - Line of Best Fit]\n";
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "y = " << line_coeffs[0] << " + " << line_coeffs[1] << "x\n";

        // Fit a parabola (degree 2)
        numcomp::Vector quad_coeffs = numcomp::polynomialFit(xs, ys, 2);
        
        std::cout << "\n[Degree 2 - Best Fit Parabola]\n";
        std::cout << "y = " << quad_coeffs[0] 
                  << " + " << quad_coeffs[1] << "x"
                  << " + " << quad_coeffs[2] << "x^2\n";

    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
    }

    return 0;
}