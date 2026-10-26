#include <numerical-computing/newton_interpolation.hpp>
#include <iostream>
#include <vector>
#include <iomanip>

int main() {
    std::cout << "--- Newton Interpolation Example ---\n\n";

    // Dataset from the numerical methods project plan
    std::vector<double> xs = {0.0, 1.0, 2.0};
    std::vector<double> ys = {1.0, 3.0, 7.0};

    std::cout << "Data Points (x, y):\n";
    for (std::size_t i = 0; i < xs.size(); ++i) {
        std::cout << "  (" << xs[i] << ", " << ys[i] << ")\n";
    }

    try {
        // Calculate the divided differences
        std::vector<double> coeffs = numcomp::newtonCoefficients(xs, ys);
        
        std::cout << "\nDivided Differences (Coefficients):\n";
        for (std::size_t i = 0; i < coeffs.size(); ++i) {
            std::cout << "  c" << i << " = " << coeffs[i] << "\n";
        }

        // Evaluate the polynomial at x = 3
        double target_x = 3.0;
        double result = numcomp::newtonEvaluate(xs, coeffs, target_x);
        
        std::cout << "\n[Evaluation]\n";
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "P(" << target_x << ") = " << result << "\n";

    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
    }

    return 0;
}