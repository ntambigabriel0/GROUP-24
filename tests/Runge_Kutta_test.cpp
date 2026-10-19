#include <iostream>
#include <cmath>
#include <cassert>
#include “algorithm.hpp”

void run_proof_tests() {
    std::cout << “Initializing RK4 Algorithmic Proof Tests...\n”;

    // Standard Math Test Model: dy/dt = y
    // Analytical correct solution curve: y(t) = e^t. Target checking: y(2) = e^2 ≈ 7.389056
    auto exponential_ode = [](double t, double y) { return y; };

    double t0 = 0.0, y0 = 1.0;
    double t_end = 2.0;
    double step_size = 0.1;

    // Instantiate object using class patterns
    RK4Solver solver(step_size);
    solver.solve(exponential_ode, t0, y0, t_end);

    // Pull tracking limits out of the arrays
    double final_calculated_y = solver.y_history.back();
    double exact_analytical_y = std::exp(2.0);
    double absolute_error = std::abs(final_calculated_y - exact_analytical_y);

    std::cout << “-> Computed RK4 Value at t=2: “ << final_calculated_y << “\n”;
    std::cout << “-> Exact Analytical Value:     “ << exact_analytical_y << “\n”;
    std::cout << “-> Absolute Deviation Error:    “ << absolute_error << “\n”;

    // accuracy verification bound check for Runge-Kutta 4 Order
    double tolerance = 1e-4;
    assert(absolute_error < tolerance && “RK4 failed accuracy verification test!”);

    // Boundary Edge Case Test
    RK4Solver boundary_solver(0.1);
    boundary_solver.solve(exponential_ode, 0.0, 5.0, 0.0);
    assert(boundary_solver.y_history.size() == 1 && boundary_solver.y_history.front() == 5.0);

    std::cout << “\u2705 All automated verification and proof tests passed successfully!\n\n”;
}

int main() {
    run_proof_tests();
    return 0;
}