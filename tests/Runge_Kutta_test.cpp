#include <iostream>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <numerical-computing/Runge_Kutta.hpp>

// Always-on check (unlike assert, it is not removed by -DNDEBUG)
static int failures = 0;
#define CHECK(cond, msg)                                              \
    do {                                                              \
        if (!(cond)) {                                                \
            std::cerr << "FAIL: " << (msg) << "  [" #cond "]\n";      \
            ++failures;                                               \
        }                                                             \
    } while (0)

void run_proof_tests() {
    std::cout << "Initializing RK4 Algorithmic Proof Tests...\n";

    // dy/dt = y, exact solution y(t) = e^t
    auto exponential_ode = [](double t, double y) { (void)t; return y; };

    double t0 = 0.0, y0 = 1.0;
    double t_end = 2.0;
    double step_size = 0.1;

    numcomp::RK4Solver solver(step_size);
    solver.solve(exponential_ode, t0, y0, t_end);

    double final_calculated_y = solver.y_history.back();
    double exact_analytical_y = std::exp(2.0);
    double absolute_error = std::abs(final_calculated_y - exact_analytical_y);

    std::cout << "-> Computed RK4 Value at t=2: " << final_calculated_y << "\n";
    std::cout << "-> Exact Analytical Value:     " << exact_analytical_y << "\n";
    std::cout << "-> Absolute Deviation Error:    " << absolute_error << "\n";

    // Accuracy check
    CHECK(absolute_error < 1e-4, "RK4 failed accuracy verification test!");

    // Endpoint lands exactly on t_end
    CHECK(solver.t_history.back() == t_end, "final t must equal t_end exactly");
    CHECK(solver.t_history.size() == solver.y_history.size(), "history sizes must match");

    // Step count: 0 -> 1 at h = 0.1 must give 11 points, not 12
    numcomp::RK4Solver count_solver(0.1);
    count_solver.solve(exponential_ode, 0.0, 1.0, 1.0);
    CHECK(count_solver.t_history.size() == 11, "0->1 at h=0.1 must have exactly 11 points");
    CHECK(count_solver.t_history.back() == 1.0, "0->1 must end exactly at t=1");

    // Non-divisible span: last step is shortened, endpoint still exact
    numcomp::RK4Solver partial_solver(0.3);
    partial_solver.solve(exponential_ode, 0.0, 1.0, 1.0);
    CHECK(partial_solver.t_history.size() == 5, "0->1 at h=0.3 must have 5 points");
    CHECK(partial_solver.t_history.back() == 1.0, "partial last step must land on t_end");

    // Zero-length interval
    numcomp::RK4Solver boundary_solver(0.1);
    boundary_solver.solve(exponential_ode, 0.0, 5.0, 0.0);
    CHECK(boundary_solver.y_history.size() == 1, "zero-length interval must have 1 point");
    CHECK(boundary_solver.y_history.front() == 5.0, "zero-length interval must return y0");

    // Tiny positive span still takes one step and lands on t_end
    numcomp::RK4Solver tiny_solver(0.1);
    tiny_solver.solve(exponential_ode, 0.0, 1.0, 1e-12);
    CHECK(tiny_solver.t_history.size() == 2, "tiny span must take exactly one step");
    CHECK(tiny_solver.t_history.back() == 1e-12, "tiny span must land on t_end");

    // Invalid step sizes must throw instead of hanging
    bool threw = false;
    try { numcomp::RK4Solver bad(0.0); } catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw, "h = 0 must throw std::invalid_argument");

    threw = false;
    try { numcomp::RK4Solver bad(-0.1); } catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw, "negative h must throw std::invalid_argument");

    // Reversed interval must throw
    threw = false;
    try { numcomp::RK4Solver s(0.1); s.solve(exponential_ode, 1.0, 1.0, 0.0); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw, "t_end < t0 must throw std::invalid_argument");

    if (failures == 0) {
        std::cout << "\u2705 All automated verification and proof tests passed successfully!\n\n";
    } else {
        std::cout << "\u274C " << failures << " test check(s) failed.\n\n";
    }
}

int main() {
    run_proof_tests();
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
