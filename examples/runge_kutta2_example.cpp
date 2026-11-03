#include <cmath>
#include <iomanip>
#include <iostream>

#include <numerical-computing/runge_kutta2.hpp>

int main() {
    auto f = [](double /*t*/, double y) { return y; };  // y' = y, y(0) = 1

    const double y0 = 1.0;
    const double t0 = 0.0;
    const double t_end = 1.0;
    const double h = 0.1;

    try {
        double result = rungeKutta2(f, y0, t0, t_end, h);
        double exact = std::exp(t_end);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "=== RK2 Midpoint Example ===\n";
        std::cout << "Final t           : " << t_end << "\n";
        std::cout << "Calculated y(1.0) : " << result << "\n";
        std::cout << "Exact Value       : " << exact << "\n";
        std::cout << "Absolute Error    : " << std::fabs(result - exact) << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
