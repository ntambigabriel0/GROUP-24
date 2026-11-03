#include <cmath>
#include <iomanip>
#include <iostream>
#include <utility>
#include <vector>
#include "runge_kutta2.hpp"

int main() {
    auto f = [](double, double y) { return y; };

    const double y0 = 1.0, t0 = 0.0, t_end = 1.0, h = 0.1;

    try {
    const double result = rungeKutta2(f, y0, t0, t_end, h);
        const double exact = std::exp(t_end);

        std::cout << std::fixed << std::setprecision(6);
        std::cout << "=== RK2 Midpoint Example ===\n";
        std::cout << "Final t           : " << t_end << "\n";
        std::cout << "Calculated y(1.0) : " << result << "\n";
        std::cout << "Exact Value       : " << exact << "\n";
        std::cout << "Absolute Error    : " << std::abs(result - exact) << "\n\n";

        const double k1 = f(t0, y0);
        const double k2 = f(t0 + h / 2, y0 + h / 2 * k1);
        std::cout << "First step (h = " << h << ")\n";
        std::cout << "Start slope       : " << k1 << "\n";
        std::cout << "Midpoint slope    : " << k2 << "\n";
        std::cout << "New value         : " << rk2_step(f, t0, y0, h)
                  << "  (true " << std::exp(t0 + h) << ")\n\n";

        const std::vector<RK2Point> path = rk2_trajectory(f, y0, t0, t_end, h);

        std::cout << "     t         y        exact\n";
        for (const auto& p : path) {
            std::cout << std::setw(6) << p.first << "  " << std::setw(8) << p.second
                      << "  " << std::setw(8) << std::exp(p.first) << "\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
