#include <iostream>
#include <iomanip>
#include <cmath>
#include "rk2.hpp"

int main() {
    
    auto f = [](double t, double y) { return y; };

    double y0 = 1.0;
    double t0 = 0.0;
    double t_end = 1.0;
    double h = 0.1;

    double result = rungeKutta2(f, y0, t0, t_end, h);
    double exact = std::exp(t_end);

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=== RK2 Midpoint Example ===\n";
    std::cout << "Final t           : " << t_end << "\n";
    std::cout << "Calculated y(1.0) : " << result << "\n";
    std::cout << "Exact Value       : " << exact << "\n";
    std::cout << "Absolute Error    : " << std::abs(result - exact) << "\n";

    return 0;
}