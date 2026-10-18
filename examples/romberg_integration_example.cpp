#include <iostream>
#include <iomanip>
#include <cmath>
#include "romberg.hpp"

int main() {
    // Define the function: f(x) = sin(x)
    auto f = [](double x) { return std::sin(x); };
    
    double a = 0.0;
    double b = M_PI;
    double exactValue = 2.0; 

    int levels = 5;
    double result = rombergIntegration(f, a, b, levels);

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "=== Romberg Integration Example ===\n";
    std::cout << "Function          : f(x) = sin(x)\n";
    std::cout << "Interval          : [0, Pi]\n";
    std::cout << "Levels            : " << levels << "\n";
    std::cout << "Computed Integral : " << result << "\n";
    std::cout << "Exact Value       : " << exactValue << "\n";
    std::cout << "Absolute Error    : " << std::abs(result - exactValue) << "\n";

    return 0;
}