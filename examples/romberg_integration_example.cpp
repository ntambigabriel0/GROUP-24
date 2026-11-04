#include <iostream>
#include <iomanip>
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#include <numerical-computing/romberg_integration.hpp>

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
