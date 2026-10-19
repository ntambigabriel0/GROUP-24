#include <iostream>
#include <cmath>
#include <cassert>
#include "rk2.hpp"

void testExponentialGrowth() {
    
    auto f = [](double t, double y) { return y; };
    double result = rungeKutta2(f, 1.0, 0.0, 1.0, 0.01);
    double exact = std::exp(1.0);
    
    assert(std::abs(result - exact) < 1e-3);
    std::cout << "[PASS] Test Exponential Growth: y(1) = " << result << "\n";
}

void testLinearDecay() {
    
    auto f = [](double t, double y) { return -2.0 * y; };
    double result = rungeKutta2(f, 1.0, 0.0, 1.0, 0.01);
    double exact = std::exp(-2.0);
    
    assert(std::abs(result - exact) < 1e-3);
    std::cout << "[PASS] Test Linear Decay: y(1) = " << result << "\n";
}

int main() {
    std::cout << "Running RK2 Midpoint Tests...\n";
    testExponentialGrowth();
    testLinearDecay();
    std::cout << "All tests passed successfully!\n";
    return 0;
}