#include <iostream>
#include <cmath>
#include <cassert>
#include "romberg.hpp"

void testPolynomial() {
    // Integral of x^3 from 0 to 2 = [x^4 / 4] = 4.0
    auto f = [](double x) { return x * x * x; };
    double result = rombergIntegration(f, 0.0, 2.0, 4);
    assert(std::abs(result - 4.0) < 1e-6);
    std::cout << "[PASS] Test Polynomial: int(x^3, 0..2) = " << result << "\n";
}

void testTrigonometric() {
    // Integral of cos(x) from 0 to pi/2 = 1.0
    auto f = [](double x) { return std::cos(x); };
    double result = rombergIntegration(f, 0.0, M_PI_2, 4);
    assert(std::abs(result - 1.0) < 1e-6);
    std::cout << "[PASS] Test Trigonometric: int(cos(x), 0..pi/2) = " << result << "\n";
}

void testExponential() {
    // Integral of e^x from 0 to 1 = e - 1 ≈ 1.718281828
    auto f = [](double x) { return std::exp(x); };
    double result = rombergIntegration(f, 0.0, 1.0, 5);
    double exact = std::exp(1.0) - 1.0;
    assert(std::abs(result - exact) < 1e-7);
    std::cout << "[PASS] Test Exponential: int(exp(x), 0..1) = " << result << "\n";
}

int main() {
    std::cout << "Running Romberg Integration Tests...\n";
    testPolynomial();
    testTrigonometric();
    testExponential();
    std::cout << "All tests passed successfully!\n";
    return 0;
}