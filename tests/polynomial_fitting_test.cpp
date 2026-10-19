#include <numerical-computing/polynomial_fitting.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

using numcomp::Vector;

static int failures = 0;

static void check(bool condition, const char* name) {
    if (condition) {
        std::cout << "PASS: " << name << "\n";
    } else {
        std::cout << "FAIL: " << name << "\n";
        ++failures;
    }
}

static bool near(double a, double b, double eps = 1e-6) {
    return std::fabs(a - b) <= eps;
}

template <typename F>
static bool throwsInvalidArgument(F f) {
    try {
        f();
    } catch (const std::invalid_argument&) {
        return true;
    }
    return false;
}

int main() {
    // Test 1: Project Plan Example (Straight Line) 
    {
        std::vector<double> xs = {0, 1, 2, 3};
        std::vector<double> ys = {1, 3, 5, 8};
        
        // Fit a degree 1 polynomial (a straight line)
        Vector coeffs = numcomp::polynomialFit(xs, ys, 1);
        
        check(coeffs.size() == 2, "Degree 1 fit returns 2 coefficients");
        check(near(coeffs[0], 0.8) && near(coeffs[1], 2.3), 
              "Best straight line matches project plan (y = 0.8 + 2.3x)");
    }

    // Test 2: Exact Quadratic Fit
    {
        // y = 1 + 2x + 3x^2
        std::vector<double> xs = {0, 1, 2, 3, 4};
        std::vector<double> ys = {1, 6, 17, 34, 57};
        
        Vector coeffs = numcomp::polynomialFit(xs, ys, 2);
        
        check(coeffs.size() == 3, "Degree 2 fit returns 3 coefficients");
        check(near(coeffs[0], 1.0) && near(coeffs[1], 2.0) && near(coeffs[2], 3.0), 
              "Exact quadratic data yields correct coefficients");
    }

    //  Test 3: Invalid Inputs 
    check(throwsInvalidArgument([] {
        numcomp::polynomialFit({}, {}, 1);
    }), "Empty data vectors throw std::invalid_argument");

    check(throwsInvalidArgument([] {
        numcomp::polynomialFit({1, 2}, {1, 2, 3}, 1);
    }), "Mismatched data vector sizes throw std::invalid_argument");

    check(throwsInvalidArgument([] {
        numcomp::polynomialFit({1, 2}, {3, 4}, -1);
    }), "Negative degree throws std::invalid_argument");

    check(throwsInvalidArgument([] {
        numcomp::polynomialFit({1, 2}, {3, 4}, 2);
    }), "Too few data points for requested degree throws std::invalid_argument");

    std::cout << (failures == 0 ? "All tests passed\n" : "Some tests failed\n");
    return failures == 0 ? 0 : 1;
}