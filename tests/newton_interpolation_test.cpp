#include <numerical-computing/newton_interpolation.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

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
    try { f(); } 
    catch (const std::invalid_argument&) { return true; }
    return false;
}

int main() {
    //  Test 1: Project Plan Example 
    {
        std::vector<double> xs = {0, 1, 2};
        std::vector<double> ys = {1, 3, 7};
        
        std::vector<double> coeffs = numcomp::newtonCoefficients(xs, ys);
        
        check(coeffs.size() == 3, "Calculates correct number of coefficients");
        check(near(coeffs[0], 1.0) && near(coeffs[1], 2.0) && near(coeffs[2], 1.0), 
              "Calculates correct divided differences (slopes 2 and 4, slope-of-slopes 1)");
              
        double result = numcomp::newtonEvaluate(xs, coeffs, 3.0);
        check(near(result, 13.0), "Evaluates correctly at x = 3 to get y = 13");
    }

    // Test 2: Invalid Inputs 
    check(throwsInvalidArgument([] {
        numcomp::newtonCoefficients({}, {});
    }), "Empty vectors throw std::invalid_argument");

    check(throwsInvalidArgument([] {
        numcomp::newtonCoefficients({1, 2}, {1, 2, 3});
    }), "Mismatched sizes throw std::invalid_argument");

    check(throwsInvalidArgument([] {
        numcomp::newtonCoefficients({1, 1, 2}, {3, 4, 5});
    }), "Duplicate x values throw std::invalid_argument");

    std::cout << (failures == 0 ? "All tests passed\n" : "Some tests failed\n");
    return failures == 0 ? 0 : 1;
}