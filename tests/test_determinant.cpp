#include <numerical-computing/determinant.hpp>
#include <numerical-computing/lu.h>

#include <cmath>
#include <iostream>
#include <stdexcept>

using numcomp::Matrix;

static int failures = 0;

// Prints PASS or FAIL and counts failed tests.
static void check(bool condition, const char* name) {
    if (condition) {
        std::cout << "PASS: " << name << "\n";
    } else {
        std::cout << "FAIL: " << name << "\n";
        ++failures;
    }
}

// Compares two floating-point numbers within a tolerance.
static bool near(double a, double b, double eps = 1e-9) {
    return std::fabs(a - b) <= eps;
}

int main() {
    // Test 1: A 2x2 matrix.
    Matrix A = {
        {2, 1},
        {1, 3}
    };
    check(near(numcomp::determinant(A), 5.0),
          "2x2 matrix determinant");

    // Test 2: A 3x3 matrix.
    Matrix B = {
        {1, 2, 3},
        {0, 1, 4},
        {5, 6, 0}
    };
    check(near(numcomp::determinant(B), 1.0),
          "3x3 matrix determinant");

    // Test 3: A matrix requiring a row swap.
    Matrix C = {
        {0, 1},
        {2, 3}
    };
    check(near(numcomp::determinant(C), -2.0),
          "determinant with row swap");

    // Test 4: A singular matrix.
    Matrix D = {
        {1, 2},
        {2, 4}
    };
    check(near(numcomp::determinant(D), 0.0),
          "singular matrix determinant");

    // Test 5: Identity matrix.
    Matrix E = {
        {1, 0},
        {0, 1}
    };
    check(near(numcomp::determinant(E), 1.0),
          "identity matrix determinant");

    // Test 6: A 1x1 matrix.
    Matrix F = {
        {7}
    };
    check(near(numcomp::determinant(F), 7.0),
          "1x1 matrix determinant");

    // Test 7: Reject a non-square matrix.
    bool exceptionCaught = false;

    try {
        Matrix G = {
            {1, 2, 3},
            {4, 5, 6}
        };
        numcomp::determinant(G);
    } catch (const std::invalid_argument&) {
        exceptionCaught = true;
    }

    check(exceptionCaught, "non-square matrix rejected");

    // Print the final test summary.
    std::cout << "\n";
    if (failures == 0) {
        std::cout << "All determinant tests passed!\n";
    } else {
        std::cout << failures << " test(s) failed.\n";
    }

    return failures == 0 ? 0 : 1;
}
