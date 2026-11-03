#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "numerical-computing/gaussian.hpp"

// Check whether two values are approximately equal
bool approximatelyEqual(double a, double b) {
    return std::fabs(a - b) < 1e-6;
}

// Test whether the algorithm produces the expected solution
void testSolution(
    const std::string& testName,
    Matrix A,
    std::vector<double> b,
    std::vector<double> expected
) {
    try {
        std::vector<double> result = gaussianElimination(A, b);

        bool passed = result.size() == expected.size();

        if (passed) {
            for (std::size_t i = 0; i < expected.size(); i++) {
                if (!approximatelyEqual(result[i], expected[i])) {
                    passed = false;
                    break;
                }
            }
        }

        if (passed) {
            std::cout << testName << ": PASSED" << std::endl;
        } else {
            std::cout << testName << ": FAILED" << std::endl;
            std::cout << "Expected: ";

            for (double value : expected) {
                std::cout << value << " ";
            }

            std::cout << "\nActual: ";

            for (double value : result) {
                std::cout << value << " ";
            }

            std::cout << std::endl;
        }
    }
    catch (const std::exception& e) {
        std::cout << testName << ": FAILED - "
                  << e.what() << std::endl;
    }
}

// Test whether the algorithm correctly detects a singular matrix
void testSingularMatrix() {
    Matrix A(2, 2);

    A.at(0, 0) = 1;
    A.at(0, 1) = 2;
    A.at(1, 0) = 2;
    A.at(1, 1) = 4;

    std::vector<double> b = {3, 6};

    try {
        gaussianElimination(A, b);
        std::cout << "Singular matrix test: FAILED" << std::endl;
    }
    catch (const std::runtime_error&) {
        std::cout << "Singular matrix test: PASSED" << std::endl;
    }
    catch (...) {
        std::cout << "Singular matrix test: FAILED - "
                  << "Unexpected exception." << std::endl;
    }
}

// Test whether an empty matrix is handled correctly
void testEmptyMatrix() {
    Matrix A(0, 0);
    std::vector<double> b;

    try {
        gaussianElimination(A, b);
        std::cout << "Empty matrix test: FAILED" << std::endl;
    }
    catch (const std::invalid_argument&) {
        std::cout << "Empty matrix test: PASSED" << std::endl;
    }
    catch (...) {
        std::cout << "Empty matrix test: FAILED - "
                  << "Unexpected exception." << std::endl;
    }
}

int main() {

    std::cout << "--- GAUSSIAN ELIMINATION TESTS ---"
              << std::endl << std::endl;

    // TEST 1: Basic 2x2 system
    Matrix A1(2, 2);

    A1.at(0, 0) = 2;
    A1.at(0, 1) = 1;
    A1.at(1, 0) = 1;
    A1.at(1, 1) = 3;

    std::vector<double> b1 = {3, 5};

    testSolution(
        "Test 1: Basic 2x2 system",
        A1, b1, {0.8, 1.4}
    );

    // TEST 2: Standard 3x3 system
    Matrix A2(3, 3);

    A2.at(0, 0) = 3;
    A2.at(0, 1) = 2;
    A2.at(0, 2) = -1;

    A2.at(1, 0) = 2;
    A2.at(1, 1) = -2;
    A2.at(1, 2) = 4;

    A2.at(2, 0) = -1;
    A2.at(2, 1) = 0.5;
    A2.at(2, 2) = -1;

    std::vector<double> b2 = {1, -2, 0};

    testSolution(
        "Test 2: Standard 3x3 system",
        A2, b2, {1, -2, -2}
    );

    // TEST 3: System requiring partial pivoting
    Matrix A3(2, 2);

    A3.at(0, 0) = 0;
    A3.at(0, 1) = 2;
    A3.at(1, 0) = 1;
    A3.at(1, 1) = 3;

    std::vector<double> b3 = {4, 7};

    testSolution(
        "Test 3: Partial pivoting",
        A3, b3, {1, 2}
    );

    // TEST 4: Singular matrix
    testSingularMatrix();

    // TEST 5: Empty matrix
    testEmptyMatrix();

    std::cout << std::endl
              << "All tests completed." << std::endl;

    return 0;
}
