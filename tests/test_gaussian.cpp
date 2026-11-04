#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include <numerical-computing/gaussian.hpp>

using numcomp::Matrix;
using numcomp::Vector;
using numcomp::gaussianElimination;

static int failures = 0;

static void check(bool ok, const std::string& name) {
    std::cout << (ok ? "  PASS: " : "  FAIL: ") << name << std::endl;
    if (!ok) ++failures;
}

static bool approximatelyEqual(double a, double b) {
    return std::fabs(a - b) < 1e-6;
}

// Solve A x = b and compare with the expected solution
static void testSolution(const std::string& name, Matrix A,
                         Vector b, Vector expected) {
    try {
        Vector result = gaussianElimination(A, b);
        bool ok = result.size() == expected.size();
        for (std::size_t i = 0; ok && i < expected.size(); i++) {
            if (!approximatelyEqual(result[i], expected[i])) ok = false;
        }
        check(ok, name);
    } catch (const std::exception& e) {
        check(false, name + " (threw: " + e.what() + ")");
    }
}

// Expect gaussianElimination to throw exception type E
template <typename E>
static void testThrows(const std::string& name, Matrix A, Vector b) {
    bool caught = false;
    try {
        gaussianElimination(A, b);
    } catch (const E&) {
        caught = true;
    } catch (...) {
    }
    check(caught, name);
}

int main() {
    std::cout << "--- GAUSSIAN ELIMINATION TESTS ---\n" << std::endl;

    // Basic 2x2 system
    testSolution("Basic 2x2 system",
                 {{2, 1},
                  {1, 3}},
                 {3, 5}, {0.8, 1.4});

    // Standard 3x3 system
    testSolution("Standard 3x3 system",
                 {{3, 2, -1},
                  {2, -2, 4},
                  {-1, 0.5, -1}},
                 {1, -2, 0}, {1, -2, -2});

    // System requiring partial pivoting (zero in the first pivot position)
    testSolution("Partial pivoting",
                 {{0, 2},
                  {1, 3}},
                 {4, 7}, {1, 2});

    // 1x1 system
    testSolution("1x1 system", {{4}}, {8}, {2});

    // Tiny leading pivot is handled by pivoting
    testSolution("Tiny leading pivot",
                 {{1e-14, 1},
                  {1, 1}},
                 {1, 2}, {1, 1});

    // Singular 2x2 (second row is twice the first)
    testThrows<std::runtime_error>("Singular 2x2 throws runtime_error",
                                   {{1, 2},
                                    {2, 4}},
                                   {3, 6});

    // Singular 3x3 (zero only appears at the last pivot)
    testThrows<std::runtime_error>("Singular 3x3 throws runtime_error",
                                   {{1, 2, 3},
                                    {2, 3, 4},
                                    {3, 4, 5}},
                                   {1, 2, 3});

    // 1x1 zero matrix
    testThrows<std::runtime_error>("1x1 zero matrix throws runtime_error",
                                   {{0}}, {1});

    // Invalid input
    testThrows<std::invalid_argument>("Empty matrix throws invalid_argument",
                                      Matrix{}, Vector{});
    testThrows<std::invalid_argument>("Non-square matrix throws invalid_argument",
                                      {{1, 2, 3},
                                       {4, 5, 6}},
                                      {1, 2});
    testThrows<std::invalid_argument>("Ragged matrix throws invalid_argument",
                                      {{1, 2},
                                       {3}},
                                      {1, 2});
    testThrows<std::invalid_argument>("Wrong-sized b throws invalid_argument",
                                      {{1, 0},
                                       {0, 1}},
                                      {1});

    std::cout << "\nAll tests finished running. Failures: " << failures << std::endl;
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
