#include <numerical-computing/least_squares.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace numcomp;

namespace {

int passed = 0;
int failed = 0;

bool closeTo(double actual, double expected,
             double tol = 1e-6) {
    return std::abs(actual - expected) <= tol;
}

void check(bool condition, const std::string& name) {
    if (condition) {
        ++passed;
        std::cout << "PASS: " << name << '\n';
    } else {
        ++failed;
        std::cout << "FAIL: " << name << '\n';
    }
}

template <typename Function>
void checkThrows(Function function, const std::string& name) {
    try {
        function();
        check(false, name);
    } catch (const std::invalid_argument&) {
        check(true, name);
    } catch (...) {
        check(false, name);
    }
}

} // namespace

int main() {
    // Test 1: Fit the three points from the example.
    {
        std::vector<Vector> A = {
            {1.0, 1.0, 1.0},
            {0.0, 1.0, 2.0}
        };
        Vector b = {1.0, 2.0, 2.0};

        Vector x = leastSquares(A, b);

        check(
            x.size() == 2 &&
            closeTo(x[0], 7.0 / 6.0) &&
            closeTo(x[1], 0.5),
            "Best-fit line"
        );
    }

    // Test 2: Exact solution for a square system.
    {
        std::vector<Vector> A = {
            {1.0, 0.0},
            {0.0, 1.0}
        };
        Vector b = {3.0, 4.0};

        Vector x = leastSquares(A, b);

        check(
            closeTo(x[0], 3.0) &&
            closeTo(x[1], 4.0),
            "Identity matrix"
        );
    }

    // Test 3: Reject an empty matrix.
    checkThrows(
        [] {
            leastSquares({}, {1.0, 2.0});
        },
        "Empty matrix"
    );

    // Test 4: Reject a mismatched column length.
    checkThrows(
        [] {
            std::vector<Vector> A = {
                {1.0, 2.0},
                {3.0}
            };
            leastSquares(A, {1.0, 2.0});
        },
        "Inconsistent column sizes"
    );

    // Test 5: Reject fewer equations than unknowns.
    checkThrows(
        [] {
            std::vector<Vector> A = {
                {1.0},
                {2.0}
            };
            leastSquares(A, {3.0});
        },
        "Underdetermined system"
    );

    // Test 6: Reject a negative tolerance.
    checkThrows(
        [] {
            std::vector<Vector> A = {
                {1.0, 1.0},
                {0.0, 1.0}
            };
            leastSquares(A, {1.0, 2.0}, -1.0);
        },
        "Negative tolerance"
    );

    std::cout << "\nPassed: " << passed
              << "\nFailed: " << failed << '\n';

    return failed == 0 ? 0 : 1;
}
