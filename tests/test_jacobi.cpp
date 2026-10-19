#include <numerical-computing/jacobi.hpp>
#include <numerical-computing/lu.h>  // used only to cross-check the answers

#include <cmath>
#include <iostream>
#include <stdexcept>

using numcomp::Matrix;
using numcomp::Vector;

static int failures = 0;  // counts how many checks failed

// Prints PASS or FAIL for one check and counts failures
static void check(bool condition, const char* name) {
    if (condition) {
        std::cout << "PASS: " << name << "\n";
    } else {
        std::cout << "FAIL: " << name << "\n";
        ++failures;
    }
}

// True if two numbers are equal within a small tolerance (never compare decimals with ==)
static bool near(double a, double b, double eps = 1e-8) {
    return std::fabs(a - b) <= eps;
}

// Largest |A*x - b|: how far x is from solving the system
static double residual(const Matrix& A, const Vector& x, const Vector& b) {
    double worst = 0.0;
    for (std::size_t i = 0; i < A.size(); ++i) {
        double r = -b[i];
        for (std::size_t j = 0; j < A.size(); ++j) r += A[i][j] * x[j];
        worst = std::fmax(worst, std::fabs(r));
    }
    return worst;
}

// True if running f() throws std::invalid_argument
template <typename F>
static bool throwsInvalidArgument(F f) {
    try {
        f();
    } catch (const std::invalid_argument&) {
        return true;  // the error we wanted
    }
    return false;  // no error was thrown
}

int main() {
    //  Normal use: inputs where we know what to expect 
    const Matrix plan = {{4, 1}, {2, 5}};  // 4x + y = 6 and 2x + 5y = 12
    const Vector planB = {6, 12};
    {
        // The plan's rounds, from (0,0): (1.5, 2.4), (0.9, 1.8), (1.05, 2.04)
        auto r1 = numcomp::jacobi(plan, planB, {0, 0}, 1e-15, 1);
        auto r2 = numcomp::jacobi(plan, planB, {0, 0}, 1e-15, 2);
        auto r3 = numcomp::jacobi(plan, planB, {0, 0}, 1e-15, 3);
        check(near(r1.x[0], 1.5) && near(r1.x[1], 2.4), "plan example round 1 is (1.5, 2.4)");
        check(near(r2.x[0], 0.9) && near(r2.x[1], 1.8), "plan example round 2 is (0.9, 1.8)");
        check(near(r3.x[0], 1.05) && near(r3.x[1], 2.04), "plan example round 3 is (1.05, 2.04)");
    }
    {
        auto r = numcomp::jacobi(plan, planB);
        check(near(r.x[0], 1) && near(r.x[1], 2) && r.converged, "plan example settles on (1, 2)");
    }
    {
        // Each round must use only the PREVIOUS round's values, not this round's.
        auto r2 = numcomp::jacobi(plan, planB, {0, 0}, 1e-15, 2);
        check(near(r2.x[0], 0.9), "each round uses only the previous round's values");
    }
    {
        // 3x3 diagonally dominant: the answer must satisfy the equations and match LU
        Matrix A = {{10, -1, 2}, {-1, 11, -1}, {2, -1, 10}};
        Vector b = {6, 25, -11};
        auto r = numcomp::jacobi(A, b);
        check(r.converged && residual(A, r.x, b) < 1e-8, "3x3 solution satisfies A*x = b");

        Vector viaLU = numcomp::luSolve(numcomp::luDecompose(A), b);
        bool same = true;
        for (std::size_t i = 0; i < 3; ++i) same = same && near(r.x[i], viaLU[i], 1e-7);
        check(same, "3x3 answer matches the direct LU solution");
    }
    {
        // Starting from the exact answer, the very first round changes nothing
        auto r = numcomp::jacobi(plan, planB, {1, 2});
        check(r.converged && r.iterations == 1 && near(r.change, 0.0),
              "starting at the exact answer converges in one round");
    }
    {
        // A good starting guess needs fewer rounds than the default zeros
        auto fromZero = numcomp::jacobi(plan, planB);
        auto fromClose = numcomp::jacobi(plan, planB, {0.99, 1.99});
        check(fromClose.iterations < fromZero.iterations, "a closer starting guess needs fewer rounds");
    }
    {
        // The inputs must not be modified
        Matrix A = plan;
        Vector b = planB;
        numcomp::jacobi(A, b);
        check(A == plan && b == planB, "inputs are not modified");
    }

    //  Edge cases: valid input that is awkward 
    {
        Matrix A = {{4}};
        auto r = numcomp::jacobi(A, {8});
        check(near(r.x[0], 2) && r.converged, "1x1 system 4x = 8");
    }
    {
        // Diagonal matrix: each equation has one unknown, so it is solved at once
        Matrix A = {{2, 0, 0}, {0, 4, 0}, {0, 0, 5}};
        auto r = numcomp::jacobi(A, {2, 8, 10});
        check(near(r.x[0], 1) && near(r.x[1], 2) && near(r.x[2], 2) && r.converged,
              "diagonal system");
    }
    {
        // Not diagonally dominant: the values grow, so it must report failure
        Matrix A = {{1, 2}, {3, 1}};
        auto r = numcomp::jacobi(A, {3, 4}, {}, 1e-10, 50);
        check(!r.converged && r.iterations == 50, "non-dominant matrix reports converged = false");
    }
    {
        // Left running long enough the values overflow: it must stop, not loop on garbage
        Matrix A = {{1, 2}, {3, 1}};
        auto r = numcomp::jacobi(A, {3, 4}, {}, 1e-10, 100000);
        check(!r.converged && r.iterations < 100000, "overflowing values stop early");
    }
    {
        // Too few rounds allowed: it must say so rather than pretend
        auto r = numcomp::jacobi(plan, planB, {}, 1e-12, 3);
        check(!r.converged && r.iterations == 3, "running out of rounds reports converged = false");
    }
    {
        // Negative entries and negative answers
        Matrix A = {{5, -1}, {-1, 5}};
        Vector b = {-4, 4};
        auto r = numcomp::jacobi(A, b);
        check(r.converged && residual(A, r.x, b) < 1e-8, "negative entries and negative answers");
    }

    //  Invalid input: each of these must throw 
    check(throwsInvalidArgument([] { numcomp::jacobi({}, {}); }),
          "empty matrix throws");
    check(throwsInvalidArgument([] { numcomp::jacobi({{1, 2, 3}, {4, 5, 6}}, {1, 2}); }),
          "non-square matrix throws");
    check(throwsInvalidArgument([] { numcomp::jacobi({{1, 2}, {3}}, {1, 2}); }),
          "ragged rows throw");
    check(throwsInvalidArgument([] { numcomp::jacobi({{4, 1}, {2, 5}}, {1, 2, 3}); }),
          "wrong-sized b throws");
    check(throwsInvalidArgument([] { numcomp::jacobi({{4, 1}, {2, 5}}, {6, 12}, {1, 2, 3}); }),
          "wrong-sized starting guess throws");
    check(throwsInvalidArgument([] { numcomp::jacobi({{0, 1}, {1, 0}}, {1, 2}); }),
          "zero on the diagonal throws");
    check(throwsInvalidArgument([] { numcomp::jacobi({{4, 1}, {2, 5}}, {6, 12}, {}, 0.0); }),
          "zero tolerance throws");
    check(throwsInvalidArgument([] { numcomp::jacobi({{4, 1}, {2, 5}}, {6, 12}, {}, -1.0); }),
          "negative tolerance throws");
    check(throwsInvalidArgument([] { numcomp::jacobi({{4, 1}, {2, 5}}, {6, 12}, {}, 1e-10, 0); }),
          "maxIterations of 0 throws");

    std::cout << (failures == 0 ? "All tests passed\n" : "Some tests failed\n");
    return failures == 0 ? 0 : 1;  // exit code 0 = success
}
