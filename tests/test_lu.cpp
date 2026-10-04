#include <numerical-computing/lu.h>

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
static bool near(double a, double b, double eps = 1e-9) {
    return std::fabs(a - b) <= eps;
}

// True if P*A equals L*U, i.e. row i of L*U is row perm[i] of A
static bool factorsMatch(const Matrix& A, const numcomp::LUResult& f, double eps = 1e-9) {
    const std::size_t n = A.size();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double lu = 0.0;
            for (std::size_t k = 0; k < n; ++k) lu += f.L[i][k] * f.U[k][j];
            if (!near(lu, A[static_cast<std::size_t>(f.perm[i])][j], eps)) return false;
        }
    }
    return true;
}

// True if L has 1s on the diagonal and zeros above it
static bool isUnitLower(const Matrix& L) {
    for (std::size_t i = 0; i < L.size(); ++i) {
        if (!near(L[i][i], 1.0)) return false;
        for (std::size_t j = i + 1; j < L.size(); ++j) {
            if (!near(L[i][j], 0.0)) return false;
        }
    }
    return true;
}

// True if U has zeros below the diagonal
static bool isUpper(const Matrix& U) {
    for (std::size_t i = 0; i < U.size(); ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            if (!near(U[i][j], 0.0)) return false;
        }
    }
    return true;
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

// True if running f() throws SingularMatrixError
template <typename F>
static bool throwsSingular(F f) {
    try {
        f();
    } catch (const numcomp::SingularMatrixError&) {
        return true;
    }
    return false;
}

int main() {
    // ---- Normal use: inputs where we know what to expect ----
    {
        // The example from the plan: the only multiplier is 0.5
        Matrix A = {{2, 1}, {1, 3}};
        auto f = numcomp::luDecompose(A);
        check(near(f.L[1][0], 0.5) && near(f.U[0][0], 2) && near(f.U[0][1], 1) &&
                  near(f.U[1][0], 0) && near(f.U[1][1], 2.5) && f.swaps == 0,
              "plan example gives multiplier 0.5 and U = [[2,1],[0,2.5]]");
        check(factorsMatch(A, f), "plan example: L*U gives back A");
        check(isUnitLower(f.L) && isUpper(f.U), "plan example: L and U have the right shape");

        auto x = numcomp::luSolve(f, {3, 5});
        check(near(x[0], 0.8) && near(x[1], 1.4), "plan example solves to x = 0.8, y = 1.4");
    }
    {
        // 3x3 that forces row swaps
        Matrix A = {{0, 2, 1}, {1, 1, 1}, {4, 1, 0}};
        auto f = numcomp::luDecompose(A);
        check(f.swaps == 2, "3x3 needing pivoting reports 2 swaps");
        check(factorsMatch(A, f), "3x3 needing pivoting: P*A equals L*U");
        check(isUnitLower(f.L) && isUpper(f.U), "3x3 needing pivoting: L and U have the right shape");
        Vector b = {3, 3, 5};
        check(residual(A, numcomp::luSolve(f, b), b) < 1e-9, "3x3 needing pivoting solves A*x = b");
    }
    {
        // One factorization, several right-hand sides
        Matrix A = {{4, 3, 2}, {2, 5, 1}, {1, 2, 6}};
        auto f = numcomp::luDecompose(A);
        bool ok = true;
        for (const Vector& b : std::vector<Vector>{{1, 0, 0}, {0, 1, 0}, {7, -2, 3}}) {
            ok = ok && residual(A, numcomp::luSolve(f, b), b) < 1e-9;
        }
        check(ok, "one factorization solves several right-hand sides");
    }
    {
        // Identity matrix: L and U are both the identity, no swaps
        Matrix I = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        auto f = numcomp::luDecompose(I);
        check(f.L == I && f.U == I && f.swaps == 0, "identity matrix is unchanged");
    }
    {
        // The input must not be modified
        Matrix A = {{2, 1}, {1, 3}};
        Matrix copy = A;
        numcomp::luDecompose(A);
        check(A == copy, "input matrix is not modified");
    }

    // ---- Edge cases: valid input that is awkward ----
    {
        Matrix A = {{5}};
        auto f = numcomp::luDecompose(A);
        check(near(f.L[0][0], 1) && near(f.U[0][0], 5) && f.swaps == 0, "1x1 matrix");
        check(near(numcomp::luSolve(f, {10})[0], 2), "1x1 matrix solves 5x = 10");
    }
    {
        // Zero in the top-left corner: a swap is unavoidable
        Matrix A = {{0, 1}, {1, 0}};
        auto f = numcomp::luDecompose(A);
        check(f.swaps == 1 && f.perm[0] == 1 && f.perm[1] == 0, "zero in the corner forces one swap");
        check(factorsMatch(A, f), "zero in the corner: P*A equals L*U");
    }
    {
        // Tiny pivot: swapping in the larger entry keeps the answer accurate
        Matrix A = {{1e-20, 1}, {1, 1}};
        auto f = numcomp::luDecompose(A);
        auto x = numcomp::luSolve(f, {1, 2});
        check(near(x[0], 1, 1e-9) && near(x[1], 1, 1e-9), "tiny pivot stays accurate thanks to pivoting");
    }

    // ---- Invalid input: each of these must throw ----
    check(throwsInvalidArgument([] { numcomp::luDecompose({}); }),
          "empty matrix throws");
    check(throwsInvalidArgument([] { numcomp::luDecompose({{1, 2, 3}, {4, 5, 6}}); }),
          "non-square matrix throws");
    check(throwsInvalidArgument([] { numcomp::luDecompose({{1, 2}, {3}}); }),
          "ragged rows throw");
    check(throwsSingular([] { numcomp::luDecompose({{1, 2}, {2, 4}}); }),
          "singular matrix throws SingularMatrixError");
    check(throwsInvalidArgument([] { numcomp::luDecompose({{1, 2}, {2, 4}}); }),
          "SingularMatrixError can also be caught as std::invalid_argument");
    check(throwsSingular([] { numcomp::luDecompose({{0, 0}, {0, 0}}); }),
          "all-zero matrix throws SingularMatrixError");
    check(throwsSingular([] { numcomp::luDecompose({{1, 0, 2}, {2, 0, 4}, {3, 0, 1}}); }),
          "matrix with a zero column throws SingularMatrixError");
    check(throwsInvalidArgument([] {
              auto f = numcomp::luDecompose({{2, 1}, {1, 3}});
              numcomp::luSolve(f, {1, 2, 3});
          }),
          "luSolve with wrong-sized b throws");

    std::cout << (failures == 0 ? "All tests passed\n" : "Some tests failed\n");
    return failures == 0 ? 0 : 1;  // exit code 0 = success, so ctest can tell
}
