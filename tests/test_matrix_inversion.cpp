#include <numerical-computing/matrix_inversion.hpp>
#include <numerical-computing/lu.h>

#include <cmath>
#include <iostream>
#include <stdexcept>

using numcomp::Vector;
using numcomp::Matrix;

static int failures = 0;

static void check(bool condition, const char* name) {
    if (condition) {
        std::cout << "PASS: " << name << "\n";
    } else {
        std::cout << "FAIL: " << name << "\n";
        ++failures;
    }
}

// Helper to multiply matrices: returns A * B
static Matrix multiply(const Matrix& A, const Matrix& B) {
    const std::size_t n = A.size();
    Matrix C(n, Vector(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

// Helper to check if a matrix is the Identity matrix (within tolerance)
static bool isNearIdentity(const Matrix& M, double eps = 1e-9) {
    const std::size_t n = M.size();
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            double expected = (i == j) ? 1.0 : 0.0;
            if (std::fabs(M[i][j] - expected) > eps) return false;
        }
    }
    return true;
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
    // ---- Normal Use ----
    {
        Matrix A = {{2, 1}, {1, 3}};
        Matrix invA = numcomp::invertMatrix(A);
        check(isNearIdentity(multiply(A, invA)), "2x2 matrix A * A^-1 is Identity");
    }
    {
        Matrix A = {{4, 3, 2}, {2, 5, 1}, {1, 2, 6}};
        Matrix invA = numcomp::invertMatrix(A);
        check(isNearIdentity(multiply(A, invA)), "3x3 matrix A * A^-1 is Identity");
    }
    {
        Matrix A = {{0, 2, 1}, {1, 1, 1}, {4, 1, 0}};
        Matrix invA = numcomp::invertMatrix(A);
        check(isNearIdentity(multiply(A, invA)), "Matrix requiring pivoting inverts correctly");
    }

    // ---- Invalid Input ----
    check(throwsInvalidArgument([] { numcomp::invertMatrix({}); }), "Empty matrix throws");
    check(throwsInvalidArgument([] { numcomp::invertMatrix({Vector{}}); }), "Empty row throws");
    check(throwsInvalidArgument([] { numcomp::invertMatrix({{1, 2}, {1}}); }), "Mismatched sizes throw");
    check(throwsInvalidArgument([] { numcomp::invertMatrix({{1, 2, 3}, {4, 5, 6}}); }), "Non-square throws");
    check(throwsSingular([] { numcomp::invertMatrix({{1, 2}, {2, 4}}); }), "Singular matrix throws SingularMatrixError");

    std::cout << (failures == 0 ? "All tests passed\n" : "Some tests failed\n");
    return failures == 0 ? 0 : 1;
}