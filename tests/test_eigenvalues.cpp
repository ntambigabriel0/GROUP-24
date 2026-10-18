#include <numerical-computing/eigenvalues.hpp>

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

// True if the computed eigenvalues match the expected ones (both largest first)
static bool sameValues(const Vector& got, const Vector& expected, double eps = 1e-8) {
    if (got.size() != expected.size()) return false;
    for (std::size_t i = 0; i < got.size(); ++i) {
        if (!near(got[i], expected[i], eps)) return false;
    }
    return true;
}

// A quick power iteration to find the largest eigenvalue.

static double powerIterationLargest(const Matrix& A, int rounds = 2000) {
    const std::size_t n = A.size();
    Vector v(n, 1.0);
    double lambda = 0.0;
    for (int r = 0; r < rounds; ++r) {
        double len = 0.0;
        for (double x : v) len += x * x;
        len = std::sqrt(len);
        for (double& x : v) x /= len;                1
        Vector w(n, 0.0);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) w[i] += A[j][i] * v[j];
        lambda = 0.0;
        for (std::size_t i = 0; i < n; ++i) lambda += v[i] * w[i];  
        v = w;
    }
    return lambda;
}

// True if running f() throws std::invalid_argument
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
    //  Normal use: inputs where we know what to expect 
    {
        // The example 
        Matrix A = {{2, 1}, {1, 2}};
        auto r = numcomp::eigenvalues(A);
        check(sameValues(r.values, {3, 1}) && r.converged, "plan example gives 3 and 1");
    }
    {
        // Symmetric 3x3 with known eigenvalues 
        Matrix A = {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}};
        auto r = numcomp::eigenvalues(A);
        check(sameValues(r.values, {2 + std::sqrt(2.0), 2, 2 - std::sqrt(2.0)}) && r.converged,
              "symmetric 3x3 gives 2+sqrt(2), 2, 2-sqrt(2)");
    }
    {
        // Non-symmetric but with real eigenvalues
        Matrix A = {{4, 2}, {1, 3}};
        auto r = numcomp::eigenvalues(A);
        check(sameValues(r.values, {5, 2}) && r.converged,
              "non-symmetric matrix with real eigenvalues gives 5 and 2");
    }
    {
        // Property check: the eigenvalues must add up to the trace of the matrix
        Matrix A = {{4, 1, 2, 0}, {1, 3, 0, 1}, {2, 0, 5, 1}, {0, 1, 1, 2}};  
        auto r = numcomp::eigenvalues(A);
        double sum = 0.0;
        for (double v : r.values) sum += v;
        check(r.converged && near(sum, 4 + 3 + 5 + 2), "4x4 symmetric: eigenvalues add up to the trace");
    }
    {
        // The plan asks to compare the largest value with Power Iteration
        Matrix A = {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}};
        auto r = numcomp::eigenvalues(A);
        check(near(r.values[0], powerIterationLargest(A), 1e-6),
              "largest eigenvalue agrees with power iteration");
    }
    {
        // The input must not be modified
        Matrix A = {{2, 1}, {1, 2}};
        Matrix copy = A;
        numcomp::eigenvalues(A);
        check(A == copy, "input matrix is not modified");
    }

    //  Edge cases: valid input that is awkward 
    {
        Matrix A = {{7}};
        auto r = numcomp::eigenvalues(A);
        check(sameValues(r.values, {7}) && r.converged && r.iterations == 0, "1x1 matrix");
    }
    {
        // Already diagonal: nothing to do, answer is read straight off
        Matrix A = {{5, 0, 0}, {0, 3, 0}, {0, 0, 1}};
        auto r = numcomp::eigenvalues(A);
        check(sameValues(r.values, {5, 3, 1}) && r.iterations == 0, "diagonal matrix needs no rounds");
    }
    {
        // Eigenvalues come back largest first even if the diagonal is in another order
        Matrix A = {{1, 0}, {0, 4}};
        check(sameValues(numcomp::eigenvalues(A).values, {4, 1}), "values are sorted largest first");
    }
    {
        // Singular matrix (eigenvalue 0): the QR step cannot factor it, so a shift is used
        Matrix A = {{1, 1}, {1, 1}};
        auto r = numcomp::eigenvalues(A);
        check(sameValues(r.values, {2, 0}) && r.converged, "singular matrix gives 2 and 0");
    }
    {
        Matrix A = {{0, 0}, {0, 0}};
        auto r = numcomp::eigenvalues(A);
        check(sameValues(r.values, {0, 0}) && r.converged, "all-zero matrix gives 0 and 0");
    }
    {
        // Complex eigenvalues (90-degree rotation): it must say it did not converge
        Matrix A = {{0, 1}, {-1, 0}};
        auto r = numcomp::eigenvalues(A, 1e-10, 50);
        check(!r.converged && r.iterations == 50, "complex eigenvalues report converged = false");
    }
    {
        // Too few rounds allowed: it must say so rather than pretend
        Matrix A = {{2, 1}, {1, 2}};
        auto r = numcomp::eigenvalues(A, 1e-12, 1);
        check(!r.converged && r.iterations == 1, "running out of rounds reports converged = false");
    }

    //  Invalid input: each of these must throw 
    check(throwsInvalidArgument([] { numcomp::eigenvalues({}); }),
          "empty matrix throws");
    check(throwsInvalidArgument([] { numcomp::eigenvalues({{1, 2, 3}, {4, 5, 6}}); }),
          "non-square matrix throws");
    check(throwsInvalidArgument([] { numcomp::eigenvalues({{1, 2}, {3}}); }),
          "ragged columns throw");
    check(throwsInvalidArgument([] { numcomp::eigenvalues({{2, 1}, {1, 2}}, 0.0); }),
          "zero tolerance throws");
    check(throwsInvalidArgument([] { numcomp::eigenvalues({{2, 1}, {1, 2}}, -1.0); }),
          "negative tolerance throws");
    check(throwsInvalidArgument([] { numcomp::eigenvalues({{2, 1}, {1, 2}}, 1e-10, -5); }),
          "negative maxIterations throws");

    std::cout << (failures == 0 ? "All tests passed\n" : "Some tests failed\n");
    return failures == 0 ? 0 : 1;  // exit code 0 = success
}
