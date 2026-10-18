#include <numerical-computing/gram_schmidt.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>

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

// True if every vector has length 1 and every pair is at 90 degrees.
// (Decimals are compared with a small tolerance, never with ==)
static bool isOrthonormal(const std::vector<Vector>& q, double eps = 1e-9) {
    for (std::size_t i = 0; i < q.size(); ++i) {
        for (std::size_t j = 0; j < q.size(); ++j) {
            const double expected = (i == j) ? 1.0 : 0.0;  // 1 with itself, 0 with others
            if (std::fabs(numcomp::dot(q[i], q[j]) - expected) > eps) {
                return false;
            }
        }
    }
    return true;
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
    // ---- Normal use: inputs where we know what to expect ----
    {
        std::vector<Vector> in = {{1, 1, 0}, {1, 0, 1}, {0, 1, 1}};
        check(isOrthonormal(numcomp::gramSchmidt(in)), "3 vectors in R^3 become orthonormal");
    }
    {
        std::vector<Vector> in = {{3, 4}};
        auto q = numcomp::gramSchmidt(in);
        check(std::fabs(q[0][0] - 0.6) < 1e-12 && std::fabs(q[0][1] - 0.8) < 1e-12,
              "single vector is just normalized");
    }
    {
        std::vector<Vector> in = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
        auto q = numcomp::gramSchmidt(in);
        check(q == in, "already orthonormal input is unchanged");
    }
    {
        std::vector<Vector> in = {{1, 2, 3, 4}, {2, 1, 0, 1}};
        check(isOrthonormal(numcomp::gramSchmidt(in)), "fewer vectors than dimensions");
    }

    // ---- Edge case: valid input that is awkward ----
    {
        // nearly parallel vectors: modified Gram-Schmidt should stay accurate
        std::vector<Vector> in = {{1, 1e-8, 0}, {1, 0, 0}, {1, 0, 1e-8}};
        check(isOrthonormal(numcomp::gramSchmidt(in, 1e-12), 1e-6),
              "nearly dependent vectors stay orthogonal");
    }

    // ---- Invalid input: each of these must throw ----
    check(throwsInvalidArgument([] { numcomp::gramSchmidt({}); }),
          "empty list throws");
    check(throwsInvalidArgument([] { numcomp::gramSchmidt({Vector{}}); }),
          "empty vector throws");
    check(throwsInvalidArgument([] { numcomp::gramSchmidt({{1, 2}, {1, 2, 3}}); }),
          "mismatched sizes throw");
    check(throwsInvalidArgument([] { numcomp::gramSchmidt({{1, 0}, {0, 1}, {1, 1}}); }),
          "more vectors than dimensions throws");
    check(throwsInvalidArgument([] { numcomp::gramSchmidt({{1, 2, 3}, {2, 4, 6}}); }),
          "dependent vectors throw");
    check(throwsInvalidArgument([] { numcomp::gramSchmidt({{0, 0, 0}}); }),
          "zero vector throws");

    std::cout << (failures == 0 ? "All tests passed\n" : "Some tests failed\n");
    return failures == 0 ? 0 : 1;  // exit code 0 = success, so ctest can tell
}
