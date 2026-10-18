#include <numerical-computing/lagrange.hpp>
#include <numerical-computing/polynomial_interpolation.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace numcomp;

static int failures = 0;

static void check(bool ok, const char* name) {
    std::cout << (ok ? "PASS: " : "FAIL: ") << name << "\n";
    if (!ok) ++failures;
}

static bool near(double a, double b, double eps = 1e-9) { return std::fabs(a - b) < eps; }

static bool vecNear(const Vector& a, const Vector& b, double eps = 1e-9) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!near(a[i], b[i], eps)) return false;
    return true;
}

template <class Ex, class F>
static bool throwsType(F f) {
    try { f(); } catch (const Ex&) { return true; } catch (...) { return false; }
    return false;
}

int main() {
    // Points (0,1), (1,3), (2,7) lie on x^2 + x + 1, so P(3) = 13 (the plan's example)
    {
        Vector c = polynomialInterpolation({0, 1, 2}, {1, 3, 7});
        check(vecNear(c, {1, 1, 1}), "plan example: coefficients are 1, 1, 1");
        check(near(evaluatePolynomial(c, 3), 13), "plan example: P(3) = 13");
    }
    // Two points give a straight line: (1,2), (3,6) -> y = 2x
    check(vecNear(polynomialInterpolation({1, 3}, {2, 6}), {0, 2}), "two points give a line");
    // One point gives a constant
    check(vecNear(polynomialInterpolation({5}, {7}), {7}), "one point gives a constant");
    // Cubic x^3 - 2x + 1 sampled at -1, 0, 1, 2
    check(vecNear(polynomialInterpolation({-1, 0, 1, 2}, {2, 1, 0, 5}), {1, -2, 0, 1}),
          "cubic is recovered exactly");
    // Unsorted x values work the same
    check(vecNear(polynomialInterpolation({2, 0, 1}, {7, 1, 3}), {1, 1, 1}), "unsorted x values");

    // The polynomial passes through every data point
    {
        Vector xs = {-2, -1, 0.5, 1, 3};
        Vector ys = {4, -1, 2, 0, 7};
        Vector c = polynomialInterpolation(xs, ys);
        bool all = true;
        for (std::size_t i = 0; i < xs.size(); ++i) all = all && near(evaluatePolynomial(c, xs[i]), ys[i]);
        check(all, "passes through all data points");
    }
    // Cross-check against Lagrange (the plan says they must agree)
    {
        Vector xs = {0, 1, 2, 4, 5};
        Vector ys = {1, 3, 2, 8, 6};
        Vector c = polynomialInterpolation(xs, ys);
        bool same = true;
        for (double q : {-1.0, 0.5, 1.5, 3.0, 6.0})
            same = same && near(evaluatePolynomial(c, q), lagrangeInterpolate(xs, ys, q), 1e-7);
        check(same, "agrees with Lagrange interpolation");
    }
    // evaluatePolynomial on its own
    check(near(evaluatePolynomial({1, 2, 3}, 2), 17), "evaluatePolynomial 1 + 2x + 3x^2 at 2");
    check(near(evaluatePolynomial({}, 5), 0), "evaluatePolynomial of empty list is 0");

    // Errors
    check(throwsType<std::invalid_argument>([] { polynomialInterpolation({}, {}); }),
          "empty input throws invalid_argument");
    check(throwsType<std::invalid_argument>([] { polynomialInterpolation({0, 1}, {1}); }),
          "size mismatch throws invalid_argument");
    check(throwsType<std::invalid_argument>([] { polynomialInterpolation({1, 1}, {2, 3}); }),
          "duplicate x values throw invalid_argument");
    check(throwsType<std::invalid_argument>(
              [] { polynomialInterpolation({0, std::numeric_limits<double>::quiet_NaN()}, {1, 2}); }),
          "NaN input throws invalid_argument");
    check(throwsType<std::invalid_argument>(
              [] { polynomialInterpolation({0, 1}, {1, std::numeric_limits<double>::infinity()}); }),
          "infinite input throws invalid_argument");

    if (failures == 0) {
        std::cout << "All tests passed\n";
        return 0;
    }
    std::cout << failures << " test(s) FAILED\n";
    return 1;
}
