#include <numerical-computing/lagrange.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using numcomp::lagrangeInterpolate;
using Vec = std::vector<double>;

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

// Decimals are compared with a small tolerance, never with ==
static bool near(double a, double b, double eps = 1e-9) {
    return std::fabs(a - b) < eps;
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
    // Normal use: answers we can work out by hand 
    {
       
        Vec xs = {0, 1, 2}, ys = {1, 3, 7};
        check(near(lagrangeInterpolate(xs, ys, 3), 13.0), "quadratic: extrapolates to x=3 -> 13");
        check(near(lagrangeInterpolate(xs, ys, 0.5), 1.75), "quadratic: x=0.5 -> 1.75");
        check(near(lagrangeInterpolate(xs, ys, -1), 1.0), "quadratic: x=-1 -> 1");
    }
    {
        
        Vec xs = {-2, 0, 1, 4, 5}, ys = {3, -1, 2, 8, -6};
        bool ok = true;
        for (std::size_t i = 0; i < xs.size(); ++i) {
            ok = ok && near(lagrangeInterpolate(xs, ys, xs[i]), ys[i]);
        }
        check(ok, "returns the original y at every given x");
    }
    {
   
        check(near(lagrangeInterpolate({0, 1}, {1, 3}, 10), 21.0), "two points give a line");
    }
    {
        
        check(near(lagrangeInterpolate({5}, {7}, -100), 7.0), "single point gives a constant");
    }
    {
       
        auto f = [](double x) { return x * x * x - 2 * x; };
        Vec xs = {-1, 0, 2, 3}, ys;
        for (double x : xs) ys.push_back(f(x));
        check(near(lagrangeInterpolate(xs, ys, 1.5), f(1.5)) &&
                  near(lagrangeInterpolate(xs, ys, -3.2), f(-3.2), 1e-8),
              "cubic is reproduced exactly");
    }
    {
        
        Vec a = lagrangeInterpolate({0, 1, 2}, {1, 3, 7}, Vec{3});
        Vec b = lagrangeInterpolate({2, 0, 1}, {7, 1, 3}, Vec{3});
        check(near(a[0], b[0]), "unsorted points give the same answer");
    }
    {
       
        Vec xs = {0, 0.5, 1.0, 1.5, 2.0}, ys;
        for (double x : xs) ys.push_back(std::sin(x));
        check(std::fabs(lagrangeInterpolate(xs, ys, 0.75) - std::sin(0.75)) < 1e-3,
              "approximates sin(0.75) closely");
    }
    {
        Vec out = lagrangeInterpolate({0, 1, 2}, {1, 3, 7}, Vec{0, 1, 2, 3});
        check(out.size() == 4 && near(out[0], 1) && near(out[1], 3) &&
                  near(out[2], 7) && near(out[3], 13),
              "vector version evaluates many x values");
    }
    {
        check(lagrangeInterpolate({0, 1}, {1, 3}, Vec{}).empty(), "empty query list -> empty result");
    }

    {
        check(near(lagrangeInterpolate({1e-15, 2e-15}, {1, 2}, 1.5e-15), 1.5), "tiny-scale x values accepted");
        check(near(lagrangeInterpolate({1e8, 2e8}, {1, 2}, 1.5e8), 1.5), "large-scale x values accepted");
    }

    //     Bad input: must throw std::invalid_argument
    check(throwsInvalidArgument([] { lagrangeInterpolate({}, {}, 1.0); }), "empty data throws");
    check(throwsInvalidArgument([] { lagrangeInterpolate({1, 2}, {1}, 1.0); }), "size mismatch throws");
    check(throwsInvalidArgument([] { lagrangeInterpolate({1, 1}, {2, 3}, 1.0); }), "duplicate x throws");
    check(throwsInvalidArgument([] {
              lagrangeInterpolate({0, std::numeric_limits<double>::quiet_NaN()}, {1, 2}, 1.0);
          }), "NaN input throws");
    check(throwsInvalidArgument([] {
              lagrangeInterpolate({0, 1}, {1, std::numeric_limits<double>::infinity()}, 1.0);
          }), "infinite input throws");
    check(throwsInvalidArgument([] { lagrangeInterpolate({1, 1}, {2, 3}, Vec{1.0}); }),
          "vector version also validates");

    std::cout << (failures == 0 ? "\nAll tests passed.\n" : "\nSome tests FAILED.\n");
    return failures == 0 ? 0 : 1;  // non-zero exit tells CMake the test failed
}
