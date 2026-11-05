#include <numerical-computing/lagrange.hpp>
#include "test_util.hpp"

#include <limits>

using numcomp::lagrangeInterpolate;
using Vec = std::vector<double>;
using namespace testutil;

int main() {
    // Answers we can work out by hand (points 0,1,2 with values 1,3,7)
    Vec xs = {0, 1, 2}, ys = {1, 3, 7};
    check(near(lagrangeInterpolate(xs, ys, 3), 13.0) && near(lagrangeInterpolate(xs, ys, 0.5), 1.75) &&
              near(lagrangeInterpolate(xs, ys, -1), 1.0),
          "quadratic: x=3 -> 13, x=0.5 -> 1.75, x=-1 -> 1");
    check(near(lagrangeInterpolate({0, 1}, {1, 3}, 10), 21.0), "two points give a line");
    check(near(lagrangeInterpolate({5}, {7}, -100), 7.0), "single point gives a constant");

    // Passes through every given point, in any order
    Vec px = {-2, 0, 1, 4, 5}, py = {3, -1, 2, 8, -6};
    bool hits = true;
    for (std::size_t i = 0; i < px.size(); ++i) hits = hits && near(lagrangeInterpolate(px, py, px[i]), py[i]);
    check(hits, "returns the original y at every given x");
    check(near(lagrangeInterpolate({2, 0, 1}, {7, 1, 3}, 3.0), 13.0), "unsorted points give the same answer");

    // A cubic through 4 points is reproduced exactly
    auto f = [](double x) { return x * x * x - 2 * x; };
    Vec cx = {-1, 0, 2, 3}, cy;
    for (double x : cx) cy.push_back(f(x));
    check(near(lagrangeInterpolate(cx, cy, 1.5), f(1.5)) && near(lagrangeInterpolate(cx, cy, -3.2), f(-3.2), 1e-8),
          "cubic is reproduced exactly");

    // Vector version, and very small / very large x values
    Vec out = lagrangeInterpolate(xs, ys, Vec{0, 1, 2, 3});
    check(out.size() == 4 && near(out[0], 1) && near(out[1], 3) && near(out[2], 7) && near(out[3], 13),
          "vector version evaluates many x values");
    check(near(lagrangeInterpolate({1e-15, 2e-15}, {1, 2}, 1.5e-15), 1.5) &&
              near(lagrangeInterpolate({1e8, 2e8}, {1, 2}, 1.5e8), 1.5),
          "tiny and large x values accepted");

    // Bad input must throw
    const double nan = std::numeric_limits<double>::quiet_NaN(), inf = std::numeric_limits<double>::infinity();
    check(throwsInvalid([] { lagrangeInterpolate({}, {}, 1.0); }), "empty data throws");
    check(throwsInvalid([] { lagrangeInterpolate({1, 2}, {1}, 1.0); }), "size mismatch throws");
    check(throwsInvalid([] { lagrangeInterpolate({1, 1}, {2, 3}, 1.0); }), "duplicate x throws");
    check(throwsInvalid([=] { lagrangeInterpolate({0, nan}, {1, 2}, 1.0); }), "NaN input throws");
    check(throwsInvalid([=] { lagrangeInterpolate({0, 1}, {1, inf}, 1.0); }), "infinite input throws");
    check(throwsInvalid([] { lagrangeInterpolate({1, 1}, {2, 3}, Vec{1.0}); }), "vector version also validates");

    return finish();
}
