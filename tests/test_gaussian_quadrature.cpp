#include <numerical-computing/gaussian_quadrature.hpp>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

int passed = 0;
int failed = 0;

void check(bool ok, const std::string& name) {
    if (ok) { ++passed; std::cout << "PASS: " << name << "\n"; }
    else    { ++failed; std::cout << "FAIL: " << name << "\n"; }
}

bool near(double a, double b, double tol = 1e-12) { return std::fabs(a - b) <= tol; }

template <typename F>
bool throwsInvalid(F call) {
    try { call(); } catch (const std::invalid_argument&) { return true; } catch (...) { return false; }
    return false;
}

}  // namespace

int main() {
    using namespace numcomp;

    // ---- Rules: known values ----
    {
        QuadratureRule r1 = gaussLegendreRule(1);
        check(near(r1.nodes[0], 0.0) && near(r1.weights[0], 2.0), "1-point rule: node 0, weight 2");

        QuadratureRule r2 = gaussLegendreRule(2);
        const double s = 1.0 / std::sqrt(3.0);
        check(near(r2.nodes[0], -s) && near(r2.nodes[1], s), "2-point rule: nodes -+1/sqrt(3)");
        check(near(r2.weights[0], 1.0) && near(r2.weights[1], 1.0), "2-point rule: weights 1, 1");

        QuadratureRule r3 = gaussLegendreRule(3);
        check(near(r3.nodes[0], -std::sqrt(0.6)) && near(r3.nodes[1], 0.0) && near(r3.nodes[2], std::sqrt(0.6)),
              "3-point rule: nodes -+sqrt(3/5) and 0");
        check(near(r3.weights[0], 5.0 / 9.0) && near(r3.weights[1], 8.0 / 9.0) && near(r3.weights[2], 5.0 / 9.0),
              "3-point rule: weights 5/9, 8/9, 5/9");
    }

    // ---- Rules: structure for many sizes ----
    {
        bool sumOk = true, symOk = true, orderOk = true, insideOk = true;
        for (int n = 1; n <= 100; ++n) {
            QuadratureRule r = gaussLegendreRule(n);
            double sum = 0.0;
            for (int i = 0; i < n; ++i) {
                sum += r.weights[i];
                if (r.weights[i] <= 0.0) sumOk = false;
                if (!near(r.nodes[i], -r.nodes[n - 1 - i], 1e-13)) symOk = false;
                if (!near(r.weights[i], r.weights[n - 1 - i], 1e-13)) symOk = false;
                if (i > 0 && !(r.nodes[i] > r.nodes[i - 1])) orderOk = false;
                if (std::fabs(r.nodes[i]) >= 1.0) insideOk = false;
            }
            if (!near(sum, 2.0, 1e-12)) sumOk = false;
        }
        check(sumOk, "n = 1..100: weights are positive and add up to 2");
        check(symOk, "n = 1..100: nodes and weights are symmetric");
        check(orderOk, "n = 1..100: nodes are in increasing order");
        check(insideOk, "n = 1..100: nodes lie strictly inside (-1, 1)");
    }

    // ---- Exactness: n points integrate degree 2n-1 exactly ----
    {
        bool allExact = true;
        for (int n = 1; n <= 12; ++n) {
            const int deg = 2 * n - 1;
            // integral of x^deg over [0, 2] is 2^(deg+1) / (deg+1)
            double exact = std::pow(2.0, deg + 1) / (deg + 1);
            double got = gaussLegendre([deg](double x) { return std::pow(x, deg); }, 0.0, 2.0, n);
            if (std::fabs(got - exact) > 1e-11 * std::fabs(exact)) allExact = false;
        }
        check(allExact, "n = 1..12: exact for x^(2n-1) on [0, 2]");

        // One degree higher should NOT be exact.
        double got = gaussLegendre([](double x) { return x * x * x * x; }, 0.0, 1.0, 2);  // degree 4 > 3
        check(!near(got, 0.2, 1e-6), "2 points: x^4 is not exact (degree above 2n-1)");
    }

    // ---- Standard integrals ----
    check(near(gaussLegendre([](double x) { return std::sin(x); }, 0.0, std::acos(-1.0), 10), 2.0, 1e-12),
          "integral of sin on [0, pi] = 2");
    check(near(gaussLegendre([](double x) { return std::exp(x); }, 0.0, 1.0, 8), std::exp(1.0) - 1.0, 1e-13),
          "integral of e^x on [0, 1] = e - 1");
    check(near(gaussLegendre([](double x) { return 1.0 / x; }, 1.0, 2.0, 10), std::log(2.0), 1e-10),
          "integral of 1/x on [1, 2] = ln 2");
    check(near(gaussLegendre([](double) { return 3.0; }, -1.0, 4.0, 1), 15.0), "constant function: 3 * length");
    check(near(gaussLegendre([](double x) { return x; }, -2.0, 2.0, 3), 0.0), "odd function over symmetric interval = 0");

    // ---- Interval handling ----
    {
        auto f = [](double x) { return x * x; };
        double fwd = gaussLegendre(f, 0.0, 3.0, 3);
        double rev = gaussLegendre(f, 3.0, 0.0, 3);
        check(near(fwd, 9.0) && near(rev, -9.0), "reversed limits change the sign");
        check(gaussLegendre(f, 2.0, 2.0, 3) == 0.0, "zero-width interval gives 0");
        check(near(gaussLegendre(f, -5.0, -2.0, 3), (-8.0 + 125.0) / 3.0), "interval with negative limits");
    }

    // ---- Accuracy improves with more points ----
    {
        auto f = [](double x) { return std::exp(-x * x); };
        const double exact = 0.746824132812427;  // integral of e^(-x^2) on [0, 1]
        double e2 = std::fabs(gaussLegendre(f, 0.0, 1.0, 2) - exact);
        double e4 = std::fabs(gaussLegendre(f, 0.0, 1.0, 4) - exact);
        double e6 = std::fabs(gaussLegendre(f, 0.0, 1.0, 6) - exact);
        check(e2 > e4 && e4 > e6 && e6 < 1e-9, "error shrinks as points go 2 -> 4 -> 6");
    }

    // ---- Composite rule ----
    {
        auto f = [](double x) { return std::sin(x); };
        const double pi = std::acos(-1.0);
        double one = std::fabs(compositeGaussLegendre(f, 0.0, 10.0, 1, 3) - (1.0 - std::cos(10.0)));
        double many = std::fabs(compositeGaussLegendre(f, 0.0, 10.0, 20, 3) - (1.0 - std::cos(10.0)));
        check(many < one && many < 1e-7, "composite rule with more panels is more accurate");
        check(near(compositeGaussLegendre(f, 0.0, pi, 1, 10), gaussLegendre(f, 0.0, pi, 10)),
              "composite with one panel equals the single rule");
        check(near(compositeGaussLegendre([](double x) { return x * x * x; }, 0.0, 2.0, 4, 2), 4.0),
              "composite is exact for polynomials the rule integrates exactly");
        check(compositeGaussLegendre(f, 1.0, 1.0, 5, 3) == 0.0, "composite: zero-width interval gives 0");
        check(near(compositeGaussLegendre([](double x) { return x; }, 2.0, 0.0, 3, 2), -2.0),
              "composite: reversed limits change the sign");
        // |x| has a corner at 0: with a panel edge at 0 each panel is a straight line, so it is exact.
        check(near(compositeGaussLegendre([](double x) { return std::fabs(x); }, -1.0, 1.0, 2, 2), 1.0),
              "composite handles a corner when it falls on a panel edge");
    }

    // ---- Invalid input ----
    {
        auto f = [](double x) { return x; };
        check(throwsInvalid([] { gaussLegendreRule(0); }), "0 points throws");
        check(throwsInvalid([] { gaussLegendreRule(-3); }), "negative points throws");
        check(throwsInvalid([] { gaussLegendreRule(101); }), "101 points throws");
        check(throwsInvalid([&] { gaussLegendre(f, 0.0, 1.0, 0); }), "gaussLegendre with 0 points throws");
        check(throwsInvalid([&] { gaussLegendre(std::function<double(double)>(), 0.0, 1.0); }), "empty function throws");
        check(throwsInvalid([&] { gaussLegendre(f, NAN, 1.0); }), "NaN limit throws");
        check(throwsInvalid([&] { gaussLegendre(f, 0.0, INFINITY); }), "infinite limit throws");
        check(throwsInvalid([&] { compositeGaussLegendre(f, 0.0, 1.0, 0); }), "composite with 0 panels throws");
        check(throwsInvalid([&] { compositeGaussLegendre(f, 0.0, 1.0, 4, 0); }), "composite with 0 points throws");
    }

    std::cout << "\nPassed: " << passed << "  Failed: " << failed << "\n";
    if (failed == 0) { std::cout << "All tests passed\n"; return 0; }
    return 1;
}
