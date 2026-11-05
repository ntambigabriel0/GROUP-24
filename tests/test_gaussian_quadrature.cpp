#include <numerical-computing/gaussian_quadrature.hpp>
#include "test_util.hpp"

#include <cmath>

using namespace numcomp;
using namespace testutil;

int main() {
    const double pi = std::acos(-1.0);

    // Rules: known values
    QuadratureRule r2 = gaussLegendreRule(2), r3 = gaussLegendreRule(3);
    const double s = 1.0 / std::sqrt(3.0);
    check(near(r2.nodes[0], -s, 1e-12) && near(r2.nodes[1], s, 1e-12) && near(r2.weights[0], 1.0, 1e-12),
          "2-point rule: nodes -+1/sqrt(3), weights 1, 1");
    check(near(r3.nodes[0], -std::sqrt(0.6), 1e-12) && near(r3.nodes[1], 0.0, 1e-12) &&
              near(r3.weights[0], 5.0 / 9.0, 1e-12) && near(r3.weights[1], 8.0 / 9.0, 1e-12),
          "3-point rule: nodes -+sqrt(3/5), 0 and weights 5/9, 8/9");

    // Rules: structure for every size 1..100
    bool sumOk = true, shapeOk = true;
    for (int n = 1; n <= 100; ++n) {
        QuadratureRule r = gaussLegendreRule(n);
        double sum = 0.0;
        for (int i = 0; i < n; ++i) {
            sum += r.weights[i];
            sumOk = sumOk && r.weights[i] > 0.0;
            shapeOk = shapeOk && near(r.nodes[i], -r.nodes[n - 1 - i], 1e-13) &&   // symmetric
                      near(r.weights[i], r.weights[n - 1 - i], 1e-13) &&
                      std::fabs(r.nodes[i]) < 1.0 &&                                // inside (-1, 1)
                      (i == 0 || r.nodes[i] > r.nodes[i - 1]);                      // increasing
        }
        sumOk = sumOk && near(sum, 2.0, 1e-12);
    }
    check(sumOk, "n = 1..100: weights are positive and add up to 2");
    check(shapeOk, "n = 1..100: nodes symmetric, increasing and inside (-1, 1)");

    // n points integrate degree 2n-1 exactly, but not degree 2n
    bool exact = true;
    for (int n = 1; n <= 12; ++n) {
        const int deg = 2 * n - 1;
        double want = std::pow(2.0, deg + 1) / (deg + 1);  // integral of x^deg on [0, 2]
        exact = exact && near(gaussLegendre([deg](double x) { return std::pow(x, deg); }, 0.0, 2.0, n), want, 1e-11 * want);
    }
    check(exact, "n = 1..12: exact for x^(2n-1) on [0, 2]");
    check(!near(gaussLegendre([](double x) { return x * x * x * x; }, 0.0, 1.0, 2), 0.2, 1e-6),
          "2 points: x^4 is not exact (degree above 2n-1)");

    // Standard integrals and interval handling
    check(near(gaussLegendre([](double x) { return std::sin(x); }, 0.0, pi, 10), 2.0, 1e-12), "integral of sin on [0, pi] = 2");
    check(near(gaussLegendre([](double x) { return std::exp(x); }, 0.0, 1.0, 8), std::exp(1.0) - 1.0, 1e-13), "integral of e^x on [0, 1] = e - 1");
    auto sq = [](double x) { return x * x; };
    check(near(gaussLegendre(sq, 0.0, 3.0, 3), 9.0, 1e-12) && near(gaussLegendre(sq, 3.0, 0.0, 3), -9.0, 1e-12),
          "reversed limits change the sign");
    check(gaussLegendre(sq, 2.0, 2.0, 3) == 0.0, "zero-width interval gives 0");

    // Accuracy improves with more points
    auto g = [](double x) { return std::exp(-x * x); };
    const double ref = 0.746824132812427;  // integral of e^(-x^2) on [0, 1]
    double e2 = std::fabs(gaussLegendre(g, 0.0, 1.0, 2) - ref), e4 = std::fabs(gaussLegendre(g, 0.0, 1.0, 4) - ref),
           e6 = std::fabs(gaussLegendre(g, 0.0, 1.0, 6) - ref);
    check(e2 > e4 && e4 > e6 && e6 < 1e-9, "error shrinks as points go 2 -> 4 -> 6");

    // Composite rule
    auto sn = [](double x) { return std::sin(x); };
    double truth = 1.0 - std::cos(10.0);
    check(std::fabs(compositeGaussLegendre(sn, 0.0, 10.0, 20, 3) - truth) < std::fabs(compositeGaussLegendre(sn, 0.0, 10.0, 1, 3) - truth) &&
              near(compositeGaussLegendre(sn, 0.0, 10.0, 20, 3), truth, 1e-7),
          "composite: more panels is more accurate");
    check(near(compositeGaussLegendre(sn, 0.0, pi, 1, 10), gaussLegendre(sn, 0.0, pi, 10), 1e-12), "composite with one panel equals the single rule");
    check(near(compositeGaussLegendre([](double x) { return x * x * x; }, 0.0, 2.0, 4, 2), 4.0, 1e-12), "composite is exact for cubics with 2 points");
    check(compositeGaussLegendre(sn, 1.0, 1.0, 5, 3) == 0.0 &&
              near(compositeGaussLegendre([](double x) { return x; }, 2.0, 0.0, 3, 2), -2.0, 1e-12),
          "composite: zero-width gives 0, reversed limits change the sign");

    // Invalid input must throw
    auto id = [](double x) { return x; };
    check(throwsInvalid([] { gaussLegendreRule(0); }) && throwsInvalid([] { gaussLegendreRule(-3); }) &&
              throwsInvalid([] { gaussLegendreRule(101); }),
          "rule size 0, negative and 101 throw");
    check(throwsInvalid([&] { gaussLegendre(id, 0.0, 1.0, 0); }), "gaussLegendre with 0 points throws");
    check(throwsInvalid([&] { gaussLegendre(std::function<double(double)>(), 0.0, 1.0); }), "empty function throws");
    check(throwsInvalid([&] { gaussLegendre(id, NAN, 1.0); }) && throwsInvalid([&] { gaussLegendre(id, 0.0, INFINITY); }),
          "NaN and infinite limits throw");
    check(throwsInvalid([&] { compositeGaussLegendre(id, 0.0, 1.0, 0); }) && throwsInvalid([&] { compositeGaussLegendre(id, 0.0, 1.0, 4, 0); }),
          "composite with 0 panels or 0 points throws");

    return finish();
}
