#include <numerical-computing/gaussian_quadrature.hpp>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main() {
    using namespace numcomp;
    const double pi = std::acos(-1.0);
    std::cout << std::setprecision(12);

    // 1. The rule itself: nodes and weights on [-1, 1].
    std::cout << "--- 1. The 3-point Gauss-Legendre rule on [-1, 1] ---\n";
    QuadratureRule rule = gaussLegendreRule(3);
    for (std::size_t i = 0; i < rule.nodes.size(); ++i)
        std::cout << "  node " << rule.nodes[i] << "   weight " << rule.weights[i] << "\n";

    // 2. A standard integral.
    std::cout << "\n--- 2. Integral of sin(x) from 0 to pi (exact value 2) ---\n";
    for (int n : {2, 3, 5, 8}) {
        double v = gaussLegendre([](double x) { return std::sin(x); }, 0.0, pi, n);
        std::cout << "  n = " << n << "  result = " << v << "  error = " << std::fabs(v - 2.0) << "\n";
    }

    // 3. Exactness: n points integrate polynomials up to degree 2n - 1 exactly.
    std::cout << "\n--- 3. Exactness: 3 points on x^k over [0, 1] (exact value 1/(k+1)) ---\n";
    for (int k : {3, 4, 5, 6}) {
        double v = gaussLegendre([k](double x) { return std::pow(x, k); }, 0.0, 1.0, 3);
        std::cout << "  x^" << k << ": " << v << "  exact " << 1.0 / (k + 1)
                  << (k <= 5 ? "  (degree <= 5, exact)" : "  (degree 6 > 5, not exact)") << "\n";
    }

    // 4. Convergence on a smooth function.
    const double expExact = std::exp(1.0) - 1.0;
    std::cout << "\n--- 4. Integral of e^x from 0 to 1 (exact value e - 1 = " << expExact << ") ---\n";
    for (int n = 1; n <= 5; ++n) {
        double v = gaussLegendre([](double x) { return std::exp(x); }, 0.0, 1.0, n);
        std::cout << "  n = " << n << "  error = " << std::fabs(v - expExact) << "\n";
    }

    // 5. Applied problem: probability that a standard normal value lies in [-1, 1].
    std::cout << "\n--- 5. Normal distribution: P(-1 <= Z <= 1) (exact value 0.682689492137) ---\n";
    auto pdf = [pi](double x) { return std::exp(-0.5 * x * x) / std::sqrt(2.0 * pi); };
    std::cout << "  5 points:  " << gaussLegendre(pdf, -1.0, 1.0, 5) << "\n";
    std::cout << "  10 points: " << gaussLegendre(pdf, -1.0, 1.0, 10) << "\n";

    // 6. Applied problem: work done by a spring, force F(x) = 200 x newtons, stretched 0 to 0.3 m.
    std::cout << "\n--- 6. Work done stretching a spring, F(x) = 200x from 0 to 0.3 m (exact 9 J) ---\n";
    std::cout << "  2 points: " << gaussLegendre([](double x) { return 200.0 * x; }, 0.0, 0.3, 2) << " J\n";

    // 7. Composite rule for a function with many oscillations.
    const double exact = 1.0 - std::cos(20.0);
    auto f = [](double x) { return std::sin(x); };
    std::cout << "\n--- 7. Integral of sin(x) from 0 to 20 (exact value " << exact << ") ---\n";
    std::cout << "  one 5-point rule:         error = " << std::fabs(gaussLegendre(f, 0.0, 20.0, 5) - exact) << "\n";
    std::cout << "  10 panels, 5 points each: error = " << std::fabs(compositeGaussLegendre(f, 0.0, 20.0, 10, 5) - exact) << "\n";

    // 8. Reversed limits.
    std::cout << "\n--- 8. Reversed limits ---\n";
    std::cout << "  integral of x^2 from 0 to 3: " << gaussLegendre([](double x) { return x * x; }, 0.0, 3.0, 3) << "\n";
    std::cout << "  integral of x^2 from 3 to 0: " << gaussLegendre([](double x) { return x * x; }, 3.0, 0.0, 3) << "\n";

    // 9. Invalid input is reported with exceptions.
    std::cout << "\n--- 9. Invalid input ---\n";
    try {
        gaussLegendreRule(0);
    } catch (const std::invalid_argument& e) {
        std::cout << "  Error caught: " << e.what() << "\n";
    }
    try {
        compositeGaussLegendre([](double x) { return x; }, 0.0, 1.0, 0);
    } catch (const std::invalid_argument& e) {
        std::cout << "  Error caught: " << e.what() << "\n";
    }
    return 0;
}
