#include <numerical-computing/gaussian_quadrature.hpp>

#include <cmath>
#include <stdexcept>

namespace numcomp {

namespace {

const double kPi = 3.14159265358979323846;

void checkPointCount(int n) {
    if (n < 1 || n > kMaxQuadraturePoints)
        throw std::invalid_argument("gaussLegendre: number of points must be between 1 and 100");
}

void checkInputs(const std::function<double(double)>& f, double a, double b) {
    if (!f) throw std::invalid_argument("gaussLegendre: function is empty");
    if (!std::isfinite(a) || !std::isfinite(b))
        throw std::invalid_argument("gaussLegendre: limits must be finite numbers");
}

// Evaluates P_n(x) and P_{n-1}(x) with the three-term recurrence
// (k+1) P_{k+1} = (2k+1) x P_k - k P_{k-1}.
void legendre(int n, double x, double& pn, double& pnMinus1) {
    double p0 = 1.0;
    double p1 = x;
    for (int k = 1; k < n; ++k) {
        double p2 = ((2.0 * k + 1.0) * x * p1 - k * p0) / (k + 1.0);
        p0 = p1;
        p1 = p2;
    }
    pn = p1;
    pnMinus1 = p0;
}

}  // namespace

QuadratureRule gaussLegendreRule(int n) {
    checkPointCount(n);

    QuadratureRule rule;
    rule.nodes.assign(n, 0.0);
    rule.weights.assign(n, 0.0);

    const int half = (n + 1) / 2;  // roots come in +/- pairs, so find half of them
    for (int i = 0; i < half; ++i) {
        // Good starting guess for the i-th largest root of P_n.
        double x = std::cos(kPi * (i + 0.75) / (n + 0.5));
        double pn = 0.0, pn1 = 0.0;

        for (int iter = 0; iter < 100; ++iter) {
            legendre(n, x, pn, pn1);
            double dp = n * (x * pn - pn1) / (x * x - 1.0);  // derivative of P_n at x
            double dx = pn / dp;                             // Newton step
            x -= dx;
            if (std::fabs(dx) < 1e-15) break;
        }

        // Derivative at the final root gives the weight: w = 2 / ((1 - x^2) P_n'(x)^2).
        legendre(n, x, pn, pn1);
        double dp = n * (x * pn - pn1) / (x * x - 1.0);
        double w = 2.0 / ((1.0 - x * x) * dp * dp);

        if (n % 2 == 1 && i == half - 1) x = 0.0;  // the middle root of an odd rule is exactly 0

        rule.nodes[i] = -x;
        rule.nodes[n - 1 - i] = x;
        rule.weights[i] = w;
        rule.weights[n - 1 - i] = w;
    }
    return rule;
}

double gaussLegendre(const std::function<double(double)>& f, double a, double b, int n) {
    checkInputs(f, a, b);
    checkPointCount(n);
    if (a == b) return 0.0;

    const QuadratureRule rule = gaussLegendreRule(n);
    const double mid = 0.5 * (a + b);
    const double halfLength = 0.5 * (b - a);

    // Map each node from [-1, 1] onto [a, b], then scale the sum by half the interval length.
    double sum = 0.0;
    for (int i = 0; i < n; ++i)
        sum += rule.weights[i] * f(mid + halfLength * rule.nodes[i]);
    return halfLength * sum;
}

double compositeGaussLegendre(const std::function<double(double)>& f,
                              double a, double b, int panels, int n) {
    checkInputs(f, a, b);
    checkPointCount(n);
    if (panels < 1)
        throw std::invalid_argument("compositeGaussLegendre: panels must be at least 1");
    if (a == b) return 0.0;

    const QuadratureRule rule = gaussLegendreRule(n);  // build the rule once, reuse it on every panel
    const double width = (b - a) / panels;
    const double halfLength = 0.5 * width;

    double total = 0.0;
    for (int p = 0; p < panels; ++p) {
        const double mid = a + (p + 0.5) * width;
        double sum = 0.0;
        for (int i = 0; i < n; ++i)
            sum += rule.weights[i] * f(mid + halfLength * rule.nodes[i]);
        total += halfLength * sum;
    }
    return total;
}

}  // namespace numcomp
