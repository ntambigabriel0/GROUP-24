
#include "runge_kutta2.hpp"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

namespace {

std::size_t stepCount(const RK2Function& f, double y0, double t0, double t_end, double h) {
    if (!f) {
        throw std::invalid_argument("rungeKutta2: f must be callable");
    }
    if (!std::isfinite(h) || h <= 0.0) {
        throw std::invalid_argument("rungeKutta2: step size h must be positive and finite");
    }
    if (!std::isfinite(t0) || !std::isfinite(t_end) || !std::isfinite(y0)) {
        throw std::invalid_argument("rungeKutta2: y0, t0 and t_end must be finite");
    }
    if (t_end < t0) {
        throw std::invalid_argument("rungeKutta2: t_end must be >= t0");
    }

    const double ratio = (t_end - t0) / h;
    if (!(ratio < static_cast<double>(std::numeric_limits<std::size_t>::max()))) {
        throw std::invalid_argument("rungeKutta2: step count is too large");
    }
    const double tol = std::sqrt(std::numeric_limits<double>::epsilon()) * (ratio > 1.0 ? ratio : 1.0);
    return static_cast<std::size_t>(std::ceil(ratio - tol));
}

template <typename Observer>
double integrate(const RK2Function& f, double y0, double t0, double t_end, double h, Observer observe) {
    const std::size_t n = stepCount(f, y0, t0, t_end, h);

    double y = y0;
    observe(t0, y);
    for (std::size_t i = 0; i < n; ++i) {
        const double t = t0 + static_cast<double>(i) * h;
        const bool last = (i + 1 == n);
        const double step = last ? (t_end - t) : h;
        y = rk2_step(f, t, y, step);
        observe(last ? t_end : t + step, y);
    }
    return y;
}

}

double rk2_step(const RK2Function& f, double t, double y, double h) {
    const double k1 = f(t, y);
    const double k2 = f(t + h / 2, y + h / 2 * k1);
    return y + h * k2;
}

std::vector<RK2Point> rk2_trajectory(const RK2Function& f, double y0, double t0, double t_end, double h) {
    std::vector<RK2Point> path;
    path.reserve(stepCount(f, y0, t0, t_end, h) + 1);
    integrate(f, y0, t0, t_end, h, [&path](double t, double y) { path.emplace_back(t, y); });
    return path;
}

double rungeKutta2(const RK2Function& f, double y0, double t0, double t_end, double h) {
    return integrate(f, y0, t0, t_end, h, [](double, double) {});
}
