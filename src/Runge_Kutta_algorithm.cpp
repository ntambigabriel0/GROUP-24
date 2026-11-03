#include <numerical-computing/Runge_Kutta.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace numcomp {

// Constructor: reject h <= 0 (h = 0 would never advance t and hang solve)
RK4Solver::RK4Solver(double step_size) : h(step_size) {
    if (!std::isfinite(step_size) || step_size <= 0.0) {
        throw std::invalid_argument("RK4Solver: step size must be finite and > 0");
    }
}

// Private stepping calculation
double RK4Solver::rk4_step(const std::function<double(double, double)>& f, double t, double y, double step_size) {
    double k1 = f(t, y);
    double k2 = f(t + step_size / 2.0, y + step_size * k1 / 2.0);
    double k3 = f(t + step_size / 2.0, y + step_size * k2 / 2.0);
    double k4 = f(t + step_size, y + step_size * k3);

    return y + (step_size / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

// Public interface runner loop
void RK4Solver::solve(const std::function<double(double, double)>& f, double t0, double y0, double t_end) {
    if (!std::isfinite(t0) || !std::isfinite(y0) || !std::isfinite(t_end)) {
        throw std::invalid_argument("RK4Solver::solve: inputs must be finite");
    }
    if (t_end < t0) {
        throw std::invalid_argument("RK4Solver::solve: t_end must be >= t0");
    }

    t_history.clear();
    y_history.clear();

    t_history.push_back(t0);
    y_history.push_back(y0);

    // Fixed step count avoids floating-point drift from repeated t += h.
    // The small tolerance stops rounding noise in span/h adding a sliver step.
    const double span = t_end - t0;
    std::size_t n = static_cast<std::size_t>(std::ceil(span / h - 1e-9));
    if (span > 0.0) {
        n = std::max<std::size_t>(n, 1);  // a tiny positive span still takes one step
    }

    double t = t0;
    double y = y0;

    for (std::size_t i = 1; i <= n; ++i) {
        double step = h;
        double t_next = t0 + static_cast<double>(i) * h;

        // Last step lands exactly on t_end (no overshoot, no sliver step)
        if (i == n) {
            t_next = t_end;
            step = t_end - t;
        }

        y = rk4_step(f, t, y, step);
        t = t_next;

        t_history.push_back(t);
        y_history.push_back(y);
    }
}

}  // namespace numcomp
