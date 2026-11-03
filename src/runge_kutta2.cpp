
#include <numerical-computing/runge_kutta2.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

double rungeKutta2(const std::function<double(double, double)>& f,
                   double y0, double t0, double t_end, double h) {
    if (!f) {
        throw std::invalid_argument("rungeKutta2: f must be callable.");
    }
    if (!std::isfinite(h) || h <= 0.0) {
        throw std::invalid_argument("rungeKutta2: step size h must be finite and > 0.");
    }
    if (!std::isfinite(t0) || !std::isfinite(t_end) || !std::isfinite(y0)) {
        throw std::invalid_argument("rungeKutta2: t0, t_end and y0 must be finite.");
    }
    if (t_end < t0) {
        throw std::invalid_argument("rungeKutta2: t_end must be >= t0.");
    }

    const double span = t_end - t0;
    if (span == 0.0) {
        return y0;
    }

    const double ratio = span / h;
    const double steps_d = std::ceil(ratio - 1e-9 * ratio);
    if (!(steps_d < static_cast<double>(std::numeric_limits<long long>::max()))) {
        throw std::invalid_argument("rungeKutta2: step count too large.");
    }
    const long long steps = static_cast<long long>(steps_d < 1.0 ? 1.0 : steps_d);

    double y = y0;
    for (long long i = 0; i < steps; ++i) {
        const double t = t0 + static_cast<double>(i) * h;  
        const double step = (i == steps - 1) ? (t_end - t) : h;

        const double k1 = f(t, y);
        const double k2 = f(t + 0.5 * step, y + 0.5 * step * k1);
        y += step * k2;
    }
    return y;
}
