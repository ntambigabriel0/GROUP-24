#include “algorithm.hpp”

// Constructor 
RK4Solver::RK4Solver(double step_size) : h(step_size) {}

// Private stepping calculation RK4Solver::rk4_step(std::functional<double(double, double)> f, double t, double y, double step_size) {
    double k1 = f(t, y);
    double k2 = f(t + step_size / 2.0, y + step_size * k1 / 2.0);
    double k3 = f(t + step_size / 2.0, y + step_size * k2 / 2.0);
    double k4 = f(t + step_size, y + step_size * k3);
    
    return y + (step_size / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

// Public interface runner loop
void RK4Solver::solve(std::functional<double(double, double)> f, double t0, double y0, double t_end) {
    // Clear array
    t_history.clear();
    y_history.clear();

    double t = t0;
    double y = y0;
    double current_h = h;

    t_history.push_back(t);
    y_history.push_back(y);

    while (t < t_end) {
        // Prevent path overshoot 
        if (t + current_h > t_end) {
            current_h = t_end - t;
        }

        y = rk4_step(f, t, y, current_h);
        t += current_h;

        t_history.push_back(t);
        y_history.push_back(y);
    }
}