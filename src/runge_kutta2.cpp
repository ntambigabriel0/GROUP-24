#include "rk2.hpp"

double rungeKutta2(std::function<double(double, double)> f, double y0, double t0, double t_end, double h) {
    double t = t0;
    double y = y0;
    double current_h = h;

    while (t < t_end) {
        // Adjust final step if it overshoots t_end
        if (t + current_h > t_end) {
            current_h = t_end - t;
        }

        // RK2 (Midpoint Method) Logic:
        double k1 = f(t, y);
        double t_mid = t + 0.5 * current_h;
        double y_mid = y + 0.5 * current_h * k1;
        double k2 = f(t_mid, y_mid);
        
        y += current_h * k2;
        t += current_h;
    }

    return y;
}