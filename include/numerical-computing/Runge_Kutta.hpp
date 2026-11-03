#ifndef NUMERICAL_COMPUTING_RUNGE_KUTTA_HPP
#define NUMERICAL_COMPUTING_RUNGE_KUTTA_HPP

#include <functional>
#include <vector>

namespace numcomp {

class RK4Solver {
private:
    // Private Data Member
    double h;

    // Private Member Function
    double rk4_step(const std::function<double(double, double)>& f, double t, double y, double step_size);

public:
    // Public Data Members
    std::vector<double> t_history;
    std::vector<double> y_history;

    // Constructor (throws std::invalid_argument if step_size is not finite and > 0)
    explicit RK4Solver(double step_size);

    // Public Member Function
    // Throws std::invalid_argument if any input is not finite or t_end < t0.
    void solve(const std::function<double(double, double)>& f, double t0, double y0, double t_end);
};

}  // namespace numcomp

#endif
