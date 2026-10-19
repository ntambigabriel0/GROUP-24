#ifndef ALGORITHM_HPP
#define ALGORITHM_HPP

#include <vector>
#include <functional>

class RK4Solver {
private:
    // Private Data Member 
    double h; 

    // Private Member Function 
    double rk4_step(std::functional<double(double, double)> f, double t, double y, double step_size);

public:
    // Public Data Members 
    std::vector<double> t_history;
    std::vector<double> y_history;

    // Constructor 
    RK4Solver(double step_size);

    // Public Member Function 
    void solve(std::functional<double(double, double)> f, double t0, double y0, double t_end);
};

#endif