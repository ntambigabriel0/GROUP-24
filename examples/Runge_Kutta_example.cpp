#include <iostream>
#include <iomanip>
#include “algorithm.hpp”

int main() {
    std::cout << “--- Dynamic Physical System Simulation: Newton’s Law of Cooling ---\n”;
    
    // Mathematical Framework: dT/dt = -k * (T - T_env)
    const double k = 0.05;
    const double T_env = 20.0;
    
    auto cooling_system_ode = [k, T_env](double t, double T) {
        return -k * (T - T_env);
    };

    double t0 = 0.0;          // Start timer tracking at 0 minutes
    double T0 = 90.0;          // Hot fluid tracking initial unit in Celsius
    double t_end = 30.0;       // Compute data out to 30 minutes
    double step_size = 2.0;    // Snapshot increments every 2 minutes

    std::cout << “Starting parameters: Initial Temp=” << T0 << “C, Ambient Room=” << T_env 
              << “C, Cooling Constant=” << k << “\n”;
    
    // Configure independent instance space
    RK4Solver cooling_solver(step_size);
    cooling_solver.solve(cooling_system_ode, t0, T0, t_end);

    // Format neat array tables directly to terminate output
    std::cout << “\nSimulation Log:\n”;
    std::cout << std::left << std::setw(15) << “Time (min)” << std::setw(20) << “Temperature (C)” << “\n”;
    std::cout << std::string(32, ‘-‘) << “\n”;

    for (size_t i = 0; i < cooling_solver.t_history.size(); ++i) {
        double time = cooling_solver.t_history[i];
        
        // Print compressed logging tracks spaced every 6 minutes 
        if (static_cast<int>(time) % 6 == 0) {
            std::cout << std::left << std::setw(15) << std::fixed << std::setprecision(1) << time 
                      << std::setw(20) << std::setprecision(2) << cooling_solver.y_history[i] << “\n”;
        }
    }

    std::cout << “\n\u25C6 Target reached: At “ << cooling_solver.t_history.back() 
              << “ mins, system output registers “ << cooling_solver.y_history.back() << “C.\n”;

    return 0;
}