#ifndef NUMERICAL_COMPUTING_RUNGE_KUTTA2_HPP
#define NUMERICAL_COMPUTING_RUNGE_KUTTA2_HPP

#include <functional>

/**
 * @brief Solves the ODE dy/dt = f(t, y) with the 2nd-order Runge-Kutta
 *        (midpoint) method.
 *
 * The number of steps is computed up front as ceil((t_end - t0) / h), so no
 * floating-point drift accumulates in t. All steps have size h except the
 * last, which is shortened to land exactly on t_end.
 *
 * @param f     Derivative function, dy/dt = f(t, y)
 * @param y0    Initial value of y at t0
 * @param t0    Initial time
 * @param t_end Final time (must be >= t0)
 * @param h     Step size (must be finite and > 0)
 * @return double Approximation of y at t_end
 *
 * @throws std::invalid_argument if f is empty, h <= 0 or non-finite,
 *         any of t0/t_end/y0 is non-finite, t_end < t0, or the step count
 *         is too large to represent.
 */
double rungeKutta2(const std::function<double(double, double)>& f,
                   double y0, double t0, double t_end, double h);

#endif  // NUMERICAL_COMPUTING_RUNGE_KUTTA2_HPP
