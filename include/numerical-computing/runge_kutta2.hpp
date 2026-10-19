#ifndef RK2_HPP
#define RK2_HPP

#include <functional>

/**
 * @brief Solves an Ordinary Differential Equation (ODE) using 
 *        the 2nd Order Runge-Kutta (Midpoint) method.
 * 
 * @param f Function representing the derivative dy/dt = f(t, y)
 * @param y0 Initial value of y at t0
 * @param t0 Initial time
 * @param t_end Final time
 * @param h Step size
 * @return double Final computed value of y at t_end
 */
double rungeKutta2(std::function<double(double, double)> f, double y0, double t0, double t_end, double h);

#endif