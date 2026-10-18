#ifndef ROMBERG_HPP
#define ROMBERG_HPP

#include <functional>
#include <vector>

/**
 * @brief Computes the definite integral of f(x) from a to b using Romberg Integration.
 * 
 * @param f Function to integrate (std::function<double(double)>)
 * @param a Lower integration limit
 * @param b Upper integration limit
 * @param n Number of levels/rows in the Romberg table (default: 5)
 * @return double Refined estimate of the definite integral
 */
double rombergIntegration(std::function<double(double)> f, double a, double b, int n = 5);

#endif