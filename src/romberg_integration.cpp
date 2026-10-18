#include "romberg.hpp"
#include <cmath>

double rombergIntegration(std::function<double(double)> f, double a, double b, int n) {
    if (n <= 0) return 0.0;

    std::vector<std::vector<double>> R(n, std::vector<double>(n, 0.0));
    double h = b - a;
    
    
    R[0][0] = 0.5 * h * (f(a) + f(b));
    
    for (int i = 1; i < n; ++i) {
        h *= 0.5; // Halve the step size
        double sum = 0.0;
        int points = 1 << (i - 1); 
        
        for (int k = 1; k <= points; ++k) {
            sum += f(a + (2 * k - 1) * h);
        }
        
        
        R[i][0] = 0.5 * R[i - 1][0] + h * sum;
        
        
        for (int j = 1; j <= i; ++j) {
            double factor = std::pow(4.0, j);
            R[i][j] = (factor * R[i][j - 1] - R[i - 1][j - 1]) / (factor - 1.0);
        }
    }
    
    return R[n - 1][n - 1]; 
}