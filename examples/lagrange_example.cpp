#include <numerical-computing/lagrange.hpp>

#include <iostream>
#include <stdexcept>
#include <vector>

int main() {
    // Three known points, all on the curve y = x^2 + x + 1
    std::vector<double> xs = {0, 1, 2};
    std::vector<double> ys = {1, 3, 7};

    // Ask: what is y when x = 3? (a value we never measured)
    double y = numcomp::lagrangeInterpolate(xs, ys, 3.0);
    std::cout << "y at x = 3 is " << y << "  (expected 13)\n";

    // Evaluate at several x values at once
    std::vector<double> queries = {0.5, 1.5, 2.5};
    std::vector<double> results = numcomp::lagrangeInterpolate(xs, ys, queries);
    for (std::size_t i = 0; i < queries.size(); ++i) {
        std::cout << "x = " << queries[i] << " -> y = " << results[i] << "\n";
    }

    // Bad input: two points with the same x but different y
    try {
        numcomp::lagrangeInterpolate({1, 1}, {2, 3}, 1.0);
    } catch (const std::invalid_argument& e) {
        std::cout << "Error caught: " << e.what() << "\n";
    }
    return 0;
}
