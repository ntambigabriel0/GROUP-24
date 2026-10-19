#include <numerical-computing/jacobi.hpp>

#include <iostream>
#include <stdexcept>

int main() {
    // The example: 4x + y = 6 and 2x + 5y = 12
    numcomp::Matrix A = {{4, 1}, {2, 5}};
    numcomp::Vector b = {6, 12};

    // Watch the first three rounds, starting from (0,0)
    std::cout << "First three rounds from (0,0):\n";
    for (int rounds = 1; rounds <= 3; ++rounds) {
        numcomp::JacobiResult r = numcomp::jacobi(A, b, {0, 0}, 1e-15, rounds);
        std::cout << "  round " << rounds << ": (" << r.x[0] << ", " << r.x[1] << ")\n";
    }

    // Now let it run until the values stop changing
    numcomp::JacobiResult r = numcomp::jacobi(A, b);
    std::cout << "Solution: x = " << r.x[0] << ", y = " << r.x[1]
              << " after " << r.iterations << " rounds"
              << " (converged: " << (r.converged ? "yes" : "no") << ")\n";

    // A system that is not diagonally dominant: Jacobi blows up instead of settling
    numcomp::Matrix bad = {{1, 2}, {3, 1}};
    numcomp::JacobiResult d = numcomp::jacobi(bad, {3, 4}, {}, 1e-10, 50);
    std::cout << "Not diagonally dominant: converged: " << (d.converged ? "yes" : "no") << "\n";

    // Bad input: a zero on the diagonal cannot be divided by
    try {
        numcomp::jacobi({{0, 1}, {1, 0}}, {1, 2});
    } catch (const std::invalid_argument& e) {
        std::cout << "Error caught: " << e.what() << "\n";
    }
    return 0;
}
