#include <numerical-computing/eigenvalues.hpp>

#include <iostream>
#include <stdexcept>

// Prints the answer of one eigenvalue computation
static void show(const char* label, const numcomp::EigenResult& r) {
    std::cout << label << "\n  eigenvalues:";
    for (double v : r.values) std::cout << " " << v;
    std::cout << "\n  rounds: " << r.iterations
              << ", converged: " << (r.converged ? "yes" : "no") << "\n";
}

int main() {
    numcomp::Matrix A = {{2, 1}, {1, 2}};
    show("[[2,1],[1,2]] (expect 3 and 1):", numcomp::eigenvalues(A));

    // A singular matrix has an eigenvalue of 0: the shift makes it work
    numcomp::Matrix S = {{1, 1}, {1, 1}};
    show("[[1,1],[1,1]] (expect 2 and 0):", numcomp::eigenvalues(S));

    // Complex eigenvalues (a 90-degree rotation)
    numcomp::Matrix R = {{0, 1}, {-1, 0}};
    show("rotation matrix (complex eigenvalues, expect converged: no):",
         numcomp::eigenvalues(R, 1e-10, 50));

    // Bad input:
    try {
        numcomp::eigenvalues({{1, 2, 3}, {4, 5, 6}});
    } catch (const std::invalid_argument& e) {
        std::cout << "Error caught: " << e.what() << "\n";
    }
    return 0;
}
