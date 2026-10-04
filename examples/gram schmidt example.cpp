#include <numerical-computing/gram_schmidt.hpp>

#include <iostream>
#include <stdexcept>

int main() {
    // Three vectors that point in different directions
    std::vector<numcomp::Vector> input = {{1, 1, 0}, {1, 0, 1}, {0, 1, 1}};

    // Turn them into three perpendicular vectors of length 1
    auto q = numcomp::gramSchmidt(input);

    std::cout << "Orthonormal vectors:\n";
    for (const auto& v : q) {
        for (double x : v) std::cout << x << " ";
        std::cout << "\n";
    }

    // Proof: perpendicular vectors have a dot product of (about) 0
    std::cout << "q0 . q1 = " << numcomp::dot(q[0], q[1]) << "\n";

    // Bad input: the second vector is just twice the first, so it adds nothing new
    try {
        numcomp::gramSchmidt({{1, 2, 3}, {2, 4, 6}});
    } catch (const std::invalid_argument& e) {
        std::cout << "Error caught: " << e.what() << "\n";
    }
    return 0;
}
