#include <numerical-computing/gram_schmidt.hpp>
#include "test_util.hpp"

using numcomp::Vector;
using namespace testutil;

// True if every vector has length 1 and every pair is at 90 degrees
static bool isOrthonormal(const std::vector<Vector>& q, double eps = 1e-9) {
    for (std::size_t i = 0; i < q.size(); ++i)
        for (std::size_t j = 0; j < q.size(); ++j)
            if (!near(numcomp::dot(q[i], q[j]), i == j ? 1.0 : 0.0, eps)) return false;
    return true;
}

int main() {
    // Normal use
    check(isOrthonormal(numcomp::gramSchmidt({{1, 1, 0}, {1, 0, 1}, {0, 1, 1}})), "3 vectors in R^3 become orthonormal");
    check(isOrthonormal(numcomp::gramSchmidt({{1, 2, 3, 4}, {2, 1, 0, 1}})), "fewer vectors than dimensions");
    auto one = numcomp::gramSchmidt({{3, 4}});
    check(near(one[0][0], 0.6, 1e-12) && near(one[0][1], 0.8, 1e-12), "single vector is just normalized");

    // Awkward but valid: nearly parallel vectors
    check(isOrthonormal(numcomp::gramSchmidt({{1, 1e-8, 0}, {1, 0, 0}, {1, 0, 1e-8}}, 1e-12), 1e-6),
          "nearly dependent vectors stay orthogonal");

    // Invalid input: each must throw
    struct Bad { const char* name; std::vector<Vector> in; };
    for (const Bad& b : std::vector<Bad>{
             {"empty list throws", {}},
             {"empty vector throws", {Vector{}}},
             {"mismatched sizes throw", {{1, 2}, {1, 2, 3}}},
             {"more vectors than dimensions throws", {{1, 0}, {0, 1}, {1, 1}}},
             {"dependent vectors throw", {{1, 2, 3}, {2, 4, 6}}},
             {"zero vector throws", {{0, 0, 0}}}})
        check(throwsInvalid([&] { numcomp::gramSchmidt(b.in); }), b.name);

    return finish();
}
