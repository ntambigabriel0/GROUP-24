#include <numerical-computing/gram_schmidt.hpp>

#include <cmath>      // std::sqrt
#include <stdexcept>  // std::invalid_argument

namespace numcomp {

// dot(a, b) = a0*b0 + a1*b1 + ...
double dot(const Vector& a, const Vector& b) {
    if (a.size() != b.size()) {  // can only multiply matching positions
        throw std::invalid_argument("dot: vectors must have the same size");
    }
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];  // multiply matching entries and add them up
    }
    return sum;
}

// Length of a vector: sqrt of the dot product with itself
double norm(const Vector& v) {
    return std::sqrt(dot(v, v));
}

std::vector<Vector> gramSchmidt(const std::vector<Vector>& vectors,
                                double tol) {
    // ---- Step 1: check the input before doing any maths ----
    if (vectors.empty()) {
        throw std::invalid_argument("gramSchmidt: no input vectors");
    }
    const std::size_t dim = vectors[0].size();  // dimension = size of the first vector
    if (dim == 0) {
        throw std::invalid_argument("gramSchmidt: vectors must not be empty");
    }
    for (const Vector& v : vectors) {
        if (v.size() != dim) {  // every vector must match the first
            throw std::invalid_argument(
                "gramSchmidt: all vectors must have the same size");
        }
    }
    if (vectors.size() > dim) {  // e.g. 4 vectors in 3D: one must be redundant
        throw std::invalid_argument(
            "gramSchmidt: more vectors than dimensions, so they must be dependent");
    }

    // ---- Step 2: build the orthonormal list one vector at a time ----
    std::vector<Vector> q;  // the finished (orthonormal) vectors
    q.reserve(vectors.size());

    for (const Vector& original : vectors) {  // outer loop: pick the next input vector
        Vector v = original;                  // work on a copy, so the input is untouched
        const double originalLength = norm(original);  // remembered for the dependency check

        // Inner loop: remove the shadow of v on each finished vector.
        // "Modified" Gram-Schmidt: dot with the updated v each time (more stable).
        for (const Vector& basis : q) {
            const double coeff = dot(basis, v);  // shadow length (basis has length 1)
            for (std::size_t i = 0; i < dim; ++i) {
                v[i] -= coeff * basis[i];  // subtract the shadow
            }
        }

        // What is left is perpendicular to every finished vector.
        const double remaining = norm(v);
        if (remaining <= tol * originalLength) {  // nothing left: v was already covered
            throw std::invalid_argument(
                "gramSchmidt: vectors are linearly dependent");
        }

        for (double& x : v) {
            x /= remaining;  // scale to length 1
        }
        q.push_back(std::move(v));  // v is finished: add it to the list
    }
    return q;
}

}  // namespace numcomp
