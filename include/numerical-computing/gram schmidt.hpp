#ifndef NUMERICAL_COMPUTING_GRAM_SCHMIDT_HPP  // include guard: stops the file being included twice
#define NUMERICAL_COMPUTING_GRAM_SCHMIDT_HPP

#include <vector>

namespace numcomp {  // everything in the library lives in this namespace

using Vector = std::vector<double>;  // nickname: one vector = a list of numbers

/// Dot product of two vectors of the same size.
/// Throws std::invalid_argument if the sizes differ.
double dot(const Vector& a, const Vector& b);

/// Euclidean (2-norm) length of a vector.
double norm(const Vector& v);

/// Gram-Schmidt orthonormalization (modified variant).
///
/// Takes a list of vectors and returns an orthonormal list spanning the
/// same space. Output vector k is built from input vector k, so the
/// first k+1 outputs span the same space as the first k+1 inputs.
///
/// Throws std::invalid_argument if:
///   - the list is empty,
///   - any vector is empty,
///   - the vectors do not all have the same size,
///   - there are more vectors than dimensions,
///   - the vectors are linearly dependent (within the tolerance).
///
/// `tol` is a relative tolerance: a vector is treated as dependent when
/// what remains after removing earlier components is smaller than
/// tol times its original length.
std::vector<Vector> gramSchmidt(const std::vector<Vector>& vectors,
                                double tol = 1e-10);  // default tolerance if the user gives none

}  // namespace numcomp

#endif  // NUMERICAL_COMPUTING_GRAM_SCHMIDT_HPP
