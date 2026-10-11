#include <numerical-computing/determinant.hpp>
#include <numerical-computing/lu.h>

#include <cstddef>
#include <stdexcept>

namespace numcomp {

double determinant(const Matrix& A, double tol) {
    // Step 1: Validate the tolerance.
    if (tol < 0.0) {
        throw std::invalid_argument(
            "determinant: tolerance cannot be negative");
    }

    // Step 2: Calculate the LU decomposition.
    LUResult result;

    try {
        result = luDecompose(A, tol);
    } catch (const SingularMatrixError&) {
        // A singular matrix has determinant zero.
        return 0.0;
    }

    // Step 3: Determine the sign from the row swaps.
    double det = (result.swaps % 2 == 0) ? 1.0 : -1.0;

    // Step 4: Multiply the diagonal elements of U.
    const std::size_t n = result.U.size();

    for (std::size_t i = 0; i < n; ++i) {
        det *= result.U[i][i];
    }

    // Step 5: Return the determinant.
    return det;
}

}  // namespace numcomp
