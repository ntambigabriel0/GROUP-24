#include <numerical-computing/jacobi.hpp>

#include <algorithm>  // std::max
#include <cmath>      // std::fabs, std::isfinite
#include <stdexcept>  // std::invalid_argument
#include <utility>    // std::move

namespace numcomp {

JacobiResult jacobi(const Matrix& A,
                    const Vector& b,
                    const Vector& x0,
                    double tol,
                    int maxIterations) {
    //  Step 1: check the input before doing any maths 
    const std::size_t n = A.size();
    if (n == 0) {
        throw std::invalid_argument("jacobi: matrix must not be empty");
    }
    for (const Vector& row : A) {
        if (row.size() != n) {
            throw std::invalid_argument("jacobi: matrix must be square");
        }
    }
    if (b.size() != n) {
        throw std::invalid_argument("jacobi: b must have one entry per row of the matrix");
    }
    if (!x0.empty() && x0.size() != n) {
        throw std::invalid_argument("jacobi: x0 must have one entry per row of the matrix");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (A[i][i] == 0.0) {  // exact zero: we are about to divide by this entry
            throw std::invalid_argument(
                "jacobi: a diagonal entry is zero, reorder the rows so it is not");
        }
    }
    if (!(tol > 0.0)) {
        throw std::invalid_argument("jacobi: tol must be positive");
    }
    if (maxIterations < 1) {
        throw std::invalid_argument("jacobi: maxIterations must be at least 1");
    }

    //  Step 2: set up the old and new copies 
    Vector oldX = x0.empty() ? Vector(n, 0.0) : x0;  // values from the previous round
    Vector newX(n, 0.0);                             // values being worked out this round
    int rounds = 0;
    double change = 0.0;
    bool converged = false;

    //  Step 3: repeat rounds until the values stop changing 
    while (rounds < maxIterations) {
        for (std::size_t i = 0; i < n; ++i) {
            double sum = b[i];
            for (std::size_t j = 0; j < n; ++j) {
                if (j != i) sum -= A[i][j] * oldX[j];  
            }
            newX[i] = sum / A[i][i];                   
        }
        ++rounds;

        // How far did the biggest value move this round?
        change = 0.0;
        bool finite = true;
        for (std::size_t i = 0; i < n; ++i) {
            if (!std::isfinite(newX[i])) finite = false;
            change = std::max(change, std::fabs(newX[i] - oldX[i]));
        }

        oldX = newX;  // this round's values become the previous ones
        if (!finite) break;       // the values blew up, it is diverging
        if (change <= tol) {      // nothing moves any more
            converged = true;
            break;
        }
    }

    return JacobiResult{std::move(oldX), rounds, converged, change};
}

}  // namespace numcomp
