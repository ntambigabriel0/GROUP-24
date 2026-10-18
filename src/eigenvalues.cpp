#include <numerical-computing/eigenvalues.hpp>
#include <numerical-computing/qr_decomposition.hpp>

#include <algorithm>  // std::sort
#include <cmath>      // std::fabs
#include <functional> // std::greater
#include <stdexcept>  // std::invalid_argument
#include <utility>    // std::move

namespace numcomp {

namespace {  


constexpr double kQrTol = 1e-10;


double maxBelowDiagonal(const Matrix& A) {
    double worst = 0.0;
    for (std::size_t j = 0; j < A.size(); ++j) {
        for (std::size_t i = j + 1; i < A.size(); ++i) {  
            worst = std::max(worst, std::fabs(A[j][i]));
        }
    }
    return worst;
}

Matrix nextRound(const Matrix& A) {
    const std::size_t n = A.size();

    const auto [Q, R] = qrDecomposition(A, kQrTol);

  
    Matrix next(n, Vector(n, 0.0));
    for (std::size_t j = 0; j < n; ++j) {
        for (std::size_t i = 0; i < n; ++i) {
            double sum = 0.0;
            for (std::size_t k = 0; k < n; ++k) {
                sum += R[k][i] * Q[j][k];
            }
            next[j][i] = sum;
        }
    }
    return next;
}

bool runRounds(Matrix& A, double tol, int maxIterations, int& rounds) {
    rounds = 0;
    while (true) {
        if (maxBelowDiagonal(A) <= tol) return true;   
        if (rounds >= maxIterations) return false;     
        A = nextRound(A);
        ++rounds;
    }
}

double rowSumNorm(const Matrix& A) {
    double worst = 0.0;
    for (std::size_t i = 0; i < A.size(); ++i) {
        double sum = 0.0;
        for (std::size_t j = 0; j < A.size(); ++j) sum += std::fabs(A[j][i]);
        worst = std::max(worst, sum);
    }
    return worst;
}

} 
EigenResult eigenvalues(const Matrix& A, double tol, int maxIterations) {
    //  Step 1: check the input before doing any maths
    const std::size_t n = A.size();
    if (n == 0) {
        throw std::invalid_argument("eigenvalues: matrix must not be empty");
    }
    for (const Vector& column : A) {
        if (column.size() != n) {
            throw std::invalid_argument("eigenvalues: matrix must be square");
        }
    }
    if (!(tol > 0.0)) {
        throw std::invalid_argument("eigenvalues: tol must be positive");
    }
    if (maxIterations < 0) {
        throw std::invalid_argument("eigenvalues: maxIterations must not be negative");
    }

    //  Step 2: run the QR rounds 
    Matrix work = A;     
    double shift = 0.0;  
    int rounds = 0;
    bool converged = false;
    try {
        converged = runRounds(work, tol, maxIterations, rounds);
    } catch (const std::invalid_argument&) {
        shift = rowSumNorm(A) + 1.0;
        work = A;
        for (std::size_t i = 0; i < n; ++i) work[i][i] += shift;
        converged = runRounds(work, tol, maxIterations, rounds);
    }

    //  Step 3: the diagonal holds the eigenvalues 
    Vector values(n);
    for (std::size_t i = 0; i < n; ++i) {
        values[i] = work[i][i] - shift;
    }
    std::sort(values.begin(), values.end(), std::greater<double>());  // largest first

    return EigenResult{std::move(values), rounds, converged};
}

}  // namespace numcomp
