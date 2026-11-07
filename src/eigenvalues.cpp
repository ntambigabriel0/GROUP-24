#include <numerical-computing/eigenvalues.hpp>
#include <numerical-computing/qr_decomposition.hpp>

#include <algorithm>   // std::sort, std::max
#include <cmath>       // std::fabs, std::sqrt
#include <functional>  // std::greater
#include <stdexcept>   // std::invalid_argument
#include <utility>     // std::move

namespace numcomp {

namespace {  // helpers used only inside this file

// Tolerance handed to qrDecomposition (its own "are these columns dependent?" test).
// This is NOT the stopping tolerance of the eigenvalue loop: they mean different things.
// It is deliberately loose (1e-6): Gram-Schmidt loses accuracy on nearly singular
// matrices, so QR is asked to refuse those, and the shift is nudged instead.
constexpr double kQrTol = 1e-6;

// How many times to nudge the shift when QR refuses a (nearly) singular matrix.
constexpr int kMaxNudges = 6;

// The biggest absolute value in the last row of the m x m block, left of the diagonal.
// The matrix is a list of columns, so row m-1, column j is B[j][m-1].
// When these are all tiny, the last diagonal entry is an eigenvalue.
double lastRowLeftOfDiagonal(const Matrix& B, std::size_t m) {
    double worst = 0.0;
    for (std::size_t j = 0; j + 1 < m; ++j) {
        worst = std::max(worst, std::fabs(B[j][m - 1]));
    }
    return worst;
}

// The biggest absolute value anywhere in the matrix.
double maxAbsEntry(const Matrix& B) {
    double worst = 0.0;
    for (const Vector& column : B) {
        for (double x : column) worst = std::max(worst, std::fabs(x));
    }
    return worst;
}

// Wilkinson shift: the eigenvalue of the bottom-right 2x2 block that is closer to
// the last diagonal entry. Shifting by it makes the last eigenvalue appear very fast.
// If that 2x2 block has complex eigenvalues, the last diagonal entry is used instead.
double chooseShift(const Matrix& B, std::size_t m) {
    const double a = B[m - 2][m - 2];  // row m-2, column m-2
    const double b = B[m - 1][m - 2];  // row m-2, column m-1
    const double c = B[m - 2][m - 1];  // row m-1, column m-2
    const double d = B[m - 1][m - 1];  // row m-1, column m-1 (the last diagonal entry)

    const double p = (a - d) / 2.0;
    const double disc = p * p + b * c;  // negative means complex eigenvalues
    if (disc < 0.0) return d;
    const double r = std::sqrt(disc);

    // The 2x2 eigenvalues are d + p + r and d + p - r. Take the one closer to d.
    // (Written as bc/(...) to avoid subtracting two nearly equal numbers.)
    if (p >= 0.0) {
        return (p + r == 0.0) ? d : d - (b * c) / (p + r);
    }
    return d + (b * c) / (r - p);
}

// One shifted round on the m x m block B: factor (B - mu*I) = Q*R, then B = R*Q + mu*I.
// The eigenvalues do not change. Returns false if QR refuses the matrix.
bool shiftedRound(Matrix& B, double mu) {
    const std::size_t m = B.size();

    Matrix shifted = B;
    for (std::size_t i = 0; i < m; ++i) shifted[i][i] -= mu;

    try {
        // Binds the two fields of the result by position (Q first, R second),
        // so this works whatever the fields are named. Both are lists of columns.
        const auto [Q, R] = qrDecomposition(shifted, kQrTol);

        // Entry (row i, column j) of R*Q is the sum over k of R(i,k) * Q(k,j).
        // In column storage, R(i,k) is R[k][i] and Q(k,j) is Q[j][k].
        Matrix next(m, Vector(m, 0.0));
        for (std::size_t j = 0; j < m; ++j) {
            for (std::size_t i = 0; i < m; ++i) {
                double sum = 0.0;
                for (std::size_t k = 0; k < m; ++k) sum += R[k][i] * Q[j][k];
                next[j][i] = sum;
            }
            next[j][j] += mu;  // put the shift back
        }
        B = std::move(next);
        return true;
    } catch (const std::invalid_argument&) {
        return false;  // the shifted matrix is (nearly) singular: columns are dependent
    }
}

// A good shift is close to an eigenvalue, which makes B - mu*I nearly singular,
// and QR refuses those. So if QR refuses, nudge the shift slightly and try again.
bool roundWithNudging(Matrix& B, double mu) {
    const double nudge = 1e-4 * maxAbsEntry(B);  // relative to the size of the matrix
    double offset = 0.0;
    for (int attempt = 0; attempt < kMaxNudges; ++attempt) {
        if (shiftedRound(B, mu + offset)) return true;
        offset = (offset == 0.0) ? nudge : offset * 10.0;
    }
    return false;
}

}  // namespace

EigenResult eigenvalues(const Matrix& A, double tol, int maxIterations) {
    // ---- Step 1: check the input before doing any maths ----
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

    // ---- Step 2: find the eigenvalues one at a time, from the bottom up ----
    // The tolerance is relative: entries count as tiny compared with the biggest entry of A.
    const double threshold = tol * maxAbsEntry(A);
    Matrix B = A;       // the part still being worked on (A itself is never changed)
    std::size_t m = n;  // its size: it shrinks each time an eigenvalue is found
    Vector found;       // eigenvalues found so far
    int rounds = 0;     // QR rounds done in total
    bool converged = true;

    while (m > 0) {
        // Found one: nothing is left of the last row except its diagonal entry
        if (m == 1 || lastRowLeftOfDiagonal(B, m) <= threshold) {
            found.push_back(B[m - 1][m - 1]);
            B.pop_back();                        // drop the last column
            for (Vector& column : B) column.pop_back();  // drop the last row
            --m;
            continue;
        }
        // Out of rounds, or QR could not continue: give up honestly
        if (rounds >= maxIterations || !roundWithNudging(B, chooseShift(B, m))) {
            converged = false;
            break;
        }
        ++rounds;
    }

    // ---- Step 3: if it gave up, the diagonal of what is left is the best guess ----
    for (std::size_t i = 0; i < m; ++i) found.push_back(B[i][i]);
    std::sort(found.begin(), found.end(), std::greater<double>());  // largest first

    return EigenResult{std::move(found), rounds, converged};
}

}  // namespace numcomp
