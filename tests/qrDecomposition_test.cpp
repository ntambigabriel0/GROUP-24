#include <numerical-computing/gram_schmidt.hpp>
#include <numerical-computing/qr_decomposition.hpp>

#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using numcomp::QRResult;
using numcomp::Vector;

namespace {

int failures = 0;

// Print PASSED/FAILED for one check and count the failures.
void report(const std::string& name, bool ok) {
    std::cout << name << ": " << (ok ? "PASSED" : "FAILED") << "\n";
    if (!ok) ++failures;
}

bool close(double a, double b, double eps = 1e-9) {
    return std::fabs(a - b) <= eps;
}

// True if Q has orthonormal columns: Q_i . Q_j is 1 when i == j, otherwise 0.
bool isOrthonormal(const std::vector<Vector>& Q, double eps = 1e-9) {
    for (std::size_t i = 0; i < Q.size(); ++i) {
        for (std::size_t j = 0; j < Q.size(); ++j) {
            const double expected = (i == j) ? 1.0 : 0.0;
            if (!close(numcomp::dot(Q[i], Q[j]), expected, eps)) return false;
        }
    }
    return true;
}

// True if R is upper triangular: R[j][i] (row i of column j) is 0 when i > j.
bool isUpperTriangular(const std::vector<Vector>& R, double eps = 1e-9) {
    for (std::size_t j = 0; j < R.size(); ++j) {
        for (std::size_t i = j + 1; i < R[j].size(); ++i) {
            if (!close(R[j][i], 0.0, eps)) return false;
        }
    }
    return true;
}

// True if Q * R gives back A. Column j of A = sum over i of R[j][i] * Q[i].
bool reconstructs(const std::vector<Vector>& A, const QRResult& r, double eps = 1e-9) {
    for (std::size_t j = 0; j < A.size(); ++j) {
        for (std::size_t row = 0; row < A[j].size(); ++row) {
            double sum = 0.0;
            for (std::size_t i = 0; i <= j; ++i) sum += r.R[j][i] * r.Q[i][row];
            if (!close(sum, A[j][row], eps)) return false;
        }
    }
    return true;
}

// Run all the structural checks (shape, orthonormal Q, triangular R, Q*R = A).
void checkFactorization(const std::string& label, const std::vector<Vector>& A) {
    QRResult r = numcomp::qrDecomposition(A);
    bool shapeOk = r.Q.size() == A.size() && r.R.size() == A.size();
    for (const Vector& col : r.Q) shapeOk = shapeOk && col.size() == A[0].size();
    for (const Vector& col : r.R) shapeOk = shapeOk && col.size() == A.size();
    report(label + " - sizes of Q and R", shapeOk);
    if (!shapeOk) return;  // the other checks would read out of range
    report(label + " - Q columns are orthonormal", isOrthonormal(r.Q));
    report(label + " - R is upper triangular", isUpperTriangular(r.R));
    report(label + " - Q * R reproduces A", reconstructs(A, r));
}

// True if calling qrDecomposition(A) throws std::invalid_argument.
bool throwsInvalid(const std::vector<Vector>& A) {
    try {
        numcomp::qrDecomposition(A);
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

}  // namespace

int main() {
    std::cout << "--- QR DECOMPOSITION TESTS ---\n\n";

    // Test 1: the original 2x2 example, with exact expected values.
    // A columns {1,0} and {2,3}: Q is the identity, R columns are {1,0} and {2,3}.
    {
        std::vector<Vector> A = {{1.0, 0.0}, {2.0, 3.0}};
        QRResult r = numcomp::qrDecomposition(A);
        bool ok = close(r.Q[0][0], 1.0) && close(r.Q[0][1], 0.0) &&
                  close(r.Q[1][0], 0.0) && close(r.Q[1][1], 1.0) &&
                  close(r.R[0][0], 1.0) && close(r.R[0][1], 0.0) &&
                  close(r.R[1][0], 2.0) && close(r.R[1][1], 3.0);
        report("Test 1: 2x2 known values", ok);
        checkFactorization("Test 1", A);
    }

    // Test 2: 2x2 where Q is not the identity.
    // A columns {3,4} and {1,0}: Q0 = {0.6, 0.8}, Q1 = {0.8, -0.6}, R = [[5, 0.6], [0, 0.8]].
    {
        std::vector<Vector> A = {{3.0, 4.0}, {1.0, 0.0}};
        QRResult r = numcomp::qrDecomposition(A);
        bool ok = close(r.Q[0][0], 0.6) && close(r.Q[0][1], 0.8) &&
                  close(r.Q[1][0], 0.8) && close(r.Q[1][1], -0.6) &&
                  close(r.R[0][0], 5.0) && close(r.R[0][1], 0.0) &&
                  close(r.R[1][0], 0.6) && close(r.R[1][1], 0.8);
        report("Test 2: 2x2 known values (rotated)", ok);
        checkFactorization("Test 2", A);
    }

    // Test 3: 3x3 matrix, structure checks only.
    checkFactorization("Test 3: 3x3", {{12.0, 6.0, -4.0}, {-51.0, 167.0, 24.0}, {4.0, -68.0, -41.0}});

    // Test 4: tall matrix (4 rows, 2 columns).
    checkFactorization("Test 4: 4x2 tall", {{1.0, 2.0, 3.0, 4.0}, {1.0, 0.0, -1.0, 2.0}});

    // Test 5: identity matrix gives Q = I and R = I.
    {
        std::vector<Vector> A = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
        QRResult r = numcomp::qrDecomposition(A);
        bool ok = true;
        for (std::size_t i = 0; i < 3; ++i) {
            for (std::size_t j = 0; j < 3; ++j) {
                const double expected = (i == j) ? 1.0 : 0.0;
                ok = ok && close(r.Q[i][j], expected) && close(r.R[i][j], expected);
            }
        }
        report("Test 5: identity gives Q = I and R = I", ok);
    }

    // Test 6: the diagonal of R is the length of the remaining part of each column,
    // so it must be positive for independent columns.
    {
        QRResult r = numcomp::qrDecomposition({{2.0, -1.0, 0.0}, {-1.0, 2.0, -1.0}, {0.0, -1.0, 2.0}});
        bool ok = true;
        for (std::size_t j = 0; j < 3; ++j) ok = ok && r.R[j][j] > 0.0;
        report("Test 6: diagonal of R is positive", ok);
    }

    // Test 7: invalid input must throw std::invalid_argument.
    report("Test 7a: empty matrix throws", throwsInvalid({}));
    report("Test 7b: dependent columns throw", throwsInvalid({{1.0, 2.0}, {2.0, 4.0}}));
    report("Test 7c: more columns than rows throws", throwsInvalid({{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}}));
    report("Test 7d: columns of different sizes throw", throwsInvalid({{1.0, 0.0}, {0.0, 1.0, 2.0}}));

    std::cout << "\nAll tests finished running. Failures: " << failures << "\n";
    return failures == 0 ? 0 : 1;
}
