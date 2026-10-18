#include "numerical-computing/cholesky.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

using numcomp::Matrix;

bool approximatelyEqual(
    double a,
    double b,
    double tol = 1e-9
) {
    return std::abs(a - b) <= tol;
}

void testKnownMatrix() {
    Matrix A = {
        {  4,  12, -16 },
        { 12,  37, -43 },
        { -16, -43,  98 }
    };

    Matrix L = numcomp::choleskyDecomposition(A);

    Matrix expected = {
        {  2, 0, 0 },
        {  6, 1, 0 },
        { -8, 5, 3 }
    };

    for (std::size_t i = 0; i < A.size(); ++i) {
        for (std::size_t j = 0; j < A.size(); ++j) {
            assert(approximatelyEqual(L[i][j], expected[i][j]));
        }
    }

    // Verify A = L * L^T.
    for (std::size_t i = 0; i < A.size(); ++i) {
        for (std::size_t j = 0; j < A.size(); ++j) {
            double sum = 0.0;

            for (std::size_t k = 0; k < A.size(); ++k) {
                sum += L[i][k] * L[j][k];
            }

            assert(approximatelyEqual(sum, A[i][j]));
        }
    }
}

void testOneByOneMatrix() {
    Matrix A = {{25.0}};

    Matrix L = numcomp::choleskyDecomposition(A);

    assert(approximatelyEqual(L[0][0], 5.0));
}

void testNonSquareMatrix() {
    Matrix A = {
        {1.0, 2.0},
        {3.0, 4.0},
        {5.0, 6.0}
    };

    bool threw = false;

    try {
        numcomp::choleskyDecomposition(A);
    } catch (const std::invalid_argument&) {
        threw = true;
    }

    assert(threw);
}

void testNonSymmetricMatrix() {
    Matrix A = {
        {4.0, 1.0},
        {2.0, 3.0}
    };

    bool threw = false;

    try {
        numcomp::choleskyDecomposition(A);
    } catch (const std::invalid_argument&) {
        threw = true;
    }

    assert(threw);
}

void testNotPositiveDefinite() {
    Matrix A = {
        {1.0, 2.0},
        {2.0, 1.0}
    };

    bool threw = false;

    try {
        numcomp::choleskyDecomposition(A);
    } catch (const std::invalid_argument&) {
        threw = true;
    }

    assert(threw);
}

void testEmptyMatrix() {
    Matrix A;

    bool threw = false;

    try {
        numcomp::choleskyDecomposition(A);
    } catch (const std::invalid_argument&) {
        threw = true;
    }

    assert(threw);
}

int main() {
    testKnownMatrix();
    testOneByOneMatrix();
    testNonSquareMatrix();
    testNonSymmetricMatrix();
    testNotPositiveDefinite();
    testEmptyMatrix();

    std::cout << "All Cholesky tests passed.\n";
    return 0;
}
