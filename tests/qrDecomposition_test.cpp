#include <numerical-computing/gram_schmidt.hpp>
#include <numerical-computing/qr_decomposition.hpp>
#include <iostream>
#include <vector>

int main() {
    std::cout << "--- Testing Tier 2: QR Decomposition Matrix Factorization ---\n";

    // Create a simple test matrix A via column vectors
    numcomp::Vector col1 = {1.0, 0.0};
    numcomp::Vector col2 = {2.0, 3.0};
    std::vector<numcomp::Vector> A = {col1, col2};

    try {
        // Run Tier 2 Factorization
        numcomp::QRResult result = numcomp::qrDecomposition(A);

        // Display published Q columns
        std::cout << "\n[Published Matrix Q Columns]:\n";
        for (std::size_t j = 0; j < result.Q.size(); ++j) {
            std::cout << "Q_col[" << j << "]: ";
            for (double val : result.Q[j]) std::cout << val << " ";
            std::cout << "\n";
        }

        // Display published R columns (Upper Triangular)
        std::cout << "\n[Published Matrix R Columns]:\n";
        for (std::size_t j = 0; j < result.R.size(); ++j) {
            std::cout << "R_col[" << j << "]: ";
            for (double val : result.R[j]) std::cout << val << " ";
            std::cout << "\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Pipeline Execution Failed: " << e.what() << std::endl;
    }

    return 0;
}