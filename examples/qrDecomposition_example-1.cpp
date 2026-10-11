#include <numerical-computing/qr_decomposition.hpp>

#include <exception>
#include <iostream>
#include <vector>

int main() {
    // Define a simple 3x3 matrix A using distinct column vectors:
    // col1 = [1, 0, 0]^T
    // col2 = [2, 3, 0]^T
    // col3 = [4, 5, 6]^T
    numcomp::Vector col1 = {1.0, 0.0, 0.0};
    numcomp::Vector col2 = {2.0, 3.0, 0.0};
    numcomp::Vector col3 = {4.0, 5.0, 6.0};
    
    // Group the vectors together into a standard matrix format
    std::vector<numcomp::Vector> A = {col1, col2, col3};

    try {
        // Run the Tier 2 process 
        numcomp::QRResult factorization = numcomp::qrDecomposition(A);

        // Display published Q columns
        std::cout << "[Published Matrix Q (Orthonormal Columns)]:\n";
        for (std::size_t j = 0; j < factorization.Q.size(); ++j) {
            std::cout << "  Q_col[" << j << "]: ";
            for (double val : factorization.Q[j]) {
                std::cout << val << "  ";
            }
            std::cout << "\n";
        }

        // Display published R columns
        std::cout << "\n[Published Matrix R (Upper Triangular Columns)]:\n";
        for (std::size_t j = 0; j < factorization.R.size(); ++j) {
            std::cout << "  R_col[" << j << "]: ";
            for (double val : factorization.R[j]) {
                std::cout << val << "  ";
            }
            std::cout << "\n";
        }

        std::cout << "\nStatus: Success. Data formatted for Tier 3 consumption.\n";

    } catch (const std::exception& e) {
        std::cerr << "\nPipeline error caught: " << e.what() << std::endl;
    }

    return 0;
}
