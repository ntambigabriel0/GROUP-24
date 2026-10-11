
#include <iostream>
#include <numerical-computing/determinant.hpp>

using numcomp::Matrix;

int main() { 
    // Step 1: Create a 2x2 matrix
    Matrix A = {
        {2, 1},
        {1, 3}
    };

    // Step 2: Calculate the determinant
    double result = numcomp::determinant(A);

    // Step 3: Display the matrix
    std::cout << "Matrix A:" << std::endl;
    std::cout << "[ 2  1 ]" << std::endl;
    std::cout << "[ 1  3 ]" << std::endl;

    // Step 4: Display the determinant
    std::cout << "Determinant = " << result << std::endl;

    return 0;
}
