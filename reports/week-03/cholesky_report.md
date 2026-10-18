Week 3 Progress Report: Cholesky Decomposition

Module / Tier: Cholesky Decomposition

Project: GROUP-24 Numerical Computing Library

Files:

- "include/numerical-computing/cholesky.hpp"
- "src/cholesky.cpp"
- "examples/example_cholesky.cpp"
- "tests/test_cholesky.cpp"

Completed

- Implemented the Cholesky Decomposition module in the "numcomp" namespace.
- Declared the "choleskyDecomposition()" function in the header file.
- Implemented the decomposition algorithm in the source file.
- Added input validation for matrix dimensions, symmetry, finite values, and positive definiteness.
- Created an example program demonstrating how to use the module.
- Successfully compiled and executed the example program in VS Code.
- Verified that the example produced the expected lower triangular matrix for the selected input.

How the Cholesky Code Was Built

1. Defined the problem. Cholesky decomposition is a numerical method used to factorize a symmetric positive-definite matrix into the product of a lower triangular matrix and its transpose.
   
   The decomposition is expressed as:
   
   [
   A=LL^T
   ]
   
   Where:
   
   - (A) is the original matrix.
   - (L) is the lower triangular matrix.
   - (L^T) is the transpose of (L).

2. Designed the interface. The function declaration was placed in "cholesky.hpp", while the implementation was placed in "cholesky.cpp". This separates the public interface from the internal calculations.

3. Represented the matrix. The module uses the "Matrix" type alias, represented by "std::vector<std::vector<double>>", to store the input matrix and the resulting factor.

4. Validated the input. The implementation checks that the matrix is non-empty, square, symmetric within a specified tolerance, and contains finite numerical values.

5. Calculated the diagonal elements. 
6. Calculated the off-diagonal elements. 
7. Handled errors. The implementation checks the calculated diagonal values to identify matrices that are not sufficiently positive definite for the chosen tolerance. Invalid inputs are reported using exceptions.

8. Tested the example. The example program was compiled and executed successfully. The output was compared with the expected Cholesky factor for the selected matrix.

Program Organization (Classes)

1. What classes did you create?

No classes were required for this module. It is organized as a free function inside the "numcomp" namespace, together with a "Matrix" type alias.

using Matrix = std::vector<std::vector<double>>;

2. Why did you organize it this way?

Cholesky decomposition is a mathematical operation that accepts a matrix and returns another matrix. It does not require an object to maintain state between function calls. A free function keeps the implementation simple, while the namespace helps prevent naming conflicts with other modules.

3. What data does it hold?

Item| Type| Public / Private| Purpose
"Matrix"| Alias for "std::vector<std::vector<double>>"| Public| Represents a matrix
"A"| "const Matrix&"| Input parameter| Stores the original matrix
"L"| "Matrix"| Local variable| Stores the lower triangular factor
"tol"| "double"| Input parameter| Controls numerical validation

There are no persistent data members. The matrices are passed to or created within the decomposition function.

4. What functions operate on that data?

Function| Public / Private| Purpose
"choleskyDecomposition(A, tol)"| Public| Calculates the lower triangular factor
Input validation logic| Internal implementation| Checks dimensions, symmetry, and numerical validity
Diagonal calculation logic| Internal implementation| Computes diagonal entries using square roots
Off-diagonal calculation logic| Internal implementation| Computes the entries below the diagonal

5. Why are some members public and others hidden?

The "choleskyDecomposition()" function is declared in the header file, allowing other modules to call it. The calculation and validation details remain in the source file, making it possible to modify the internal implementation without unnecessarily changing the public interface.

6. How does the rest of the program use it?

The example program supplies a matrix to "choleskyDecomposition()" and displays the resulting lower triangular matrix. The test program is intended to verify the correctness of the factorization and the handling of invalid inputs.



How This Module Connects to Others

- Uses: C++ standard library components for matrix storage, numerical calculations, and error handling.
- Potentially used by: Numerical computing modules that require Cholesky factorization.
- Related applications: Solving systems of linear equations involving symmetric positive-definite matrices and other numerical methods that use matrix factorization.

Issues Faced

- Header file inclusion: The compiler initially could not locate "cholesky.hpp" because the include path was not correctly specified in the compilation command.
- Compilation setup: The correct include directory and source file had to be supplied when compiling the example.
- Matrix validation: The implementation needed checks to ensure that the input satisfies the conditions required for Cholesky decomposition.
- Numerical precision: Floating-point calculations require tolerances when checking symmetry and numerical results.
- Project integration: The module needed to follow the naming conventions and directory structure of the GROUP-24 repository.

The header inclusion issue was resolved, after which the example program compiled and ran successfully.

In Progress

- Running the complete Cholesky test suite.
- Verifying the behavior of the implementation for invalid matrices.
- Confirming that all module files are integrated correctly into the project.
- Reviewing the implementation before submitting it for integration into the shared repository.

Challenges / Blockers

- Ensuring that the implementation handles unsuitable matrices correctly.
- Confirming that all test cases pass.
- Maintaining consistent file paths and include statements.
- Integrating the module with the existing build configuration and GitHub workflow.
