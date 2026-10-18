WEEK 3 PROGRESS REPORT

1. Completed

- Implemented the Power Iteration algorithm in C++ using Visual Studio Code.
- Used the "<iostream>" and "<cmath>" libraries for input/output and mathematical operations.
- Added input validation to ensure that the matrix size is between 1 and 20 and that the initial vector is not zero.
- Implemented matrix-vector multiplication, vector normalization, eigenvalue estimation, and convergence checking.
- Set the convergence tolerance to "0.0001" and the maximum number of iterations to 100.
- Tested the code using the matrix (\begin{bmatrix}2&1\1&2\end{bmatrix}) with the initial vector ([1,0]). The dominant eigenvalue converged to approximately 3.
- Prepared a draft header file, "power_iteration.hpp", for the function declaration.

2. How the Power Iteration Code Was Built

Power Iteration is a numerical method used to estimate the dominant eigenvalue and its corresponding eigenvector of a square matrix.

The implementation was developed through the following steps:

1. Defining the inputs: Declared arrays to store the matrix, initial vector, and resulting vector.
2. Validating the input: Checked that the matrix order is between 1 and 20 and that the initial vector is not zero.
3. Normalizing the vector: Divided each component of the initial vector by its largest absolute component.
4. Matrix-vector multiplication: Multiplied the matrix by the current vector to obtain a new vector.
5. Estimating the eigenvalue: Used the Rayleigh quotient to calculate an approximation of the dominant eigenvalue.
6. Checking convergence: Compared successive eigenvalue estimates using the tolerance of "0.0001".
7. Displaying the results: Printed the estimated dominant eigenvalue and eigenvector.

3. Issues and Challenges Faced

- Understanding how repeated matrix-vector multiplication leads to the dominant eigenvalue.
- Implementing vector normalization correctly to prevent excessively large or small values.
- Determining an appropriate convergence tolerance and stopping condition.
- Organizing the standalone program into a reusable function suitable for the Numerical Methods Library.

4. Use of AI

- Used AI to understand the Power Iteration algorithm and how to translate the mathematical steps into C++.
- Used AI assistance to review the code, explain errors and possible improvements, and prepare the header-file declaration.
- Used AI to identify a suitable test matrix with a known dominant eigenvalue to check the program's output.

5. Work in Progress

- Converting the standalone implementation into a reusable function that matches the header file.
- Integrating the Power Iteration module into the Numerical Methods Library.
- Preparing a separate test file and testing additional matrices, invalid inputs, and convergence failure cases.
