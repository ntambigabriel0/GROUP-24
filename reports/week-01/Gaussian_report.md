WEEK 1 PROGRESS REPORT

1. Completed

- Implemented the Gaussian Elimination with Partial Pivoting module in "include/numerical-computing/gaussian_elimination.hpp" and "src/gaussian_elimination.cpp" within the "numcomp" namespace.
- Added input validation to ensure that matrix (A) is square and that the size of vector (b) matches the dimensions of (A).
- Implemented forward elimination with partial pivoting, row swapping, and back substitution.
- Wrote tests in "tests/test_gaussian_elimination.cpp" covering normal 2×2 and 3×3 systems, a system requiring partial pivoting and row swapping, and a singular matrix to test error detection.
- Created a usage example demonstrating a system of linear equations, such as (2x+y=3) and (x+3y=5), yielding (x=0.8) and (y=1.4)

2. How the Gaussian Elimination Code Was Built

Gaussian Elimination is a numerical method used to solve systems of simultaneous linear equations. It transforms the original system into an upper triangular form through forward elimination, followed by back substitution to determine the unknown variables.

The implementation was developed through the following steps:

1. Defining the problem: Established Gaussian Elimination as the method for solving systems of linear equations by transforming the coefficient matrix into an upper triangular matrix.
2. Validating the input: Checked that matrix (A) is square and that vector (b) has the same number of entries as the number of rows in (A).
3. Selecting the pivot: For each column (k), identified the row containing the largest absolute value in that column to select the pivot row.
4. Checking for singularity: Checked whether the pivot was approximately zero during elimination and at the final pivot. This helps detect singular matrices or systems without a unique solution.
5. Performing partial pivoting and row swapping: When the selected pivot row differed from the current row, swapped the corresponding rows of matrix (A) and the matching entries of vector (b).
6. Executing forward elimination: Calculated elimination factors for the rows below each pivot and updated the entries of matrix (A) and vector (b) to obtain an upper triangular system.
7. Applying back substitution: Worked from the last row upwards, calculating each unknown variable and storing the results in the solution vector (x).
8. Testing the implementation: Verified the implementation using different systems, including 2×2 and 3×3 systems, a system requiring partial pivoting and row swapping, and a singular matrix.

3. Issues I Faced

- Handling singular matrices: Ensuring that the algorithm correctly identifies approximately zero pivots and reports when a matrix is singular or the system does not have a unique solution.
- Implementing partial pivoting: Ensuring that row swaps are correctly applied to both matrix (A) and vector (b) while maintaining numerical stability and accuracy.
- Back-substitution ordering: Ensuring that the indices are handled correctly so that the unknown variables are calculated from the bottom row upwards.

4. In Progress

  -Continuing to study the implementation of numerical methods to improve my understanding of algorithms.
Planning the implementation of next numerical computing module using an iterative approach.

5. Challenges

- Understanding the implementation process and ensuring that the Gaussian Elimination algorithm produces accurate results.
- Handling partial pivoting, row swapping, and singular matrix detection correctly.
- Ensuring that the implementation is properly tested and integrated into the existing project structure.

6. Work Planned for Next Week

- Continue developing the numerical computing library using an iterative approach, implementing the next module step by step.
- Write and run tests for each iteration to identify errors and ensure correctness.
- Develop usage examples for the newly implemented module.
- Review and improve the implementation based on test results and feedback.

7. AI Use

- Tool: Claude
- Purpose: Assistance with formatting and structuring the Gaussian Elimination progress report.
- Reason: To maintain consistency with the formatting style of the work.
