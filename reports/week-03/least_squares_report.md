Week 3 Progress Report: Least Squares Solution

Module / Tier: Least Squares Solution (Tier 3)
Files: "include/numerical-computing/least_squares.hpp", "src/least_squares.cpp", "tests/test_least_squares.cpp", "examples/least_squares_example.cpp"

Completed

- Worked on the Least Squares Solution module within the GROUP-24 numerical computing project.
- Identified the purpose of the least squares method in solving overdetermined systems of linear equations.
- Planned the module interface and separation of declarations in the header file from implementation in the source file.
- Identified the need to validate input matrices and vectors before performing calculations.
- Identified test cases for checking the accuracy of the least squares solution.
- Considered how the module would interact with other numerical computing modules in the project.

Note: Confirm the completed items against your actual code and progress before submitting the report.

How the Least Squares Code Was Built

1. Defined the problem. The least squares method finds an approximate solution to a system of linear equations (Ax=b), particularly when the system has more equations than unknowns or when an exact solution does not exist.

2. Designed the interface. The module uses a header file to declare the functions and a source file to contain their implementations. This separates the public interface from the internal calculations.

3. Represented the input data. The method uses a matrix (A) containing the coefficients of the system and a vector (b) containing the known values.

4. Selected the solution approach. A suitable implementation can use QR decomposition to solve the least squares problem. If (A=QR), the solution is obtained by solving (Rx=Q^Tb).

5. Performed the calculations. The implementation uses the selected numerical method to calculate the approximate solution that minimizes the sum of squared residuals:


6. Considered input validation. The program should check that the dimensions of the matrix and vector are compatible and that the input is suitable for the chosen solution method.

7. Planned error handling. The implementation should handle cases where the system cannot be solved reliably, such as when the matrix does not have sufficient rank.

8. Tested the solution. The results should be checked using known systems, overdetermined systems, and cases that test invalid input or numerical limitations.

Program Organization (Classes)

1. What classes did you create?

The module can be implemented using free functions within the "numcomp" namespace, without creating classes. Matrices and vectors can be represented using suitable C++ data structures.

The exact organization should reflect the final implementation in the repository.

2. Why did you organize it this way?

The least squares method is a numerical calculation that takes a matrix and a vector as input and produces an approximate solution. Free functions provide a straightforward way to perform the calculation without requiring objects to maintain state between function calls.

Using the "numcomp" namespace also helps prevent naming conflicts with functions in other modules.

3. What data does it hold?

Item| Type| Public / Private| Purpose
Matrix (A)| Matrix representation| Input| Stores the coefficients of the system
Vector (b)| Vector representation| Input| Stores the known values
Vector (x)| Vector representation| Output| Stores the approximate solution

These describe the mathematical data used by the method. The actual C++ types and visibility should match the code.

4. What functions operate on that data?

Function| Public / Private| What it does
Least squares solver| Public, if declared in the header| Calculates an approximate solution to (Ax=b)
Matrix and vector helper functions| Depends on implementation| Support the numerical calculations
Input validation helpers| Depends on implementation| Check dimensions and other input requirements

5. Why are some members public and others hidden?

Functions declared in the header file form the public interface and can be called by other modules. The implementation details remain in the source file, allowing the internal calculations to be modified without unnecessarily changing the interface used by other parts of the project.

6. How does the rest of the program use it?

The test file can call the least squares solver to verify the calculated solution, while the example program can demonstrate how to supply a matrix and a vector and obtain an approximate solution.

The module may also depend on QR decomposition if that is the approach selected for solving the system.


How This Module Connects to Others

- Uses: C++ standard library components and any numerical computing modules required by the chosen implementation.
- Used by: Other modules or example programs that require an approximate solution to a system of linear equations.
- Related module: QR Decomposition, if the least squares solver uses QR to obtain the solution.

Issues Faced

- Matrix dimensions: Ensuring that the number of rows in the matrix matches the length of the right-hand-side vector.
- Numerical accuracy: Floating-point calculations can introduce small errors, so results should be compared using a suitable tolerance.
- Rank deficiency: A matrix without sufficient rank can make the solution non-unique or prevent the selected method from working reliably.
- Integration: Ensuring that the function declarations and implementation are consistent with the project's structure and naming conventions.
- Git/GitHub workflow: Coordinating code changes, reviewing contributions, and integrating the module into the shared repository.


In Progress

- Completing or reviewing the implementation of the least squares solver.
- Verifying the interface and file organization.
- Developing or refining tests for different types of linear systems.
- Integrating the module with QR decomposition, if required.
- Reviewing code through pull requests.

Challenges / Blockers

- Ensuring that the numerical method produces accurate and reliable results.
- Handling invalid input and systems that do not have a unique solution.
- Confirming that the module builds correctly within the shared project.
- Coordinating code reviews and resolving integration issues.

Next Week

- Complete any outstanding implementation work.
- Run tests on known systems and compare the computed solutions with expected results.
- Verify input validation and error handling.
- Confirm compatibility with QR decomposition and other dependent modules.
- Perform a clean build using the project's build instructions.
- Review and merge the module's pull request once the code has been tested.
