Progress Report: Polynomial Fitting
by Muwonge Emmanuel
Module / Tier: Polynomial Fitting (Tier 4)

Files:polynomial_fitting.hpp, src/polynomial_fitting.cpp, polynomial_fitting_test.cpp, polynomial_fitting_example.cpp"

Completed

* Wrote the Tier 4 least-squares polynomial fitting module in C++.
* Reviewed Fortunate's Tier 3 Least-Squares code to make sure I understood exactly how to format the data for her solver.
* Wrote the logic to take standard x coordinates and turn them into a Vandermonde matrix (where rows are [1, x, x^2, ...]).
* Hooked up my matrix generator directly to Fortunate's Tier 3 solver so it can find the best-fit polynomial coefficients.
* Wrote up this progress report to document how it all fits together.

How the Polynomial Fitting Code Was Built

1. Defined the problem. I needed to fit a polynomial curve of a specific degree to scattered, noisy data points using least-squares.
2. Checked the existing code. I looked at the Tier 3 module so my code could just build the right matrix and pass the actual heavy math down the chain.
3. Built the matrix. I wrote a loop that takes the x values and generates columns for increasing powers of x.
4. Calculated the fit. I fed that newly built matrix and the y values into Fortunate's solver to get the final coefficients.
5. Tested it. I ran the exact example from the project plan—passing points (0,1), (1,3), (2,5), (3,8)—and confirmed it returned the correct straight line of y = 0.8 + 2.3x.
6. Reviewed the work. Checked the repo to make sure everything integrates properly.

Program Organization (Classes)

1. What classes did you create?
I didn't create any new classes for this. It is implemented as a free function.

2. Why did you organize it this way?
Since the project already uses shared vector and matrix types, creating a whole new class just for curve fitting was unnecessary. A simple function keeps things lightweight and feeds data directly into the Tier 3 solver.

3. What data does it hold?

| Item | Type | Purpose |
| --- | --- | --- |
| x_vals | std::vector | The input x coordinates. |
| y_vals | std::vector | The input y coordinates. |
| Design Matrix | std::vector[numcomp::Vector] | The matrix we build (powers of x) to feed the Tier 3 solver. |
| Coefficients | numcomp::Vector | The final polynomial equation numbers. |

4. What functions operate on that data?

| Function | Purpose |
| --- | --- |
| Polynomial Fitting | Builds the matrix and returns the final equation. |
| Least-Squares | Fortunate's Tier 3 solver that does the heavy lifting. |
| Tests | Checks if the curve fits the data correctly. |

5. Why are some members public and others hidden?
The main polynomial fitting function is declared public in the header so teammates can call it. The matrix construction loop is kept in the source file to keep the interface clean.

6. How does the rest of the program use it?
You include the header, pass in your x and y coordinate vectors, and specify the desired polynomial degree. The function handles the rest.

How This Module Connects to Others

* Uses: Fortunate's Tier 3 Least-Squares module.
* Used by: Anyone who needs to fit a curve to noisy data.
* Dependency: This module is **not standalone**. It sits at the top of the orthogonal chain (Tier 4) and depends directly on Fortunate's Tier 3 Least-Squares solver, which in turn depends on Musiime's Tier 2 QR decomposition and Ntambi's Tier 1 Gram-Schmidt module.

Issues Faced

* Syncing inputs: I had to make sure my matrix column formatting exactly matched what the Tier 3 solver expected before I could test anything.
*Keeping track of updates to the least_squares.hpp header online so my local tests wouldn't break