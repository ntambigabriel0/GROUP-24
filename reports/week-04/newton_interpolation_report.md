Progress Report: Newton Interpolation
by Muwonge Emmanuel
Module / Tier: Newton Interpolation (Standalone)

Files:newton_interpolation.hpp, newton_interpolation.cpp, newton_interpolation_test.cpp, newton_example.cpp

Completed

* Developed the Newton Interpolation algorithm in C++ using divided differences.


* Created the logic to compute the table of slopes and slope-of-slopes for any set of input data points.


* Built an evaluation function to quickly estimate y values using Horner's method.
* Added test cases confirming our outputs match the (0,1), (1,3), (2,7) example from the project plan.



How the Code Was Built

1. Defined the problem. I needed a standalone interpolator that is cheap to update when new data points show up.


2. Computed divided differences. Wrote a loop to build the triangular slope table and extracted the top edge as coefficients.


3. Designed the evaluator. Used Horner's method to calculate the final polynomial value cleanly without messy power expansions.
4. Tested it. Ran the points (0,1), (1,3), (2,7) through it, getting slopes of 2 and 4, and a slope-of-slopes of 1, which correctly gives 13 at x=3.



Program Organization (Classes)

1. What classes did you create?
I didn't create any new classes. It is written entirely using standard free functions.


2. Why did you organize it this way?
It keeps the module lightweight and matches how Ntambi structured his Lagrange module, accepting standard C++ vectors directly.


3. What data does it hold?

| Item | Type | Purpose |
| x_vals | std::vector | The known x coordinates. |
| y_vals | std::vector | The known y coordinates. |
| Divided Differences | std::vector | The calculated slopes acting as polynomial coefficients|
| Result | double | The final estimated y value. |

4. What functions operate on that data?

| Function | Purpose |
| Coefficient Generator - Builds the divided difference table and returns coefficients. |
| Evaluator - Evaluates the polynomial at a target x value. |
| Tests - Cross-checks results against project plan expectations.



5. Why are some members public and others hidden?
The main coefficient and evaluation functions are public in the header so teammates can use them. The internal table generation logic is hidden away in the source file.


6. How does the rest of the program use it?
You include the header, pass your coordinate vectors and a target x value, and it returns the interpolated result.

How This Module Connects to Others

* Uses: Just the standard C++ library since it's a standalone tool.


* Used by: Anyone needing fast, easily expandable interpolation.


* Dependency: None. However, the plan notes it should be cross-checked against Kirya's Polynomial Interpolation and Ntambi's Lagrange Interpolation to ensure mathematical consistency.