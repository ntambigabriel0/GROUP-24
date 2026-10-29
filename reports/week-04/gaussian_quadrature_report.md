# Week 4 Progress Report: Gaussian Quadrature

**Module:** Gaussian Quadrature (Gauss-Legendre)
**Files:** `include/numerical-computing/gaussian_quadrature.hpp`, `src/gaussian_quadrature.cpp`, `tests/test_gaussian_quadrature.cpp`, `examples/gaussian_quadrature_example.cpp`

## Completed
- Implemented the Gaussian Quadrature module in the `numcomp` namespace: `gaussLegendreRule` (builds the nodes and weights for n points), `gaussLegendre` (integral of a function over `[a, b]` with one rule), and `compositeGaussLegendre` (splits `[a, b]` into equal panels and applies a rule on each one)
- Added input validation: number of points outside 1 to 100, an empty function, NaN or infinite limits, and fewer than 1 panel
- Handled the interval cases: reversed limits (`a > b`) change the sign of the result, and `a == b` returns 0
- Wrote 35 checks covering known rules, rule structure for every n from 1 to 100, exactness, standard integrals, interval handling, accuracy, the composite rule, and invalid inputs. All pass with no compiler warnings (`-Wall -Wextra -Wpedantic`)
- Wrote a usage example with nine demonstrations, including two applied problems (a normal-distribution probability and the work done stretching a spring)

## How the Gaussian Quadrature Code Was Built
1. **Defined the problem.** Approximate the integral of f over `[a, b]` as a weighted sum of f at a few chosen points. Gauss-Legendre picks the points (nodes) and weights so that an n-point rule integrates every polynomial of degree 2n - 1 or lower exactly. That is about twice what n equally spaced points can reach.
2. **Designed the interface first.** Put the declarations in the header and kept the implementation in `src/`. The header holds one small struct, `QuadratureRule`, and three functions. `n` defaults to 5 for the single rule and 3 per panel for the composite rule.
3. **Validated input before any maths.** Reject bad point counts, an empty function, non-finite limits, and zero panels, using `std::invalid_argument`.
4. **Found the nodes.** The nodes on `[-1, 1]` are the roots of the Legendre polynomial `P_n`. `P_n(x)` is computed with the three-term recurrence `(k+1) P_{k+1} = (2k+1) x P_k - k P_{k-1}`. Each root is found with Newton's method from the starting guess `cos(pi (i + 0.75) / (n + 0.5))`, stopping when the step is below `1e-15`.
5. **Computed the weights.** Each weight is `2 / ((1 - x^2) P_n'(x)^2)`, using the derivative at the final root. Roots come in plus and minus pairs, so only half are computed and the other half are mirrored. For odd n the middle node is set to exactly 0.
6. **Mapped the rule onto `[a, b]`.** Each node `t` on `[-1, 1]` becomes `mid + halfLength * t`, and the weighted sum is multiplied by `halfLength`. A reversed interval gives a negative `halfLength`, which is why the sign flips with no extra code.
7. **Added the composite version.** It splits `[a, b]` into equal panels, builds the rule once, and reuses it on every panel, so the roots are not recomputed for each panel.
8. **Tested as we went.** Checked the 1, 2 and 3 point rules against values known by hand (for example weights 5/9, 8/9, 5/9), confirmed the weights are positive, add up to 2, and the nodes are symmetric, ordered and inside `(-1, 1)` for every n from 1 to 100, confirmed `x^(2n-1)` is integrated exactly for n up to 12 while one degree higher is not, then compared against integrals with known answers (`sin`, `e^x`, `1/x`), and finally tested the invalid inputs.

Results from the example program (exact answers in brackets):

| Integral | Points | Error |
|----------|--------|-------|
| sin(x) on [0, pi] (exact 2) | 2 / 3 / 5 / 8 | 6.4e-2 / 1.4e-3 / 1.1e-7 / 4.2e-15 |
| e^x on [0, 1] (exact e - 1) | 1 / 3 / 5 | 7.0e-2 / 8.2e-7 / 6.5e-13 |
| Normal probability P(-1 <= Z <= 1) (exact 0.682689492137) | 10 | matches to 12 digits |
| Spring work, F(x) = 200x on [0, 0.3] (exact 9 J) | 2 | 9 J, exact |

## Program Organization (Classes)

### 1. What classes did you create?
No classes with behavior. The module has one plain struct, `QuadratureRule`, and the rest is free functions inside the `numcomp` namespace.

### 2. Why did you organize it this way?
Integration is a calculation: a function and limits go in and a number comes out. Nothing has to be remembered between calls, so there is no state for a class to manage. The one struct exists only because a rule is two lists that belong together (nodes and weights), and returning them as a pair is clearer than returning two separate vectors. The namespace stops our names from clashing with other code.

### 3. What data does it hold?
`QuadratureRule` holds two `std::vector<double>` fields: `nodes` (the sample points, smallest first) and `weights` (the matching weights, which add up to 2). Everything else lives only inside the functions while they run. The function to integrate is passed in as a `const std::function<double(double)>&`, so lambdas and ordinary functions both work.

### 4. What functions operate on that data?
| Function | Public / Private | What it does |
|----------|------------------|--------------|
| `gaussLegendreRule(n)` | Public | Builds the n-point rule on `[-1, 1]` |
| `gaussLegendre(f, a, b, n)` | Public | Integral of `f` over `[a, b]` with one n-point rule |
| `compositeGaussLegendre(f, a, b, panels, n)` | Public | Integral of `f` over `[a, b]` using equal panels, one rule per panel |
| `checkPointCount(n)` | Private (anonymous namespace in `gaussian_quadrature.cpp`) | Throws `std::invalid_argument` if n is not between 1 and 100 |
| `checkInputs(f, a, b)` | Private (anonymous namespace) | Throws `std::invalid_argument` for an empty function or non-finite limits |
| `legendre(n, x, pn, pnMinus1)` | Private (anonymous namespace) | Evaluates `P_n(x)` and `P_{n-1}(x)` with the recurrence |

### 5. Why are some members public and others hidden?
The three integration functions and `QuadratureRule` are declared in the header, so other modules can call them. The checks and the Legendre evaluation sit in an anonymous namespace in `gaussian_quadrature.cpp`, so they are visible only inside that file. Users only see the function signatures and their documented behavior, so the root-finding can change later without breaking their code.

## How This Module Connects to Others
- **Uses:** nothing from the library; only the C++ standard library.
- **Used by:** nothing yet. Romberg Integration solves the same kind of problem (a definite integral of a function over `[a, b]`), so the two can be compared on the same integrals as a cross-check once the Romberg module builds. Gauss-Legendre is also closely related to Lagrange Interpolation: the weights are what you get by integrating the Lagrange polynomial through the nodes, which is why the rule is exact up to degree 2n - 1.

## Issues Faced
- **A test threshold that was too strict.** The composite-rule accuracy test originally required an error below `1e-9` using 20 panels of 3 points each on `sin(x)` over `[0, 10]`. The real error is about `1.4e-8`, so one test failed (34 of 35 passed). The code was correct: the error shrinks by about 64 times each time the panels double, which is what a 3-point rule should do. The threshold was loosened to `1e-7`, and all 35 now pass.
- **Newton's method needs a good starting guess.** A poor guess can land on the wrong root or the same root twice. The cosine starting guess fixes this, and the structure tests (ordered, symmetric, inside `(-1, 1)`) confirm it for every n from 1 to 100.
- **Floating-point comparison.** Computed nodes and weights are rarely exactly equal to the known values, so the tests compare with a small tolerance instead of `==`. The middle node of an odd rule is set to exactly 0 so symmetry holds exactly.
- **NaN and infinite limits.** These would silently spread through the sum and give a meaningless answer, so they are rejected at the start.
- **One rule is not enough for wiggly functions.** A single 5-point rule on `sin(x)` over `[0, 20]` has an error of about 4.5, while 10 panels of 5 points each give about `2.8e-10`. This is why the composite function exists.
- **Known limitations.** The rule works best when the function is smooth. A corner or jump inside a panel reduces the accuracy, although putting a panel edge on the corner fixes it (tested with `|x|`). Functions that blow up at an endpoint are not handled. `gaussLegendre` rebuilds the rule on every call, which is fine for one integral but wasteful in a tight loop; the composite version avoids this by building it once. The number of points is capped at 100.

## In Progress
- Adding the Gaussian Quadrature files to the repository and opening a pull request
- Registering the module in `CMakeLists.txt` and updating the README (the module status table, the "What the library does" table, and the note that Gaussian quadrature has no source file yet)
- Reviewing each other's code through pull requests

## Challenges / Blockers
- The Romberg cross-check is blocked until the Romberg module builds: its source, example and test include `romberg.hpp`, but the header is named `romberg_integration.hpp`
- `CMakeLists.txt` currently registers only Gram-Schmidt, so each module's files still have to be added by hand
- Making sure every member has visible commits
