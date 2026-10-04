# Week 1 Progress Report

## Completed
- Set up the GitHub repository (GROUP-24) with the recommended project structure
- Added branch protection on `main` (pull request and 1 approving review required)
- Implemented the Gram-Schmidt module (`include/numerical-computing/gram_schmidt.hpp`, `src/gram_schmidt.cpp`) in the `numcomp` namespace: `dot`, `norm`, and `gramSchmidt`
- Added input validation: empty list, empty vectors, mismatched sizes, more vectors than dimensions, linearly dependent vectors
- Wrote tests (`tests/test_gram_schmidt.cpp`) covering normal use, an edge case (nearly dependent vectors), and invalid inputs
- Wrote a usage example (`examples/gram_schmidt_example.cpp`)
- Created `CMakeLists.txt` to build the library, tests, and example; tests run with `ctest`

## How the Gram-Schmidt code was built
1. **Defined the problem.** Gram-Schmidt takes a list of linearly independent vectors and returns a list of vectors of length 1 that are all perpendicular to each other and span the same space.
2. **Designed the interface first.** Put declarations in the header (`Vector` alias for `std::vector<double>`, `dot`, `norm`, `gramSchmidt`) and kept the implementation in `src/`, so the public API is separate from the logic.
3. **Built the helpers.** `dot(a, b)` multiplies matching entries and sums them. `norm(v)` is the square root of `dot(v, v)`. These are reused inside `gramSchmidt` to avoid duplicated code.
4. **Validated input before any maths.** Reject an empty list, empty vectors, vectors of different sizes, and more vectors than dimensions (they must be dependent), using `std::invalid_argument`.
5. **Built the orthonormal list one vector at a time.** For each input vector, copy it, subtract its projection onto every finished vector, then divide by what is left to get length 1.
6. **Used the modified variant.** Each projection uses the already-updated vector instead of the original, which is more numerically stable than classical Gram-Schmidt.
7. **Detected dependent vectors.** If the leftover length is smaller than `tol` times the original length, the vector adds nothing new and an exception is thrown.
8. **Tested as we went.** Checked results by confirming every vector has length 1 and every pair has dot product 0, then tested the invalid inputs.

## Issues faced
- **Floating-point comparison.** Decimal results are never exactly 0 or 1, so the tests compare with a small tolerance instead of `==`.
- **Deciding when a vector is "dependent".** The leftover after subtraction is rarely exactly 0 for dependent input, so a relative tolerance (`tol`) was needed instead of checking for zero.
- **Loss of accuracy with nearly parallel vectors.** Classical Gram-Schmidt loses orthogonality here, so we switched to the modified version and added a test for it.
- **More vectors than dimensions.** Such input can never be independent, so it is rejected up front with a clear error message.
- **Git/GitHub workflow.** Getting used to pull requests and the review requirement from the branch protection rules, and creating folders such as `reports/week-01/` through the GitHub web editor.

## In Progress
- README.md (hand-typed)
- Reviewing each other's code through pull requests

## Challenges/Blockers
- Learning the pull request workflow under the branch protection rules
- Making sure every member has visible commits

## Next Week
- Add the next library module(s) with matching tests and examples
- Verify a clean build from scratch following the README steps
- Merge pull requests so every member has commits in the history

## AI Use
- Tool: Claude
- Purpose: Help with GitHub steps (creating folders, `.gitkeep`) and drafting this report
- Reason: To save time on repository setup and report formatting
