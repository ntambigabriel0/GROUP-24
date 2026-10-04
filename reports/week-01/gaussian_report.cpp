/*
GAUSSIAN ELIMINATION ALGORITHM REPORT

Algorithm: Gaussian Elimination with Partial Pivoting


OVERVIEW:
    Matrix A
    Vector b


OUTPUT:
    Solution vector x


STEPS:


1. Check that A is a square matrix.
   Check that the size of b matches A.

2. FOR each column k from 0 to n-2:

   a. Find the row containing the largest
      absolute value in column k.

   b. Set this row as the pivot row.

   c. If the pivot is approximately zero:
          Report that the matrix is singular
          or has no unique solution.

   d. If the pivot row is different from k:
          Swap the two rows of A.
          Swap the corresponding values of b.

   e. FOR each row i below the pivot:

          Calculate:

              factor = A[i][k] / A[k][k]

          FOR each column j:

              A[i][j] = A[i][j]
                        - factor * A[k][j]

          Update:

              b[i] = b[i] - factor * b[k]


3. Check the final pivot.

   If the final pivot is approximately zero:

       Report that the matrix is singular
       or has no unique solution.


4. BACK SUBSTITUTION:

   FOR i from n-1 down to 0:

       sum = b[i]

       FOR j from i+1 to n-1:

           sum = sum - A[i][j] * x[j]

       x[i] = sum / A[i][i]


5. Return x.

END

 DESCRIPTION
Gaussian Elimination is a numerical method used to
solve a system of simultaneous linear equations.

The algorithm changes the original system into an
upper triangular form using forward elimination.

After the elimination process, back substitution is
used to calculate the unknown variables.

The implementation uses partial pivoting. Partial
pivoting selects the largest available pivot in the
current column and swaps rows when necessary.

The algorithm also checks for invalid matrix dimensions
and singular matrices.
STEPS
1. Matrix validation
2. Pivot selection
3. Partial pivoting
4. Row swapping
5. Forward elimination
6. Upper triangular matrix formation
7. Back substitution
8. Solution output
EXAMPLE
Consider the system:

        2x + y = 3

        x + 3y = 5


Matrix A:

        [ 2   1 ]
        [ 1   3 ]

Vector b:

        [ 3 ]
        [ 5 ]

After applying Gaussian Elimination and
back substitution:
        x = 0.8
        y = 1.4
Therefore:
 Solution = (0.8, 1.4)

 TESTING

The implementation is tested using different systems.

Test 1:
    A normal 2 x 2 system.

Test 2:
    A 3 x 3 system.

Test 3:
    A system requiring partial pivoting and
    row swapping.

Test 4:
    A singular matrix to test error detection.


Expected testing result:

    Test 1 -> PASS
    Test 2 -> PASS
    Test 3 -> PASS
    Test 4 -> PASS
CONCLUSION

The Gaussian Elimination algorithm provides a direct
method for solving systems of linear equations.

The implementation uses partial pivoting, forward
elimination and back substitution.

The testing program verifies the correctness of the
algorithm using different types of systems, including
a case requiring row swapping and a singular matrix.

*/