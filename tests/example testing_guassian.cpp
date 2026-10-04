#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace std;

// Solves Ax = b using Gaussian Elimination with partial pivoting
vector<double> gaussianElimination(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    // Check if the system is valid
    if (A[0].size() != n || b.size() != n) {
        throw invalid_argument("Matrix dimensions must be square and match the vector b.");
    }

    // 1. Forward Elimination Phase
    for (int k = 0; k < n - 1; k++) {
        
        // Partial Pivoting: find the best row below row k
        int maxRow = k;
        for (int i = k + 1; i < n; i++) {
            if (abs(A[i][k]) > abs(A[maxRow][k])) {
                maxRow = i;
            }
        }

        // Check if the best pivot is basically zero
        if (abs(A[maxRow][k]) < 1e-12) {
            throw runtime_error("Matrix is singular (no unique solution exists).");
        }

        // Swap the current row with the max row if needed
        if (maxRow != k) {
            swap(A[k], A[maxRow]);
            swap(b[k], b[maxRow]);
        }

        // Eliminate entries underneath the pivot element
        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            
            // subtract row k from row i
            for (int j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Check the final element on the main diagonal
    if (abs(A[n - 1][n - 1]) < 1e-12) {
        throw runtime_error("Matrix is singular (no unique solution exists).");
    }

    // 2. Back Substitution Phase
    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        double currentSum = b[i];
        
        for (int j = i + 1; j < n; j++) {
            currentSum -= A[i][j] * x[j];
        }
        
        x[i] = currentSum / A[i][i];
    }

    return x;
}

int main() {
    cout << "--- GAUSSIAN ELIMINATION TESTS ---\n" << endl;

    // TEST 1: Basic 2x2 Matrix system
    vector<vector<double>> A1 = {{2, 1}, {1, 3}};
    vector<double> b1 = {3, 5};
    
    try {
        cout << "[Test 1] Running 2x2 System..." << endl;
        vector<double> x1 = gaussianElimination(A1, b1);
        cout << "Solution x: ";
        for (double val : x1) cout << val << " ";
        cout << "\n(Expected roughly: 0.8 1.4)\n" << endl;
    } catch (exception& e) {
        cout << "Error: " << e.what() << endl;
    }

    //  TEST 2: Standard 3x3  MatrixSystem 
    vector<vector<double>> A2 = {
        {3, 2, -1},
        {2, -2, 4},
        {-1, 0.5, -1}
    };
    vector<double> b2 = {1, -2, 0};

    try {
        cout << "[Test 2] Running 3x3 System..." << endl;
        vector<double> x2 = gaussianElimination(A2, b2);
        cout << "Solution x: ";
        for (double val : x2) cout << val << " ";
        cout << "\n(Expected roughly: 1 -2 -2)\n" << endl;
    } catch (exception& e) {
        cout << "Error: " << e.what() << endl;
    }

    // TEST 3: System requiring pivoting
    // Top-left is 0, so it must pivot or it fails by dividing by zero
    vector<vector<double>> A3 = {{0, 2}, {1, 3}};
    vector<double> b3 = {4, 7};

    try {
        cout << "[Test 3] Running Partial Pivoting Check..." << endl;
        vector<double> x3 = gaussianElimination(A3, b3);
        cout << "Solution x: ";
        for (double val : x3) cout << val << " ";
        cout << "\n(Expected roughly: 1 2)\n" << endl;
    } catch (exception& e) {
        cout << "Error: " << e.what() << endl;
    }

    //  TEST 4: Singular Matrix Test 
    vector<vector<double>> A4 = {{1, 2}, {2, 4}};
    vector<double> b4 = {3, 6};

    try {
        cout << "[Test 4] Running Singular Matrix..." << endl;
        vector<double> x4 = gaussianElimination(A4, b4);
        cout << "Oops, failed to catch singular matrix!" << endl;
    } catch (exception& e) {
        cout << "Caught expected error: " << e.what() << endl;
        cout << "Test passed." << endl;
    }

    cout << "\nAll tests finished running." << endl;
    return 0;
}