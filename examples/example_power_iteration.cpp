#ifndef POWER_ITERATION_EXAMPLE_H
#define POWER_ITERATION_EXAMPLE_H

#include <string>
using namespace std;

struct PowerIterationTest {
    string name;
    int n;
    double A[2][2];
    double initial[2];
    double expectedEigenvalue;
    double tolerance;
};

const PowerIterationTest example = {
    "Power Iteration Example 1",
    2,
    {{2, 1}, {1, 2}},
    {1, 0},
    3.0,
    0.001
};

#endif
