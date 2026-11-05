// Small helpers shared by the tests (include with: #include "test_util.hpp")
#pragma once
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace testutil {

inline int failures = 0;  // counts how many checks failed

// Prints PASS or FAIL for one check and counts failures
inline void check(bool ok, const char* name) {
    std::cout << (ok ? "PASS: " : "FAIL: ") << name << "\n";
    if (!ok) ++failures;
}

// Decimals are compared with a small tolerance, never with ==
inline bool near(double a, double b, double eps = 1e-9) { return std::fabs(a - b) <= eps; }

// True if running f() throws std::invalid_argument
template <typename F>
bool throwsInvalid(F f) {
    try { f(); } catch (const std::invalid_argument&) { return true; } catch (...) {}
    return false;
}

// Call at the end of main(): prints the summary and gives the exit code for ctest
inline int finish() {
    std::cout << (failures == 0 ? "All tests passed\n" : "Some tests failed\n");
    return failures == 0 ? 0 : 1;
}

}  // namespace testutil
