#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

#include <numerical-computing/runge_kutta2.hpp>

static int failures = 0;

static void check(bool ok, const char* name) {
    std::cout << (ok ? "  PASS: " : "  FAIL: ") << name << std::endl;
    if (!ok) ++failures;
}

template <typename F>
static bool throwsInvalid(F fn) {
    try {
        fn();
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

int main() {
    std::cout << "--- RUNGE-KUTTA 2 TESTS ---\n" << std::endl;

    auto grow = [](double, double y) { return y; };
    auto decay = [](double, double y) { return -2.0 * y; };

    std::cout << "[Test 1] Exponential growth y' = y" << std::endl;
    {
        double r = rungeKutta2(grow, 1.0, 0.0, 1.0, 0.01);
        check(std::fabs(r - std::exp(1.0)) < 1e-3, "y(1) close to e");
    }

    std::cout << "[Test 2] Exponential decay y' = -2y" << std::endl;
    {
        double r = rungeKutta2(decay, 1.0, 0.0, 1.0, 0.01);
        check(std::fabs(r - std::exp(-2.0)) < 1e-3, "y(1) close to e^-2");
    }

    std::cout << "[Test 3] Exact step count (no floating-point drift)" << std::endl;
    {
        int calls = 0;
        auto counted = [&](double, double y) { ++calls; return y; };
        rungeKutta2(counted, 1.0, 0.0, 1.0, 0.1);
        check(calls == 20, "h=0.1 over [0,1] -> 10 steps (20 evaluations)");

        calls = 0;
        rungeKutta2(counted, 1.0, 0.0, 0.3, 0.1);
        check(calls == 6, "h=0.1 over [0,0.3] -> 3 steps");

        calls = 0;
        rungeKutta2(counted, 1.0, 0.0, 1.0, 0.001);
        check(calls == 2000, "h=0.001 over [0,1] -> 1000 steps");

        calls = 0;
        rungeKutta2(counted, 1.0, 0.0, 1.0, 0.3);
        check(calls == 8, "h=0.3 over [0,1] -> 3 full steps + 1 short step");
    }

    std::cout << "[Test 4] Last step lands exactly on t_end" << std::endl;
    {
        // y' = 1 is integrated exactly by RK2, so y(t_end) = y0 + (t_end - t0).
        auto one = [](double, double) { return 1.0; };
        double r = rungeKutta2(one, 0.0, 0.0, 1.0, 0.3);
        check(std::fabs(r - 1.0) < 1e-12, "y' = 1 gives y(1) = 1 with uneven last step");
    }

    std::cout << "[Test 5] Step larger than interval" << std::endl;
    {
        auto one = [](double, double) { return 1.0; };
        double r = rungeKutta2(one, 0.0, 0.0, 0.5, 2.0);
        check(std::fabs(r - 0.5) < 1e-12, "single short step covers [0, 0.5]");
    }

    std::cout << "[Test 6] Zero-length interval" << std::endl;
    {
        double r = rungeKutta2(grow, 3.5, 2.0, 2.0, 0.1);
        check(r == 3.5, "returns y0 when t_end == t0");
    }

    std::cout << "[Test 7] Non-zero start time" << std::endl;
    {
        // y' = t, y(1) = 0.5  ->  y(t) = t^2/2 ; y(3) = 4.5 (RK2 exact for linear f)
        auto lin = [](double t, double) { return t; };
        double r = rungeKutta2(lin, 0.5, 1.0, 3.0, 0.1);
        check(std::fabs(r - 4.5) < 1e-9, "y' = t from t=1 to t=3");
    }

    std::cout << "[Test 8] Second-order convergence" << std::endl;
    {
        double e1 = std::fabs(rungeKutta2(grow, 1.0, 0.0, 1.0, 0.1) - std::exp(1.0));
        double e2 = std::fabs(rungeKutta2(grow, 1.0, 0.0, 1.0, 0.05) - std::exp(1.0));
        double ratio = e1 / e2;
        check(ratio > 3.6 && ratio < 4.4, "halving h cuts the error by about 4x");
    }

    std::cout << "[Test 9] Invalid input throws invalid_argument" << std::endl;
    {
        check(throwsInvalid([&] { rungeKutta2(grow, 1.0, 0.0, 1.0, 0.0); }), "h = 0");
        check(throwsInvalid([&] { rungeKutta2(grow, 1.0, 0.0, 1.0, -0.1); }), "h < 0");
        check(throwsInvalid([&] {
            rungeKutta2(grow, 1.0, 0.0, 1.0, std::numeric_limits<double>::quiet_NaN());
        }), "h = NaN");
        check(throwsInvalid([&] {
            rungeKutta2(grow, 1.0, 0.0, 1.0, std::numeric_limits<double>::infinity());
        }), "h = inf");
        check(throwsInvalid([&] { rungeKutta2(grow, 1.0, 1.0, 0.0, 0.1); }), "t_end < t0");
        check(throwsInvalid([&] {
            rungeKutta2(grow, 1.0, 0.0, std::numeric_limits<double>::infinity(), 0.1);
        }), "t_end = inf");
        check(throwsInvalid([&] {
            rungeKutta2(std::function<double(double, double)>(), 1.0, 0.0, 1.0, 0.1);
        }), "empty function");
        check(throwsInvalid([&] { rungeKutta2(grow, 1.0, 0.0, 1e300, 1e-300); }),
              "absurd step count");
    }

    std::cout << "\nAll tests finished running. Failures: " << failures << std::endl;
    return failures == 0 ? 0 : 1;
}
