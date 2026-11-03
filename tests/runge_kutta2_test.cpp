#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>
#include "runge_kutta2.hpp"

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { \
    std::cerr << "FAIL line " << __LINE__ << ": " #cond "\n"; ++failures; } } while (0)

template <typename F>
static bool throwsInvalidArgument(F fn) {
    try { fn(); } catch (const std::invalid_argument&) { return true; } catch (...) { return false; }
    return false;
}

int main() {
    auto f = [](double, double y) { return y; };
    auto one = [](double, double) { return 1.0; };

    CHECK(std::abs(rungeKutta2(f, 1.0, 0.0, 1.0, 0.1) - std::exp(1.0)) < 1e-2);
    CHECK(std::abs(rungeKutta2(f, 1.0, 0.0, 1.0, 0.001) - std::exp(1.0)) < 1e-6);

    double e1 = std::abs(rungeKutta2(f, 1.0, 0.0, 1.0, 0.1)  - std::exp(1.0));
    double e2 = std::abs(rungeKutta2(f, 1.0, 0.0, 1.0, 0.05) - std::exp(1.0));
    CHECK(e1 / e2 > 3.5 && e1 / e2 < 4.5);

    CHECK(std::abs(rk2_step(f, 0.0, 1.0, 0.1) - 1.105) < 1e-12);

    CHECK(std::abs(rungeKutta2(one, 0.0, 0.0, 1.1, 0.1) - 1.1) < 1e-12);  
    CHECK(std::abs(rungeKutta2(one, 0.0, 0.0, 0.7, 0.1) - 0.7) < 1e-12);
    CHECK(std::abs(rungeKutta2(one, 0.0, 0.0, 1.0, 0.3) - 1.0) < 1e-12);   

    
    CHECK(rungeKutta2(f, 2.0, 1.0, 1.0, 0.1) == 2.0);
    CHECK(std::abs(rungeKutta2(one, 0.0, 0.0, 0.5, 2.0) - 0.5) < 1e-12);

    {
        auto path = rk2_trajectory(f, 1.0, 0.0, 1.0, 0.1);
        CHECK(path.size() == 11);
        CHECK(path.front().first == 0.0 && path.front().second == 1.0);
        CHECK(path.back().first == 1.0);
        CHECK(path.back().second == rungeKutta2(f, 1.0, 0.0, 1.0, 0.1));

        path = rk2_trajectory(f, 1.0, 0.0, 1.1, 0.1);
        CHECK(path.size() == 12 && path.back().first == 1.1);   
        path = rk2_trajectory(f, 1.0, 0.0, 1.0, 0.3);
        CHECK(path.size() == 5 && path.back().first == 1.0);    }

    CHECK(throwsInvalidArgument([&] { rungeKutta2(f, 1.0, 0.0, 1.0, 0.0); }));
    CHECK(throwsInvalidArgument([&] { rungeKutta2(f, 1.0, 0.0, 1.0, -0.1); }));
    CHECK(throwsInvalidArgument([&] { rungeKutta2(f, 1.0, 0.0, 1.0, std::nan("")); }));
    CHECK(throwsInvalidArgument([&] { rungeKutta2(f, 1.0, 1.0, 0.0, 0.1); }));
    CHECK(throwsInvalidArgument([&] { rk2_trajectory(f, 1.0, 0.0, 1.0, 0.0); }));

    if (failures == 0) std::cout << "All tests passed\n";
    return failures == 0 ? 0 : 1;
}
