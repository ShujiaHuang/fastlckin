/**
 * @file test_optimizer.cpp
 * @brief Unit tests for optimizer.h (Nelder-Mead)
 */

#include "test_harness.h"
#include "optimizer.h"
#include <cmath>

TEST_CASE(nm_rosenbrock) {
    // Rosenbrock function: f(x,y) = (1-x)^2 + 100*(y-x^2)^2
    // Minimum at (1, 1)
    auto rosenbrock = [](const std::vector<double>& x) -> double {
        double a = 1.0 - x[0];
        double b = x[1] - x[0] * x[0];
        return a * a + 100.0 * b * b;
    };

    fastlckin::NelderMeadConfig cfg;
    cfg.xtol = 1e-6;
    cfg.ftol = 1e-10;
    cfg.max_iter = 50000;
    cfg.max_fun_evals = 50000;

    auto result = fastlckin::nelder_mead(rosenbrock, {-1.0, 1.0}, cfg);
    CHECK_NEAR(result.x[0], 1.0, 0.01);
    CHECK_NEAR(result.x[1], 1.0, 0.01);
    CHECK(result.fval < 1e-4);
}

TEST_CASE(nm_quadratic) {
    // Simple quadratic: f(x,y) = (x-3)^2 + (y+2)^2
    auto quad = [](const std::vector<double>& x) -> double {
        double a = x[0] - 3.0;
        double b = x[1] + 2.0;
        return a * a + b * b;
    };

    auto result = fastlckin::nelder_mead(quad, {0.0, 0.0});
    CHECK_NEAR(result.x[0], 3.0, 0.05);
    CHECK_NEAR(result.x[1], -2.0, 0.05);
    CHECK(result.fval < 0.01);
}

TEST_CASE(nm_constant_penalty) {
    // Function with constraint-like penalty
    auto func = [](const std::vector<double>& x) -> double {
        if (x[0] < 0 || x[1] < 0) return 1e10;
        return (x[0] - 0.5) * (x[0] - 0.5) + (x[1] - 0.25) * (x[1] - 0.25);
    };

    auto result = fastlckin::nelder_mead(func, {0.1, 0.1});
    CHECK_NEAR(result.x[0], 0.5, 0.05);
    CHECK_NEAR(result.x[1], 0.25, 0.05);
}

TEST_CASE(nm_1d) {
    // 1D function: f(x) = (x - 7)^2
    auto func = [](const std::vector<double>& x) -> double {
        return (x[0] - 7.0) * (x[0] - 7.0);
    };

    auto result = fastlckin::nelder_mead(func, {0.0});
    CHECK_NEAR(result.x[0], 7.0, 0.05);
}

TEST_CASE(nm_converged_flag) {
    auto func = [](const std::vector<double>& x) -> double {
        return x[0] * x[0] + x[1] * x[1];
    };

    fastlckin::NelderMeadConfig cfg;
    cfg.xtol = 1e-4;
    cfg.ftol = 1e-8;
    cfg.max_iter = 50000;

    auto result = fastlckin::nelder_mead(func, {5.0, 5.0}, cfg);
    CHECK(result.converged);
    CHECK_NEAR(result.fval, 0.0, 1e-4);
}

int main() {
    return RUN_ALL_TESTS();
}
