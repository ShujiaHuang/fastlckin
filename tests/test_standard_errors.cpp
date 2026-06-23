/**
 * @file test_standard_errors.cpp
 * @brief Unit tests for Fisher information matrix and standard error computation
 */

#include "test_harness.h"
#include "optimizer.h"
#include <cmath>
#include <vector>

using namespace fastlckin;

TEST_CASE(standard_errors_quadratic_known) {
    // f(x,y) = (x-3)^2 + (y+2)^2
    // Hessian = [[2, 0], [0, 2]] (constant)
    // Covariance = H^{-1} = [[0.5, 0], [0, 0.5]]
    // SE = [sqrt(0.5), sqrt(0.5)] = [0.707, 0.707]
    
    auto quad = [](const std::vector<double>& x) -> double {
        double a = x[0] - 3.0;
        double b = x[1] + 2.0;
        return a * a + b * b;
    };
    
    std::vector<double> x_opt = {3.0, -2.0};
    auto se = compute_standard_errors(quad, x_opt);
    
    CHECK(se.size() == 2);
    CHECK_NEAR(se[0], std::sqrt(0.5), 0.01);
    CHECK_NEAR(se[1], std::sqrt(0.5), 0.01);
}

TEST_CASE(standard_errors_asymmetric_quadratic) {
    // f(x,y) = 2*(x-1)^2 + 8*(y-2)^2
    // Hessian = [[4, 0], [0, 16]]
    // Covariance = [[0.25, 0], [0, 0.0625]]
    // SE = [0.5, 0.25]
    
    auto func = [](const std::vector<double>& x) -> double {
        return 2.0 * (x[0] - 1.0) * (x[0] - 1.0) + 
               8.0 * (x[1] - 2.0) * (x[1] - 2.0);
    };
    
    std::vector<double> x_opt = {1.0, 2.0};
    auto se = compute_standard_errors(func, x_opt);
    
    CHECK(se.size() == 2);
    CHECK_NEAR(se[0], 0.5, 0.01);
    CHECK_NEAR(se[1], 0.25, 0.01);
}

TEST_CASE(standard_errors_with_correlation) {
    // f(x,y) = (x-1)^2 + (x-1)*(y-2) + (y-2)^2
    // Hessian = [[2, 1], [1, 2]]
    // det = 4 - 1 = 3
    // Covariance = (1/3) * [[2, -1], [-1, 2]]
    // SE = [sqrt(2/3), sqrt(2/3)] = [0.816, 0.816]
    
    auto func = [](const std::vector<double>& x) -> double {
        double dx = x[0] - 1.0;
        double dy = x[1] - 2.0;
        return dx * dx + dx * dy + dy * dy;
    };
    
    std::vector<double> x_opt = {1.0, 2.0};
    auto se = compute_standard_errors(func, x_opt);
    
    CHECK(se.size() == 2);
    CHECK_NEAR(se[0], std::sqrt(2.0/3.0), 0.01);
    CHECK_NEAR(se[1], std::sqrt(2.0/3.0), 0.01);
}

TEST_CASE(standard_errors_singular_hessian) {
    // f(x,y) = (x+y-3)^2
    // Hessian = [[2, 2], [2, 2]]
    // det = 0 (singular)
    // Should return empty
    
    auto func = [](const std::vector<double>& x) -> double {
        double s = x[0] + x[1] - 3.0;
        return s * s;
    };
    
    std::vector<double> x_opt = {1.5, 1.5};  // x+y=3
    auto se = compute_standard_errors(func, x_opt);
    
    CHECK(se.empty());
}

TEST_CASE(standard_errors_after_optimization) {
    // Test integration: optimize then compute SE
    auto quad = [](const std::vector<double>& x) -> double {
        return (x[0] - 5.0) * (x[0] - 5.0) + 
               4.0 * (x[1] + 3.0) * (x[1] + 3.0);
    };
    
    NelderMeadConfig cfg;
    cfg.xtol = 1e-6;
    cfg.ftol = 1e-8;
    cfg.max_iter = 5000;
    
    auto result = nelder_mead(quad, {0.0, 0.0}, cfg);
    
    CHECK(result.converged);
    CHECK_NEAR(result.x[0], 5.0, 0.01);
    CHECK_NEAR(result.x[1], -3.0, 0.01);
    
    // Compute SE at optimum
    // Hessian = [[2, 0], [0, 8]]
    // Covariance = [[0.5, 0], [0, 0.125]]
    // SE = [0.707, 0.354]
    auto se = compute_standard_errors(quad, result.x);
    
    CHECK(se.size() == 2);
    CHECK_NEAR(se[0], std::sqrt(0.5), 0.02);
    CHECK_NEAR(se[1], std::sqrt(0.125), 0.02);
}

TEST_CASE(standard_errors_1d) {
    // 1D case: f(x) = 3*(x-2)^2
    // Hessian = [6]
    // Covariance = [1/6]
    // SE = sqrt(1/6) = 0.408
    
    auto func = [](const std::vector<double>& x) -> double {
        return 3.0 * (x[0] - 2.0) * (x[0] - 2.0);
    };
    
    std::vector<double> x_opt = {2.0};
    auto se = compute_standard_errors(func, x_opt);
    
    // Currently only supports 2D, so should return empty
    CHECK(se.empty());
}

TEST_CASE(standard_errors_step_size_sensitivity) {
    // Test that different step sizes give similar results
    auto quad = [](const std::vector<double>& x) -> double {
        return (x[0] - 1.0) * (x[0] - 1.0) + 
               (x[1] - 1.0) * (x[1] - 1.0);
    };
    
    std::vector<double> x_opt = {1.0, 1.0};
    
    auto se1 = compute_standard_errors(quad, x_opt, 1e-4);
    auto se2 = compute_standard_errors(quad, x_opt, 1e-5);
    auto se3 = compute_standard_errors(quad, x_opt, 1e-6);
    
    // All should be close to sqrt(0.5) = 0.707
    CHECK(se1.size() == 2);
    CHECK(se2.size() == 2);
    CHECK(se3.size() == 2);
    
    // Check that results are consistent across step sizes
    CHECK_NEAR(se1[0], se2[0], 0.01);
    CHECK_NEAR(se2[0], se3[0], 0.01);
    CHECK_NEAR(se1[1], se2[1], 0.01);
    CHECK_NEAR(se2[1], se3[1], 0.01);
}

int main() {
    return RUN_ALL_TESTS();
}
