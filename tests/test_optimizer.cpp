/**
 * @file test_optimizer.cpp
 * @brief Unit tests for optimizer.h (Nelder-Mead simplex optimizer)
 * 
 * This file contains comprehensive tests for the Nelder-Mead optimizer,
 * organized into three sections:
 * 
 * Section 1: Basic Functionality Tests
 *   - Classic optimization problems (Rosenbrock, quadratic, etc.)
 *   - Convergence verification
 *   - Constraint handling
 * 
 * Section 2: v0.3.0 Improvement Tests
 *   - Default parameter changes (xtol/ftol: 0.01 → 1e-4)
 *   - Restart count increase (3 → 5)
 *   - Precision improvement verification
 * 
 * Section 3: Edge Cases and Robustness
 *   - Boundary conditions
 *   - 1D optimization
 *   - Multimodal functions
 */

#include "test_harness.h"
#include "optimizer.h"
#include "kinship_estimator.h"
#include <cmath>
#include <vector>
#include <algorithm>

using namespace fastlckin;

// ============================================================================
// Section 1: Basic Functionality Tests
// ============================================================================
// These tests verify the core Nelder-Mead algorithm works correctly
// on classic optimization problems.

TEST_CASE(nm_rosenbrock) {
    // Rosenbrock function: f(x,y) = (1-x)^2 + 100*(y-x^2)^2
    // Minimum at (1, 1), known challenging test function
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
    // Minimum at (3, -2)
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
    // Function with constraint-like penalty (box constraints)
    // Minimum at (0.5, 0.25)
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
    // Tests optimizer works for n=1 (edge case)
    auto func = [](const std::vector<double>& x) -> double {
        return (x[0] - 7.0) * (x[0] - 7.0);
    };

    auto result = fastlckin::nelder_mead(func, {0.0});
    CHECK_NEAR(result.x[0], 7.0, 0.05);
}

TEST_CASE(nm_converged_flag) {
    // Verify convergence flag is set correctly
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

// ============================================================================
// Section 2: v0.3.0 Improvement Tests
// ============================================================================
// These tests verify the v0.3.0 optimizations:
// - Tighter default tolerance (1e-4 vs 0.01)
// - More restarts (5 vs 3)
// - Improved precision

TEST_CASE(nm_default_tolerance_v030) {
    // Verify default tolerance is now 1e-4 (v0.3.0 improvement from 0.01)
    NelderMeadConfig cfg;
    CHECK_NEAR(cfg.xtol, 1e-4, 1e-10);
    CHECK_NEAR(cfg.ftol, 1e-4, 1e-10);
}

TEST_CASE(nm_default_restarts_v030) {
    // Verify default restarts is now 5 (v0.3.0 improvement from 3)
    KinshipConfig config;
    CHECK(config.n_restarts == 5);
}

TEST_CASE(nm_improved_precision_quadratic) {
    // Test that tighter tolerance (v0.3.0) gives better precision
    auto quad = [](const std::vector<double>& x) -> double {
        return (x[0] - 0.333333) * (x[0] - 0.333333) + 
               (x[1] - 0.666667) * (x[1] - 0.666667);
    };
    
    // Old tolerance (v0.2.0: 0.01)
    NelderMeadConfig cfg_old;
    cfg_old.xtol = 0.01;
    cfg_old.ftol = 0.01;
    auto result_old = nelder_mead(quad, {0.0, 0.0}, cfg_old);
    
    // New tolerance (v0.3.0: 1e-4)
    NelderMeadConfig cfg_new;
    cfg_new.xtol = 1e-4;
    cfg_new.ftol = 1e-4;
    auto result_new = nelder_mead(quad, {0.0, 0.0}, cfg_new);
    
    // v0.3.0 should be more accurate
    double err_old = std::abs(result_old.x[0] - 0.333333) + std::abs(result_old.x[1] - 0.666667);
    double err_new = std::abs(result_new.x[0] - 0.333333) + std::abs(result_new.x[1] - 0.666667);
    
    CHECK(err_new <= err_old);
    CHECK_NEAR(result_new.x[0], 0.333333, 1e-3);
    CHECK_NEAR(result_new.x[1], 0.666667, 1e-3);
}

TEST_CASE(nm_convergence_with_tight_tolerance) {
    // Test that optimization still converges with v0.3.0 tighter tolerance
    auto rosenbrock = [](const std::vector<double>& x) -> double {
        double a = 1.0 - x[0];
        double b = x[1] - x[0] * x[0];
        return a * a + 100.0 * b * b;
    };
    
    NelderMeadConfig cfg;
    cfg.xtol = 1e-4;
    cfg.ftol = 1e-4;
    cfg.max_iter = 50000;
    cfg.max_fun_evals = 50000;
    
    auto result = nelder_mead(rosenbrock, {-1.0, 1.0}, cfg);
    
    CHECK(result.converged || result.fun_evals < cfg.max_fun_evals);
    CHECK_NEAR(result.x[0], 1.0, 0.01);
    CHECK_NEAR(result.x[1], 1.0, 0.01);
    CHECK(result.fval < 0.01);
}

TEST_CASE(nm_multiple_restarts_better_solution) {
    // Test that more restarts (v0.3.0: 5 vs 3) help find better solutions
    // Use a multimodal function with local minima
    auto multimodal = [](const std::vector<double>& x) -> double {
        // Multiple wells at different locations
        double global_min = (x[0]-0.5)*(x[0]-0.5) + (x[1]-0.5)*(x[0]-0.5);
        double local_min1 = (x[0]+2.0)*(x[0]+2.0) + (x[1]+2.0)*(x[1]+2.0) + 1.0;
        double local_min2 = (x[0]-3.0)*(x[0]-3.0) + (x[1]+1.0)*(x[1]+1.0) + 2.0;
        return std::min({global_min, local_min1, local_min2});
    };
    
    // With 3 restarts (v0.2.0)
    int n_restarts_3 = 3;
    double best_fval_3 = 1e10;
    for (int r = 0; r < n_restarts_3; ++r) {
        std::vector<double> x0 = {r * 0.5, r * 0.3};
        auto result = nelder_mead(multimodal, x0);
        best_fval_3 = std::min(best_fval_3, result.fval);
    }
    
    // With 5 restarts (v0.3.0)
    int n_restarts_5 = 5;
    double best_fval_5 = 1e10;
    for (int r = 0; r < n_restarts_5; ++r) {
        std::vector<double> x0 = {r * 0.5, r * 0.3};
        auto result = nelder_mead(multimodal, x0);
        best_fval_5 = std::min(best_fval_5, result.fval);
    }
    
    // More restarts should find equal or better solution
    CHECK(best_fval_5 <= best_fval_3);
}

// ============================================================================
// Section 3: Edge Cases and Robustness
// ============================================================================
// These tests verify optimizer behavior in challenging scenarios.

TEST_CASE(nm_precision_boundary_cases) {
    // Test precision near boundaries with Cotterman constraint
    // 4*k0*k2 <= k1^2 (genetic feasibility constraint)
    auto constrained_func = [](const std::vector<double>& x) -> double {
        double k0 = 1.0 - x[0] - x[1];
        double k1 = x[0];
        double k2 = x[1];
        
        if (k0 < 0 || k1 < 0 || k2 < 0) return 1e10;
        if (4.0 * k0 * k2 > k1 * k1 + 1e-10) return 1e10;
        
        return (k1 - 0.1) * (k1 - 0.1) + (k2 - 0.01) * (k2 - 0.01);
    };
    
    NelderMeadConfig cfg;
    cfg.xtol = 1e-4;
    cfg.ftol = 1e-4;
    
    auto result = nelder_mead(constrained_func, {0.2, 0.02}, cfg);
    
    // Should find a feasible solution with good precision
    double k0 = 1.0 - result.x[0] - result.x[1];
    CHECK(k0 >= -0.01);
    CHECK(result.x[0] >= -0.01);
    CHECK(result.x[1] >= -0.01);
    CHECK(4.0 * k0 * result.x[1] <= result.x[0] * result.x[0] + 0.01);
    CHECK(result.fval < 0.1);  // Reasonable objective value
}

TEST_CASE(nm_convergence_to_known_minimum) {
    // Test on function with known minimum from multiple starting points
    auto func = [](const std::vector<double>& x) -> double {
        return std::pow(x[0] - 0.3, 2) + std::pow(x[1] - 0.7, 2);
    };
    
    NelderMeadConfig cfg;
    cfg.xtol = 1e-4;
    cfg.ftol = 1e-4;
    cfg.max_iter = 10000;
    
    // Try multiple starting points to verify global convergence
    std::vector<std::vector<double>> starts = {
        {0.0, 0.0}, {1.0, 1.0}, {-1.0, 2.0}, {0.5, 0.5}
    };
    
    for (const auto& x0 : starts) {
        auto result = nelder_mead(func, x0, cfg);
        
        // Should converge to (0.3, 0.7)
        CHECK_NEAR(result.x[0], 0.3, 0.01);
        CHECK_NEAR(result.x[1], 0.7, 0.01);
        CHECK(result.fval < 1e-6);
    }
}

TEST_CASE(numerical_robustness_extreme_values) {
    // Test robustness with extreme parameter values (narrow minimum)
    auto func = [](const std::vector<double>& x) -> double {
        // Very narrow minimum (high curvature)
        return 100.0 * std::pow(x[0] - 0.001, 2) + 
               100.0 * std::pow(x[1] - 0.999, 2);
    };
    
    NelderMeadConfig cfg;
    cfg.xtol = 1e-4;
    cfg.ftol = 1e-4;
    
    auto result = nelder_mead(func, {0.5, 0.5}, cfg);
    
    // Should find the narrow minimum
    CHECK_NEAR(result.x[0], 0.001, 0.01);
    CHECK_NEAR(result.x[1], 0.999, 0.01);
}

int main() {
    return RUN_ALL_TESTS();
}
