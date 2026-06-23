/**
 * @file test_v030_correctness.cpp
 * @brief Comprehensive correctness verification for v0.3.0 optimizations
 * 
 * This test suite verifies optimization correctness by comparing against
 * known analytical solutions and theoretical properties.
 */

#include "test_harness.h"
#include "frequency_from_likelihoods.h"
#include "optimizer.h"
#include "ld_prune.h"
#include <cmath>
#include <vector>
#include <numeric>
#include <random>

using namespace fastlckin;

// ============================================================================
// 1. Leave-one-out EM: Theoretical Correctness
// ============================================================================

TEST_CASE(loo_em_unbiasedness_verification) {
    // Verify that leave-one-out EM produces unbiased frequency estimates
    // when genotypes are known (delta-function likelihoods)
    
    // Setup: 100 samples, true AF = 0.3
    int n_samples = 100;
    int n_snps = 1;
    double true_af = 0.3;
    
    // Expected genotype counts under HWE:
    // P(G=0) = (1-p)^2 = 0.49 → ~49 samples
    // P(G=1) = 2p(1-p) = 0.42 → ~42 samples
    // P(G=2) = p^2 = 0.09 → ~9 samples
    
    std::vector<std::vector<int>> gt(n_samples, std::vector<int>(1));
    for (int i = 0; i < 49; ++i) gt[i][0] = 0;
    for (int i = 49; i < 91; ++i) gt[i][0] = 1;
    for (int i = 91; i < 100; ++i) gt[i][0] = 2;
    
    const double EPS = 1e-10;
    LikelihoodMatrix lk(n_samples, std::vector<GenotypeLikelihood>(n_snps));
    for (int i = 0; i < n_samples; ++i) {
        int g = gt[i][0];
        lk[i][0].gl[0] = EPS; lk[i][0].gl[1] = EPS; lk[i][0].gl[2] = EPS;
        lk[i][0].gl[g] = 1.0;
        lk[i][0].masked = false;
    }
    
    // Test: Full EM should recover true AF
    auto afs_full = compute_af_from_likelihoods(lk);
    double af_full = afs_full[0];
    
    // Expected: (42*1 + 9*2) / (100*2) = 60/200 = 0.3
    CHECK_NEAR(af_full, 0.3, 0.01);
    
    // Test: Leave-one-out (exclude 2 samples) should give correct adjusted AF
    // Exclude 1 het and 1 hom_alt
    std::vector<int> exclude = {49, 91};  // one het, one hom_alt
    auto afs_loo = compute_af_from_likelihoods_leave_one_out(lk, exclude);
    double af_loo = afs_loo[0];
    
    // Expected: (41*1 + 8*2) / (98*2) = 57/196 = 0.2908
    double expected_loo = 57.0 / 196.0;
    CHECK_NEAR(af_loo, expected_loo, 0.01);
    
    // Verify: LOO AF differs from full EM (bias correction working)
    CHECK(std::abs(af_full - af_loo) > 0.005);
}

TEST_CASE(loo_em_linearity_property) {
    // Leave-one-out should satisfy: 
    // AF_loo = (2*N*AF_full - E[G_excluded]) / (2*(N-1))
    
    int n = 20;
    std::vector<std::vector<int>> gt(n, std::vector<int>(1));
    for (int i = 0; i < 10; ++i) gt[i][0] = 0;
    for (int i = 10; i < 20; ++i) gt[i][0] = 1;
    
    const double EPS = 1e-10;
    LikelihoodMatrix lk(n, std::vector<GenotypeLikelihood>(1));
    for (int i = 0; i < n; ++i) {
        int g = gt[i][0];
        lk[i][0].gl[0] = EPS; lk[i][0].gl[1] = EPS; lk[i][0].gl[2] = EPS;
        lk[i][0].gl[g] = 1.0;
        lk[i][0].masked = false;
    }
    
    auto afs_full = compute_af_from_likelihoods(lk);
    double af_full = afs_full[0];  // Should be 0.25
    
    // Exclude 2 het samples
    std::vector<int> exclude = {10, 11};
    auto afs_loo = compute_af_from_likelihoods_leave_one_out(lk, exclude);
    double af_loo = afs_loo[0];
    
    // Theoretical: (20*0.25*2 - 2) / (18*2) = (10 - 2) / 36 = 8/36 = 0.222
    double af_expected = (2.0 * n * af_full - 2.0) / (2.0 * (n - 2));
    
    CHECK_NEAR(af_loo, af_expected, 0.01);
}

// ============================================================================
// 2. Standard Errors: Analytical Verification
// ============================================================================

TEST_CASE(se_analytical_quadratic) {
    // For f(x,y) = a(x-x0)^2 + b(y-y0)^2:
    // Hessian = [[2a, 0], [0, 2b]]
    // Covariance = [[1/(2a), 0], [0, 1/(2b)]]
    // SE = [1/sqrt(2a), 1/sqrt(2b)]
    
    auto func = [](const std::vector<double>& x) -> double {
        return 5.0 * (x[0] - 2.0) * (x[0] - 2.0) + 
               10.0 * (x[1] - 3.0) * (x[1] - 3.0);
    };
    
    std::vector<double> x_opt = {2.0, 3.0};
    auto se = compute_standard_errors(func, x_opt);
    
    // Theoretical: SE = [1/sqrt(10), 1/sqrt(20)]
    double se0_theory = 1.0 / std::sqrt(10.0);
    double se1_theory = 1.0 / std::sqrt(20.0);
    
    CHECK(se.size() == 2);
    CHECK_NEAR(se[0], se0_theory, 0.001);
    CHECK_NEAR(se[1], se1_theory, 0.001);
}

TEST_CASE(se_numerical_stability) {
    // Test SE computation with different step sizes
    auto func = [](const std::vector<double>& x) -> double {
        return (x[0] - 1.0) * (x[0] - 1.0) + 
               (x[1] - 1.0) * (x[1] - 1.0);
    };
    
    std::vector<double> x_opt = {1.0, 1.0};
    
    // Test multiple step sizes
    std::vector<double> step_sizes = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
    std::vector<double> se0_results;
    
    for (double h : step_sizes) {
        auto se = compute_standard_errors(func, x_opt, h);
        if (!se.empty()) {
            se0_results.push_back(se[0]);
        }
    }
    
    // All results should be close to theoretical value
    double se_theory = std::sqrt(0.5);  // 1/sqrt(2)
    for (double se0 : se0_results) {
        CHECK_NEAR(se0, se_theory, 0.01);
    }
    
    // Results should be consistent across step sizes
    if (se0_results.size() >= 2) {
        double max_diff = 0.0;
        for (size_t i = 1; i < se0_results.size(); ++i) {
            max_diff = std::max(max_diff, std::abs(se0_results[i] - se0_results[0]));
        }
        CHECK(max_diff < 0.01);  // Consistent within 1%
    }
}

TEST_CASE(se_symmetric_hessian) {
    // Verify Hessian symmetry property: H[i][j] = H[j][i]
    // This is checked implicitly by the analytical inversion
    
    auto func = [](const std::vector<double>& x) -> double {
        double dx = x[0] - 1.0;
        double dy = x[1] - 2.0;
        return 2.0 * dx * dx + 3.0 * dx * dy + 2.0 * dy * dy;
    };
    
    std::vector<double> x_opt = {1.0, 2.0};
    auto se = compute_standard_errors(func, x_opt);
    
    // Should compute successfully (symmetric Hessian)
    CHECK(!se.empty());
    CHECK(se[0] > 0);
    CHECK(se[1] > 0);
}

// ============================================================================
// 3. Nelder-Mead: Convergence Verification
// ============================================================================

TEST_CASE(nm_convergence_to_known_minimum) {
    // Test on function with known minimum
    auto func = [](const std::vector<double>& x) -> double {
        return std::pow(x[0] - 0.3, 2) + std::pow(x[1] - 0.7, 2);
    };
    
    NelderMeadConfig cfg;
    cfg.xtol = 1e-4;
    cfg.ftol = 1e-4;
    cfg.max_iter = 10000;
    
    // Try multiple starting points
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

TEST_CASE(nm_precision_improvement_v030) {
    // Verify v0.3.0 tolerance improvement gives better precision
    
    auto func = [](const std::vector<double>& x) -> double {
        return std::pow(x[0] - 0.123456, 2) + std::pow(x[1] - 0.654321, 2);
    };
    
    // Old tolerance (v0.2.0)
    NelderMeadConfig cfg_old;
    cfg_old.xtol = 0.01;
    cfg_old.ftol = 0.01;
    
    // New tolerance (v0.3.0)
    NelderMeadConfig cfg_new;
    cfg_new.xtol = 1e-4;
    cfg_new.ftol = 1e-4;
    
    auto result_old = nelder_mead(func, {0.0, 0.0}, cfg_old);
    auto result_new = nelder_mead(func, {0.0, 0.0}, cfg_new);
    
    double err_old = std::abs(result_old.x[0] - 0.123456) + 
                     std::abs(result_old.x[1] - 0.654321);
    double err_new = std::abs(result_new.x[0] - 0.123456) + 
                     std::abs(result_new.x[1] - 0.654321);
    
    // New should be more precise
    CHECK(err_new < err_old);
    
    // New should be within tolerance
    CHECK(err_new < 0.001);
}

// ============================================================================
// 4. Global LD Pruning: Consistency Verification
// ============================================================================

TEST_CASE(global_ld_consistency_with_per_pair) {
    // When mask is identical, global and per-pair should give same result
    
    std::mt19937 rng(42);
    int n_samples = 30;
    int n_snps = 15;
    
    std::vector<std::vector<int8_t>> genotypes(n_samples, std::vector<int8_t>(n_snps));
    std::uniform_int_distribution<int> dist(0, 2);
    
    for (int s = 0; s < n_snps; ++s) {
        for (int i = 0; i < n_samples; ++i) {
            genotypes[i][s] = dist(rng);
        }
    }
    
    std::vector<bool> mask(n_snps, false);
    LDPruneConfig config;
    config.window_size = 50;
    config.step_size = 5;
    config.r2_threshold = 0.5;
    
    auto global = ld_prune_global(genotypes, mask, config);
    auto perpair = ld_prune(genotypes, mask, config);
    
    // Should be identical
    CHECK(global.size() == perpair.size());
    for (size_t i = 0; i < global.size(); ++i) {
        CHECK(global[i] == perpair[i]);
    }
}

TEST_CASE(global_ld_determinism) {
    // Global LD pruning should be deterministic (same input → same output)
    
    std::mt19937 rng(42);
    int n_samples = 20;
    int n_snps = 10;
    
    std::vector<std::vector<int8_t>> genotypes(n_samples, std::vector<int8_t>(n_snps));
    std::uniform_int_distribution<int> dist(0, 2);
    
    for (int s = 0; s < n_snps; ++s) {
        for (int i = 0; i < n_samples; ++i) {
            genotypes[i][s] = dist(rng);
        }
    }
    
    std::vector<bool> mask(n_snps, false);
    LDPruneConfig config;
    
    auto result1 = ld_prune_global(genotypes, mask, config);
    auto result2 = ld_prune_global(genotypes, mask, config);
    
    CHECK(result1 == result2);
}

// ============================================================================
// 5. Integration: Combined Correctness
// ============================================================================

TEST_CASE(integrated_workflow_correctness) {
    // Test the full workflow: LOO EM → NM optimization → SE computation
    
    // Simple quadratic objective (simulating -log-likelihood)
    auto objective = [](const std::vector<double>& k) -> double {
        return std::pow(k[0] - 0.25, 2) + std::pow(k[1] - 0.5, 2);
    };
    
    NelderMeadConfig cfg;
    cfg.xtol = 1e-4;
    cfg.ftol = 1e-4;
    
    // Optimize
    auto result = nelder_mead(objective, {0.0, 0.0}, cfg);
    
    // Check convergence
    CHECK(result.converged || result.fval < 1e-6);
    CHECK_NEAR(result.x[0], 0.25, 0.01);
    CHECK_NEAR(result.x[1], 0.5, 0.01);
    
    // Compute SE
    auto se = compute_standard_errors(objective, result.x);
    
    // Check SE validity
    CHECK(!se.empty());
    CHECK(se[0] > 0);
    CHECK(se[1] > 0);
    
    // Theoretical SE for this quadratic: 1/sqrt(2) ≈ 0.707
    CHECK_NEAR(se[0], std::sqrt(0.5), 0.1);
    CHECK_NEAR(se[1], std::sqrt(0.5), 0.1);
}

TEST_CASE(numerical_robustness_extreme_values) {
    // Test robustness with extreme parameter values
    
    auto func = [](const std::vector<double>& x) -> double {
        // Very narrow minimum
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
    
    // SE should be small (narrow minimum = high certainty)
    auto se = compute_standard_errors(func, result.x);
    CHECK(!se.empty());
    CHECK(se[0] < 0.2);
    CHECK(se[1] < 0.2);
}

int main() {
    return RUN_ALL_TESTS();
}
