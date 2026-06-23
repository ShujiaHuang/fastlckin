/**
 * @file test_v030_kinship_accuracy.cpp
 * @brief Simplified verification for v0.3.0 kinship estimation improvements
 * 
 * This test suite provides focused verification that v0.3.0 improvements
 * are correctly integrated into the kinship estimation workflow.
 * 
 * Note: Full end-to-end testing requires complex data generation.
 * This file focuses on validating key integration points.
 */

#include "test_harness.h"
#include "frequency_from_likelihoods.h"
#include "optimizer.h"
#include "kinship_estimator.h"
#include <cmath>
#include <vector>
#include <random>

using namespace fastlckin;

// ============================================================================
// 1. Leave-one-out Integration Validation
// ============================================================================

TEST_CASE(loo_integration_workflow) {
    // Verify that LOO frequency estimation can be used in kinship workflow
    
    // Generate simple test data: 10 samples, 5 SNPs
    int n_samples = 10;
    int n_snps = 5;
    
    const double EPS = 1e-10;
    LikelihoodMatrix lk(n_samples, std::vector<GenotypeLikelihood>(n_snps));
    
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> gt_dist(0, 2);
    
    for (int i = 0; i < n_samples; ++i) {
        for (int s = 0; s < n_snps; ++s) {
            int g = gt_dist(rng);
            lk[i][s].gl[0] = EPS;
            lk[i][s].gl[1] = EPS;
            lk[i][s].gl[2] = EPS;
            lk[i][s].gl[g] = 1.0;
            lk[i][s].masked = false;
        }
    }
    
    // Compute full EM frequencies
    auto afs_full = compute_af_from_likelihoods(lk);
    
    // Compute LOO frequencies for pair (0, 1)
    std::vector<int> exclude = {0, 1};
    auto afs_loo = compute_af_from_likelihoods_leave_one_out(lk, exclude);
    
    // Verify both produce valid frequencies
    CHECK(afs_full.size() == n_snps);
    CHECK(afs_loo.size() == n_snps);
    
    for (int s = 0; s < n_snps; ++s) {
        CHECK(afs_full[s] > 0 && afs_full[s] < 1);
        CHECK(afs_loo[s] > 0 && afs_loo[s] < 1);
    }
    
    // LOO should differ from full EM (proving bias correction is applied)
    double total_diff = 0.0;
    for (int s = 0; s < n_snps; ++s) {
        total_diff += std::abs(afs_full[s] - afs_loo[s]);
    }
    double avg_diff = total_diff / n_snps;
    
    // Average difference should be non-zero (LOO is working)
    CHECK(avg_diff > 1e-6);
}

TEST_CASE(loo_handles_edge_cases) {
    // Verify LOO handles edge cases gracefully
    
    // Edge case 1: All samples identical (extreme scenario)
    int n = 5;
    const double EPS = 1e-10;
    LikelihoodMatrix lk(n, std::vector<GenotypeLikelihood>(1));
    
    for (int i = 0; i < n; ++i) {
        lk[i][0].gl[0] = EPS;
        lk[i][0].gl[1] = EPS;
        lk[i][0].gl[2] = 1.0;  // All hom_alt
        lk[i][0].masked = false;
    }
    
    auto afs_full = compute_af_from_likelihoods(lk);
    std::vector<int> exclude = {0, 1};
    auto afs_loo = compute_af_from_likelihoods_leave_one_out(lk, exclude);
    
    // Both should return valid frequencies
    CHECK(afs_full.size() == 1);
    CHECK(afs_loo.size() == 1);
    CHECK(afs_full[0] > 0.99);  // Nearly all alt alleles
    CHECK(afs_loo[0] > 0.99);
}

// ============================================================================
// 2. Standard Error Integration Validation
// ============================================================================

TEST_CASE(se_api_usability) {
    // Verify SE computation is accessible from optimization result
    
    auto func = [](const std::vector<double>& x) -> double {
        return (x[0] - 0.3) * (x[0] - 0.3) + (x[1] - 0.7) * (x[1] - 0.7);
    };
    
    NelderMeadConfig cfg;
    cfg.xtol = 1e-4;
    cfg.ftol = 1e-4;
    
    auto result = nelder_mead(func, {0.0, 0.0}, cfg);
    
    // Verify optimization succeeded
    CHECK(result.converged || result.fval < 1e-6);
    CHECK_NEAR(result.x[0], 0.3, 0.01);
    CHECK_NEAR(result.x[1], 0.7, 0.01);
    
    // Compute SE at the optimum
    auto se = compute_standard_errors(func, result.x);
    
    // SE should be valid
    CHECK(se.size() == 2);
    CHECK(se[0] > 0);
    CHECK(se[1] > 0);
    
    // Theoretical SE for this quadratic: 1/sqrt(2) ≈ 0.707
    CHECK_NEAR(se[0], std::sqrt(0.5), 0.1);
    CHECK_NEAR(se[1], std::sqrt(0.5), 0.1);
}

TEST_CASE(se_reflects_uncertainty) {
    // Verify SE actually reflects estimation uncertainty
    
    // Scenario 1: Tight minimum (high certainty) → Small SE
    auto tight_func = [](const std::vector<double>& x) -> double {
        return 100.0 * (x[0] - 0.5) * (x[0] - 0.5) + 
               100.0 * (x[1] - 0.5) * (x[1] - 0.5);
    };
    
    // Scenario 2: Broad minimum (low certainty) → Large SE
    auto broad_func = [](const std::vector<double>& x) -> double {
        return 1.0 * (x[0] - 0.5) * (x[0] - 0.5) + 
               1.0 * (x[1] - 0.5) * (x[1] - 0.5);
    };
    
    std::vector<double> x_opt = {0.5, 0.5};
    
    auto se_tight = compute_standard_errors(tight_func, x_opt);
    auto se_broad = compute_standard_errors(broad_func, x_opt);
    
    // Tight minimum should have smaller SE
    CHECK(se_tight.size() == 2);
    CHECK(se_broad.size() == 2);
    CHECK(se_tight[0] < se_broad[0]);
    CHECK(se_tight[1] < se_broad[1]);
}

// ============================================================================
// 3. Global LD Pruning Integration Validation
// ============================================================================

TEST_CASE(global_ld_prune_api) {
    // Verify global LD pruning API is accessible
    
    int n_samples = 10;
    int n_snps = 5;
    
    std::vector<std::vector<int8_t>> genotypes(n_samples, std::vector<int8_t>(n_snps));
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> dist(0, 2);
    
    for (int s = 0; s < n_snps; ++s) {
        for (int i = 0; i < n_samples; ++i) {
            genotypes[i][s] = dist(rng);
        }
    }
    
    std::vector<bool> mask(n_snps, false);
    LDPruneConfig config;
    config.r2_threshold = 0.5;
    
    // Call global LD pruning
    auto global_result = ld_prune_global(genotypes, mask, config);
    auto per_pair_result = ld_prune(genotypes, mask, config);
    
    // Should produce identical results
    CHECK(global_result.size() == per_pair_result.size());
    for (size_t i = 0; i < global_result.size(); ++i) {
        CHECK(global_result[i] == per_pair_result[i]);
    }
}

TEST_CASE(global_ld_prune_with_expected_genotypes) {
    // Verify global LD pruning works with expected genotypes
    
    int n_samples = 10;
    int n_snps = 5;
    
    std::vector<std::vector<double>> expected_g(n_samples, std::vector<double>(n_snps));
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 2.0);
    
    for (int s = 0; s < n_snps; ++s) {
        for (int i = 0; i < n_samples; ++i) {
            expected_g[i][s] = dist(rng);
        }
    }
    
    std::vector<bool> mask(n_snps, false);
    LDPruneConfig config;
    
    // Should work without crashing
    auto result = ld_prune_global_expected(expected_g, mask, config);
    CHECK(result.size() > 0);
}

// ============================================================================
// 4. v0.3.0 Default Configuration Validation
// ============================================================================

TEST_CASE(v030_default_config_updated) {
    // Verify v0.3.0 default configurations are correctly set
    
    // Optimizer defaults
    NelderMeadConfig nm_cfg;
    CHECK_NEAR(nm_cfg.xtol, 1e-4, 1e-10);
    CHECK_NEAR(nm_cfg.ftol, 1e-4, 1e-10);
    
    // Kinship estimator defaults
    KinshipConfig kin_cfg;
    CHECK(kin_cfg.n_restarts == 5);
    CHECK(kin_cfg.global_ld_prune == false);  // Default off
}

TEST_CASE(v030_improved_precision_on_kinship_objective) {
    // Verify tighter tolerance improves precision on kinship-like objective
    
    // Simulate kinship objective with known minimum
    auto kinship_obj = [](const std::vector<double>& k) -> double {
        double k0 = 1.0 - k[0] - k[1];
        double k1 = k[0];
        double k2 = k[1];
        
        // Penalize invalid regions
        if (k0 < 0 || k1 < 0 || k2 < 0) return 1e10;
        if (4.0 * k0 * k2 > k1 * k1 + 1e-10) return 1e10;
        
        // Quadratic around (k1=0.5, k2=0.25) = full siblings
        return std::pow(k1 - 0.5, 2) + std::pow(k2 - 0.25, 2);
    };
    
    // Old tolerance (v0.2.0)
    NelderMeadConfig cfg_old;
    cfg_old.xtol = 0.01;
    cfg_old.ftol = 0.01;
    auto result_old = nelder_mead(kinship_obj, {0.5, 0.2}, cfg_old);
    
    // New tolerance (v0.3.0)
    NelderMeadConfig cfg_new;
    cfg_new.xtol = 1e-4;
    cfg_new.ftol = 1e-4;
    auto result_new = nelder_mead(kinship_obj, {0.5, 0.2}, cfg_new);
    
    // New should be more precise
    double err_old = std::abs(result_old.x[0] - 0.5) + std::abs(result_old.x[1] - 0.25);
    double err_new = std::abs(result_new.x[0] - 0.5) + std::abs(result_new.x[1] - 0.25);
    
    CHECK(err_new < err_old);
    CHECK(err_new < 0.001);  // High precision
}

int main() {
    return RUN_ALL_TESTS();
}
