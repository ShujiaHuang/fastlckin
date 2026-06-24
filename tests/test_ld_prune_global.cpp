/**
 * @file test_ld_prune_global.cpp
 * @brief Unit tests for global LD pruning (v0.3.0 feature)
 */

#include "test_harness.h"
#include "ld_prune.h"
#include <cmath>
#include <vector>
#include <random>

using namespace fastlckin;

// Helper: create genotype matrix with known LD structure
static std::vector<std::vector<int8_t>> make_test_genotypes_ld(
    int n_samples, int n_snps, std::mt19937& rng)
{
    std::vector<std::vector<int8_t>> genotypes(n_samples, std::vector<int8_t>(n_snps));
    
    std::uniform_int_distribution<int> allele_dist(0, 1);
    
    // SNP 0-9: independent
    for (int s = 0; s < 10; ++s) {
        for (int i = 0; i < n_samples; ++i) {
            genotypes[i][s] = allele_dist(rng) + allele_dist(rng);
        }
    }
    
    // SNP 10-19: highly correlated with SNP 0-9 (r² ≈ 1)
    for (int s = 0; s < 10; ++s) {
        for (int i = 0; i < n_samples; ++i) {
            genotypes[i][s + 10] = genotypes[i][s];  // Perfect LD
        }
    }
    
    return genotypes;
}

TEST_CASE(global_ld_prune_basic) {
    // Test that global LD pruning removes redundant SNPs
    std::mt19937 rng(42);
    int n_samples = 50;
    int n_snps = 20;  // 10 independent + 10 in perfect LD
    
    auto genotypes = make_test_genotypes_ld(n_samples, n_snps, rng);
    
    // No initial mask
    std::vector<bool> mask(n_snps, false);
    
    LDPruneConfig config;
    config.window_size = 50;
    config.step_size = 5;
    config.r2_threshold = 0.5;
    
    auto retained = ld_prune_global(genotypes, mask, config);
    
    // Should remove ~10 SNPs (the correlated ones)
    CHECK(retained.size() < n_snps);
    CHECK(retained.size() >= 10);  // At least the independent ones
}

TEST_CASE(global_ld_prune_vs_per_pair) {
    // Global and per-pair LD pruning should give same results when mask is same
    std::mt19937 rng(42);
    int n_samples = 50;
    int n_snps = 20;
    
    auto genotypes = make_test_genotypes_ld(n_samples, n_snps, rng);
    std::vector<bool> mask(n_snps, false);
    
    LDPruneConfig config;
    config.window_size = 50;
    config.step_size = 5;
    config.r2_threshold = 0.5;
    
    auto global_retained = ld_prune_global(genotypes, mask, config);
    auto perpair_retained = ld_prune(genotypes, mask, config);
    
    // Should be identical
    CHECK(global_retained.size() == perpair_retained.size());
    for (size_t i = 0; i < global_retained.size(); ++i) {
        CHECK(global_retained[i] == perpair_retained[i]);
    }
}

TEST_CASE(global_ld_prune_with_mask) {
    // Test that global LD pruning respects the mask
    std::mt19937 rng(42);
    int n_samples = 50;
    int n_snps = 20;
    
    auto genotypes = make_test_genotypes_ld(n_samples, n_snps, rng);
    
    // Mask out first 5 SNPs
    std::vector<bool> mask(n_snps, false);
    for (int s = 0; s < 5; ++s) mask[s] = true;
    
    LDPruneConfig config;
    config.window_size = 50;
    config.step_size = 5;
    config.r2_threshold = 0.5;
    
    auto retained = ld_prune_global(genotypes, mask, config);
    
    // Masked SNPs should not be in retained
    for (int s : retained) {
        CHECK(s >= 5);
    }
}

TEST_CASE(global_ld_prune_expected_genotypes) {
    // Test global LD pruning with expected (continuous) genotypes
    std::mt19937 rng(42);
    int n_samples = 50;
    int n_snps = 20;
    
    std::vector<std::vector<double>> expected_g(n_samples, std::vector<double>(n_snps));
    std::uniform_real_distribution<double> dist(0.0, 2.0);
    
    // SNP 0-9: independent
    for (int s = 0; s < 10; ++s) {
        for (int i = 0; i < n_samples; ++i) {
            expected_g[i][s] = dist(rng);
        }
    }
    
    // SNP 10-19: correlated with SNP 0-9
    for (int s = 0; s < 10; ++s) {
        for (int i = 0; i < n_samples; ++i) {
            expected_g[i][s + 10] = expected_g[i][s] + (dist(rng) - 1.0) * 0.1;
        }
    }
    
    std::vector<bool> mask(n_snps, false);
    
    LDPruneConfig config;
    config.window_size = 50;
    config.step_size = 5;
    config.r2_threshold = 0.5;
    
    auto retained = ld_prune_global_expected(expected_g, mask, config);
    
    // Should remove correlated SNPs
    CHECK(retained.size() < n_snps);
    CHECK(retained.size() >= 10);
}

TEST_CASE(global_ld_prune_efficiency) {
    // Global LD pruning should be faster than per-pair for many samples
    // This is more of a conceptual test - we verify the API works
    
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
    
    // Should work without errors
    auto global = ld_prune_global(genotypes, mask, config);
    CHECK(!global.empty() || n_snps == 0);
}

TEST_CASE(global_ld_prune_empty_input) {
    // Test edge cases
    std::vector<std::vector<int8_t>> empty_genotypes;
    std::vector<bool> empty_mask;
    LDPruneConfig config;
    
    auto retained = ld_prune_global(empty_genotypes, empty_mask, config);
    CHECK(retained.empty());
}

TEST_CASE(global_ld_prune_all_masked) {
    // All SNPs masked
    int n_samples = 10;
    int n_snps = 5;
    
    std::vector<std::vector<int8_t>> genotypes(n_samples, std::vector<int8_t>(n_snps, 1));
    std::vector<bool> mask(n_snps, true);  // All masked
    
    LDPruneConfig config;
    auto retained = ld_prune_global(genotypes, mask, config);
    
    CHECK(retained.empty());
}

TEST_CASE(global_ld_prune_determinism) {
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

int main() {
    return RUN_ALL_TESTS();
}
