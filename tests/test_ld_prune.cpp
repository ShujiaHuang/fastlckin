/**
 * @file test_ld_prune.cpp
 * @brief Unit tests for ld_prune.h (LD pruning and r² computation)
 */

#include "test_harness.h"
#include "ld_prune.h"
#include <cmath>
#include <random>

TEST_CASE(r2_identical_snps) {
    // Two identical SNPs should have r²=1
    std::vector<int8_t> g = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0};
    CHECK_NEAR(fastlckin::compute_r2(g, g), 1.0, 1e-10);
}

TEST_CASE(r2_uncorrelated) {
    // Constructed to have ~0 correlation
    std::vector<int8_t> g1 = {0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 0, 0, 1, 1, 2, 2, 0, 0};
    std::vector<int8_t> g2 = {0, 2, 0, 2, 0, 2, 1, 1, 1, 1, 1, 1, 0, 2, 0, 2, 0, 2, 1, 1};
    double r2 = fastlckin::compute_r2(g1, g2);
    CHECK(r2 < 0.3);  // Should be low-ish
}

TEST_CASE(r2_with_missing) {
    // Missing values should be skipped
    std::vector<int8_t> g1 = {0, -1, 1, 2, -1, 0, 1, 2, 0, 1};
    std::vector<int8_t> g2 = {0,  1, 1, 2,  0, 0, 1, 2, 0, 1};
    double r2 = fastlckin::compute_r2(g1, g2);
    CHECK(r2 >= 0.0);  // Should still compute
}

TEST_CASE(r2_too_few_samples) {
    std::vector<int8_t> g1 = {0, 1, 2};
    std::vector<int8_t> g2 = {0, 1, 2};
    CHECK_NEAR(fastlckin::compute_r2(g1, g2), 0.0, 1e-15);
}

TEST_CASE(r2_monomorphic) {
    // All same value → variance=0 → r²=0
    std::vector<int8_t> g1(20, 1);
    std::vector<int8_t> g2(20, 0);
    CHECK_NEAR(fastlckin::compute_r2(g1, g2), 0.0, 1e-15);
}

TEST_CASE(ld_prune_no_pruning_needed) {
    // Independent SNPs: r²=0 between all pairs
    int n_samples = 30;
    int n_snps = 10;

    // Create independent random genotypes
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> dist(0, 2);
    std::vector<std::vector<int8_t>> genotypes(n_samples, std::vector<int8_t>(n_snps));
    for (int i = 0; i < n_samples; ++i)
        for (int j = 0; j < n_snps; ++j)
            genotypes[i][j] = static_cast<int8_t>(dist(rng));

    std::vector<bool> mask(n_snps, false);
    fastlckin::LDPruneConfig config;
    config.window_size = 50;
    config.r2_threshold = 0.8;

    auto retained = fastlckin::ld_prune(genotypes, mask, config);
    // All should be retained since they're independent
    CHECK(static_cast<int>(retained.size()) == n_snps);
}

TEST_CASE(ld_prune_removes_high_ld) {
    // Create two identical SNPs (r²=1) plus one independent
    int n_samples = 30;
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> dist(0, 2);

    std::vector<std::vector<int8_t>> genotypes(n_samples, std::vector<int8_t>(3));
    for (int i = 0; i < n_samples; ++i) {
        int8_t g = static_cast<int8_t>(dist(rng));
        genotypes[i][0] = g;       // SNP 0
        genotypes[i][1] = g;       // SNP 1 = copy of SNP 0 (r²=1)
        genotypes[i][2] = static_cast<int8_t>(dist(rng)); // independent SNP 2
    }

    std::vector<bool> mask(3, false);
    fastlckin::LDPruneConfig config;
    config.window_size = 50;
    config.r2_threshold = 0.8;

    auto retained = fastlckin::ld_prune(genotypes, mask, config);
    // One of SNP 0 or 1 should be removed
    CHECK(retained.size() == 2);
}

TEST_CASE(ld_prune_all_masked) {
    std::vector<std::vector<int8_t>> genotypes(10, std::vector<int8_t>(5, 0));
    std::vector<bool> mask(5, true);

    fastlckin::LDPruneConfig config;
    auto retained = fastlckin::ld_prune(genotypes, mask, config);
    CHECK(retained.empty());
}

TEST_CASE(ld_prune_empty) {
    std::vector<std::vector<int8_t>> genotypes;
    std::vector<bool> mask;

    fastlckin::LDPruneConfig config;
    auto retained = fastlckin::ld_prune(genotypes, mask, config);
    CHECK(retained.empty());
}

int main() {
    return RUN_ALL_TESTS();
}
