/**
 * @file test_ld_prune_gl.cpp
 * @brief Unit tests for GL-based LD pruning (expected genotype r²)
 */

#include "test_harness.h"
#include "ld_prune.h"
#include <cmath>
#include <random>
#include <vector>

TEST_CASE(r2_expected_identical) {
    // Identical expected genotypes → r²=1
    std::vector<double> g = {0.0, 0.5, 1.0, 1.5, 2.0, 0.0, 0.5, 1.0, 1.5, 2.0};
    CHECK_NEAR(fastlckin::compute_r2_expected(g, g), 1.0, 1e-10);
}

TEST_CASE(r2_expected_uncorrelated) {
    // Constructed to have low correlation
    std::vector<double> g1 = {0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 0.0, 0.0, 1.0, 1.0,
                              2.0, 2.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0, 0.0, 0.0};
    std::vector<double> g2 = {0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 1.0, 1.0, 1.0, 1.0,
                              1.0, 1.0, 0.0, 2.0, 0.0, 2.0, 0.0, 2.0, 1.0, 1.0};
    double r2 = fastlckin::compute_r2_expected(g1, g2);
    CHECK(r2 < 0.3);
}

TEST_CASE(r2_expected_with_missing) {
    // Missing values (-1.0) should be skipped
    std::vector<double> g1 = {0.0, -1.0, 1.0, 2.0, -1.0, 0.0, 1.0, 2.0, 0.0, 1.0};
    std::vector<double> g2 = {0.0,  1.0, 1.0, 2.0,  0.0, 0.0, 1.0, 2.0, 0.0, 1.0};
    double r2 = fastlckin::compute_r2_expected(g1, g2);
    CHECK(r2 >= 0.0);
}

TEST_CASE(r2_expected_too_few) {
    std::vector<double> g1 = {0.0, 1.0, 2.0};
    std::vector<double> g2 = {0.0, 1.0, 2.0};
    CHECK_NEAR(fastlckin::compute_r2_expected(g1, g2), 0.0, 1e-15);
}

TEST_CASE(r2_expected_monomorphic) {
    // All same value → variance=0 → r²=0
    std::vector<double> g1(20, 1.0);
    std::vector<double> g2(20, 0.5);
    CHECK_NEAR(fastlckin::compute_r2_expected(g1, g2), 0.0, 1e-15);
}

TEST_CASE(r2_expected_matches_hard_high_coverage) {
    // With integer-valued expected genotypes (simulating high-coverage),
    // r²_expected should match r² from int8_t version
    int n = 30;
    std::mt19937 rng(42);
    std::uniform_int_distribution<int> dist(0, 2);

    std::vector<int8_t> gi(n);
    std::vector<double> gd(n);
    for (int i = 0; i < n; ++i) {
        gi[i] = static_cast<int8_t>(dist(rng));
        gd[i] = static_cast<double>(gi[i]);
    }

    std::vector<int8_t> gi2(n);
    std::vector<double> gd2(n);
    for (int i = 0; i < n; ++i) {
        gi2[i] = static_cast<int8_t>(dist(rng));
        gd2[i] = static_cast<double>(gi2[i]);
    }

    double r2_hard = fastlckin::compute_r2(gi, gi2);
    double r2_soft = fastlckin::compute_r2_expected(gd, gd2);

    CHECK_NEAR(r2_hard, r2_soft, 1e-6);
}

TEST_CASE(ld_prune_from_gl_independent) {
    // Independent SNPs should all be retained
    int n_samples = 30;
    int n_snps = 10;

    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 2.0);
    std::vector<std::vector<double>> expected_g(n_samples, std::vector<double>(n_snps));
    for (int i = 0; i < n_samples; ++i)
        for (int j = 0; j < n_snps; ++j)
            expected_g[i][j] = dist(rng);

    std::vector<bool> mask(n_snps, false);
    fastlckin::LDPruneConfig config;
    config.window_size = 50;
    config.r2_threshold = 0.8;

    auto retained = fastlckin::ld_prune_from_gl(expected_g, mask, config);
    CHECK(static_cast<int>(retained.size()) == n_snps);
}

TEST_CASE(ld_prune_from_gl_removes_high_ld) {
    // Two identical SNPs + one independent
    int n_samples = 30;
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0.0, 2.0);

    std::vector<std::vector<double>> expected_g(n_samples, std::vector<double>(3));
    for (int i = 0; i < n_samples; ++i) {
        double g = dist(rng);
        expected_g[i][0] = g;           // SNP 0
        expected_g[i][1] = g;           // SNP 1 = copy (r²=1)
        expected_g[i][2] = dist(rng);   // independent SNP 2
    }

    std::vector<bool> mask(3, false);
    fastlckin::LDPruneConfig config;
    config.window_size = 50;
    config.r2_threshold = 0.8;

    auto retained = fastlckin::ld_prune_from_gl(expected_g, mask, config);
    CHECK(retained.size() == 2);
}

TEST_CASE(ld_prune_from_gl_all_masked) {
    std::vector<std::vector<double>> expected_g(10, std::vector<double>(5, 0.5));
    std::vector<bool> mask(5, true);

    fastlckin::LDPruneConfig config;
    auto retained = fastlckin::ld_prune_from_gl(expected_g, mask, config);
    CHECK(retained.empty());
}

TEST_CASE(ld_prune_from_gl_empty) {
    std::vector<std::vector<double>> expected_g;
    std::vector<bool> mask;

    fastlckin::LDPruneConfig config;
    auto retained = fastlckin::ld_prune_from_gl(expected_g, mask, config);
    CHECK(retained.empty());
}

int main() {
    return RUN_ALL_TESTS();
}
