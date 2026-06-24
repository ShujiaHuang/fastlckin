/**
 * @file test_frequency_leave_one_out.cpp
 * @brief Unit tests for leave-one-out EM allele frequency estimation
 */

#include "test_harness.h"
#include "frequency_from_likelihoods.h"
#include <cmath>
#include <vector>

using namespace fastlckin;

// Helper: create a LikelihoodMatrix with known hard genotypes (delta-function)
static LikelihoodMatrix make_hard_lk_matrix(int n_samples, int n_snps,
                                             const std::vector<std::vector<int>>& genotypes)
{
    const double EPS = 1e-10;
    LikelihoodMatrix lk(n_samples, std::vector<GenotypeLikelihood>(n_snps));
    for (int i = 0; i < n_samples; ++i) {
        for (int s = 0; s < n_snps; ++s) {
            int g = genotypes[i][s];
            if (g < 0) {
                lk[i][s].masked = true;
                lk[i][s].gl[0] = EPS;
                lk[i][s].gl[1] = EPS;
                lk[i][s].gl[2] = EPS;
            } else {
                lk[i][s].gl[0] = EPS;
                lk[i][s].gl[1] = EPS;
                lk[i][s].gl[2] = EPS;
                lk[i][s].gl[g] = 1.0;
                lk[i][s].masked = false;
            }
        }
    }
    return lk;
}

TEST_CASE(leave_one_out_basic) {
    // 20 samples, 1 SNP
    // True AF (all samples) = 0.3 → 12 alt alleles out of 40
    // Genotypes: 8 hom_ref(0), 12 het(1)
    int n = 20;
    std::vector<std::vector<int>> gt(n, std::vector<int>(1));
    for (int i = 0; i < 8; ++i) gt[i][0] = 0;
    for (int i = 8; i < 20; ++i) gt[i][0] = 1;

    auto lk = make_hard_lk_matrix(n, 1, gt);

    // Exclude samples 0 and 1 (both hom_ref)
    // Remaining: 6 hom_ref + 12 het → 12 alt alleles out of 36
    // AF = 12/36 = 0.333...
    std::vector<int> exclude = {0, 1};
    auto afs = compute_af_from_likelihoods_leave_one_out(lk, exclude);

    CHECK_NEAR(afs[0], 12.0/36.0, 0.01);
}

TEST_CASE(leave_one_out_exclude_het_samples) {
    // Same setup: 8 hom_ref + 12 het
    int n = 20;
    std::vector<std::vector<int>> gt(n, std::vector<int>(1));
    for (int i = 0; i < 8; ++i) gt[i][0] = 0;
    for (int i = 8; i < 20; ++i) gt[i][0] = 1;

    auto lk = make_hard_lk_matrix(n, 1, gt);

    // Exclude samples 8 and 9 (both het)
    // Remaining: 8 hom_ref + 10 het → 10 alt alleles out of 36
    // AF = 10/36 = 0.277...
    std::vector<int> exclude = {8, 9};
    auto afs = compute_af_from_likelihoods_leave_one_out(lk, exclude);

    CHECK_NEAR(afs[0], 10.0/36.0, 0.01);
}

TEST_CASE(leave_one_out_vs_full_comparison) {
    // Demonstrate the bias reduction: comparing full EM vs leave-one-out
    int n = 10;
    std::vector<std::vector<int>> gt(n, std::vector<int>(1));
    // 5 hom_ref, 5 het → AF = 5/20 = 0.25
    for (int i = 0; i < 5; ++i) gt[i][0] = 0;
    for (int i = 5; i < 10; ++i) gt[i][0] = 1;

    auto lk = make_hard_lk_matrix(n, 1, gt);

    // Full EM
    auto afs_full = compute_af_from_likelihoods(lk);
    CHECK_NEAR(afs_full[0], 0.25, 0.01);

    // Leave-one-out (exclude 2 het samples)
    // Remaining: 5 hom_ref + 3 het → 3 alt alleles out of 16
    // AF = 3/16 = 0.1875
    std::vector<int> exclude = {5, 6};
    auto afs_loo = compute_af_from_likelihoods_leave_one_out(lk, exclude);

    CHECK_NEAR(afs_loo[0], 0.1875, 0.01);

    // Verify they're different (leave-one-out corrects the bias)
    CHECK(std::abs(afs_full[0] - afs_loo[0]) > 0.05);
}

TEST_CASE(leave_one_out_exclude_all) {
    // Edge case: exclude all samples
    int n = 5;
    std::vector<std::vector<int>> gt(n, std::vector<int>(1, 0));

    auto lk = make_hard_lk_matrix(n, 1, gt);
    std::vector<int> exclude = {0, 1, 2, 3, 4};
    auto afs = compute_af_from_likelihoods_leave_one_out(lk, exclude);

    // With no valid samples, AF should be 0
    CHECK_NEAR(afs[0], 0.0, 1e-10);
}

TEST_CASE(leave_one_out_multiple_snps) {
    // 10 samples, 3 SNPs
    int n = 10;
    std::vector<std::vector<int>> gt(n, std::vector<int>(3));

    // SNP 0: 5 hom_ref, 5 het → AF_full = 0.25
    for (int i = 0; i < 5; ++i) gt[i][0] = 0;
    for (int i = 5; i < n; ++i) gt[i][0] = 1;

    // SNP 1: 8 hom_ref, 2 het → AF_full = 0.1
    for (int i = 0; i < 8; ++i) gt[i][1] = 0;
    for (int i = 8; i < n; ++i) gt[i][1] = 1;

    // SNP 2: all hom_alt → AF_full = 1.0
    for (int i = 0; i < n; ++i) gt[i][2] = 2;

    auto lk = make_hard_lk_matrix(n, 3, gt);

    // Exclude samples 0 and 5
    std::vector<int> exclude = {0, 5};
    auto afs = compute_af_from_likelihoods_leave_one_out(lk, exclude);

    CHECK(afs.size() == 3);

    // SNP 0: exclude 1 hom_ref + 1 het → 4 hom_ref + 4 het → 4/16 = 0.25
    CHECK_NEAR(afs[0], 0.25, 0.01);

    // SNP 1: exclude 1 hom_ref → 7 hom_ref + 2 het → 2/16 = 0.125
    CHECK_NEAR(afs[1], 0.125, 0.01);

    // SNP 2: exclude 2 hom_alt → 8 hom_alt → 16/16 = 1.0
    CHECK_NEAR(afs[2], 1.0, 0.01);
}

TEST_CASE(leave_one_out_with_missing) {
    // Some samples already masked + additional exclusion
    int n = 10;
    std::vector<std::vector<int>> gt(n, std::vector<int>(1));
    for (int i = 0; i < 5; ++i) gt[i][0] = 0;
    for (int i = 5; i < 8; ++i) gt[i][0] = 1;
    for (int i = 8; i < 10; ++i) gt[i][0] = -1;  // missing

    auto lk = make_hard_lk_matrix(n, 1, gt);

    // Exclude sample 0 (hom_ref)
    // Valid after exclusion: 4 hom_ref (samples 1-4) + 3 het (samples 5-7)
    // AF = 3 / (2 * 7) = 3/14 = 0.2143
    std::vector<int> exclude = {0};
    auto afs = compute_af_from_likelihoods_leave_one_out(lk, exclude);

    CHECK_NEAR(afs[0], 3.0/14.0, 0.01);
}

TEST_CASE(leave_one_out_empty_matrix) {
    LikelihoodMatrix empty;
    std::vector<int> exclude = {0, 1};
    auto afs = compute_af_from_likelihoods_leave_one_out(empty, exclude);
    CHECK(afs.empty());
}

TEST_CASE(leave_one_out_invalid_indices) {
    // Invalid indices should be ignored gracefully
    int n = 5;
    std::vector<std::vector<int>> gt(n, std::vector<int>(1, 0));

    auto lk = make_hard_lk_matrix(n, 1, gt);
    std::vector<int> exclude = {-1, 100, 0};  // invalid + valid
    auto afs = compute_af_from_likelihoods_leave_one_out(lk, exclude);

    // Only sample 0 is excluded (valid index)
    // Remaining: 4 hom_ref → AF = 0
    CHECK(afs[0] < 0.01);
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

TEST_CASE(loo_handles_edge_cases) {
    // All samples identical (hom_alt) — LOO should still work
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

    CHECK(afs_full.size() == 1);
    CHECK(afs_loo.size() == 1);
    CHECK(afs_full[0] > 0.99);
    CHECK(afs_loo[0] > 0.99);
}

int main() {
    return RUN_ALL_TESTS();
}
