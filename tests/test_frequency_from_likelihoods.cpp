/**
 * @file test_frequency_from_likelihoods.cpp
 * @brief Unit tests for EM-based allele frequency estimation from genotype likelihoods
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

TEST_CASE(af_from_likelihoods_known_frequency) {
    // 20 samples, 1 SNP with known genotypes
    // True AF = 0.3 → 12 alt alleles out of 40
    // Genotypes: 8 hom_ref(0), 12 het(1), 0 hom_alt(2)
    // That gives AF = 12/40 = 0.3
    int n = 20;
    std::vector<std::vector<int>> gt(n, std::vector<int>(1));
    for (int i = 0; i < 8; ++i) gt[i][0] = 0;   // hom ref
    for (int i = 8; i < 20; ++i) gt[i][0] = 1;   // het

    auto lk = make_hard_lk_matrix(n, 1, gt);
    auto afs = compute_af_from_likelihoods(lk);

    // With hard genotypes (delta-function), EM should converge to the exact count
    CHECK_NEAR(afs[0], 0.3, 0.01);
}

TEST_CASE(af_from_likelihoods_monomorphic_ref) {
    // All samples are hom_ref → AF should be ~0
    int n = 10;
    std::vector<std::vector<int>> gt(n, std::vector<int>(1, 0));

    auto lk = make_hard_lk_matrix(n, 1, gt);
    auto afs = compute_af_from_likelihoods(lk);

    CHECK(afs[0] < 0.01);
}

TEST_CASE(af_from_likelihoods_monomorphic_alt) {
    // All samples are hom_alt → AF should be ~1
    int n = 10;
    std::vector<std::vector<int>> gt(n, std::vector<int>(1, 2));

    auto lk = make_hard_lk_matrix(n, 1, gt);
    auto afs = compute_af_from_likelihoods(lk);

    CHECK(afs[0] > 0.99);
}

TEST_CASE(af_from_likelihoods_with_missing) {
    // Some samples masked → should ignore them
    int n = 10;
    std::vector<std::vector<int>> gt(n, std::vector<int>(1));
    for (int i = 0; i < 5; ++i) gt[i][0] = 0;   // hom ref
    for (int i = 5; i < 8; ++i) gt[i][0] = 1;   // het
    for (int i = 8; i < 10; ++i) gt[i][0] = -1;  // missing

    auto lk = make_hard_lk_matrix(n, 1, gt);
    auto afs = compute_af_from_likelihoods(lk);

    // 8 valid samples: 5 hom_ref + 3 het → 3 alt alleles out of 16
    // AF = 3/16 = 0.1875
    CHECK_NEAR(afs[0], 0.1875, 0.01);
}

TEST_CASE(af_from_likelihoods_all_masked) {
    // All samples masked → AF should be 0
    int n = 5;
    std::vector<std::vector<int>> gt(n, std::vector<int>(1, -1));

    auto lk = make_hard_lk_matrix(n, 1, gt);
    auto afs = compute_af_from_likelihoods(lk);

    CHECK_NEAR(afs[0], 0.0, 1e-10);
}

TEST_CASE(af_from_likelihoods_multiple_snps) {
    // 10 samples, 3 SNPs with different frequencies
    int n = 10;
    std::vector<std::vector<int>> gt(n, std::vector<int>(3));

    // SNP 0: all hom_ref → AF ≈ 0
    for (int i = 0; i < n; ++i) gt[i][0] = 0;
    // SNP 1: 5 het, 5 hom_ref → AF = 5/20 = 0.25
    for (int i = 0; i < 5; ++i) gt[i][1] = 1;
    for (int i = 5; i < n; ++i) gt[i][1] = 0;
    // SNP 2: all hom_alt → AF ≈ 1
    for (int i = 0; i < n; ++i) gt[i][2] = 2;

    auto lk = make_hard_lk_matrix(n, 3, gt);
    auto afs = compute_af_from_likelihoods(lk);

    CHECK(afs.size() == 3);
    CHECK_NEAR(afs[0], 0.0, 0.01);
    CHECK_NEAR(afs[1], 0.25, 0.01);
    CHECK_NEAR(afs[2], 1.0, 0.01);
}

TEST_CASE(af_from_likelihoods_empty) {
    LikelihoodMatrix empty;
    auto afs = compute_af_from_likelihoods(empty);
    CHECK(afs.empty());
}

TEST_CASE(expected_genotypes_hard_lk) {
    // With delta-function likelihoods, expected genotypes should match hard genotypes
    int n = 5;
    std::vector<std::vector<int>> gt(n, std::vector<int>(1));
    gt[0][0] = 0; gt[1][0] = 1; gt[2][0] = 2; gt[3][0] = 0; gt[4][0] = 1;

    auto lk = make_hard_lk_matrix(n, 1, gt);
    std::vector<double> afs = {0.3};

    auto eg = compute_expected_genotypes(lk, afs);

    // With delta-function likelihoods, posterior is dominated by likelihood
    // E[G] should be close to the hard genotype
    CHECK_NEAR(eg[0][0], 0.0, 0.01);
    CHECK_NEAR(eg[1][0], 1.0, 0.01);
    CHECK_NEAR(eg[2][0], 2.0, 0.01);
}

TEST_CASE(expected_genotypes_masked) {
    // Masked samples should get -1.0
    int n = 3;
    std::vector<std::vector<int>> gt = {{0}, {-1}, {1}};

    auto lk = make_hard_lk_matrix(n, 1, gt);
    std::vector<double> afs = {0.3};

    auto eg = compute_expected_genotypes(lk, afs);

    CHECK_NEAR(eg[0][0], 0.0, 0.01);
    CHECK(eg[1][0] < -0.5);  // masked → -1.0
    CHECK_NEAR(eg[2][0], 1.0, 0.01);
}

int main() {
    return RUN_ALL_TESTS();
}
