/**
 * @file test_ibd_model.cpp
 * @brief Comprehensive unit tests for ibd_model.h (IBS|IBD probabilities)
 */

#include "test_harness.h"
#include "ibd_model.h"
#include <cmath>
#include <vector>

// ── Mij function tests ──────────────────────────────────────────

TEST_CASE(Mij_basic) {
    CHECK_NEAR(fastlckin::Mij(0.3, 0.0, 0), 0.3, 1e-12);
    CHECK_NEAR(fastlckin::Mij(0.3, 0.0, 1), 0.3, 1e-12);
    CHECK_NEAR(fastlckin::Mij(0.3, 0.0, 2), 0.3, 1e-12);
}

TEST_CASE(Mij_fst_correction) {
    double p = 0.4, fst = 0.1;
    CHECK_NEAR(fastlckin::Mij(p, fst, 0), (1.0 - fst) * p, 1e-12);
    CHECK_NEAR(fastlckin::Mij(p, fst, 1), (1.0 - fst) * p + fst, 1e-12);
    CHECK_NEAR(fastlckin::Mij(p, fst, 2), (1.0 - fst) * p + 2.0 * fst, 1e-12);
}

TEST_CASE(Mij_boundary_values) {
    CHECK_NEAR(fastlckin::Mij(0.0, 0.05, 0), 0.0, 1e-15);
    CHECK_NEAR(fastlckin::Mij(0.0, 0.05, 1), 0.05, 1e-15);
    CHECK_NEAR(fastlckin::Mij(1.0, 0.05, 0), 0.95, 1e-15);
    CHECK_NEAR(fastlckin::Mij(1.0, 0.05, 1), 1.0, 1e-15);
}

// ── Normalization: sum_c B[c][z] == 1 for each z ──────────────

TEST_CASE(ibd_probs_sum_to_one_fst0) {
    auto B = fastlckin::compute_ibs_ibd_probs(0.3, 0.0);
    for (int z = 0; z < 3; ++z) {
        double sum = 0.0;
        for (int c = 0; c < 9; ++c) sum += B[c][z];
        CHECK_NEAR(sum, 1.0, 1e-10);
    }
}

TEST_CASE(ibd_probs_sum_to_one_fst_nonzero) {
    auto B = fastlckin::compute_ibs_ibd_probs(0.25, 0.01);
    for (int z = 0; z < 3; ++z) {
        double sum = 0.0;
        for (int c = 0; c < 9; ++c) sum += B[c][z];
        CHECK_NEAR(sum, 1.0, 1e-10);
    }
}

TEST_CASE(ibd_probs_sum_to_one_multiple_fst) {
    std::vector<double> fst_values = {0.0, 0.001, 0.01, 0.05, 0.1, 0.2};
    std::vector<double> af_values  = {0.05, 0.1, 0.25, 0.3, 0.5, 0.7, 0.9, 0.95};
    for (double fst : fst_values) {
        for (double af : af_values) {
            auto B = fastlckin::compute_ibs_ibd_probs(af, fst);
            for (int z = 0; z < 3; ++z) {
                double sum = 0.0;
                for (int c = 0; c < 9; ++c) sum += B[c][z];
                CHECK_NEAR(sum, 1.0, 1e-10);
            }
        }
    }
}

// ── Structural zeros ────────────────────────────────────────────

TEST_CASE(ibd_probs_ibs0_zero_at_ibd2) {
    auto B = fastlckin::compute_ibs_ibd_probs(0.4, 0.0);
    CHECK_NEAR(B[fastlckin::PPQQ][2], 0.0, 1e-15);
    CHECK_NEAR(B[fastlckin::QQPP][2], 0.0, 1e-15);
    auto B2 = fastlckin::compute_ibs_ibd_probs(0.3, 0.05);
    CHECK_NEAR(B2[fastlckin::PPQQ][2], 0.0, 1e-15);
    CHECK_NEAR(B2[fastlckin::QQPP][2], 0.0, 1e-15);
}

TEST_CASE(ibd_probs_ibs0_zero_at_ibd1) {
    auto B = fastlckin::compute_ibs_ibd_probs(0.3, 0.0);
    CHECK_NEAR(B[fastlckin::PPQQ][1], 0.0, 1e-15);
    CHECK_NEAR(B[fastlckin::QQPP][1], 0.0, 1e-15);
    auto B2 = fastlckin::compute_ibs_ibd_probs(0.4, 0.05);
    CHECK_NEAR(B2[fastlckin::PPQQ][1], 0.0, 1e-15);
    CHECK_NEAR(B2[fastlckin::QQPP][1], 0.0, 1e-15);
}

TEST_CASE(ibd_probs_ibs1_zero_at_ibd2) {
    auto B = fastlckin::compute_ibs_ibd_probs(0.3, 0.0);
    CHECK_NEAR(B[fastlckin::PPPQ][2], 0.0, 1e-15);
    CHECK_NEAR(B[fastlckin::PQPP][2], 0.0, 1e-15);
    CHECK_NEAR(B[fastlckin::PQQQ][2], 0.0, 1e-15);
    CHECK_NEAR(B[fastlckin::QQPQ][2], 0.0, 1e-15);
    auto B2 = fastlckin::compute_ibs_ibd_probs(0.25, 0.05);
    CHECK_NEAR(B2[fastlckin::PPPQ][2], 0.0, 1e-15);
    CHECK_NEAR(B2[fastlckin::PQPP][2], 0.0, 1e-15);
    CHECK_NEAR(B2[fastlckin::PQQQ][2], 0.0, 1e-15);
    CHECK_NEAR(B2[fastlckin::QQPQ][2], 0.0, 1e-15);
}

// ── Symmetry ────────────────────────────────────────────────────

TEST_CASE(ibd_probs_symmetry_fst0) {
    auto B = fastlckin::compute_ibs_ibd_probs(0.3, 0.0);
    auto B2 = fastlckin::compute_ibs_ibd_probs(0.7, 0.0);
    CHECK_NEAR(B[fastlckin::PPQQ][0], B2[fastlckin::QQPP][0], 1e-10);
}

TEST_CASE(ibd_probs_symmetry_swap_individuals) {
    for (double af : {0.2, 0.3, 0.5}) {
        auto B = fastlckin::compute_ibs_ibd_probs(af, 0.01);
        for (int z = 0; z < 3; ++z) {
            CHECK_NEAR(B[fastlckin::PPPQ][z], B[fastlckin::PQPP][z], 1e-12);
            CHECK_NEAR(B[fastlckin::PQQQ][z], B[fastlckin::QQPQ][z], 1e-12);
            CHECK_NEAR(B[fastlckin::PPQQ][z], B[fastlckin::QQPP][z], 1e-12);
        }
    }
}

TEST_CASE(ibd_probs_pq_symmetry) {
    auto B = fastlckin::compute_ibs_ibd_probs(0.5, 0.01);
    for (int z = 0; z < 3; ++z) {
        CHECK_NEAR(B[fastlckin::PPPP][z], B[fastlckin::QQQQ][z], 1e-12);
        CHECK_NEAR(B[fastlckin::PPPQ][z], B[fastlckin::PQQQ][z], 1e-12);
    }
}

// ── Non-negativity ──────────────────────────────────────────────

TEST_CASE(ibd_probs_nonnegative) {
    for (double af : {0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99}) {
        for (double fst : {0.0, 0.001, 0.01, 0.05, 0.1, 0.3}) {
            auto B = fastlckin::compute_ibs_ibd_probs(af, fst);
            for (int c = 0; c < 9; ++c)
                for (int z = 0; z < 3; ++z)
                    CHECK(B[c][z] >= -1e-15);
        }
    }
}

// ── Exact values at FST=0 ───────────────────────────────────────

TEST_CASE(ibd_probs_fst0_exact_values) {
    double p = 0.3, q = 0.7;
    auto B = fastlckin::compute_ibs_ibd_probs(p, 0.0);
    // IBD=2: genotype probabilities
    CHECK_NEAR(B[fastlckin::PPPP][2], p * p, 1e-10);
    CHECK_NEAR(B[fastlckin::PQPQ][2], 2.0 * p * q, 1e-10);
    CHECK_NEAR(B[fastlckin::QQQQ][2], q * q, 1e-10);
    // IBD=0: product of genotype probabilities
    CHECK_NEAR(B[fastlckin::PPQQ][0], p * p * q * q, 1e-10);
    CHECK_NEAR(B[fastlckin::PPPP][0], p * p * p * p, 1e-10);
    CHECK_NEAR(B[fastlckin::PQPQ][0], 4.0 * p * p * q * q, 1e-10);
    CHECK_NEAR(B[fastlckin::PPPQ][0], 2.0 * p * p * p * q, 1e-10);
}

// ── IBD=2 structure ─────────────────────────────────────────────

TEST_CASE(ibd2_only_matching_genotypes) {
    for (double af : {0.1, 0.3, 0.5, 0.7}) {
        auto B = fastlckin::compute_ibs_ibd_probs(af, 0.01);
        CHECK_NEAR(B[fastlckin::PPQQ][2], 0.0, 1e-15);
        CHECK_NEAR(B[fastlckin::QQPP][2], 0.0, 1e-15);
        CHECK_NEAR(B[fastlckin::PPPQ][2], 0.0, 1e-15);
        CHECK_NEAR(B[fastlckin::PQPP][2], 0.0, 1e-15);
        CHECK_NEAR(B[fastlckin::PQQQ][2], 0.0, 1e-15);
        CHECK_NEAR(B[fastlckin::QQPQ][2], 0.0, 1e-15);
        CHECK(B[fastlckin::PPPP][2] > 0.0);
        CHECK(B[fastlckin::PQPQ][2] > 0.0);
        CHECK(B[fastlckin::QQQQ][2] > 0.0);
    }
}

// ── FST monotonicity ────────────────────────────────────────────

TEST_CASE(ibd_probs_fst_increases_homozygosity) {
    auto B_low  = fastlckin::compute_ibs_ibd_probs(0.3, 0.0);
    auto B_high = fastlckin::compute_ibs_ibd_probs(0.3, 0.1);
    CHECK(B_high[fastlckin::PPPP][0] > B_low[fastlckin::PPPP][0]);
    CHECK(B_high[fastlckin::PQPQ][0] < B_low[fastlckin::PQPQ][0]);
}

// ── Precompute matrix ───────────────────────────────────────────

TEST_CASE(precompute_matrix_size) {
    std::vector<fastlckin::SNPInfo> snps(10);
    for (int i = 0; i < 10; ++i) snps[i].af = 0.1 * (i + 1);
    auto mat = fastlckin::precompute_ibs_ibd_matrix(snps, 0.01);
    CHECK(mat.size() == 9);
    CHECK(mat[0].size() == 10);
}

TEST_CASE(precompute_matrix_normalization) {
    std::vector<fastlckin::SNPInfo> snps(5);
    for (int i = 0; i < 5; ++i) snps[i].af = 0.1 * (i + 1);
    auto mat = fastlckin::precompute_ibs_ibd_matrix(snps, 0.02);
    for (size_t s = 0; s < snps.size(); ++s) {
        for (int z = 0; z < 3; ++z) {
            double sum = 0.0;
            for (int c = 0; c < 9; ++c) sum += mat[c][s][z];
            CHECK_NEAR(sum, 1.0, 1e-10);
        }
    }
}

int main() {
    return RUN_ALL_TESTS();
}
