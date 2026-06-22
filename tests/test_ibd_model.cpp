/**
 * @file test_ibd_model.cpp
 * @brief Unit tests for ibd_model.h (Anderson & Weir probabilities)
 */

#include "test_harness.h"
#include "ibd_model.h"
#include <cmath>

TEST_CASE(Mij_basic) {
    // When FST=0, M(p, 0, i) = p regardless of i
    CHECK_NEAR(fastlckin::Mij(0.3, 0.0, 0), 0.3, 1e-12);
    CHECK_NEAR(fastlckin::Mij(0.3, 0.0, 1), 0.3, 1e-12);
    CHECK_NEAR(fastlckin::Mij(0.3, 0.0, 2), 0.3, 1e-12);
}

TEST_CASE(Mij_fst_correction) {
    // M(p, FST, i) = (1-FST)*p + i*FST
    double p = 0.4, fst = 0.1;
    CHECK_NEAR(fastlckin::Mij(p, fst, 0), (1.0 - fst) * p, 1e-12);
    CHECK_NEAR(fastlckin::Mij(p, fst, 1), (1.0 - fst) * p + fst, 1e-12);
    CHECK_NEAR(fastlckin::Mij(p, fst, 2), (1.0 - fst) * p + 2.0 * fst, 1e-12);
}

TEST_CASE(ibd_probs_sum_to_one_fst0) {
    // With FST=0, for each IBD state z, sum over all 9 combos should be ~1
    auto B = fastlckin::compute_ibs_ibd_probs(0.3, 0.0);

    for (int z = 0; z < 3; ++z) {
        double sum = 0.0;
        for (int c = 0; c < 9; ++c) {
            sum += B[c][z];
        }
        CHECK_NEAR(sum, 1.0, 1e-10);
    }
}

TEST_CASE(ibd_probs_sum_to_one_fst_nonzero) {
    // With FST=0.01, same check
    auto B = fastlckin::compute_ibs_ibd_probs(0.25, 0.01);

    for (int z = 0; z < 3; ++z) {
        double sum = 0.0;
        for (int c = 0; c < 9; ++c) {
            sum += B[c][z];
        }
        CHECK_NEAR(sum, 1.0, 1e-10);
    }
}

TEST_CASE(ibd_probs_ibs0_zero_at_ibd2) {
    // IBS=0 combos (PPQQ, QQPP) must have P(IBS=0|IBD=2) = 0
    auto B = fastlckin::compute_ibs_ibd_probs(0.4, 0.0);
    CHECK_NEAR(B[fastlckin::PPQQ][2], 0.0, 1e-15);
    CHECK_NEAR(B[fastlckin::QQPP][2], 0.0, 1e-15);
}

TEST_CASE(ibd_probs_ibs0_zero_at_ibd1) {
    // IBS=0 combos must have P(IBS=0|IBD=1) = 0
    auto B = fastlckin::compute_ibs_ibd_probs(0.3, 0.0);
    CHECK_NEAR(B[fastlckin::PPQQ][1], 0.0, 1e-15);
    CHECK_NEAR(B[fastlckin::QQPP][1], 0.0, 1e-15);
}

TEST_CASE(ibd_probs_ibs1_zero_at_ibd2) {
    // IBS=1 combos must have P(IBS=1|IBD=2) = 0
    auto B = fastlckin::compute_ibs_ibd_probs(0.3, 0.0);
    CHECK_NEAR(B[fastlckin::PPPQ][2], 0.0, 1e-15);
    CHECK_NEAR(B[fastlckin::PQPP][2], 0.0, 1e-15);
    CHECK_NEAR(B[fastlckin::PQQQ][2], 0.0, 1e-15);
    CHECK_NEAR(B[fastlckin::QQPQ][2], 0.0, 1e-15);
}

TEST_CASE(ibd_probs_symmetry_fst0) {
    // At FST=0, PPQQ and QQPP with swapped p/q should match
    auto B = fastlckin::compute_ibs_ibd_probs(0.3, 0.0);
    // PPQQ with p=0.3 should equal QQPP with same p (since QQPP = pp↔qq swap)
    // Actually PPQQ(af=p) = QQPP(af=1-p)
    auto B2 = fastlckin::compute_ibs_ibd_probs(0.7, 0.0);
    CHECK_NEAR(B[fastlckin::PPQQ][0], B2[fastlckin::QQPP][0], 1e-10);
}

TEST_CASE(ibd_probs_nonnegative) {
    // All probabilities must be non-negative
    for (double af : {0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95}) {
        auto B = fastlckin::compute_ibs_ibd_probs(af, 0.01);
        for (int c = 0; c < 9; ++c) {
            for (int z = 0; z < 3; ++z) {
                CHECK(B[c][z] >= -1e-15);
            }
        }
    }
}

TEST_CASE(precompute_matrix_size) {
    std::vector<fastlckin::SNPInfo> snps(10);
    for (int i = 0; i < 10; ++i) {
        snps[i].af = 0.1 * (i + 1);
    }
    auto mat = fastlckin::precompute_ibs_ibd_matrix(snps, 0.01);
    CHECK(mat.size() == 9);
    CHECK(mat[0].size() == 10);
}

int main() {
    return RUN_ALL_TESTS();
}
