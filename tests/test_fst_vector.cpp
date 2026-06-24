/**
 * @file test_fst_vector.cpp
 * @brief Unit tests for v0.4.0 per-SNP FST vector support
 */

#include "test_harness.h"
#include "ibd_model.h"
#include "kinship_estimator.h"
#include <fstream>
#include <cstdio>
#include <cmath>

using namespace fastlckin;

static const std::string FST_FILE = "tests/data/test_fst.tsv";

static void cleanup_fst_file() {
    std::remove(FST_FILE.c_str());
}

// ── Test 1: Uniform FST vector equals global FST ──
TEST_CASE(fst_vector_uniform_equals_global) {
    std::vector<SNPInfo> snps(10);
    for (int i = 0; i < 10; ++i) {
        snps[i].af = 0.1 + i * 0.08;
    }

    double global_fst = 0.02;
    std::vector<double> fst_vector(10, global_fst);

    auto mat_global = precompute_ibs_ibd_matrix(snps, global_fst);
    auto mat_vector = precompute_ibs_ibd_matrix_fst_vector(snps, fst_vector);

    // Should be identical
    CHECK(mat_global.size() == mat_vector.size());
    for (size_t c = 0; c < 9; ++c) {
        CHECK(mat_global[c].size() == mat_vector[c].size());
        for (size_t s = 0; s < snps.size(); ++s) {
            for (int z = 0; z < 3; ++z) {
                CHECK_NEAR(mat_global[c][s][z], mat_vector[c][s][z], 1e-12);
            }
        }
    }
}

// ── Test 2: Different FST values produce different probabilities ──
TEST_CASE(fst_vector_different_values) {
    std::vector<SNPInfo> snps(3);
    snps[0].af = 0.3;
    snps[1].af = 0.3;
    snps[2].af = 0.3;

    std::vector<double> fst_vector = {0.0, 0.05, 0.1};

    auto mat = precompute_ibs_ibd_matrix_fst_vector(snps, fst_vector);

    // Same AF but different FST → different probabilities
    // IBD=0, PPPP combo should differ across SNPs
    bool any_different = false;
    for (int c = 0; c < 9; ++c) {
        for (int z = 0; z < 3; ++z) {
            if (std::abs(mat[c][0][z] - mat[c][1][z]) > 1e-10) {
                any_different = true;
            }
        }
    }
    CHECK(any_different);
}

// ── Test 3: FST file parsing ──
TEST_CASE(fst_file_parsing) {
    cleanup_fst_file();

    // Create FST file matching test data
    // Test data: chr1, positions 10000, 11000, 12000, ...
    std::ofstream ofs(FST_FILE);
    ofs << "CHR\tPOS\tFST\n";
    for (int s = 0; s < 50; ++s) {
        int pos = 10000 + s * 1000;  // 1-based positions (matching .bim)
        double fst = 0.01 + s * 0.001;
        ofs << "chr1\t" << pos << "\t" << fst << "\n";
    }
    ofs.close();

    KinshipConfig config;
    config.plink_prefix = "tests/data/test";
    config.output_path = "tests/data/test_fst_output.tsv";
    config.input_mode = InputMode::PLINK_ONLY;
    config.ld_config.skip = true;
    config.classify = false;
    config.verbose = true;
    config.threads = 1;
    config.fst_path = FST_FILE;
    config.fst = 0.0;  // Default for unmatched

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    CHECK(!results.empty());
    for (const auto& r : results) {
        CHECK(!r.failed);
    }

    cleanup_fst_file();
    std::remove(config.output_path.c_str());
}

// ── Test 4: FST vector size mismatch handled gracefully ──
TEST_CASE(fst_vector_size_mismatch) {
    std::vector<SNPInfo> snps(5);
    for (int i = 0; i < 5; ++i) snps[i].af = 0.3;

    // Shorter FST vector: only 3 values for 5 SNPs
    std::vector<double> short_fst(3, 0.05);

    // Should not throw; unmatched SNPs use FST=0
    auto mat = precompute_ibs_ibd_matrix_fst_vector(snps, short_fst);
    CHECK(mat.size() == 9);
    CHECK(mat[0].size() == 5);

    // First 3 SNPs use FST=0.05, last 2 use FST=0
    auto probs_fst05 = compute_ibs_ibd_probs(0.3, 0.05);
    auto probs_fst0  = compute_ibs_ibd_probs(0.3, 0.0);

    // SNP 0 should match FST=0.05
    for (int c = 0; c < 9; ++c)
        for (int z = 0; z < 3; ++z)
            CHECK_NEAR(mat[c][0][z], probs_fst05[c][z], 1e-12);

    // SNP 4 should match FST=0
    for (int c = 0; c < 9; ++c)
        for (int z = 0; z < 3; ++z)
            CHECK_NEAR(mat[c][4][z], probs_fst0[c][z], 1e-12);
}

// ── Test 5: FST normalization is preserved for each FST value ──
TEST_CASE(fst_vector_normalization) {
    std::vector<SNPInfo> snps(5);
    std::vector<double> fst_values = {0.0, 0.01, 0.05, 0.1, 0.2};
    for (int i = 0; i < 5; ++i) {
        snps[i].af = 0.3;
    }

    auto mat = precompute_ibs_ibd_matrix_fst_vector(snps, fst_values);

    // Each SNP should have normalized probabilities (sum_c B[c][z] = 1 for each z)
    for (size_t s = 0; s < snps.size(); ++s) {
        for (int z = 0; z < 3; ++z) {
            double sum = 0.0;
            for (int c = 0; c < 9; ++c) {
                sum += mat[c][s][z];
            }
            CHECK_NEAR(sum, 1.0, 1e-10);
        }
    }
}

// ── Test 6: End-to-end: uniform FST file produces same result as global FST ──
TEST_CASE(fst_vector_e2e_uniform) {
    cleanup_fst_file();

    // Create FST file with uniform FST = 0.01 for all SNPs
    std::ofstream ofs(FST_FILE);
    ofs << "CHR\tPOS\tFST\n";
    for (int s = 0; s < 50; ++s) {
        int pos = 10000 + s * 1000;
        ofs << "chr1\t" << pos << "\t0.01\n";
    }
    ofs.close();

    // Run with global FST=0.01
    KinshipConfig config_global;
    config_global.plink_prefix = "tests/data/test";
    config_global.output_path = "tests/data/test_fst_global.tsv";
    config_global.input_mode = InputMode::PLINK_ONLY;
    config_global.ld_config.skip = true;
    config_global.fst = 0.01;
    config_global.verbose = false;
    config_global.threads = 1;
    config_global.pairs_path = "";

    KinshipEstimator est_global(config_global);
    auto results_global = est_global.run();

    // Run with FST file (uniform 0.01)
    KinshipConfig config_vec;
    config_vec.plink_prefix = "tests/data/test";
    config_vec.output_path = "tests/data/test_fst_vec.tsv";
    config_vec.input_mode = InputMode::PLINK_ONLY;
    config_vec.ld_config.skip = true;
    config_vec.fst = 0.01;
    config_vec.fst_path = FST_FILE;
    config_vec.verbose = false;
    config_vec.threads = 1;
    config_vec.pairs_path = "";

    KinshipEstimator est_vec(config_vec);
    auto results_vec = est_vec.run();

    CHECK(results_global.size() == results_vec.size());
    for (size_t i = 0; i < results_global.size(); ++i) {
        CHECK_NEAR(results_global[i].pi_hat, results_vec[i].pi_hat, 1e-4);
        CHECK_NEAR(results_global[i].k0, results_vec[i].k0, 1e-4);
        CHECK_NEAR(results_global[i].k1, results_vec[i].k1, 1e-4);
        CHECK_NEAR(results_global[i].k2, results_vec[i].k2, 1e-4);
    }

    cleanup_fst_file();
    std::remove(config_global.output_path.c_str());
    std::remove(config_vec.output_path.c_str());
}

int main() {
    return RUN_ALL_TESTS();
}
