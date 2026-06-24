/**
 * @file test_v040_integration.cpp
 * @brief Integration tests for v0.4.0 features
 *
 * Verifies all new features work together and maintain backward compatibility.
 */

#include "test_harness.h"
#include "kinship_estimator.h"
#include "ibd_model.h"
#include <fstream>
#include <cstdio>
#include <cmath>
#include <random>

using namespace fastlckin;

static const std::string DATA_DIR = "tests/data";

static KinshipConfig make_plink_config() {
    KinshipConfig config;
    config.plink_prefix = DATA_DIR + "/test";
    config.output_path = DATA_DIR + "/test_v040_integration.tsv";
    config.input_mode = InputMode::PLINK_ONLY;
    config.ld_config.skip = true;
    config.classify = false;
    config.verbose = false;
    config.threads = 1;
    config.pairs_path = "";
    return config;
}

static void write_pairs_file(const std::string& path,
    const std::vector<std::pair<std::string, std::string>>& pairs) {
    std::ofstream ofs(path);
    for (const auto& p : pairs) {
        ofs << p.first << "\t" << p.second << "\n";
    }
}

static void write_fst_file(const std::string& path, double uniform_fst) {
    std::ofstream ofs(path);
    ofs << "CHR\tPOS\tFST\n";
    for (int s = 0; s < 50; ++s) {
        int pos = 10000 + s * 1000;
        ofs << "chr1\t" << pos << "\t" << uniform_fst << "\n";
    }
}

// ── Test 1: --pairs end-to-end ──
TEST_CASE(integration_pairs_feature) {
    std::string pairs_file = DATA_DIR + "/test_int_pairs.tsv";
    write_pairs_file(pairs_file, {
        {"Sample1", "Sample2"},
        {"Sample3", "Sample4"},
        {"Sample5", "Sample10"}
    });

    KinshipConfig config = make_plink_config();
    config.pairs_path = pairs_file;
    config.classify = true;

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    CHECK(results.size() == 3);
    for (const auto& r : results) {
        CHECK(!r.failed);
        CHECK(r.n_snps > 0);
        CHECK(!r.relationship.empty());
    }

    std::remove(pairs_file.c_str());
    std::remove(config.output_path.c_str());
}

// ── Test 2: Custom classification end-to-end ──
TEST_CASE(integration_custom_classification) {
    KinshipConfig config = make_plink_config();
    config.classify = true;

    // Use lenient thresholds
    config.classify_config.duplicate_threshold = 0.3;
    config.classify_config.first_degree_threshold = 0.1;
    config.classify_config.second_degree_threshold = 0.05;
    config.classify_config.third_degree_threshold = 0.01;
    config.classify_config.duplicate_k2_threshold = 0.3;
    config.classify_config.use_custom = true;

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    CHECK(!results.empty());
    // With lenient thresholds, many pairs should be classified as related
    int n_classified = 0;
    for (const auto& r : results) {
        if (!r.failed && r.relationship != "Unrelated") {
            ++n_classified;
        }
    }
    CHECK(n_classified > 0);  // At least some should be classified

    std::remove(config.output_path.c_str());
}

// ── Test 3: FST vector end-to-end ──
TEST_CASE(integration_fst_vector) {
    std::string fst_file = DATA_DIR + "/test_int_fst.tsv";
    write_fst_file(fst_file, 0.01);  // Uniform FST = 0.01

    KinshipConfig config = make_plink_config();
    config.fst = 0.01;
    config.fst_path = fst_file;

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    CHECK(!results.empty());
    for (const auto& r : results) {
        CHECK(!r.failed);
    }

    std::remove(fst_file.c_str());
    std::remove(config.output_path.c_str());
}

// ── Test 4: KING screening end-to-end ──
TEST_CASE(integration_king_screening) {
    KinshipConfig config = make_plink_config();
    config.screening_config.enable_screening = true;
    config.screening_config.pi_hat_threshold = 0.0442;
    config.screening_config.verbose = true;

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    // Results should be a subset of full estimation
    CHECK(results.size() <= 190);

    std::remove(config.output_path.c_str());
}

// ── Test 5: Combined features (--pairs + --screen) ──
TEST_CASE(integration_combined_features) {
    std::string pairs_file = DATA_DIR + "/test_int_combo_pairs.tsv";
    // Request many pairs
    std::vector<std::pair<std::string, std::string>> pairs;
    for (int i = 1; i <= 10; ++i) {
        for (int j = i + 1; j <= 10; ++j) {
            pairs.push_back({"Sample" + std::to_string(i),
                             "Sample" + std::to_string(j)});
        }
    }
    write_pairs_file(pairs_file, pairs);

    KinshipConfig config = make_plink_config();
    config.pairs_path = pairs_file;
    config.classify = true;
    config.screening_config.enable_screening = true;
    config.screening_config.pi_hat_threshold = 0.0;

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    // Should have some results (might be filtered by screening)
    CHECK(results.size() <= 45);  // C(10,2) = 45
    for (const auto& r : results) {
        CHECK(!r.failed);
        CHECK(!r.relationship.empty());
    }

    std::remove(pairs_file.c_str());
    std::remove(config.output_path.c_str());
}

// ── Test 6: Regression test - v0.4.0 without new features matches v0.3.0 ──
TEST_CASE(integration_backward_compatibility) {
    // Run with default config (no new v0.4.0 features enabled)
    KinshipConfig config = make_plink_config();
    // Ensure all v0.4.0 features are disabled
    config.pairs_path = "";
    config.fst_path = "";
    config.classify_config.use_custom = false;
    config.screening_config.enable_screening = false;

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    // Should produce 190 pairs (C(20,2))
    CHECK(results.size() == 190);

    // All results should be valid
    int n_valid = 0;
    for (const auto& r : results) {
        if (!r.failed && r.n_snps > 0) {
            ++n_valid;
            // IBD coefficients should sum to ~1
            CHECK_NEAR(r.k0 + r.k1 + r.k2, 1.0, 0.01);
            // PI_HAT should be in valid range
            CHECK(r.pi_hat >= -0.1 && r.pi_hat <= 1.1);
        }
    }
    CHECK(n_valid >= 180);  // At least 180/190 valid

    std::remove(config.output_path.c_str());
}

// ── Test 7: IBS|IBD normalization with FST vector ──
TEST_CASE(integration_fst_vector_normalization) {
    std::vector<SNPInfo> snps(20);
    std::vector<double> fst_values(20);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> af_dist(0.05, 0.95);
    std::uniform_real_distribution<double> fst_dist(0.0, 0.1);

    for (int i = 0; i < 20; ++i) {
        snps[i].af = af_dist(rng);
        fst_values[i] = fst_dist(rng);
    }

    auto mat = precompute_ibs_ibd_matrix_fst_vector(snps, fst_values);

    // Verify normalization for each SNP
    for (int s = 0; s < 20; ++s) {
        for (int z = 0; z < 3; ++z) {
            double sum = 0.0;
            for (int c = 0; c < 9; ++c) {
                sum += mat[c][s][z];
            }
            CHECK_NEAR(sum, 1.0, 1e-10);
        }
    }
}

int main() {
    return RUN_ALL_TESTS();
}
