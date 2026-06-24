/**
 * @file test_pairs_filter.cpp
 * @brief Unit tests for v0.4.0 --pairs parameter support
 *
 * Tests the ability to estimate only specific sample pairs
 * instead of all O(n^2) pairs.
 */

#include "test_harness.h"
#include "kinship_estimator.h"
#include <fstream>
#include <cstdio>
#include <random>
#include <chrono>

using namespace fastlckin;

static const std::string DATA_DIR = "tests/data";
static const std::string PAIRS_FILE = "tests/data/test_pairs.tsv";

// Helper: clean up temp files
static void cleanup_pairs_file() {
    std::remove(PAIRS_FILE.c_str());
}

// Helper: write a pairs file
static void write_pairs_file(const std::vector<std::pair<std::string, std::string>>& pairs) {
    std::ofstream ofs(PAIRS_FILE);
    for (const auto& p : pairs) {
        ofs << p.first << "\t" << p.second << "\n";
    }
}

// Helper: create a basic PLINK-only config
static KinshipConfig make_plink_config() {
    KinshipConfig config;
    config.plink_prefix = DATA_DIR + "/test";
    config.output_path = DATA_DIR + "/test_pairs_output.tsv";
    config.input_mode = InputMode::PLINK_ONLY;
    config.ld_config.skip = true;  // Skip LD pruning for speed
    config.classify = false;
    config.verbose = false;
    config.threads = 1;
    return config;
}

// ── Test 1: Empty pairs file returns empty results ──
TEST_CASE(pairs_empty_file) {
    cleanup_pairs_file();

    // Write empty pairs file (only comments)
    {
        std::ofstream ofs(PAIRS_FILE);
        ofs << "# This is a comment\n";
        ofs << "\n";  // Empty line
    }

    KinshipConfig config = make_plink_config();
    config.pairs_path = PAIRS_FILE;

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    CHECK(results.empty());

    // Clean up
    cleanup_pairs_file();
    std::remove(config.output_path.c_str());
}

// ── Test 2: Single pair estimation ──
TEST_CASE(pairs_single_pair) {
    cleanup_pairs_file();

    // Request just one pair
    write_pairs_file({{"Sample1", "Sample2"}});

    KinshipConfig config = make_plink_config();
    config.pairs_path = PAIRS_FILE;

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    CHECK(results.size() == 1);
    if (!results.empty()) {
        // Check sample names are present (order may vary due to sorting)
        bool has_pair = (results[0].ind1 == "Sample1" && results[0].ind2 == "Sample2") ||
                        (results[0].ind1 == "Sample2" && results[0].ind2 == "Sample1");
        CHECK(has_pair);
        CHECK(!results[0].failed);
        CHECK(results[0].pi_hat >= -0.1 && results[0].pi_hat <= 1.1);
    }

    cleanup_pairs_file();
    std::remove(config.output_path.c_str());
}

// ── Test 3: Multiple pairs ──
TEST_CASE(pairs_multiple_pairs) {
    cleanup_pairs_file();

    write_pairs_file({
        {"Sample1", "Sample2"},
        {"Sample3", "Sample5"},
        {"Sample10", "Sample15"}
    });

    KinshipConfig config = make_plink_config();
    config.pairs_path = PAIRS_FILE;

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    CHECK(results.size() == 3);
    for (const auto& r : results) {
        CHECK(!r.failed);
        CHECK(r.n_snps > 0);
    }

    cleanup_pairs_file();
    std::remove(config.output_path.c_str());
}

// ── Test 4: Invalid sample names are skipped gracefully ──
TEST_CASE(pairs_invalid_sample_names) {
    cleanup_pairs_file();

    write_pairs_file({
        {"Sample1", "Sample2"},       // Valid
        {"UnknownA", "UnknownB"},      // Both unknown
        {"Sample3", "NonExistent"}    // One unknown
    });

    KinshipConfig config = make_plink_config();
    config.pairs_path = PAIRS_FILE;
    config.verbose = true;  // To see warnings

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    // Only the valid pair should be estimated
    CHECK(results.size() == 1);
    if (!results.empty()) {
        bool has_pair = (results[0].ind1 == "Sample1" && results[0].ind2 == "Sample2") ||
                        (results[0].ind1 == "Sample2" && results[0].ind2 == "Sample1");
        CHECK(has_pair);
    }

    cleanup_pairs_file();
    std::remove(config.output_path.c_str());
}

// ── Test 5: Subset vs full estimation consistency ──
TEST_CASE(pairs_subset_consistency) {
    cleanup_pairs_file();

    // First: run full estimation (no pairs file)
    KinshipConfig config_full = make_plink_config();
    config_full.output_path = DATA_DIR + "/test_full_output.tsv";

    KinshipEstimator est_full(config_full);
    auto results_full = est_full.run();

    // Now: request a subset that includes some of those pairs
    write_pairs_file({
        {"Sample1", "Sample2"},
        {"Sample1", "Sample3"}
    });

    KinshipConfig config_sub = make_plink_config();
    config_sub.pairs_path = PAIRS_FILE;
    config_sub.output_path = DATA_DIR + "/test_sub_output.tsv";

    KinshipEstimator est_sub(config_sub);
    auto results_sub = est_sub.run();

    CHECK(results_sub.size() == 2);

    // The subset results should match the corresponding full results
    for (const auto& sub_r : results_sub) {
        bool found = false;
        for (const auto& full_r : results_full) {
            if (full_r.ind1 == sub_r.ind1 && full_r.ind2 == sub_r.ind2) {
                found = true;
                CHECK_NEAR(full_r.k0, sub_r.k0, 1e-6);
                CHECK_NEAR(full_r.k1, sub_r.k1, 1e-6);
                CHECK_NEAR(full_r.k2, sub_r.k2, 1e-6);
                CHECK_NEAR(full_r.pi_hat, sub_r.pi_hat, 1e-6);
                break;
            }
        }
        CHECK(found);
    }

    cleanup_pairs_file();
    std::remove(config_full.output_path.c_str());
    std::remove(config_sub.output_path.c_str());
}

// ── Test 6: Reversed pair order (Sample2, Sample1) still works ──
TEST_CASE(pairs_reversed_order) {
    cleanup_pairs_file();

    // Write pair in reversed order
    write_pairs_file({{"Sample5", "Sample1"}});

    KinshipConfig config = make_plink_config();
    config.pairs_path = PAIRS_FILE;

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    CHECK(results.size() == 1);
    if (!results.empty()) {
        CHECK(!results[0].failed);
        // The pair should be normalized (smaller index first in internal representation)
        // but output uses sample names
        bool valid_names = (results[0].ind1 == "Sample1" && results[0].ind2 == "Sample5") ||
                           (results[0].ind1 == "Sample5" && results[0].ind2 == "Sample1");
        CHECK(valid_names);
    }

    cleanup_pairs_file();
    std::remove(config.output_path.c_str());
}

int main() {
    return RUN_ALL_TESTS();
}
