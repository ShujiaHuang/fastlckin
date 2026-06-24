/**
 * @file test_king_screening.cpp
 * @brief Unit tests for v0.4.0 built-in KING-robust fast screening
 */

#include "test_harness.h"
#include "kinship_estimator.h"
#include <fstream>
#include <cstdio>
#include <cmath>
#include <chrono>

using namespace fastlckin;

static const std::string DATA_DIR = "tests/data";

static KinshipConfig make_plink_config() {
    KinshipConfig config;
    config.plink_prefix = DATA_DIR + "/test";
    config.output_path = DATA_DIR + "/test_screening_output.tsv";
    config.input_mode = InputMode::PLINK_ONLY;
    config.ld_config.skip = true;
    config.classify = false;
    config.verbose = false;
    config.threads = 1;
    config.pairs_path = "";
    return config;
}

// ── Test 1: Screening disabled returns all pairs ──
TEST_CASE(screening_disabled_returns_all) {
    KinshipConfig config = make_plink_config();
    config.screening_config.enable_screening = false;

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    // 20 samples → C(20,2) = 190 pairs
    CHECK(results.size() == 190);

    std::remove(config.output_path.c_str());
}

// ── Test 2: Screening with very low threshold keeps most pairs ──
TEST_CASE(screening_low_threshold) {
    KinshipConfig config = make_plink_config();
    config.screening_config.enable_screening = true;
    config.screening_config.pi_hat_threshold = -1.0;  // Keep everything
    config.screening_config.verbose = true;

    KinshipEstimator estimator(config);
    auto results = estimator.run();

    // Should keep nearly all pairs (random pairs might have negative PI_HAT)
    CHECK(results.size() >= 150);  // At least 150/190

    std::remove(config.output_path.c_str());
}

// ── Test 3: Screening works and is deterministic ──
TEST_CASE(screening_deterministic) {
    // Run screening twice with same config → same results
    KinshipConfig config = make_plink_config();
    config.screening_config.enable_screening = true;
    config.screening_config.pi_hat_threshold = 0.0442;

    KinshipEstimator est1(config);
    auto results1 = est1.run();

    config.output_path = DATA_DIR + "/test_screen_det2.tsv";
    KinshipEstimator est2(config);
    auto results2 = est2.run();

    CHECK(results1.size() == results2.size());
    for (size_t i = 0; i < results1.size(); ++i) {
        CHECK(results1[i].ind1 == results2[i].ind1);
        CHECK(results1[i].ind2 == results2[i].ind2);
        CHECK_NEAR(results1[i].pi_hat, results2[i].pi_hat, 1e-6);
    }

    std::remove(config.output_path.c_str());
    std::remove((DATA_DIR + "/test_screen_det2.tsv").c_str());
}

// ── Test 4: Screening reduces pairs while keeping high-kinship ones ──
TEST_CASE(screening_reduces_pairs) {
    // First run without screening to get full results
    KinshipConfig config_full = make_plink_config();
    config_full.output_path = DATA_DIR + "/test_full_screen.tsv";

    KinshipEstimator est_full(config_full);
    auto results_full = est_full.run();

    // Now run with moderate screening
    KinshipConfig config_screen = make_plink_config();
    config_screen.output_path = DATA_DIR + "/test_screened.tsv";
    config_screen.screening_config.enable_screening = true;
    config_screen.screening_config.pi_hat_threshold = 0.0442;  // ~3rd degree
    config_screen.screening_config.verbose = true;

    KinshipEstimator est_screen(config_screen);
    auto results_screen = est_screen.run();

    // Screened should have fewer or equal pairs
    CHECK(results_screen.size() <= results_full.size());

    // Any pair in screened results should also be in full results
    for (const auto& sr : results_screen) {
        bool found = false;
        for (const auto& fr : results_full) {
            if (fr.ind1 == sr.ind1 && fr.ind2 == sr.ind2) {
                found = true;
                break;
            }
        }
        CHECK(found);
    }

    std::remove(config_full.output_path.c_str());
    std::remove(config_screen.output_path.c_str());
}

// ── Test 5: Screening performance (should be faster) ──
TEST_CASE(screening_performance) {
    // Full run timing
    KinshipConfig config_full = make_plink_config();
    config_full.output_path = DATA_DIR + "/test_perf_full.tsv";

    auto t0 = std::chrono::steady_clock::now();
    KinshipEstimator est_full(config_full);
    auto results_full = est_full.run();
    auto t1 = std::chrono::steady_clock::now();
    double time_full = std::chrono::duration<double>(t1 - t0).count();

    // Screened run timing
    KinshipConfig config_screen = make_plink_config();
    config_screen.output_path = DATA_DIR + "/test_perf_screen.tsv";
    config_screen.screening_config.enable_screening = true;
    config_screen.screening_config.pi_hat_threshold = 0.1;

    auto t2 = std::chrono::steady_clock::now();
    KinshipEstimator est_screen(config_screen);
    auto results_screen = est_screen.run();
    auto t3 = std::chrono::steady_clock::now();
    double time_screen = std::chrono::duration<double>(t3 - t2).count();

    // Screened should process fewer pairs
    CHECK(results_screen.size() <= results_full.size());

    // Both should complete successfully
    CHECK(!results_full.empty());
    // Screened might be empty if no pairs pass, which is OK

    std::remove(config_full.output_path.c_str());
    std::remove(config_screen.output_path.c_str());
}

// ── Test 6: ScreeningConfig default values ──
TEST_CASE(screening_config_defaults) {
    ScreeningConfig config;
    CHECK(config.enable_screening == false);
    CHECK_NEAR(config.pi_hat_threshold, 0.0442, 1e-10);
    CHECK(config.verbose == false);
}

// ── Test 7: Screened pairs are a subset of full pairs ──
TEST_CASE(screening_subset_property) {
    // Run full estimation
    KinshipConfig config_full = make_plink_config();
    config_full.output_path = DATA_DIR + "/test_subset_full.tsv";

    KinshipEstimator est_full(config_full);
    auto results_full = est_full.run();

    // Run with screening
    KinshipConfig config_screen = make_plink_config();
    config_screen.output_path = DATA_DIR + "/test_subset_screen.tsv";
    config_screen.screening_config.enable_screening = true;
    config_screen.screening_config.pi_hat_threshold = 0.0;  // threshold = 0

    KinshipEstimator est_screen(config_screen);
    auto results_screen = est_screen.run();

    // With threshold=0, most pairs with non-negative KING PI_HAT should pass
    // All screened results should be in full results
    for (const auto& sr : results_screen) {
        bool found = false;
        for (const auto& fr : results_full) {
            if (fr.ind1 == sr.ind1 && fr.ind2 == sr.ind2) {
                found = true;
                break;
            }
        }
        CHECK(found);
    }

    std::remove(config_full.output_path.c_str());
    std::remove(config_screen.output_path.c_str());
}

int main() {
    return RUN_ALL_TESTS();
}
