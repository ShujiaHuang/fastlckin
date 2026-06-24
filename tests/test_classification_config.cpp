/**
 * @file test_classification_config.cpp
 * @brief Unit tests for v0.4.0 user-customizable classification thresholds
 */

#include "test_harness.h"
#include "kinship_estimator.h"

using namespace fastlckin;

// ── Test 1: Default thresholds match v0.3.0 behavior ──
TEST_CASE(classify_default_thresholds) {
    // Typical Duplicate/MZ: pi_hat=0.9, k2=0.9
    CHECK(KinshipEstimator::classify_relationship(0.01, 0.09, 0.9, 0.9) == "Duplicate/MZ");
    
    // Typical Parent-Offspring: pi_hat=0.5, k0=0.0, k1=1.0
    CHECK(KinshipEstimator::classify_relationship(0.0, 1.0, 0.0, 0.5) == "Parent-Offspring");
    
    // Typical Full-Sibling: pi_hat=0.5, k0=0.25, k1=0.5, k2=0.25
    CHECK(KinshipEstimator::classify_relationship(0.25, 0.5, 0.25, 0.5) == "Full-Sibling");
    
    // Typical Second-degree: pi_hat=0.25
    CHECK(KinshipEstimator::classify_relationship(0.5, 0.5, 0.0, 0.25) == "Second-degree");
    
    // Typical Third-degree: pi_hat=0.125
    CHECK(KinshipEstimator::classify_relationship(0.75, 0.25, 0.0, 0.125) == "Third-degree");
    
    // Unrelated: pi_hat=0.01
    CHECK(KinshipEstimator::classify_relationship(0.98, 0.02, 0.0, 0.01) == "Unrelated");
}

// ── Test 2: Custom thresholds work correctly ──
TEST_CASE(classify_custom_thresholds) {
    ClassificationConfig config;
    config.duplicate_threshold = 0.8;
    config.first_degree_threshold = 0.4;
    config.second_degree_threshold = 0.2;
    config.third_degree_threshold = 0.1;
    config.duplicate_k2_threshold = 0.85;
    config.use_custom = true;
    
    // PI_HAT=0.9, k2=0.9 → Duplicate (above 0.8, k2 > 0.85)
    CHECK(KinshipEstimator::classify_relationship(0.01, 0.09, 0.9, 0.9, config) == "Duplicate/MZ");
    
    // PI_HAT=0.5 → First-degree (in [0.4, 0.8])
    CHECK(KinshipEstimator::classify_relationship(0.0, 1.0, 0.0, 0.5, config) == "Parent-Offspring");
    
    // PI_HAT=0.25 → Second-degree (in [0.2, 0.4))
    CHECK(KinshipEstimator::classify_relationship(0.5, 0.5, 0.0, 0.25, config) == "Second-degree");
    
    // PI_HAT=0.15 → Third-degree (in [0.1, 0.2))
    CHECK(KinshipEstimator::classify_relationship(0.7, 0.3, 0.0, 0.15, config) == "Third-degree");
    
    // PI_HAT=0.05 → Unrelated (below 0.1)
    CHECK(KinshipEstimator::classify_relationship(0.9, 0.1, 0.0, 0.05, config) == "Unrelated");
}

// ── Test 3: Extreme strict thresholds → all Unrelated ──
TEST_CASE(classify_extreme_strict) {
    ClassificationConfig config;
    config.duplicate_threshold = 0.99;
    config.first_degree_threshold = 0.98;
    config.second_degree_threshold = 0.97;
    config.third_degree_threshold = 0.96;
    config.duplicate_k2_threshold = 0.99;
    config.use_custom = true;
    
    // Even closely related pairs classified as Unrelated
    CHECK(KinshipEstimator::classify_relationship(0.0, 1.0, 0.0, 0.5, config) == "Unrelated");
    CHECK(KinshipEstimator::classify_relationship(0.25, 0.5, 0.25, 0.5, config) == "Unrelated");
    CHECK(KinshipEstimator::classify_relationship(0.5, 0.5, 0.0, 0.25, config) == "Unrelated");
    
    // Only extremely high pi_hat gets classified
    CHECK(KinshipEstimator::classify_relationship(0.0, 0.005, 0.995, 0.995, config) == "Duplicate/MZ");
}

// ── Test 4: Extreme lenient thresholds → most are Duplicate/MZ ──
TEST_CASE(classify_extreme_lenient) {
    ClassificationConfig config;
    config.duplicate_threshold = 0.1;
    config.first_degree_threshold = 0.05;
    config.second_degree_threshold = 0.02;
    config.third_degree_threshold = 0.01;
    config.duplicate_k2_threshold = 0.0;
    config.use_custom = true;
    
    // Most pairs now classified as Duplicate/MZ (need k2 > 0)
    CHECK(KinshipEstimator::classify_relationship(0.0, 0.5, 0.5, 0.75, config) == "Duplicate/MZ");
    CHECK(KinshipEstimator::classify_relationship(0.0, 0.5, 0.5, 0.5, config) == "Duplicate/MZ");
    CHECK(KinshipEstimator::classify_relationship(0.0, 0.5, 0.5, 0.25, config) == "Duplicate/MZ");
    
    // Very low pi_hat still Unrelated
    CHECK(KinshipEstimator::classify_relationship(0.98, 0.02, 0.0, 0.005, config) == "Unrelated");
}

// ── Test 5: Default config produces same results as old function ──
TEST_CASE(classify_default_config_equivalence) {
    ClassificationConfig defaults;
    
    // Test several points across the classification space
    struct TestPoint {
        double k0, k1, k2, pi_hat;
    };
    
    std::vector<TestPoint> points = {
        {0.01, 0.09, 0.9, 0.9},     // Duplicate
        {0.0, 1.0, 0.0, 0.5},        // Parent-Offspring
        {0.25, 0.5, 0.25, 0.5},      // Full-Sibling
        {0.5, 0.5, 0.0, 0.25},       // Second-degree
        {0.75, 0.25, 0.0, 0.125},    // Third-degree
        {0.98, 0.02, 0.0, 0.01},     // Unrelated
        {0.02, 0.08, 0.9, 0.94},     // Duplicate
        {0.01, 0.98, 0.01, 0.5},     // Parent-Offspring
    };
    
    for (const auto& tp : points) {
        std::string old_result = KinshipEstimator::classify_relationship(
            tp.k0, tp.k1, tp.k2, tp.pi_hat
        );
        std::string new_result = KinshipEstimator::classify_relationship(
            tp.k0, tp.k1, tp.k2, tp.pi_hat, defaults
        );
        CHECK(old_result == new_result);
    }
}

// ── Test 6: ClassificationConfig default values ──
TEST_CASE(classify_config_defaults) {
    ClassificationConfig config;
    
    CHECK_NEAR(config.duplicate_threshold, 0.708, 1e-10);
    CHECK_NEAR(config.first_degree_threshold, 0.354, 1e-10);
    CHECK_NEAR(config.second_degree_threshold, 0.177, 1e-10);
    CHECK_NEAR(config.third_degree_threshold, 0.0884, 1e-10);
    CHECK_NEAR(config.duplicate_k2_threshold, 0.8, 1e-10);
    CHECK(config.use_custom == false);
}

// ── Test 7: End-to-end with custom thresholds via KinshipEstimator ──
TEST_CASE(classify_end_to_end_custom) {
    KinshipConfig config;
    config.plink_prefix = "tests/data/test";
    config.output_path = "tests/data/test_classify_custom.tsv";
    config.input_mode = InputMode::PLINK_ONLY;
    config.ld_config.skip = true;
    config.classify = true;
    config.verbose = false;
    config.threads = 1;
    config.pairs_path = "";
    
    // Set very strict thresholds → everything should be Unrelated
    config.classify_config.duplicate_threshold = 0.99;
    config.classify_config.first_degree_threshold = 0.98;
    config.classify_config.second_degree_threshold = 0.97;
    config.classify_config.third_degree_threshold = 0.96;
    config.classify_config.use_custom = true;
    
    KinshipEstimator estimator(config);
    auto results = estimator.run();
    
    CHECK(!results.empty());
    for (const auto& r : results) {
        if (!r.failed) {
            CHECK(r.relationship == "Unrelated");
        }
    }
    
    std::remove(config.output_path.c_str());
}

int main() {
    return RUN_ALL_TESTS();
}
