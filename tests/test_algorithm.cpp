/**
 * @file test_algorithm.cpp
 * @brief Unit tests for algorithm.h (math utilities)
 */

#include "test_harness.h"
#include "algorithm.h"
#include <random>

TEST_CASE(safe_log_positive) {
    CHECK_NEAR(fastlckin::safe_log(1.0), 0.0, 1e-12);
    CHECK_NEAR(fastlckin::safe_log(std::exp(1.0)), 1.0, 1e-12);
}

TEST_CASE(safe_log_zero) {
    double v = fastlckin::safe_log(0.0);
    CHECK(std::isinf(v) && v < 0);
}

TEST_CASE(safe_log_negative) {
    double v = fastlckin::safe_log(-1.0);
    CHECK(std::isinf(v) && v < 0);
}

TEST_CASE(log_sum_exp_basic) {
    // log(exp(1) + exp(2) + exp(3))
    std::vector<double> vals = {1.0, 2.0, 3.0};
    double expected = std::log(std::exp(1.0) + std::exp(2.0) + std::exp(3.0));
    CHECK_NEAR(fastlckin::log_sum_exp(vals), expected, 1e-10);
}

TEST_CASE(log_sum_exp_single) {
    std::vector<double> vals = {5.0};
    CHECK_NEAR(fastlckin::log_sum_exp(vals), 5.0, 1e-12);
}

TEST_CASE(log_sum_exp_empty) {
    std::vector<double> vals;
    double v = fastlckin::log_sum_exp(vals);
    CHECK(std::isinf(v) && v < 0);
}

TEST_CASE(log_sum_exp_numerical_stability) {
    // Very large values that would overflow naive exp()
    std::vector<double> vals = {1000.0, 1001.0, 1002.0};
    double result = fastlckin::log_sum_exp(vals);
    // Should be ~1002 + log(1 + exp(-1) + exp(-2)) ≈ 1002 + 0.46
    CHECK(result > 1002.0 && result < 1003.0);
}

TEST_CASE(random_ibd_start_valid) {
    std::mt19937 rng(42);
    for (int i = 0; i < 100; ++i) {
        auto x = fastlckin::random_ibd_start(rng);
        CHECK(x.size() == 2);
        double k1 = x[0], k2 = x[1];
        double k0 = 1.0 - k1 - k2;
        CHECK(k0 >= -1e-10);
        CHECK(k0 <= 1.0 + 1e-10);
        CHECK(k1 >= -1e-10);
        CHECK(k2 >= -1e-10);
        CHECK(4.0 * k0 * k2 <= k1 * k1 + 1e-6);
    }
}

int main() {
    return RUN_ALL_TESTS();
}
