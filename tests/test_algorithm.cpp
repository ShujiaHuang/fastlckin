/**
 * @file test_algorithm.cpp
 * @brief Unit tests for algorithm.h (random_ibd_start)
 */

#include "test_harness.h"
#include "algorithm.h"
#include <random>

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
