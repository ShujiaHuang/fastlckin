/**
 * @file algorithm.cpp
 * @brief Mathematical utility functions implementation
 * @author Shujia Huang
 * @date 2026-06-23
 */

#include "algorithm.h"

namespace fastlckin {

std::vector<double> random_ibd_start(std::mt19937& rng) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (int attempt = 0; attempt < 10000; ++attempt) {
        double k1 = dist(rng);
        double k2 = dist(rng) * (1.0 - k1);
        double k0 = 1.0 - k1 - k2;

        if (k0 < 0.0 || k0 > 1.0) continue;
        if (4.0 * k0 * k2 > k1 * k1) continue;

        return {k1, k2};
    }

    // Fallback: unrelated starting point
    return {0.0, 0.0};
}

}  // namespace fastlckin
