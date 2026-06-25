#ifndef _FASTLCKIN_ALGORITHM_H_
#define _FASTLCKIN_ALGORITHM_H_

/**
 * @file algorithm.h
 * @brief Mathematical utility functions for kinship estimation
 * @author Shujia Huang
 * @date 2026-06-23
 */

#include <vector>
#include <random>

namespace fastlckin {

/// Generate a random starting point (k1, k2) satisfying IBD constraints:
///   k0 = 1 - k1 - k2 >= 0
///   4 * k0 * k2 <= k1^2
std::vector<double> random_ibd_start(std::mt19937& rng);

}  // namespace fastlckin

#endif
