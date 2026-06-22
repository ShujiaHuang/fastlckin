#ifndef _FASTLCKIN_ALGORITHM_H_
#define _FASTLCKIN_ALGORITHM_H_

/**
 * @file algorithm.h
 * @brief Mathematical utility functions for kinship estimation
 * @author Shujia Huang
 * @date 2025-06-23
 */

#include <cmath>
#include <vector>
#include <limits>
#include <random>

namespace fastlckin {

/// Safe log: avoids log(0) by returning -infinity
inline double safe_log(double x) {
    return (x > 0.0) ? std::log(x) : -std::numeric_limits<double>::infinity();
}

/// Numerically stable log-sum-exp
double log_sum_exp(const std::vector<double>& values);

/// Generate a random starting point (k1, k2) satisfying IBD constraints:
///   k0 = 1 - k1 - k2 >= 0
///   4 * k0 * k2 <= k1^2
std::vector<double> random_ibd_start(std::mt19937& rng);

}  // namespace fastlckin

#endif
