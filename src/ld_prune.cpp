/**
 * @file ld_prune.cpp
 * @brief C++ built-in LD pruning implementation
 * @author Shujia Huang
 * @date 2025-06-23
 */

#include "ld_prune.h"
#include <cmath>
#include <algorithm>
#include <numeric>
#include <unordered_set>

namespace fastlckin {

double compute_r2(const std::vector<int8_t>& g1, const std::vector<int8_t>& g2) {
    int n = static_cast<int>(g1.size());
    if (n != static_cast<int>(g2.size()) || n < 5) return 0.0;

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0;
    double sum_x2 = 0.0, sum_y2 = 0.0;
    int count = 0;

    for (int i = 0; i < n; ++i) {
        if (g1[i] < 0 || g2[i] < 0) continue;
        double x = static_cast<double>(g1[i]);
        double y = static_cast<double>(g2[i]);
        sum_x += x;
        sum_y += y;
        sum_xy += x * y;
        sum_x2 += x * x;
        sum_y2 += y * y;
        ++count;
    }

    if (count < 5) return 0.0;

    double n_d = static_cast<double>(count);
    double mean_x = sum_x / n_d;
    double mean_y = sum_y / n_d;

    double var_x = sum_x2 / n_d - mean_x * mean_x;
    double var_y = sum_y2 / n_d - mean_y * mean_y;

    if (var_x < 1e-12 || var_y < 1e-12) return 0.0;

    double cov = sum_xy / n_d - mean_x * mean_y;
    double r = cov / std::sqrt(var_x * var_y);

    return r * r;
}

// Extract column j from genotypes[sample][snp] as a vector
static std::vector<int8_t> get_snp_column(
    const std::vector<std::vector<int8_t>>& genotypes, int col)
{
    int n_samples = static_cast<int>(genotypes.size());
    std::vector<int8_t> result(n_samples);
    for (int i = 0; i < n_samples; ++i) {
        result[i] = genotypes[i][col];
    }
    return result;
}

std::vector<int> ld_prune(
    const std::vector<std::vector<int8_t>>& genotypes,
    const std::vector<bool>& mask,
    const LDPruneConfig& config)
{
    int n_snps = static_cast<int>(mask.size());
    if (n_snps == 0) return {};

    // Build list of candidate SNP indices (unmasked)
    std::vector<int> candidates;
    for (int i = 0; i < n_snps; ++i) {
        if (!mask[i]) candidates.push_back(i);
    }

    if (candidates.empty()) return {};

    // Track which candidates are excluded
    std::unordered_set<int> excluded;

    int n_cand = static_cast<int>(candidates.size());
    int window = config.window_size;
    int step = config.step_size;

    // Sliding window pruning
    for (int start = 0; start < n_cand; start += step) {
        int end = std::min(start + window, n_cand);

        // Collect active SNPs in this window
        std::vector<int> active;
        for (int i = start; i < end; ++i) {
            if (excluded.find(candidates[i]) == excluded.end()) {
                active.push_back(candidates[i]);
            }
        }

        if (active.size() < 2) continue;

        // Cache genotype columns for this window
        std::vector<std::vector<int8_t>> cols(active.size());
        for (size_t i = 0; i < active.size(); ++i) {
            cols[i] = get_snp_column(genotypes, active[i]);
        }

        // Iteratively remove high-LD SNPs
        bool changed = true;
        while (changed) {
            changed = false;
            int n_active = static_cast<int>(active.size());

            // Count high-LD partners for each SNP
            std::vector<int> ld_count(n_active, 0);
            std::vector<std::pair<int, int>> high_ld_pairs;

            for (int i = 0; i < n_active - 1; ++i) {
                if (excluded.find(active[i]) != excluded.end()) continue;
                for (int j = i + 1; j < n_active; ++j) {
                    if (excluded.find(active[j]) != excluded.end()) continue;

                    double r2 = compute_r2(cols[i], cols[j]);
                    if (r2 > config.r2_threshold) {
                        ld_count[i]++;
                        ld_count[j]++;
                        high_ld_pairs.push_back({i, j});
                    }
                }
            }

            // Remove the SNP with most high-LD partners
            if (!high_ld_pairs.empty()) {
                int worst = -1;
                int worst_count = 0;
                for (int i = 0; i < n_active; ++i) {
                    if (excluded.find(active[i]) != excluded.end()) continue;
                    if (ld_count[i] > worst_count) {
                        worst_count = ld_count[i];
                        worst = i;
                    }
                }
                if (worst >= 0) {
                    excluded.insert(active[worst]);
                    changed = true;
                }
            }
        }
    }

    // Collect retained SNPs
    std::vector<int> retained;
    for (int idx : candidates) {
        if (excluded.find(idx) == excluded.end()) {
            retained.push_back(idx);
        }
    }

    return retained;
}

double compute_r2_expected(const std::vector<double>& g1, const std::vector<double>& g2) {
    int n = static_cast<int>(g1.size());
    if (n != static_cast<int>(g2.size()) || n < 5) return 0.0;

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0;
    double sum_x2 = 0.0, sum_y2 = 0.0;
    int count = 0;

    for (int i = 0; i < n; ++i) {
        if (g1[i] < -0.5 || g2[i] < -0.5) continue;  // -1.0 = missing
        sum_x += g1[i];
        sum_y += g2[i];
        sum_xy += g1[i] * g2[i];
        sum_x2 += g1[i] * g1[i];
        sum_y2 += g2[i] * g2[i];
        ++count;
    }

    if (count < 5) return 0.0;

    double n_d = static_cast<double>(count);
    double mean_x = sum_x / n_d;
    double mean_y = sum_y / n_d;

    double var_x = sum_x2 / n_d - mean_x * mean_x;
    double var_y = sum_y2 / n_d - mean_y * mean_y;

    if (var_x < 1e-12 || var_y < 1e-12) return 0.0;

    double cov = sum_xy / n_d - mean_x * mean_y;
    double r = cov / std::sqrt(var_x * var_y);

    return r * r;
}

// Extract column j from expected genotypes[sample][snp] as a double vector
static std::vector<double> get_snp_column_double(
    const std::vector<std::vector<double>>& genotypes, int col)
{
    int n_samples = static_cast<int>(genotypes.size());
    std::vector<double> result(n_samples);
    for (int i = 0; i < n_samples; ++i) {
        result[i] = genotypes[i][col];
    }
    return result;
}

std::vector<int> ld_prune_from_gl(
    const std::vector<std::vector<double>>& expected_g,
    const std::vector<bool>& mask,
    const LDPruneConfig& config)
{
    int n_snps = static_cast<int>(mask.size());
    if (n_snps == 0) return {};

    // Build list of candidate SNP indices (unmasked)
    std::vector<int> candidates;
    for (int i = 0; i < n_snps; ++i) {
        if (!mask[i]) candidates.push_back(i);
    }

    if (candidates.empty()) return {};

    // Track which candidates are excluded
    std::unordered_set<int> excluded;

    int n_cand = static_cast<int>(candidates.size());
    int window = config.window_size;
    int step = config.step_size;

    // Sliding window pruning (same logic as ld_prune, but with expected genotypes)
    for (int start = 0; start < n_cand; start += step) {
        int end = std::min(start + window, n_cand);

        std::vector<int> active;
        for (int i = start; i < end; ++i) {
            if (excluded.find(candidates[i]) == excluded.end()) {
                active.push_back(candidates[i]);
            }
        }

        if (active.size() < 2) continue;

        // Cache genotype columns for this window
        std::vector<std::vector<double>> cols(active.size());
        for (size_t i = 0; i < active.size(); ++i) {
            cols[i] = get_snp_column_double(expected_g, active[i]);
        }

        // Iteratively remove high-LD SNPs
        bool changed = true;
        while (changed) {
            changed = false;
            int n_active = static_cast<int>(active.size());

            std::vector<int> ld_count(n_active, 0);
            std::vector<std::pair<int, int>> high_ld_pairs;

            for (int i = 0; i < n_active - 1; ++i) {
                if (excluded.find(active[i]) != excluded.end()) continue;
                for (int j = i + 1; j < n_active; ++j) {
                    if (excluded.find(active[j]) != excluded.end()) continue;

                    double r2 = compute_r2_expected(cols[i], cols[j]);
                    if (r2 > config.r2_threshold) {
                        ld_count[i]++;
                        ld_count[j]++;
                        high_ld_pairs.push_back({i, j});
                    }
                }
            }

            // Remove the SNP with most high-LD partners
            if (!high_ld_pairs.empty()) {
                int worst = -1;
                int worst_count = 0;
                for (int i = 0; i < n_active; ++i) {
                    if (excluded.find(active[i]) != excluded.end()) continue;
                    if (ld_count[i] > worst_count) {
                        worst_count = ld_count[i];
                        worst = i;
                    }
                }
                if (worst >= 0) {
                    excluded.insert(active[worst]);
                    changed = true;
                }
            }
        }
    }

    // Collect retained SNPs
    std::vector<int> retained;
    for (int idx : candidates) {
        if (excluded.find(idx) == excluded.end()) {
            retained.push_back(idx);
        }
    }

    return retained;
}

}  // namespace fastlckin
