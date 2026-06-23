#ifndef _FASTLCKIN_LD_PRUNE_H_
#define _FASTLCKIN_LD_PRUNE_H_

/**
 * @file ld_prune.h
 * @brief C++ built-in LD pruning (replaces external PLINK calls)
 * @author Shujia Huang
 * @date 2025-06-23
 */

#include <vector>
#include <cstdint>

namespace fastlckin {

struct LDPruneConfig {
    int window_size = 50;       ///< Window size in SNPs
    int step_size = 5;          ///< Step size
    double r2_threshold = 0.5;  ///< r² threshold
    bool skip = false;          ///< If true, skip LD pruning entirely (use all unmasked SNPs)
};

/// Compute r² (squared Pearson correlation) between two SNP genotype vectors
/// @param g1  SNP1 genotypes (0/1/2, -1=missing)
/// @param g2  SNP2 genotypes (0/1/2, -1=missing)
/// @return r² value; returns 0.0 if insufficient valid samples
double compute_r2(const std::vector<int8_t>& g1, const std::vector<int8_t>& g2);

/// Compute r² between two SNP vectors of continuous expected genotypes
/// @param g1  SNP1 expected genotypes (continuous, -1.0=missing)
/// @param g2  SNP2 expected genotypes (continuous, -1.0=missing)
/// @return r² value; returns 0.0 if insufficient valid samples
double compute_r2_expected(const std::vector<double>& g1, const std::vector<double>& g2);

/// Perform LD pruning on a set of SNPs
/// Equivalent to PLINK --indep-pairwise {window} {step} {r2}
/// @param genotypes  Genotype matrix [n_samples][n_snps]
/// @param mask       Initial mask (true = already excluded)
/// @param config     LD pruning parameters
/// @return Indices of retained SNPs (column indices into genotypes)
std::vector<int> ld_prune(
    const std::vector<std::vector<int8_t>>& genotypes,
    const std::vector<bool>& mask,
    const LDPruneConfig& config
);

/// Perform LD pruning using expected (continuous) genotypes from genotype
/// likelihoods (VCF-only mode). Same sliding-window iterative pruning as
/// ld_prune(), but uses Pearson r² computed from posterior expected genotypes.
/// @param expected_g  Expected genotype matrix [n_samples][n_snps] (-1.0=missing)
/// @param mask        Initial mask (true = already excluded)
/// @param config      LD pruning parameters
/// @return Indices of retained SNPs
std::vector<int> ld_prune_from_likelihoods(
    const std::vector<std::vector<double>>& expected_g,
    const std::vector<bool>& mask,
    const LDPruneConfig& config
);

/// Perform global LD pruning across all samples (v0.3.0 new feature).
///
/// This is more efficient than per-pair LD pruning: instead of running LD
/// pruning for each sample pair independently, we run it once using all
/// samples, then reuse the retained SNP set for all pairs.
///
/// @param genotypes     Genotype matrix [n_samples][n_snps] (int8_t or double)
/// @param global_mask   Initial mask (true = already excluded for all samples)
/// @param config        LD pruning parameters
/// @param mode          "hard" for int8_t genotypes, "expected" for double
/// @return Indices of globally retained SNPs
std::vector<int> ld_prune_global(
    const std::vector<std::vector<int8_t>>& genotypes,
    const std::vector<bool>& global_mask,
    const LDPruneConfig& config
);

std::vector<int> ld_prune_global_expected(
    const std::vector<std::vector<double>>& expected_g,
    const std::vector<bool>& global_mask,
    const LDPruneConfig& config
);

}  // namespace fastlckin

#endif
