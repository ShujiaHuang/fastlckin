/**
 * @file frequency_from_likelihoods.cpp
 * @brief EM algorithm for allele frequency estimation from genotype likelihoods
 */

#include "frequency_from_likelihoods.h"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace fastlckin {

std::vector<double> compute_af_from_likelihoods(
    const LikelihoodMatrix& lk_matrix,
    int max_iter,
    double tol,
    bool verbose)
{
    if (lk_matrix.empty()) return {};

    int n_samples = static_cast<int>(lk_matrix.size());
    int n_snps = static_cast<int>(lk_matrix[0].size());

    std::vector<double> afs(n_snps, 0.0);

    for (int s = 0; s < n_snps; ++s) {
        // Initialize AF at 0.5 (non-informative)
        double p = 0.5;

        for (int iter = 0; iter < max_iter; ++iter) {
            // E-step + M-step combined
            double expected_dosage_sum = 0.0;
            int n_valid = 0;

            double q = 1.0 - p;

            // Hardy-Weinberg prior
            double prior[3] = {q * q, 2.0 * p * q, p * p};

            for (int i = 0; i < n_samples; ++i) {
                const auto& lk = lk_matrix[i][s];
                if (lk.masked) continue;

                // Posterior ∝ likelihood × prior
                double post0 = lk.gl[0] * prior[0];
                double post1 = lk.gl[1] * prior[1];
                double post2 = lk.gl[2] * prior[2];
                double total = post0 + post1 + post2;

                if (total <= 0.0) continue;

                // Normalize and accumulate expected dosage
                expected_dosage_sum += (post1 + 2.0 * post2) / total;
                ++n_valid;
            }

            double p_new;
            if (n_valid == 0) {
                p_new = 0.0;
            } else {
                p_new = expected_dosage_sum / (2.0 * n_valid);
            }

            // Clamp to avoid boundary issues
            p_new = std::max(1e-10, std::min(1.0 - 1e-10, p_new));

            if (std::abs(p_new - p) < tol) {
                p = p_new;
                break;
            }

            p = p_new;
        }

        afs[s] = p;
    }

    if (verbose) {
        std::cerr << "[fastlckin] EM AF estimation: " << n_snps << " SNPs processed\n";
    }

    return afs;
}

std::vector<std::vector<double>> compute_expected_genotypes(
    const LikelihoodMatrix& lk_matrix,
    const std::vector<double>& afs)
{
    if (lk_matrix.empty()) return {};

    int n_samples = static_cast<int>(lk_matrix.size());
    int n_snps = static_cast<int>(lk_matrix[0].size());

    std::vector<std::vector<double>> expected_g(
        n_samples, std::vector<double>(n_snps, -1.0));

    for (int s = 0; s < n_snps; ++s) {
        double p = afs[s];
        double q = 1.0 - p;
        double prior[3] = {q * q, 2.0 * p * q, p * p};

        for (int i = 0; i < n_samples; ++i) {
            const auto& lk = lk_matrix[i][s];
            if (lk.masked) continue;

            double post0 = lk.gl[0] * prior[0];
            double post1 = lk.gl[1] * prior[1];
            double post2 = lk.gl[2] * prior[2];
            double total = post0 + post1 + post2;

            if (total <= 0.0) continue;

            expected_g[i][s] = (post1 + 2.0 * post2) / total;
        }
    }

    return expected_g;
}

}  // namespace fastlckin
