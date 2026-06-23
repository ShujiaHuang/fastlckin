#ifndef _FASTLCKIN_FREQUENCY_FROM_LIKELIHOODS_H_
#define _FASTLCKIN_FREQUENCY_FROM_LIKELIHOODS_H_

/**
 * @file frequency_from_likelihoods.h
 * @brief Allele frequency estimation from genotype likelihoods via EM algorithm
 *
 * For VCF-only mode (Mode 1): when no PLINK files are available, allele
 * frequencies are estimated directly from genotype likelihoods (PL/GL)
 * using an Expectation-Maximization (EM) algorithm with a Hardy-Weinberg
 * prior.
 *
 * EM Algorithm:
 *   E-step: P(G=g | D_i, p) = P(D_i|G=g) * P(G=g|p) / Z_i
 *           where P(G=0|p)=(1-p)^2, P(G=1|p)=2p(1-p), P(G=2|p)=p^2
 *   M-step: p = (1/2N) * SUM_i (P(G=1|D_i,p) + 2*P(G=2|D_i,p))
 *
 * Convergence: |p_new - p_old| < tol  (default tol=1e-6)
 */

#include <vector>
#include "genotype_likelihood.h"

namespace fastlckin {

/// Estimate allele frequencies from genotype likelihoods using EM.
///
/// @param lk_matrix    [sample][snp] GenotypeLikelihood (linear-scale P(D|G))
/// @param max_iter     Maximum EM iterations per SNP (default: 100)
/// @param tol          Convergence threshold |p_new - p_old| (default: 1e-6)
/// @param verbose      Print convergence info (default: false)
/// @return Vector of alt allele frequencies, one per SNP
std::vector<double> compute_af_from_likelihoods(
    const LikelihoodMatrix& lk_matrix,
    int max_iter = 100,
    double tol = 1e-6,
    bool verbose = false
);

/// Compute posterior expected genotypes given likelihoods and allele frequencies.
///
/// For each sample i at SNP s:
///   E[G_is] = P(G=1|D_i,p) + 2 * P(G=2|D_i,p)
///
/// where the posterior uses the Hardy-Weinberg prior with the given AF.
/// Masked samples get expected_g[i][s] = -1.0 (missing sentinel).
///
/// @param lk_matrix    [sample][snp] GenotypeLikelihood
/// @param afs          Alt allele frequency for each SNP
/// @return Expected genotype matrix [sample][snp], -1.0 for masked
std::vector<std::vector<double>> compute_expected_genotypes(
    const LikelihoodMatrix& lk_matrix,
    const std::vector<double>& afs
);

}  // namespace fastlckin

#endif
