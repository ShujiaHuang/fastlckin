#ifndef _FASTLCKIN_KINSHIP_ESTIMATOR_H_
#define _FASTLCKIN_KINSHIP_ESTIMATOR_H_

/**
 * @file kinship_estimator.h
 * @brief Main kinship estimation pipeline coordinator
 * @author Shujia Huang
 * @date 2025-06-23
 */

#include <vector>
#include <string>
#include <array>

#include "ibd_model.h"
#include "genotype_likelihood.h"
#include "ld_prune.h"
#include "optimizer.h"

namespace fastlckin {

/// Kinship estimation result for one pair
struct KinshipResult {
    std::string ind1;
    std::string ind2;
    double k0 = 0.0;          ///< IBD=0 probability
    double k1 = 0.0;          ///< IBD=1 probability
    double k2 = 0.0;          ///< IBD=2 probability
    double pi_hat = 0.0;      ///< Kinship coefficient = 0.5*k1 + k2
    int n_snps = 0;           ///< Number of SNPs used
    double log_likelihood = 0.0; ///< Best negative log-likelihood
    std::string relationship; ///< Classification label
    bool failed = false;      ///< Estimation failed flag
};

/// Configuration for the full kinship estimation pipeline
struct KinshipConfig {
    // Input files
    std::string vcf_path;
    std::string plink_prefix;
    std::string freq_path;    ///< Optional pre-computed .frq
    std::string output_path;

    // Algorithm parameters
    double fst = 0.0;
    double maf_min = 0.05;
    double maf_max = 0.95;
    int gq_min = 1;
    int n_restarts = 3;
    bool classify = false;
    bool verbose = false;

    // LD pruning parameters
    LDPruneConfig ld_config;

    // Optimizer parameters
    NelderMeadConfig nm_config;

    // Parallelism
    int threads = 1;
};

/// Main kinship estimator class
class KinshipEstimator {
public:
    explicit KinshipEstimator(const KinshipConfig& config);

    /// Run the full kinship estimation pipeline
    /// @return Results for all sample pairs
    std::vector<KinshipResult> run();

private:
    KinshipConfig _config;

    // Internal data
    std::vector<std::string> _sample_names;
    std::vector<SNPInfo> _snp_infos;          ///< VCF SNPs (filtered biallelic)
    GLMatrix _gl_matrix;                       ///< [sample][snp]
    IBS_IBD_Matrix _ibs_ibd;                   ///< [9][snp][3]
    std::vector<std::vector<int8_t>> _bed_genotypes; ///< [sample][snp] for LD

    // Mapping: VCF SNP index → .bim SNP index
    std::vector<int> _vcf_to_bim_index;

    // Internal steps
    void _load_vcf();
    void _load_frequencies();
    void _precompute_ibs_ibd();
    void _load_bed_genotypes();

    /// Estimate IBD coefficients for a single pair
    KinshipResult _estimate_pair(int ind1, int ind2);

    /// GLkin likelihood function (core objective)
    double _glkin_likelihood(
        const std::vector<double>& k,          ///< [k1, k2]
        const std::array<double, 9>& pibs,     ///< PIBS for this pair
        const std::vector<int>& snp_indices    ///< LD-pruned SNP indices
    );

    /// Classify relationship based on IBD coefficients
    static std::string classify_relationship(double k0, double k1, double k2, double pi_hat);
};

}  // namespace fastlckin

#endif
