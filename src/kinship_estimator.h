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

/// Input mode for the kinship estimation pipeline
enum class InputMode {
    VCF_ONLY,     ///< Mode 1: VCF only (AF and LD from genotype likelihoods)
    VCF_PLINK,    ///< Mode 2: VCF + PLINK (AF and LD from reference panel)
    PLINK_ONLY    ///< Mode 3: PLINK only (hard genotype mode)
};

/// Return a human-readable name for the input mode
inline const char* input_mode_name(InputMode m) {
    switch (m) {
        case InputMode::VCF_ONLY:   return "VCF-only";
        case InputMode::VCF_PLINK:  return "VCF+PLINK";
        case InputMode::PLINK_ONLY: return "PLINK-only";
    }
    return "Unknown";
}

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

    // VCF FORMAT field for Phred-scaled genotype likelihoods (default: PL)
    std::string pl_field = "PL";

    // Input mode (auto-detected from provided files)
    InputMode input_mode = InputMode::VCF_PLINK;

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
    LikelihoodMatrix _gl_matrix;                ///< [sample][snp] genotype likelihoods
    IBS_IBD_Matrix _ibs_ibd;                   ///< [9][snp][3]
    std::vector<std::vector<int8_t>> _bed_genotypes; ///< [sample][snp] for LD (Mode 2/3)
    std::vector<std::vector<double>> _expected_genotypes; ///< [sample][snp] for LD (Mode 1)

    // Mapping: VCF SNP index → .bim SNP index (Mode 2 only)
    std::vector<int> _vcf_to_bim_index;

    // Internal steps — Mode-dependent loading
    void _load_vcf();           ///< Mode 2: VCF with .bim whitelist
    void _load_vcf_only();      ///< Mode 1: VCF without .bim filter
    void _load_plink_as_gl();   ///< Mode 3: .bed → delta-function likelihoods

    // Internal steps — common
    void _load_frequencies();       ///< Mode 2: AF from .bed (or .frq)
    void _compute_af_from_likelihoods(); ///< Mode 1: AF via EM from likelihoods
    void _compute_af_from_bed();    ///< Mode 3: AF from .bed directly
    void _precompute_ibs_ibd();     ///< All modes
    void _load_bed_genotypes();     ///< Mode 2: .bed for LD pruning
    void _compute_expected_g();     ///< Mode 1: posterior E[G] for LD pruning

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
