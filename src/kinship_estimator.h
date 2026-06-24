#ifndef _FASTLCKIN_KINSHIP_ESTIMATOR_H_
#define _FASTLCKIN_KINSHIP_ESTIMATOR_H_

/**
 * @file kinship_estimator.h
 * @brief Main kinship estimation pipeline coordinator
 * @author Shujia Huang
 * @date 2026-06-23
 */

#include <vector>
#include <string>
#include <array>
#include <utility>

#include "ibd_model.h"
#include "genotype_likelihood.h"
#include "ld_prune.h"
#include "optimizer.h"
#include "relationship_likelihood.h"

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
    double k0 = 0.0;             ///< IBD=0 probability
    double k1 = 0.0;             ///< IBD=1 probability
    double k2 = 0.0;             ///< IBD=2 probability
    double pi_hat = 0.0;         ///< Kinship coefficient = 0.5*k1 + k2
    int n_snps = 0;              ///< Number of SNPs used
    double log_likelihood = 0.0; ///< Best negative log-likelihood
    std::string relationship;    ///< Classification label (threshold-based or model selection)
    bool failed = false;         ///< Estimation failed flag

    // Model selection fields (v0.7.0)
    std::array<double, NUM_RELATIONSHIP_TYPES> rel_log_likelihoods{}; ///< Per-relationship log-likelihoods
    std::string ms_relationship;  ///< Model selection classification result

    // Standard errors and 95% confidence intervals (v0.7.0)
    // -1.0 indicates NA (boundary or non-positive-definite Hessian)
    double se_k0 = -1.0;          ///< SE of k0
    double se_k1 = -1.0;          ///< SE of k1
    double se_k2 = -1.0;          ///< SE of k2
    double se_pi_hat = -1.0;      ///< SE of PI_HAT (Delta method)
    double ci_k0_lo = -1.0;       ///< 95% CI lower bound for k0
    double ci_k0_hi = -1.0;       ///< 95% CI upper bound for k0
    double ci_k1_lo = -1.0;       ///< 95% CI lower bound for k1
    double ci_k1_hi = -1.0;       ///< 95% CI upper bound for k1
    double ci_k2_lo = -1.0;       ///< 95% CI lower bound for k2
    double ci_k2_hi = -1.0;       ///< 95% CI upper bound for k2
    double ci_pi_hat_lo = -1.0;   ///< 95% CI lower bound for PI_HAT
    double ci_pi_hat_hi = -1.0;   ///< 95% CI upper bound for PI_HAT
};

/// Configuration for relationship classification thresholds (v0.4.0)
struct ClassificationConfig {
    double duplicate_threshold     = 0.708;  ///< PI_HAT > this && k2 > dup_k2 → Duplicate/MZ
    double first_degree_threshold  = 0.354;  ///< PI_HAT in [this, dup] → 1st degree
    double second_degree_threshold = 0.177;  ///< PI_HAT in [this, 1st] → 2nd degree
    double third_degree_threshold  = 0.0884; ///< PI_HAT in [this, 2nd] → 3rd degree
    double duplicate_k2_threshold  = 0.8;    ///< k2 > this required for Duplicate/MZ
    bool use_custom = false;                 ///< True if user set custom thresholds
};

/// Configuration for KING-robust pre-screening (v0.4.0)
struct ScreeningConfig {
    bool enable_screening = false;           ///< Enable fast KING-robust pre-filter
    double pi_hat_threshold = 0.0442;        ///< PI_HAT threshold (~3rd degree)
    bool verbose = false;                    ///< Verbose screening output
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
    double fst     = 0.0;
    double maf_min = 0.05;
    double maf_max = 0.95;
    int gq_min     = 1;
    int n_restarts = 5;           ///< Increased from 3 in v0.3.0 for better global optimum search
    bool classify  = false;
    bool model_selection = false;  ///< Enable model selection classification (v0.7.0)
    bool output_likelihoods = false; ///< Output per-relationship log-likelihoods (v0.7.0)
    bool verbose   = false;

    // LD pruning parameters
    LDPruneConfig ld_config;
    bool global_ld_prune = false;   ///< If true, use global LD pruning (v0.3.0)

    // Optimizer parameters
    NelderMeadConfig nm_config;

    // Parallelism
    int threads = 1;

    // v0.4.0 new features
    std::string pairs_path;       ///< Optional file listing specific pairs to estimate (TSV: ind1\tind2)
    ClassificationConfig classify_config;  ///< Custom classification thresholds
    std::string fst_path;         ///< Optional per-SNP FST file (TSV: CHR\tPOS\tFST)
    ScreeningConfig screening_config;  ///< KING-robust pre-screening config
};

/// Main kinship estimator class
class KinshipEstimator {
public:
    explicit KinshipEstimator(const KinshipConfig& config);

    /// Run the full kinship estimation pipeline
    /// @return Results for all sample pairs
    std::vector<KinshipResult> run();

    /// Classify relationship based on IBD coefficients (default thresholds)
    static std::string classify_relationship(double k0, double k1, double k2, double pi_hat);

    /// Classify with custom thresholds (v0.4.0)
    static std::string classify_relationship(
        double k0, double k1, double k2, double pi_hat,
        const ClassificationConfig& config
    );

private:
    KinshipConfig _config;

    // Internal data
    std::vector<std::string> _sample_names;
    std::vector<SNPInfo> _snp_infos;          ///< VCF SNPs (filtered biallelic)
    LikelihoodMatrix _gl_matrix;                ///< [sample][snp] genotype likelihoods
    IBS_IBD_Matrix _ibs_ibd;                   ///< [9][snp][3] (Mij model)
    RelationshipLikelihoodMatrix _rel_likelihoods; ///< Exact per-relationship probs (v0.7.0)
    std::vector<std::vector<int8_t>> _bed_genotypes; ///< [sample][snp] for LD (Mode 2/3)
    std::vector<std::vector<double>> _expected_genotypes; ///< [sample][snp] for LD (Mode 1)
    std::vector<int> _global_ld_retained;      ///< Global LD-pruned SNP indices (v0.3.0)
    std::vector<double> _fst_vector;            ///< Per-SNP FST values (v0.4.0)

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

    // v0.4.0 helper methods
    /// Load specific sample pairs from a TSV file (ind1\tind2 per line)
    std::vector<std::pair<int, int>> _load_pairs_file(const std::string& path);

    /// Load per-SNP FST values from file (v0.4.0)
    std::vector<double> _load_fst_file(const std::string& path);

    /// KING-robust PI_HAT estimate for a pair (v0.4.0)
    double _king_robust_pihat(int ind1, int ind2);

    /// Screen pairs using KING-robust, return only those above threshold (v0.4.0)
    std::vector<std::pair<int, int>> _screen_pairs(
        const std::vector<std::pair<int, int>>& all_pairs
    );

    /// Pre-compute exact relationship likelihoods (v0.7.0)
    void _precompute_relationship_likelihoods();

    /// Perform model selection for a pair (v0.7.0)
    /// @return Index of best relationship type (0-7)
    int _model_selection_pair(
        const std::vector<std::array<double, 9>>& pibs_per_snp,
        const std::vector<int>& retained,
        KinshipResult& result
    );
};

}  // namespace fastlckin

#endif
