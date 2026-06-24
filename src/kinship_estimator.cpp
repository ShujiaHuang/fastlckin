/**
 * @file kinship_estimator.cpp
 * @brief Main kinship estimation pipeline implementation
 * @author Shujia Huang
 * @date 2026-06-23
 */

#include "kinship_estimator.h"
#include "frequency.h"
#include "frequency_from_likelihoods.h"
#include "algorithm.h"
#include "version.h"
#include "external/thread_pool.h"

#include "io/vcf.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <random>
#include <chrono>
#include <map>
#include <ctime>
#include <iomanip>

namespace fastlckin {

KinshipEstimator::KinshipEstimator(const KinshipConfig& config)
    : _config(config)
{
    if (_config.output_path.empty()) {
        auto now = std::chrono::system_clock::now();
        auto t = std::chrono::system_clock::to_time_t(now);
        std::ostringstream oss;
        oss << "fastlckin_output_" << std::put_time(std::localtime(&t), "%Y%m%d_%H%M%S") << ".tsv";
        _config.output_path = oss.str();
    }
}

void KinshipEstimator::_load_vcf() {
    if (_config.verbose) std::cerr << "[fastlckin] Loading VCF: " << _config.vcf_path << "\n";

    ngslib::VCFFile vcf(_config.vcf_path, "r");
    auto& hdr = vcf.header();
    _sample_names = hdr.sample_names();
    int n_samples = static_cast<int>(_sample_names.size());

    if (_config.verbose) std::cerr << "[fastlckin]   " << n_samples << " samples found\n";

    // Build .bim SNP lookup for matching (bim uses 1-based positions)
    std::string bim_path = _config.plink_prefix + ".bim";
    auto bim_snps = load_bim_snps(bim_path);
    std::map<std::string, int> bim_lookup;
    for (int i = 0; i < static_cast<int>(bim_snps.size()); ++i) {
        std::string key = bim_snps[i].chrom + "_" + std::to_string(bim_snps[i].pos);
        bim_lookup[key] = i;
    }

    // Stream-read VCF records
    ngslib::VCFRecord rec;
    int total_records = 0;
    int kept_records = 0;

    std::vector<std::vector<GenotypeLikelihood>> snp_gls;

    while (vcf.read(rec) >= 0) {
        rec.unpack(BCF_UN_FMT);
        ++total_records;

        // Filter: biallelic, single-base
        if (rec.n_alt() != 1) continue;
        std::string ref = rec.ref();
        auto alts = rec.alt();
        if (ref.size() != 1 || alts.empty() || alts[0].size() != 1) continue;

        std::string chrom = rec.chrom(hdr);
        hts_pos_t pos = rec.pos();
        std::string alt = alts[0];

        // Match to .bim: VCF pos is 0-based, .bim pos is 1-based
        std::string key = chrom + "_" + std::to_string(pos + 1);
        auto it = bim_lookup.find(key);
        if (it == bim_lookup.end()) continue;

        auto gls = extract_genotype_likelihoods(rec, hdr, n_samples, _config.gq_min, _config.pl_field);

        SNPInfo si;
        si.chrom = chrom;
        si.pos = pos;
        si.id = key;
        si.ref = ref;
        si.alt = alt;
        si.af = 0.0;
        si.af_masked = false;

        _snp_infos.push_back(std::move(si));
        _vcf_to_bim_index.push_back(it->second);
        snp_gls.push_back(std::move(gls));
        ++kept_records;
    }

    if (_config.verbose) {
        std::cerr << "[fastlckin]   Read " << total_records << " records, kept "
                  << kept_records << " biallelic SNPs matching .bim\n";
    }

    // Transpose to [sample][snp]
    int n_snps = static_cast<int>(_snp_infos.size());
    _gl_matrix.assign(n_samples, std::vector<GenotypeLikelihood>(n_snps));
    for (int s = 0; s < n_snps; ++s) {
        for (int i = 0; i < n_samples; ++i) {
            _gl_matrix[i][s] = snp_gls[s][i];
        }
    }
}

void KinshipEstimator::_load_vcf_only() {
    if (_config.verbose) std::cerr << "[fastlckin] Loading VCF (VCF-only mode): " << _config.vcf_path << "\n";

    ngslib::VCFFile vcf(_config.vcf_path, "r");
    auto& hdr = vcf.header();
    _sample_names = hdr.sample_names();
    int n_samples = static_cast<int>(_sample_names.size());

    if (_config.verbose) std::cerr << "[fastlckin]   " << n_samples << " samples found\n";

    // Stream-read ALL biallelic VCF records (no .bim filter)
    ngslib::VCFRecord rec;
    int total_records = 0;
    int kept_records = 0;

    std::vector<std::vector<GenotypeLikelihood>> snp_gls;

    while (vcf.read(rec) >= 0) {
        rec.unpack(BCF_UN_FMT);
        ++total_records;

        // Filter: biallelic, single-base
        if (rec.n_alt() != 1) continue;
        std::string ref = rec.ref();
        auto alts = rec.alt();
        if (ref.size() != 1 || alts.empty() || alts[0].size() != 1) continue;

        std::string chrom = rec.chrom(hdr);
        hts_pos_t pos = rec.pos();
        std::string alt = alts[0];

        auto gls = extract_genotype_likelihoods(rec, hdr, n_samples, _config.gq_min, _config.pl_field);

        SNPInfo si;
        si.chrom = chrom;
        si.pos = pos;
        si.id = chrom + "_" + std::to_string(pos + 1);
        si.ref = ref;
        si.alt = alt;
        si.af = 0.0;
        si.af_masked = false;

        _snp_infos.push_back(std::move(si));
        snp_gls.push_back(std::move(gls));
        ++kept_records;
    }

    if (_config.verbose) {
        std::cerr << "[fastlckin]   Read " << total_records << " records, kept "
                  << kept_records << " biallelic SNPs\n";
    }

    // Transpose to [sample][snp]
    int n_snps = static_cast<int>(_snp_infos.size());
    _gl_matrix.assign(n_samples, std::vector<GenotypeLikelihood>(n_snps));
    for (int s = 0; s < n_snps; ++s) {
        for (int i = 0; i < n_samples; ++i) {
            _gl_matrix[i][s] = snp_gls[s][i];
        }
    }
}

void KinshipEstimator::_load_plink_as_gl() {
    if (_config.verbose) std::cerr << "[fastlckin] Loading PLINK (PLINK-only mode): " << _config.plink_prefix << "\n";

    std::string bed_path = _config.plink_prefix + ".bed";
    std::string bim_path = _config.plink_prefix + ".bim";
    std::string fam_path = _config.plink_prefix + ".fam";

    // Load sample names from .fam
    {
        std::ifstream ifs(fam_path);
        if (!ifs) throw std::runtime_error("Cannot open .fam file: " + fam_path);
        std::string line;
        while (std::getline(ifs, line)) {
            if (line.empty()) continue;
            std::istringstream iss(line);
            std::string fam_id, sample_id;
            iss >> fam_id >> sample_id;
            _sample_names.push_back(sample_id);
        }
    }
    int n_samples = static_cast<int>(_sample_names.size());
    if (_config.verbose) std::cerr << "[fastlckin]   " << n_samples << " samples found\n";

    // Load SNP metadata from .bim
    _snp_infos = load_bim_snps(bim_path);
    int n_snps = static_cast<int>(_snp_infos.size());
    if (_config.verbose) std::cerr << "[fastlckin]   " << n_snps << " SNPs in .bim\n";

    // Load hard genotypes from .bed
    std::vector<int> all_indices(n_snps);
    for (int i = 0; i < n_snps; ++i) all_indices[i] = i;
    _bed_genotypes = read_bed_genotypes(bed_path, bim_path, all_indices);

    // Convert hard genotypes to delta-function GL
    const double EPS = 1e-20;
    _gl_matrix.assign(n_samples, std::vector<GenotypeLikelihood>(n_snps));

    for (int i = 0; i < n_samples; ++i) {
        for (int s = 0; s < n_snps; ++s) {
            int8_t g = _bed_genotypes[i][s];
            if (g < 0) {
                // Missing genotype
                _gl_matrix[i][s].masked = true;
                _gl_matrix[i][s].gl[0] = EPS;
                _gl_matrix[i][s].gl[1] = EPS;
                _gl_matrix[i][s].gl[2] = EPS;
            } else {
                // Delta-function: observed genotype has likelihood 1
                _gl_matrix[i][s].gl[0] = EPS;
                _gl_matrix[i][s].gl[1] = EPS;
                _gl_matrix[i][s].gl[2] = EPS;
                _gl_matrix[i][s].gl[g] = 1.0;
                _gl_matrix[i][s].masked = false;
            }
        }
    }

    if (_config.verbose) {
        std::cerr << "[fastlckin]   Converted " << n_snps << " SNPs to delta-function GL\n";
    }
}

void KinshipEstimator::_load_frequencies() {
    if (_config.verbose) std::cerr << "[fastlckin] Loading allele frequencies...\n";

    std::string bed_path = _config.plink_prefix + ".bed";
    std::string bim_path = _config.plink_prefix + ".bim";

    std::vector<double> freqs;

    if (!_config.freq_path.empty()) {
        freqs = read_plink_freq(_config.freq_path, _snp_infos);
    } else {
        freqs = compute_allele_frequencies(bed_path, bim_path, _vcf_to_bim_index);
    }

    // Apply MAF filter and assign frequencies
    int n_masked = 0;
    for (size_t i = 0; i < _snp_infos.size(); ++i) {
        double af = (i < freqs.size()) ? freqs[i] : -1.0;
        if (af < 0.0 || af > 1.0) {
            _snp_infos[i].af_masked = true;
            ++n_masked;
            continue;
        }
        _snp_infos[i].af = af;

        double maf = std::min(af, 1.0 - af);
        if (maf < _config.maf_min || maf > _config.maf_max) {
            _snp_infos[i].af_masked = true;
            ++n_masked;
        }
    }

    if (_config.verbose) {
        std::cerr << "[fastlckin]   " << (_snp_infos.size() - n_masked) << " SNPs passed MAF filter, "
                  << n_masked << " masked\n";
    }
}

void KinshipEstimator::_compute_af_from_likelihoods() {
    if (_config.verbose) std::cerr << "[fastlckin] Computing allele frequencies from genotype likelihoods (EM algorithm)...\n";

    auto freqs = compute_af_from_likelihoods(_gl_matrix, 100, 1e-6, _config.verbose);

    // The EM algorithm converges to the frequency of the allele at genotype
    // index 2 (gl[2]).  In VCF convention, gl[0]=REF/REF and gl[2]=ALT/ALT,
    // so EM gives the ALT allele frequency.
    //
    // Our convention: snp_infos[i].af stores the ALT allele frequency,
    // consistent with the IBD model where P = gl[2] = ALT allele.
    // All three AF paths (EM, compute_allele_frequencies, read_plink_freq)
    // return ALT allele frequency for consistency.
    //
    // Apply MAF filter and assign frequencies
    int n_masked = 0;
    for (size_t i = 0; i < _snp_infos.size(); ++i) {
        double af = (i < freqs.size()) ? freqs[i] : -1.0;
        if (af < 0.0 || af > 1.0) {
            _snp_infos[i].af_masked = true;
            ++n_masked;
            continue;
        }
        _snp_infos[i].af = af;

        double maf = std::min(af, 1.0 - af);
        if (maf < _config.maf_min || maf > _config.maf_max) {
            _snp_infos[i].af_masked = true;
            ++n_masked;
        }
    }

    if (_config.verbose) {
        std::cerr << "[fastlckin]   " << (_snp_infos.size() - n_masked) << " SNPs passed MAF filter, "
                  << n_masked << " masked\n";
    }
}

void KinshipEstimator::_compute_af_from_bed() {
    if (_config.verbose) std::cerr << "[fastlckin] Computing allele frequencies from .bed...\n";

    std::string bed_path = _config.plink_prefix + ".bed";
    std::string bim_path = _config.plink_prefix + ".bim";

    int n_snps = static_cast<int>(_snp_infos.size());
    std::vector<int> all_indices(n_snps);
    for (int i = 0; i < n_snps; ++i) all_indices[i] = i;

    auto freqs = compute_allele_frequencies(bed_path, bim_path, all_indices);

    // Apply MAF filter and assign frequencies
    int n_masked = 0;
    for (size_t i = 0; i < _snp_infos.size(); ++i) {
        double af = (i < freqs.size()) ? freqs[i] : -1.0;
        if (af < 0.0 || af > 1.0) {
            _snp_infos[i].af_masked = true;
            ++n_masked;
            continue;
        }
        _snp_infos[i].af = af;

        double maf = std::min(af, 1.0 - af);
        if (maf < _config.maf_min || maf > _config.maf_max) {
            _snp_infos[i].af_masked = true;
            ++n_masked;
        }
    }

    if (_config.verbose) {
        std::cerr << "[fastlckin]   " << (_snp_infos.size() - n_masked) << " SNPs passed MAF filter, "
                  << n_masked << " masked\n";
    }
}

void KinshipEstimator::_precompute_ibs_ibd() {
    if (_config.verbose) std::cerr << "[fastlckin] Precomputing IBS|IBD probability matrix...\n";
    
    if (!_config.fst_path.empty()) {
        // v0.4.0: Load per-SNP FST from file
        _fst_vector = _load_fst_file(_config.fst_path);
        _ibs_ibd = precompute_ibs_ibd_matrix_fst_vector(_snp_infos, _fst_vector);
        if (_config.verbose) {
            std::cerr << "[fastlckin]   Using per-SNP FST values from " 
                      << _config.fst_path << " (" << _fst_vector.size() << " SNPs)\n";
        }
    } else {
        // Use global FST
        _ibs_ibd = precompute_ibs_ibd_matrix(_snp_infos, _config.fst);
        if (_config.verbose) {
            std::cerr << "[fastlckin]   Using global FST = " << _config.fst << "\n";
        }
    }

}

void KinshipEstimator::_precompute_relationship_likelihoods() {
    if (_config.verbose) std::cerr << "[fastlckin] Precomputing exact relationship likelihoods...\n";
    _rel_likelihoods = precompute_relationship_likelihoods(_snp_infos, _config.fst);
    if (_config.verbose) {
        std::cerr << "[fastlckin]   Computed likelihoods for " << NUM_RELATIONSHIP_TYPES
                  << " relationship types x " << _snp_infos.size() << " SNPs\n";
    }
}

int KinshipEstimator::_model_selection_pair(
    const std::vector<std::array<double, 9>>& pibs_per_snp,
    const std::vector<int>& retained,
    KinshipResult& result
) {
    // Step 1: Compute log-likelihood for each of the 9 relationship types
    for (int r = 0; r < NUM_RELATIONSHIP_TYPES; ++r) {
        result.rel_log_likelihoods[r] = compute_relationship_log_likelihood(
            _rel_likelihoods, pibs_per_snp, retained, static_cast<RelationshipType>(r)
        );
    }

    // Step 2: Group by IBD equivalence class and find best class
    // Due to the "random spouse" assumption, relationships with the same IBD
    // coefficients produce identical likelihoods. We select by equivalence class.
    double class_best_ll[NUM_IBD_CLASSES];
    for (int c = 0; c < NUM_IBD_CLASSES; ++c) {
        class_best_ll[c] = -1e18;
    }

    for (int r = 0; r < NUM_RELATIONSHIP_TYPES; ++r) {
        int c = static_cast<int>(relationship_to_ibd_class(static_cast<RelationshipType>(r)));
        if (result.rel_log_likelihoods[r] > class_best_ll[c]) {
            class_best_ll[c] = result.rel_log_likelihoods[r];
        }
    }

    // Step 3: Select the IBD class with highest log-likelihood
    int best_class = 0;
    double best_ll = class_best_ll[0];
    for (int c = 1; c < NUM_IBD_CLASSES; ++c) {
        if (class_best_ll[c] > best_ll) {
            best_ll = class_best_ll[c];
            best_class = c;
        }
    }

    result.ms_relationship = ibd_class_name(static_cast<IBDEquivalenceClass>(best_class));
    return best_class;
}

void KinshipEstimator::_load_bed_genotypes() {
    if (_config.verbose) std::cerr << "[fastlckin] Loading .bed genotypes for LD pruning...\n";

    std::string bed_path = _config.plink_prefix + ".bed";
    std::string bim_path = _config.plink_prefix + ".bim";

    _bed_genotypes = read_bed_genotypes(bed_path, bim_path, _vcf_to_bim_index);

    if (_config.verbose) {
        std::cerr << "[fastlckin]   " << _bed_genotypes.size() << " samples × "
                  << (_bed_genotypes.empty() ? 0 : _bed_genotypes[0].size()) << " SNPs loaded\n";
    }
}

void KinshipEstimator::_compute_expected_g() {
    if (_config.verbose) std::cerr << "[fastlckin] Computing expected genotypes for LD pruning...\n";

    std::vector<double> afs(_snp_infos.size());
    for (size_t i = 0; i < _snp_infos.size(); ++i) {
        afs[i] = _snp_infos[i].af;
    }

    _expected_genotypes = compute_expected_genotypes(_gl_matrix, afs);

    if (_config.verbose) {
        std::cerr << "[fastlckin]   " << _expected_genotypes.size() << " samples × "
                  << (_expected_genotypes.empty() ? 0 : _expected_genotypes[0].size()) << " SNPs computed\n";
    }
}

KinshipResult KinshipEstimator::_estimate_pair(int ind1, int ind2) {
    KinshipResult result;
    result.ind1 = _sample_names[ind1];
    result.ind2 = _sample_names[ind2];

    int n_snps = static_cast<int>(_snp_infos.size());

    // Build combined mask
    std::vector<bool> mask(n_snps, false);
    for (int s = 0; s < n_snps; ++s) {
        mask[s] = _snp_infos[s].af_masked ||
                  _gl_matrix[ind1][s].masked ||
                  _gl_matrix[ind2][s].masked;
    }

    // LD pruning (mode-dependent) or use all unmasked SNPs if skipped
    std::vector<int> retained;
    if (_config.ld_config.skip) {
        // Skip LD pruning: retain all unmasked SNPs
        retained.reserve(n_snps);
        for (int s = 0; s < n_snps; ++s) {
            if (!mask[s]) retained.push_back(s);
        }
    } else if (_config.global_ld_prune) {
        // Use global LD-pruned SNPs (v0.3.0)
        // Filter by pair-specific mask
        for (int s : _global_ld_retained) {
            if (!mask[s]) {
                retained.push_back(s);
            }
        }
    } else if (_config.input_mode == InputMode::VCF_ONLY) {
        retained = ld_prune_from_likelihoods(_expected_genotypes, mask, _config.ld_config);
    } else {
        retained = ld_prune(_bed_genotypes, mask, _config.ld_config);
    }
    result.n_snps = static_cast<int>(retained.size());

    if (retained.empty()) {
        result.failed = true;
        result.k0 = -9; result.k1 = -9; result.k2 = -9;
        result.pi_hat = -9; result.log_likelihood = -9;
        result.relationship = "Failed";
        return result;
    }

    // Compute PIBS for this pair: genotype likelihood products per IBS combo
    // PIBS[c] = GL1[g1] * GL2[g2] for the genotype combo c = (g1, g2)
    std::vector<std::array<double, 9>> pibs_per_snp(retained.size());
    for (size_t ri = 0; ri < retained.size(); ++ri) {
        int s = retained[ri];
        const auto& gl1 = _gl_matrix[ind1][s];
        const auto& gl2 = _gl_matrix[ind2][s];

        // P = gl[2] = ALT/ALT, Q = gl[0] = REF/REF
        // This follows the standard VCF convention where gl[2] corresponds
        // to the ALT allele, so the IBD model's af parameter = ALT frequency.
        pibs_per_snp[ri][PPQQ] = gl1.gl[2] * gl2.gl[0];
        pibs_per_snp[ri][QQPP] = gl1.gl[0] * gl2.gl[2];
        pibs_per_snp[ri][PPPQ] = gl1.gl[2] * gl2.gl[1];
        pibs_per_snp[ri][PQPP] = gl1.gl[1] * gl2.gl[2];
        pibs_per_snp[ri][PQQQ] = gl1.gl[1] * gl2.gl[0];
        pibs_per_snp[ri][QQPQ] = gl1.gl[0] * gl2.gl[1];
        pibs_per_snp[ri][PPPP] = gl1.gl[2] * gl2.gl[2];
        pibs_per_snp[ri][PQPQ] = gl1.gl[1] * gl2.gl[1];
        pibs_per_snp[ri][QQQQ] = gl1.gl[0] * gl2.gl[0];
    }

    // Objective function using the Anderson-Weir Mij model
    // Uses smooth quadratic penalty for the Franks constraint (4*k0*k2 <= k1*k1)
    auto objective = [&](const std::vector<double>& k) -> double {
        double k0 = 1.0 - k[0] - k[1];
        double k1 = k[0];
        double k2 = k[1];

        if (k0 < 0 || k0 > 1 || k1 < 0 || k1 > 1 || k2 < 0 || k2 > 1) return 1e10;

        // Smooth Franks constraint penalty: 4*k0*k2 <= k1*k1
        double franks_violation = 4.0 * k0 * k2 - k1 * k1;
        double penalty = 0.0;
        if (franks_violation > 0.0) {
            penalty = 100.0 * franks_violation * franks_violation;
        }

        double log_likelihood = 0.0;
        for (size_t ri = 0; ri < retained.size(); ++ri) {
            int s = retained[ri];
            double site_prob = 0.0;
            for (int c = 0; c < 9; ++c) {
                double ibd_sum = k0 * _ibs_ibd[c][s][0]
                               + k1 * _ibs_ibd[c][s][1]
                               + k2 * _ibs_ibd[c][s][2];
                site_prob += pibs_per_snp[ri][c] * ibd_sum;
            }
            if (site_prob <= 0.0 || std::isinf(site_prob)) return 1e10;
            log_likelihood += std::log(site_prob);
        }
        if (std::isinf(log_likelihood)) return 1e10;
        return -log_likelihood + penalty;
    };

    // ── v0.5.0: Count IBS=0 sites for PO detection ───────────────
    // Parent-offspring pairs can NEVER have opposite homozygotes (IBS=0).
    int ibs0_count = 0;
    for (size_t ri = 0; ri < retained.size(); ++ri) {
        int s = retained[ri];
        const auto& gl1 = _gl_matrix[ind1][s];
        const auto& gl2 = _gl_matrix[ind2][s];
        bool g1_hom_ref = (gl1.gl[0] > gl1.gl[1] && gl1.gl[0] > gl1.gl[2]);
        bool g1_hom_alt = (gl1.gl[2] > gl1.gl[0] && gl1.gl[2] > gl1.gl[1]);
        bool g2_hom_ref = (gl2.gl[0] > gl2.gl[1] && gl2.gl[0] > gl2.gl[2]);
        bool g2_hom_alt = (gl2.gl[2] > gl2.gl[0] && gl2.gl[2] > gl2.gl[1]);
        if ((g1_hom_ref && g2_hom_alt) || (g1_hom_alt && g2_hom_ref)) {
            ++ibs0_count;
        }
    }


    // ── v0.5.0: Strategic multi-start initialization ──────────────────
    // Deterministic starting points cover all major relationship types,
    // each slightly offset from boundaries so the Nelder-Mead simplex
    // can explore in all directions.
    static const std::vector<std::vector<double>> strategic_starts = {
        {0.001, 0.001},  // Near-unrelated:   (k0=0.998, k1=0.001, k2=0.001)
        {0.01, 0.01},    // Unrelated-like:   (k0=0.98, k1=0.01, k2=0.01)
        {0.95, 0.01},    // PO-like:          (k0=0.04, k1=0.95, k2=0.01)
        {0.50, 0.20},    // FS-like:          (k0=0.30, k1=0.50, k2=0.20)
        {0.45, 0.05},    // HS/Avunc-like:    (k0=0.50, k1=0.45, k2=0.05)
        {0.25, 0.25},    // Duplicate-like:   (k0=0.50, k1=0.25, k2=0.25)
    };

    std::mt19937 rng(ind1 * 1000 + ind2);
    NelderMeadResult best;
    best.fval = 1e10;

    // Phase 1: Strategic deterministic starts
    int n_strategic = static_cast<int>(strategic_starts.size());
    for (int restart = 0; restart < n_strategic; ++restart) {
        auto nm_result = nelder_mead(objective, strategic_starts[restart], _config.nm_config);
        if (nm_result.fval < best.fval) {
            best = nm_result;
        }
    }

    // Phase 2: Random starts for remaining restart budget
    for (int restart = n_strategic; restart < _config.n_restarts; ++restart) {
        auto x0 = random_ibd_start(rng);
        auto nm_result = nelder_mead(objective, x0, _config.nm_config);
        if (nm_result.fval < best.fval) {
            best = nm_result;
        }
    }

    // Phase 3: Boundary-near refinement (v0.5.0)
    // Always search along the PO edge (k2=0) and the Duplicate edge (k0≈0)
    // to guarantee these boundary solutions are explored.
    {
        // PO edge search: k2=0, vary k1 from 0.5 to 0.99
        // Find the best starting point on this edge, then refine with NM.
        // Also evaluate exact boundary points
        {
            std::vector<double> unrelated_pt = {0.0001, 0.0001};
            double f_unrelated = objective(unrelated_pt);
            if (f_unrelated < best.fval) {
                best.fval = f_unrelated;
                best.x = unrelated_pt;
            }
        }

        // PO edge search: scan k2=0 boundary from k1=0 to k1=1
        double best_edge_fval = 1e10;
        std::vector<double> best_edge_x0 = {0.95, 0.0};
        for (double k1 = 0.01; k1 <= 0.995; k1 += 0.05) {
            std::vector<double> pt = {k1, 0.0};
            double f = objective(pt);
            if (f < best_edge_fval) {
                best_edge_fval = f;
                best_edge_x0 = pt;
            }
        }
        // Refine from the best edge point
        auto nm_edge = nelder_mead(objective, best_edge_x0, _config.nm_config);
        if (nm_edge.fval < best.fval) {
            best = nm_edge;
        }

        // Also try near the Duplicate vertex (k1≈0, k2≈1)
        std::vector<double> dup_start = {0.001, 0.998};
        auto nm_dup = nelder_mead(objective, dup_start, _config.nm_config);
        if (nm_dup.fval < best.fval) {
            best = nm_dup;
        }

        // Near-boundary refinement from current best
        double k0_best = 1.0 - best.x[0] - best.x[1];
        double k2_best = best.x[1];
        if (k0_best < 0.15) {
            std::vector<double> refine = {0.998, std::max(0.001, k2_best)};
            auto nm_refine = nelder_mead(objective, refine, _config.nm_config);
            if (nm_refine.fval < best.fval) {
                best = nm_refine;
            }
        }
    }

    // Extract results
    double k1_opt = best.x[0];
    double k2_opt = best.x[1];
    double k0_opt = 1.0 - k1_opt - k2_opt;

    // Clamp to valid range FIRST, so SE/CI is computed at a valid point
    k0_opt = std::max(0.0, std::min(1.0, k0_opt));
    k1_opt = std::max(0.0, std::min(1.0, k1_opt));
    k2_opt = std::max(0.0, std::min(1.0, k2_opt));

    // Normalize
    double sum = k0_opt + k1_opt + k2_opt;
    if (sum > 0) {
        k0_opt /= sum;
        k1_opt /= sum;
        k2_opt /= sum;
    }

    // ── v0.5.0: PO override based on IBS=0 count ──────────────────
    bool po_override_applied = false;
    if (k2_opt < 0.02 && ibs0_count == 0 && k0_opt > 0.01) {
        k0_opt = 0.0;
        k1_opt = 1.0;
        k2_opt = 0.0;
        po_override_applied = true;
    }

    // ── v0.7.0: Standard errors via observed Fisher Information ──────────
    // Compute numerical Hessian at the CLAMPED/NORMALIZED MLE point.
    // Uses adaptive step size to maintain precision across different SNP counts.
    if (!po_override_applied) {
        auto neg_log_lik_pure = [&](double k1v, double k2v) -> double {
            double k0v = 1.0 - k1v - k2v;
            double ll = 0.0;
            for (size_t ri = 0; ri < retained.size(); ++ri) {
                int s = retained[ri];
                double site_prob = 0.0;
                for (int c = 0; c < 9; ++c) {
                    double ibd_sum = k0v * _ibs_ibd[c][s][0]
                                   + k1v * _ibs_ibd[c][s][1]
                                   + k2v * _ibs_ibd[c][s][2];
                    site_prob += pibs_per_snp[ri][c] * ibd_sum;
                }
                if (site_prob <= 0.0 || std::isinf(site_prob)) return 1e10;
                ll += std::log(site_prob);
            }
            if (std::isinf(ll)) return 1e10;
            return -ll;
        };

        // Adaptive step size: h = max(1e-6, 1e-4 * sqrt(|f00|/n_snps))
        // Balances truncation error (h⁴f''''/12) vs rounding error (ε·f/h²)
        // For N SNPs: |f| ~ N·2.3, optimal h ~ N^(1/2) · 1e-5
        double f00 = neg_log_lik_pure(k1_opt, k2_opt);
        double n_snps_d = static_cast<double>(retained.size());
        double h_eps = std::max(1e-6, 1e-4 * std::sqrt(std::abs(f00) / std::max(1.0, n_snps_d)));

        double fp1 = neg_log_lik_pure(k1_opt + h_eps, k2_opt);
        double fm1 = neg_log_lik_pure(k1_opt - h_eps, k2_opt);
        double fp2 = neg_log_lik_pure(k1_opt, k2_opt + h_eps);
        double fm2 = neg_log_lik_pure(k1_opt, k2_opt - h_eps);
        double fpp = neg_log_lik_pure(k1_opt + h_eps, k2_opt + h_eps);
        double fpm = neg_log_lik_pure(k1_opt + h_eps, k2_opt - h_eps);
        double fmp = neg_log_lik_pure(k1_opt - h_eps, k2_opt + h_eps);
        double fmm = neg_log_lik_pure(k1_opt - h_eps, k2_opt - h_eps);

        double H00 = (fp1 - 2.0 * f00 + fm1) / (h_eps * h_eps);  // d²f/dk1²
        double H11 = (fp2 - 2.0 * f00 + fm2) / (h_eps * h_eps);  // d²f/dk2²
        double H01 = (fpp - fpm - fmp + fmm) / (4.0 * h_eps * h_eps);  // d²f/dk1dk2

        double det = H00 * H11 - H01 * H01;
        // Hessian must be positive definite: H00 > 0 and det > 0
        // Also require det to be sufficiently large relative to Hessian elements
        // to avoid ill-conditioned inversion (condition number check)
        if (H00 > 0.0 && det > 0.0) {
            double cond = (H00 + H11 + std::sqrt(H01 * H01)) / std::sqrt(det);
            if (cond < 1e8) {  // Condition number threshold
                double var_k1 = H11 / det;
                double var_k2 = H00 / det;
                double cov_k12 = -H01 / det;

                // SE(k0): k0 = 1 - k1 - k2 => Var(k0) = Var(k1) + Var(k2) + 2*Cov(k1,k2)
                double var_k0 = var_k1 + var_k2 + 2.0 * cov_k12;
                double se_k0 = (var_k0 > 0.0) ? std::sqrt(var_k0) : -1.0;
                double se_k1 = (var_k1 > 0.0) ? std::sqrt(var_k1) : -1.0;
                double se_k2 = (var_k2 > 0.0) ? std::sqrt(var_k2) : -1.0;

                // SE(PI_HAT): PI_HAT = 0.5*k1 + k2 (Delta method)
                // Var(PI) = 0.25*Var(k1) + Var(k2) + 2*0.5*Cov(k1,k2)
                double var_pi = 0.25 * var_k1 + var_k2 + cov_k12;
                double se_pi = (var_pi > 0.0) ? std::sqrt(var_pi) : -1.0;

                // Boundary check: SE unreliable when MLE is on parameter boundary
                // Only check parameter bounds (k near 0), not the Cotterman constraint,
                // because the penalty function may push the MLE slightly outside the
                // feasible region.
                bool on_boundary = (k0_opt < 0.01) || (k1_opt < 0.01) || (k2_opt < 0.01);

                if (!on_boundary) {
                    result.se_k0 = se_k0;
                    result.se_k1 = se_k1;
                    result.se_k2 = se_k2;
                    result.se_pi_hat = se_pi;

                    // 95% Wald confidence intervals: estimate +/- 1.96 * SE
                    // Using the CLAMPED/NORMALIZED estimates as CI centers
                    const double z = 1.96;
                    result.ci_k0_lo = std::max(0.0, k0_opt - z * se_k0);
                    result.ci_k0_hi = std::min(1.0, k0_opt + z * se_k0);
                    result.ci_k1_lo = std::max(0.0, k1_opt - z * se_k1);
                    result.ci_k1_hi = std::min(1.0, k1_opt + z * se_k1);
                    result.ci_k2_lo = std::max(0.0, k2_opt - z * se_k2);
                    result.ci_k2_hi = std::min(1.0, k2_opt + z * se_k2);
                    result.ci_pi_hat_lo = std::max(0.0, 0.5 * k1_opt + k2_opt - z * se_pi);
                    result.ci_pi_hat_hi = std::min(1.0, 0.5 * k1_opt + k2_opt + z * se_pi);
                }
            }
        }
    }

    result.k0 = k0_opt;
    result.k1 = k1_opt;
    result.k2 = k2_opt;
    result.pi_hat = 0.5 * k1_opt + k2_opt;
    result.log_likelihood = best.fval;
    result.failed = false;

    // ── v0.7.0: Model selection ──────────────────────────────────
    if (_config.model_selection || _config.output_likelihoods) {
        int best_r = _model_selection_pair(pibs_per_snp, retained, result);
        if (_config.model_selection) {
            // Use model selection result as the primary classification
            result.relationship = result.ms_relationship;
        }
    }

    if (_config.classify && !_config.model_selection) {
        if (_config.classify_config.use_custom) {
            result.relationship = classify_relationship(
                k0_opt, k1_opt, k2_opt, result.pi_hat, _config.classify_config
            );
        } else {
            result.relationship = classify_relationship(k0_opt, k1_opt, k2_opt, result.pi_hat);
        }
    }

    return result;
}

std::string KinshipEstimator::classify_relationship(double k0, double k1, double k2, double pi_hat) {
    ClassificationConfig defaults;
    return classify_relationship(k0, k1, k2, pi_hat, defaults);
}

std::string KinshipEstimator::classify_relationship(
    double k0, double k1, double k2, double pi_hat,
    const ClassificationConfig& config
) {
    if (pi_hat > config.duplicate_threshold && k2 > config.duplicate_k2_threshold) 
        return "Duplicate/MZ";
    if (pi_hat >= config.first_degree_threshold && pi_hat <= config.duplicate_threshold) {
        return (k0 < 0.05) ? "Parent-Offspring" : "Full-Sibling";
    }
    if (pi_hat >= config.second_degree_threshold && pi_hat < config.first_degree_threshold) 
        return "Second-degree";
    if (pi_hat >= config.third_degree_threshold && pi_hat < config.second_degree_threshold) 
        return "Third-degree";
    return "Unrelated";
}

std::vector<std::pair<int, int>> KinshipEstimator::_load_pairs_file(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open pairs file: " + path);

    // Build sample name → index lookup
    std::map<std::string, int> name_to_idx;
    for (int i = 0; i < static_cast<int>(_sample_names.size()); ++i) {
        name_to_idx[_sample_names[i]] = i;
    }

    std::vector<std::pair<int, int>> pairs;
    std::string line;
    int line_num = 0;
    int n_skipped = 0;

    while (std::getline(ifs, line)) {
        ++line_num;
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string name1, name2;
        if (!(iss >> name1 >> name2)) {
            if (_config.verbose) {
                std::cerr << "[fastlckin] Warning: skipping malformed line " << line_num
                          << " in pairs file\n";
            }
            ++n_skipped;
            continue;
        }

        auto it1 = name_to_idx.find(name1);
        auto it2 = name_to_idx.find(name2);

        if (it1 == name_to_idx.end() || it2 == name_to_idx.end()) {
            if (_config.verbose) {
                std::cerr << "[fastlckin] Warning: unknown sample name on line " << line_num
                          << " (" << name1 << ", " << name2 << "), skipping\n";
            }
            ++n_skipped;
            continue;
        }

        int idx1 = it1->second;
        int idx2 = it2->second;
        // Ensure consistent ordering (smaller index first)
        if (idx1 > idx2) std::swap(idx1, idx2);
        if (idx1 != idx2) {
            pairs.push_back({idx1, idx2});
        }
    }

    if (_config.verbose && n_skipped > 0) {
        std::cerr << "[fastlckin]   Skipped " << n_skipped << " invalid/unknown lines in pairs file\n";
    }

    return pairs;
}

std::vector<double> KinshipEstimator::_load_fst_file(const std::string& path) {
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open FST file: " + path);

    // Read FST file: CHR\tPOS\tFST (with optional header)
    // Build lookup: "CHR_POS" → FST value
    std::map<std::string, double> fst_lookup;
    std::string line;
    int n_lines = 0;

    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string chr, fst_str;
        std::string pos_str;

        if (!(iss >> chr >> pos_str >> fst_str)) continue;

        // Skip header line
        if (chr == "CHR" || chr == "chr") continue;

        try {
            double fst_val = std::stod(fst_str);
            fst_lookup[chr + "_" + pos_str] = fst_val;
            ++n_lines;
        } catch (...) {
            continue;
        }
    }

    if (_config.verbose) {
        std::cerr << "[fastlckin]   Loaded " << n_lines << " FST values from file\n";
    }

    // Map FST values to _snp_infos order
    // SNPInfo stores pos as 0-based, FST file may use 1-based (PLINK .bim convention)
    std::vector<double> fst_vector(_snp_infos.size(), _config.fst);  // default: global FST
    int n_matched = 0;
    int n_unmatched = 0;

    for (size_t i = 0; i < _snp_infos.size(); ++i) {
        const auto& snp = _snp_infos[i];
        // Try both 0-based and 1-based positions
        std::string key_0based = snp.chrom + "_" + std::to_string(snp.pos);
        std::string key_1based = snp.chrom + "_" + std::to_string(snp.pos + 1);

        auto it = fst_lookup.find(key_1based);
        if (it == fst_lookup.end()) {
            it = fst_lookup.find(key_0based);
        }
        if (it == fst_lookup.end()) {
            // Try snp.id (which is CHR_POS format)
            it = fst_lookup.find(snp.id);
        }

        if (it != fst_lookup.end()) {
            fst_vector[i] = it->second;
            ++n_matched;
        } else {
            ++n_unmatched;
        }
    }

    if (_config.verbose) {
        std::cerr << "[fastlckin]   FST matching: " << n_matched << " matched, "
                  << n_unmatched << " using global FST default\n";
    }

    return fst_vector;
}

double KinshipEstimator::_king_robust_pihat(int ind1, int ind2) {
    // KING-robust formula (Manichaikul et al. 2010):
    // PI_HAT = sum_s [ (g1_s - 2*p_s) * (g2_s - 2*p_s) ] / sum_s [ 2 * p_s * (1 - p_s) ]
    // where g1_s, g2_s are expected genotypes E[G] = P(G=1) + 2*P(G=2)
    // and p_s is the allele frequency

    int n_snps = static_cast<int>(_snp_infos.size());
    double sum_num = 0.0;
    double sum_den = 0.0;

    for (int s = 0; s < n_snps; ++s) {
        if (_snp_infos[s].af_masked) continue;
        if (_gl_matrix[ind1][s].masked || _gl_matrix[ind2][s].masked) continue;

        double p = _snp_infos[s].af;
        if (p <= 0.0 || p >= 1.0) continue;

        // Expected genotype E[G] = P(G=1) + 2*P(G=2)
        double g1 = _gl_matrix[ind1][s].gl[1] + 2.0 * _gl_matrix[ind1][s].gl[2];
        double g2 = _gl_matrix[ind2][s].gl[1] + 2.0 * _gl_matrix[ind2][s].gl[2];

        // Normalize GL
        double norm1 = _gl_matrix[ind1][s].gl[0] + _gl_matrix[ind1][s].gl[1] + _gl_matrix[ind1][s].gl[2];
        double norm2 = _gl_matrix[ind2][s].gl[0] + _gl_matrix[ind2][s].gl[1] + _gl_matrix[ind2][s].gl[2];

        if (norm1 > 0 && norm2 > 0) {
            g1 /= norm1;
            g2 /= norm2;
        }

        // KING-robust formula components
        double numerator = (g1 - 2.0 * p) * (g2 - 2.0 * p);
        double denominator = 2.0 * p * (1.0 - p);

        sum_num += numerator;
        sum_den += denominator;
    }

    return (sum_den > 1e-10) ? sum_num / sum_den : 0.0;
}

std::vector<std::pair<int, int>> KinshipEstimator::_screen_pairs(
    const std::vector<std::pair<int, int>>& all_pairs
) {
    if (!_config.screening_config.enable_screening) {
        return all_pairs;  // Screening not enabled
    }

    auto t0 = std::chrono::steady_clock::now();

    std::vector<std::pair<int, int>> candidate_pairs;
    int n_screened = 0;

    for (const auto& pair : all_pairs) {
        double pihat = _king_robust_pihat(pair.first, pair.second);

        if (pihat >= _config.screening_config.pi_hat_threshold) {
            candidate_pairs.push_back(pair);
        }
        n_screened++;

        if (_config.screening_config.verbose && (n_screened % 10000 == 0)) {
            std::cerr << "[fastlckin]   Screening: " << n_screened << "/" 
                      << all_pairs.size() << " pairs, " 
                      << candidate_pairs.size() << " candidates\n";
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();

    if (_config.screening_config.verbose) {
        std::cerr << "[fastlckin]   Screening completed: " << n_screened << " pairs in "
                  << std::fixed << std::setprecision(1) << elapsed << "s\n";
        std::cerr << "[fastlckin]   " << candidate_pairs.size() << "/" << n_screened 
                  << " pairs passed threshold (" 
                  << (100.0 * candidate_pairs.size() / std::max(1, n_screened)) << "%)\n";
    }

    return candidate_pairs;
}

std::vector<KinshipResult> KinshipEstimator::run() {
    auto t0 = std::chrono::steady_clock::now();

    if (_config.verbose) {
        std::cerr << "[fastlckin] Input mode: " << input_mode_name(_config.input_mode) << "\n";
    }

    // Step 1: Load input (mode-dependent)
    switch (_config.input_mode) {
        case InputMode::VCF_ONLY:
            _load_vcf_only();
            break;
        case InputMode::VCF_PLINK:
            _load_vcf();
            if (!_config.ld_config.skip) {
                _load_bed_genotypes();  // Only needed for LD pruning
            }
            break;
        case InputMode::PLINK_ONLY:
            _load_plink_as_gl();
            break;
    }

    // Step 2: Load/compute frequencies (mode-dependent)
    switch (_config.input_mode) {
        case InputMode::VCF_ONLY:
            _compute_af_from_likelihoods();
            break;
        case InputMode::VCF_PLINK:
            _load_frequencies();
            break;
        case InputMode::PLINK_ONLY:
            _compute_af_from_bed();
            break;
    }

    // Step 3: Precompute IBS|IBD matrix (all modes)
    _precompute_ibs_ibd();

    // Step 3.5: Precompute exact relationship likelihoods (v0.7.0)
    if (_config.model_selection || _config.output_likelihoods) {
        _precompute_relationship_likelihoods();
    }

    // Step 4: Prepare LD pruning data (mode-dependent; skipped if --no-ld-prune)
    if (!_config.ld_config.skip) {
        switch (_config.input_mode) {
            case InputMode::VCF_ONLY:
                _compute_expected_g();
                break;
            case InputMode::VCF_PLINK:
                // _bed_genotypes already loaded in Step 1
                break;
            case InputMode::PLINK_ONLY:
                // _bed_genotypes already loaded in _load_plink_as_gl()
                break;
        }
        
        // Step 4.5: Global LD pruning (v0.3.0 new feature)
        if (_config.global_ld_prune) {
            if (_config.verbose) {
                std::cerr << "[fastlckin] Performing global LD pruning (once for all samples)...\n";
            }
            
            // Build global mask: SNP is masked if AF is masked
            int n_snps = static_cast<int>(_snp_infos.size());
            std::vector<bool> global_mask(n_snps, false);
            for (int s = 0; s < n_snps; ++s) {
                global_mask[s] = _snp_infos[s].af_masked;
            }
            
            // Perform global LD pruning
            if (_config.input_mode == InputMode::VCF_ONLY) {
                _global_ld_retained = ld_prune_global_expected(_expected_genotypes, global_mask, _config.ld_config);
            } else {
                _global_ld_retained = ld_prune_global(_bed_genotypes, global_mask, _config.ld_config);
            }
            
            if (_config.verbose) {
                std::cerr << "[fastlckin]   Global LD pruning: " << _global_ld_retained.size()
                          << " SNPs retained from " << n_snps << " total\n";
            }
        }
    } else {
        if (_config.verbose) {
            std::cerr << "[fastlckin] LD pruning skipped (--no-ld-prune): using all unmasked SNPs\n";
        }
    }

    // Step 5: Generate all pairs (or load from pairs file)
    int n_samples = static_cast<int>(_sample_names.size());
    std::vector<std::pair<int, int>> pairs;
    
    if (!_config.pairs_path.empty()) {
        // v0.4.0: Load specific pairs from file
        pairs = _load_pairs_file(_config.pairs_path);
        if (_config.verbose) {
            std::cerr << "[fastlckin] Loaded " << pairs.size() 
                      << " specific pairs from " << _config.pairs_path << "\n";
        }
    } else {
        // Generate all pairs
        for (int i = 0; i < n_samples; ++i) {
            for (int j = i + 1; j < n_samples; ++j) {
                pairs.push_back({i, j});
            }
        }
    }

    // Step 5.5: Apply KING-robust screening (v0.4.0)
    if (_config.screening_config.enable_screening) {
        pairs = _screen_pairs(pairs);
    }

    int n_pairs = static_cast<int>(pairs.size());
    if (_config.verbose) {
        std::cerr << "[fastlckin] Estimating " << n_pairs << " pairs with "
                  << _config.threads << " thread(s)...\n";
    }

    // Step 6: Parallel estimation
    std::vector<KinshipResult> results(n_pairs);

    if (_config.threads <= 1) {
        for (int p = 0; p < n_pairs; ++p) {
            results[p] = _estimate_pair(pairs[p].first, pairs[p].second);
            if (_config.verbose && (p + 1) % 100 == 0) {
                std::cerr << "[fastlckin]   " << (p + 1) << "/" << n_pairs << " pairs done\n";
            }
        }
    } else {
        ThreadPool pool(_config.threads);
        std::vector<std::future<KinshipResult>> futures(n_pairs);

        for (int p = 0; p < n_pairs; ++p) {
            futures[p] = pool.submit([this, &pairs, p]() {
                return this->_estimate_pair(pairs[p].first, pairs[p].second);
            });
        }

        for (int p = 0; p < n_pairs; ++p) {
            results[p] = futures[p].get();
            if (_config.verbose && (p + 1) % 100 == 0) {
                std::cerr << "[fastlckin]   " << (p + 1) << "/" << n_pairs << " pairs done\n";
            }
        }
    }

    // Step 7: Sort by (ind1, ind2) and write output
    std::sort(results.begin(), results.end(),
              [](const KinshipResult& a, const KinshipResult& b) {
                  if (a.ind1 != b.ind1) return a.ind1 < b.ind1;
                  return a.ind2 < b.ind2;
              });

    // Write output
    {
        std::ofstream ofs(_config.output_path);
        if (!ofs) throw std::runtime_error("Cannot open output: " + _config.output_path);

        // Header comments
        ofs << "# fastlckin relatedness v" << FASTLCKIN_VERSION << "\n";
        ofs << "# Mode: " << input_mode_name(_config.input_mode) << "\n";

        // Get current time
        auto now = std::chrono::system_clock::now();
        auto t = std::chrono::system_clock::to_time_t(now);
        char time_buf[64];
        std::strftime(time_buf, sizeof(time_buf), "%Y-%m-%d %H:%M:%S", std::localtime(&t));

        ofs << "# Command: fastlckin relatedness";
        if (!_config.vcf_path.empty()) ofs << " -v " << _config.vcf_path;
        if (!_config.plink_prefix.empty()) ofs << " -p " << _config.plink_prefix;
        ofs << " -t " << _config.threads << "\n";
        ofs << "# Date: " << time_buf << "\n";
        ofs << "# FST: " << _config.fst
            << "  MAF_filter: [" << _config.maf_min << ", " << _config.maf_max << "]";
        if (_config.ld_config.skip) {
            ofs << "  LD_prune: skipped (all unmasked SNPs used)\n";
        } else {
            ofs << "  LD_prune: " << _config.ld_config.window_size << "/"
                << _config.ld_config.step_size << "/" << _config.ld_config.r2_threshold << "\n";
        }

        // Column header
        ofs << "Ind1\tInd2\tk0\tk1\tk2\tPI_HAT\tN_SNPs";
        if (_config.classify || _config.model_selection) ofs << "\tRelationship";
        if (_config.model_selection) ofs << "\tMS_Relationship";
        if (_config.output_likelihoods) {
            ofs << "\tLL_UN\tLL_DUP\tLL_PO\tLL_FS\tLL_HS\tLL_AV\tLL_GP\tLL_FC\tLL_GGP";
        }
        ofs << "\tSE_k0\tSE_k1\tSE_k2\tSE_PI_HAT"
            << "\tCI_k0_lo\tCI_k0_hi\tCI_k1_lo\tCI_k1_hi"
            << "\tCI_k2_lo\tCI_k2_hi\tCI_PI_lo\tCI_PI_hi\n";

        // Helper: format double or "NA" for -1 sentinel
        auto fmt_val = [](std::ostream& os, double v) {
            if (v < -0.5) os << "NA";
            else os << std::fixed << std::setprecision(6) << v;
        };

        // Data rows
        ofs << std::fixed << std::setprecision(6);
        for (const auto& r : results) {
            ofs << r.ind1 << "\t" << r.ind2 << "\t"
                << r.k0 << "\t" << r.k1 << "\t" << r.k2 << "\t"
                << r.pi_hat << "\t" << r.n_snps;
            if (_config.classify || _config.model_selection) ofs << "\t" << r.relationship;
            if (_config.model_selection) ofs << "\t" << r.ms_relationship;
            if (_config.output_likelihoods) {
                for (int ri = 0; ri < NUM_RELATIONSHIP_TYPES; ++ri) {
                    ofs << "\t" << std::fixed << std::setprecision(4) << r.rel_log_likelihoods[ri];
                }
            }
            ofs << "\t"; fmt_val(ofs, r.se_k0); ofs << "\t";
            fmt_val(ofs, r.se_k1); ofs << "\t";
            fmt_val(ofs, r.se_k2); ofs << "\t";
            fmt_val(ofs, r.se_pi_hat);
            ofs << "\t"; fmt_val(ofs, r.ci_k0_lo); ofs << "\t";
            fmt_val(ofs, r.ci_k0_hi); ofs << "\t";
            fmt_val(ofs, r.ci_k1_lo); ofs << "\t";
            fmt_val(ofs, r.ci_k1_hi); ofs << "\t";
            fmt_val(ofs, r.ci_k2_lo); ofs << "\t";
            fmt_val(ofs, r.ci_k2_hi); ofs << "\t";
            fmt_val(ofs, r.ci_pi_hat_lo); ofs << "\t";
            fmt_val(ofs, r.ci_pi_hat_hi);
            ofs << "\n";
        }
    }

    auto t1 = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(t1 - t0).count();

    std::cerr << "[fastlckin] Done. " << n_pairs << " pairs estimated in "
              << std::fixed << std::setprecision(1) << elapsed << "s.\n";
    std::cerr << "[fastlckin] Output written to: " << _config.output_path << "\n";

    return results;
}

}  // namespace fastlckin
