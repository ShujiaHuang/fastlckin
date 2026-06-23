/**
 * @file kinship_estimator.cpp
 * @brief Main kinship estimation pipeline implementation
 * @author Shujia Huang
 * @date 2025-06-23
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
    const double EPS = 1e-10;
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
    _ibs_ibd = precompute_ibs_ibd_matrix(_snp_infos, _config.fst);
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

    // LD pruning (mode-dependent)
    std::vector<int> retained;
    if (_config.input_mode == InputMode::VCF_ONLY) {
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

    // Compute PIBS for this pair: PIBS[9] averaged across retained SNPs
    // Actually PIBS varies per SNP, so we need per-SNP PIBS
    // The likelihood function uses _ibs_ibd per SNP and PIBS per SNP per combo
    // We pass snp_indices to the likelihood function which computes per-SNP

    // Precompute PIBS[combo][snp_idx_within_retained]
    // For efficiency, we store PIBS as [9][n_retained]
    std::vector<std::array<double, 9>> pibs_per_snp(retained.size());
    for (size_t ri = 0; ri < retained.size(); ++ri) {
        int s = retained[ri];
        const auto& gl1 = _gl_matrix[ind1][s];
        const auto& gl2 = _gl_matrix[ind2][s];

        pibs_per_snp[ri][PPQQ] = gl1.gl[0] * gl2.gl[2];
        pibs_per_snp[ri][QQPP] = gl1.gl[2] * gl2.gl[0];
        pibs_per_snp[ri][PPPQ] = gl1.gl[0] * gl2.gl[1];
        pibs_per_snp[ri][PQPP] = gl1.gl[1] * gl2.gl[0];
        pibs_per_snp[ri][PQQQ] = gl1.gl[1] * gl2.gl[2];
        pibs_per_snp[ri][QQPQ] = gl1.gl[2] * gl2.gl[1];
        pibs_per_snp[ri][PPPP] = gl1.gl[0] * gl2.gl[0];
        pibs_per_snp[ri][PQPQ] = gl1.gl[1] * gl2.gl[1];
        pibs_per_snp[ri][QQQQ] = gl1.gl[2] * gl2.gl[2];
    }

    // Objective function wrapper
    auto objective = [&](const std::vector<double>& k) -> double {
        double k0 = 1.0 - k[0] - k[1];
        double k1 = k[0];
        double k2 = k[1];

        if (k0 < 0 || k0 > 1 || k1 < 0 || k1 > 1 || k2 < 0 || k2 > 1) return 1e10;
        if (4.0 * k2 * k0 > k1 * k1 + 1e-10) return 1e10;

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
        return -log_likelihood;
    };

    // Run Nelder-Mead with multiple restarts
    std::mt19937 rng(ind1 * 1000 + ind2);
    NelderMeadResult best;
    best.fval = 1e10;

    for (int restart = 0; restart < _config.n_restarts; ++restart) {
        std::vector<double> x0;
        if (restart == 0) {
            x0 = {0.0, 0.0};
        } else {
            x0 = random_ibd_start(rng);
        }

        auto nm_result = nelder_mead(objective, x0, _config.nm_config);
        if (nm_result.fval < best.fval) {
            best = nm_result;
        }
    }

    // Extract results
    double k1_opt = best.x[0];
    double k2_opt = best.x[1];
    double k0_opt = 1.0 - k1_opt - k2_opt;

    // Clamp to valid range
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

    result.k0 = k0_opt;
    result.k1 = k1_opt;
    result.k2 = k2_opt;
    result.pi_hat = 0.5 * k1_opt + k2_opt;
    result.log_likelihood = best.fval;
    result.failed = false;

    if (_config.classify) {
        result.relationship = classify_relationship(k0_opt, k1_opt, k2_opt, result.pi_hat);
    }

    return result;
}

std::string KinshipEstimator::classify_relationship(double k0, double k1, double k2, double pi_hat) {
    if (pi_hat > 0.708 && k2 > 0.8) return "Duplicate/MZ";
    if (pi_hat >= 0.354 && pi_hat <= 0.708) {
        return (k0 < 0.05) ? "Parent-Offspring" : "Full-Sibling";
    }
    if (pi_hat >= 0.177 && pi_hat < 0.354) return "Second-degree";
    if (pi_hat >= 0.0884 && pi_hat < 0.177) return "Third-degree";
    return "Unrelated";
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
            _load_bed_genotypes();
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

    // Step 4: Prepare LD pruning data (mode-dependent)
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

    // Step 5: Generate all pairs
    int n_samples = static_cast<int>(_sample_names.size());
    std::vector<std::pair<int, int>> pairs;
    for (int i = 0; i < n_samples; ++i) {
        for (int j = i + 1; j < n_samples; ++j) {
            pairs.push_back({i, j});
        }
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
            << "  MAF_filter: [" << _config.maf_min << ", " << _config.maf_max << "]"
            << "  LD_prune: " << _config.ld_config.window_size << "/"
            << _config.ld_config.step_size << "/" << _config.ld_config.r2_threshold << "\n";

        // Column header
        ofs << "Ind1\tInd2\tk0\tk1\tk2\tPI_HAT\tN_SNPs";
        if (_config.classify) ofs << "\tRelationship";
        ofs << "\n";

        // Data rows
        ofs << std::fixed << std::setprecision(6);
        for (const auto& r : results) {
            ofs << r.ind1 << "\t" << r.ind2 << "\t"
                << r.k0 << "\t" << r.k1 << "\t" << r.k2 << "\t"
                << r.pi_hat << "\t" << r.n_snps;
            if (_config.classify) ofs << "\t" << r.relationship;
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
