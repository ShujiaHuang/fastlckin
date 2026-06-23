/**
 * @file main.cpp
 * @brief fastlckin entry point: sub-command dispatcher
 * @author Shujia Huang
 * @date 2025-06-23
 */

#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <stdexcept>
#include <getopt.h>

#include "version.h"
#include "kinship_estimator.h"
#include "frequency.h"
#include <fstream>

static int usage() {
    std::cerr << "\n"
              << "fastlckin v" << FASTLCKIN_VERSION << "\n"
              << "  " << FASTLCKIN_DESCRIPTION << "\n"
              << "Author: " << FASTLCKIN_AUTHOR << "\n\n"
              << "Usage: fastlckin <command> [options]\n\n"
              << "Commands:\n"
              << "  relatedness   Estimate pairwise kinship (IBD coefficients)\n"
              << "  freq          Compute allele frequencies from PLINK .bed/.bim files\n\n"
              << "Use 'fastlckin <command> -h' for command-specific help.\n";
    return 1;
}

// ────────────────────────────────────────────────────────────────────
// Sub-command: relatedness
// ────────────────────────────────────────────────────────────────────
static int relatedness_usage(const fastlckin::KinshipConfig &config) {
    std::cerr << "\n"
              << "Usage: fastlckin relatedness [-v <VCF>] [-p <PLINK_PREFIX>] [options]\n\n"
              << "At least one of -v or -p must be provided:\n"
              << "  -v, --vcf FILE          Input VCF/BCF file (.vcf, .vcf.gz, .bcf)\n"
              << "  -p, --plink PREFIX      PLINK binary file prefix (.bed/.bim/.fam)\n\n"
              << "Input modes:\n"
              << "  -v only      VCF-only (AF and LD from genotype likelihoods)\n"
              << "  -v + -p      VCF + PLINK (AF and LD from reference panel)\n"
              << "  -p only      PLINK-only (hard genotype mode)\n\n"
              << "Optional:\n"
              << "  -F, --freq FILE         Pre-computed .frq file (default: auto-compute)\n"
              << "  -f, --fst FLOAT         Prior FST value (default: " << config.fst << ")\n"
              << "  -t, --threads INT       Number of threads (default: " << config.threads << ")\n"
              << "  -o, --output FILE       Output TSV path (default: auto-generated)\n"
              << "      --maf-min FLOAT     Min allele frequency filter (default: " << config.maf_min << ")\n"
              << "      --maf-max FLOAT     Max allele frequency filter (default: " << config.maf_max << ")\n"
              << "      --ld-window INT     LD pruning window size in SNPs (default: " << config.ld_config.window_size << ")\n"
              << "      --ld-step INT       LD pruning step size (default: " << config.ld_config.step_size << ")\n"
              << "      --ld-r2 FLOAT       LD pruning r2 threshold (default: " << config.ld_config.r2_threshold << ")\n"
              << "      --gq-min INT        Min GQ quality threshold (default: " << config.gq_min << ")\n"
              << "      --pl-field STR      VCF FORMAT field for Phred-scaled GL (default: PL)\n"
              << "      --n-restarts INT    Nelder-Mead restarts (default: " << config.n_restarts << ")\n"
              << "      --xtol FLOAT        Optimizer parameter convergence (default: " << config.nm_config.xtol << ")\n"
              << "      --ftol FLOAT        Optimizer function convergence (default: " << config.nm_config.ftol << ")\n"
              << "      --classify          Enable automatic relationship classification\n"
              << "      --verbose           Verbose logging\n"
              << "  -h, --help              Show this help message\n";
    return 1;
}

static int run_relatedness(int argc, char* argv[]) {
    fastlckin::KinshipConfig config;

    static struct option long_options[] = {
        {"vcf",        required_argument, nullptr, 'v'},
        {"plink",      required_argument, nullptr, 'p'},
        {"freq",       required_argument, nullptr, 'F'},
        {"fst",        required_argument, nullptr, 'f'},
        {"threads",    required_argument, nullptr, 't'},
        {"output",     required_argument, nullptr, 'o'},
        {"maf-min",    required_argument, nullptr, 1001},
        {"maf-max",    required_argument, nullptr, 1002},
        {"ld-window",  required_argument, nullptr, 1003},
        {"ld-step",    required_argument, nullptr, 1004},
        {"ld-r2",      required_argument, nullptr, 1005},
        {"gq-min",     required_argument, nullptr, 1006},
        {"pl-field",   required_argument, nullptr, 1012},
        {"n-restarts", required_argument, nullptr, 1007},
        {"xtol",       required_argument, nullptr, 1008},
        {"ftol",       required_argument, nullptr, 1009},
        {"classify",   no_argument,       nullptr, 1010},
        {"verbose",    no_argument,       nullptr, 1011},
        {"help",       no_argument,       nullptr, 'h'},
        {nullptr,      0,                 nullptr,  0 }
    };

    int c;
    optind = 1;
    while ((c = getopt_long(argc, argv, "v:p:F:f:t:o:h", long_options, nullptr)) != -1) {
        switch (c) {
            case 'v': config.vcf_path = optarg; break;
            case 'p': config.plink_prefix = optarg; break;
            case 'F': config.freq_path = optarg; break;
            case 'f': config.fst = std::atof(optarg); break;
            case 't': config.threads = std::atoi(optarg); break;
            case 'o': config.output_path = optarg; break;
            case 1001: config.maf_min = std::atof(optarg); break;
            case 1002: config.maf_max = std::atof(optarg); break;
            case 1003: config.ld_config.window_size = std::atoi(optarg); break;
            case 1004: config.ld_config.step_size = std::atoi(optarg); break;
            case 1005: config.ld_config.r2_threshold = std::atof(optarg); break;
            case 1006: config.gq_min = std::atoi(optarg); break;
            case 1012: config.pl_field = optarg; break;
            case 1007: config.n_restarts = std::atoi(optarg); break;
            case 1008: config.nm_config.xtol = std::atof(optarg); break;
            case 1009: config.nm_config.ftol = std::atof(optarg); break;
            case 1010: config.classify = true; break;
            case 1011: config.verbose = true; break;
            case 'h': return relatedness_usage(config);
            default:  return relatedness_usage(config);
        }
    }

    if (config.vcf_path.empty() && config.plink_prefix.empty()) {
        std::cerr << "Error: at least one of -v or -p must be provided.\n";
        return relatedness_usage(config);
    }

    // Auto-detect input mode
    if (!config.vcf_path.empty() && config.plink_prefix.empty()) {
        config.input_mode = fastlckin::InputMode::VCF_ONLY;
    } else if (!config.vcf_path.empty() && !config.plink_prefix.empty()) {
        config.input_mode = fastlckin::InputMode::VCF_PLINK;
    } else {
        config.input_mode = fastlckin::InputMode::PLINK_ONLY;
    }

    fastlckin::KinshipEstimator estimator(config);
    estimator.run();
    return 0;
}

// ────────────────────────────────────────────────────────────────────
// Sub-command: freq
// ────────────────────────────────────────────────────────────────────
static int freq_usage() {
    std::cerr << "\n"
              << "Usage: fastlckin freq -p <PLINK_PREFIX> -o <OUTPUT.frq> [options]\n\n"
              << "Required:\n"
              << "  -p, --plink PREFIX   PLINK binary file prefix (.bed/.bim/.fam)\n"
              << "  -o, --output FILE    Output .frq file path\n\n"
              << "Optional:\n"
              << "  -t, --threads INT    Number of threads (default: 1)\n"
              << "  -h, --help           Show this help message\n";
    return 1;
}

static int run_freq(int argc, char* argv[]) {
    std::string plink_prefix;
    std::string output_path;
    int threads = 1;

    static struct option long_options[] = {
        {"plink",   required_argument, nullptr, 'p'},
        {"output",  required_argument, nullptr, 'o'},
        {"threads", required_argument, nullptr, 't'},
        {"help",    no_argument,       nullptr, 'h'},
        {nullptr,   0,                 nullptr,  0 }
    };

    int c;
    optind = 1;
    while ((c = getopt_long(argc, argv, "p:o:t:h", long_options, nullptr)) != -1) {
        switch (c) {
        case 'p': plink_prefix = optarg; break;
        case 'o': output_path = optarg; break;
        case 't': threads = std::atoi(optarg); break;
        case 'h': return freq_usage();
        default:  return freq_usage();
        }
    }

    if (plink_prefix.empty() || output_path.empty()) {
        std::cerr << "Error: --plink and --output are required.\n";
        return freq_usage();
    }

    fastlckin::KinshipConfig config;
    config.plink_prefix = plink_prefix;
    config.threads = threads;

    // Build SNP info from .bim and compute frequencies
    std::string bed_path = plink_prefix + ".bed";
    std::string bim_path = plink_prefix + ".bim";

    auto snp_infos = fastlckin::load_bim_snps(bim_path);
    std::vector<int> all_indices(snp_infos.size());
    for (size_t i = 0; i < snp_infos.size(); ++i) all_indices[i] = static_cast<int>(i);

    auto freqs = fastlckin::compute_allele_frequencies(bed_path, bim_path, all_indices);

    // Write .frq file
    std::ofstream ofs(output_path);
    if (!ofs) throw std::runtime_error("Cannot open output: " + output_path);
    ofs << "CHR\tSNP\tA1\tA2\tMAF\tNCHROBS\n";
    for (size_t i = 0; i < snp_infos.size(); ++i) {
        double af = freqs[i];
        double maf = (af > 0.5) ? (1.0 - af) : af;
        std::string a1 = (af > 0.5) ? snp_infos[i].ref : snp_infos[i].alt;
        std::string a2 = (af > 0.5) ? snp_infos[i].alt : snp_infos[i].ref;
        ofs << snp_infos[i].chrom << "\t" << snp_infos[i].id << "\t"
            << a1 << "\t" << a2 << "\t" << maf << "\t" << 0 << "\n";
    }
    ofs.close();

    std::cerr << "[fastlckin] Wrote allele frequencies for " << snp_infos.size()
              << " SNPs to " << output_path << "\n";
    return 0;
}

// ────────────────────────────────────────────────────────────────────
// main
// ────────────────────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    if (argc < 2) return usage();

    std::string cmd(argv[1]);
    int ret = 0;

    try {
        if (cmd == "relatedness") {
            ret = run_relatedness(argc - 1, argv + 1);
        } else if (cmd == "freq") {
            ret = run_freq(argc - 1, argv + 1);
        } else if (cmd == "-h" || cmd == "--help") {
            ret = usage();
        } else if (cmd == "--version" || cmd == "-V") {
            std::cout << "fastlckin " << FASTLCKIN_VERSION << "\n";
        } else {
            std::cerr << "Error: unknown command '" << cmd << "'\n";
            ret = usage();
        }
    } catch (const std::exception& e) {
        std::cerr << "\n[fastlckin] Error: " << e.what() << "\n";
        return 1;
    }

    return ret;
}
