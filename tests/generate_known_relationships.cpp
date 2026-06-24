/**
 * @file generate_known_relationships.cpp
 * @brief Configurable synthetic PLINK data generator with KNOWN relationships
 *
 * Generates datasets where the ground-truth IBD coefficients (k0, k1, k2)
 * are known by construction, enabling end-to-end accuracy validation of
 * fastlckin (or any kinship estimator that reads PLINK binary files).
 *
 * Relationship types supported:
 *   Duplicate/MZ, Parent-Offspring, Full-Sibling,
 *   Half-Sibling, Avuncular (uncle/aunt–nephew/niece), Unrelated
 *
 * All counts are configurable via CLI.  The generator automatically
 * calculates the required number of independent founders.
 *
 * Build : compiled alongside other tests via CMake
 * Usage : cd <project_root>
 *         build/tests/generate_known_relationships [options]
 *
 * Example:
 *   build/tests/generate_known_relationships --founders 100 --snps 2000
 *   build/tests/generate_known_relationships --po 20 --fs 10 --hs 5
 *   build/tests/generate_known_relationships --seed 999 --prefix big_test
 */

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <sys/stat.h>

// ── CLI defaults ─────────────────────────────────────────────────────

static constexpr int    DEF_FOUNDERS   = 20;
static constexpr int    DEF_SNPS       = 500;
static constexpr int    DEF_DUP        = 1;
static constexpr int    DEF_PO         = 1;    // extra PO pairs (beyond FS/HS families)
static constexpr int    DEF_FS         = 1;    // FS families (each → 2 children)
static constexpr int    DEF_HS         = 1;    // HS pairs
static constexpr double DEF_AF_MIN     = 0.1;
static constexpr double DEF_AF_MAX     = 0.5;
static constexpr unsigned DEF_SEED     = 42;

// ── Helpers ──────────────────────────────────────────────────────────

static int transmit(int parent_gt, std::mt19937& rng) {
    std::uniform_int_distribution<int> coin(0, 1);
    if (parent_gt == 0) return 0;
    if (parent_gt == 2) return 1;
    return coin(rng);
}

static int make_child(int p1, int p2, std::mt19937& rng) {
    return transmit(p1, rng) + transmit(p2, rng);
}

static void gen_hwe(std::vector<int>& gt, int n, double p, std::mt19937& rng) {
    std::uniform_real_distribution<double> U(0.0, 1.0);
    double pAA = (1 - p) * (1 - p);
    double pAB = 2 * p * (1 - p);
    gt.resize(n);
    for (int i = 0; i < n; ++i) {
        double u = U(rng);
        gt[i] = (u < pAA) ? 0 : (u < pAA + pAB) ? 1 : 2;
    }
}

static void usage(const char* prog) {
    std::cerr <<
        "Usage: " << prog << " [options]\n\n"
        "Sample composition:\n"
        "  --founders INT    Independent (unrelated) founders     [" << DEF_FOUNDERS << "]\n"
        "                    If 0, auto-calculated from relationship counts.\n"
        "  --dup INT         Duplicate / MZ-twin pairs            [" << DEF_DUP  << "]\n"
        "  --po INT          Extra parent-offspring pairs          [" << DEF_PO   << "]\n"
        "  --fs INT          Full-sibling families (2 children ea)[" << DEF_FS   << "]\n"
        "  --hs INT          Half-sibling pairs                   [" << DEF_HS   << "]\n"
        "  (Avuncular pairs arise automatically from FS families.)\n\n"
        "SNP configuration:\n"
        "  --snps INT        Number of SNPs                       [" << DEF_SNPS << "]\n"
        "  --af-min FLOAT    Minimum allele frequency             [" << DEF_AF_MIN << "]\n"
        "  --af-max FLOAT    Maximum allele frequency             [" << DEF_AF_MAX << "]\n\n"
        "Output:\n"
        "  --prefix STR      Output file prefix                   [kinship_test]\n"
        "  --output-dir DIR  Output directory                     [tests/data]\n"
        "  --seed INT        Random seed                          [" << DEF_SEED << "]\n\n"
        "  -h, --help        Show this help message\n";
}

// ── Relationship record ──────────────────────────────────────────────

struct GT {
    std::string ind1, ind2, rel;
    double k0, k1, k2, pihat;
};

// ── Main ─────────────────────────────────────────────────────────────

int main(int argc, char* argv[]) {

    // ── Parse CLI ────────────────────────────────────────────────────
    int  n_founders_req = DEF_FOUNDERS;
    int  n_snps         = DEF_SNPS;
    int  n_dup          = DEF_DUP;
    int  n_po           = DEF_PO;
    int  n_fs           = DEF_FS;
    int  n_hs           = DEF_HS;
    double af_min       = DEF_AF_MIN;
    double af_max       = DEF_AF_MAX;
    unsigned seed       = DEF_SEED;
    std::string prefix  = "kinship_test";
    std::string out_dir = "tests/data";

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next_int = [&]() -> int { return std::atoi(argv[++i]); };
        auto next_dbl = [&]() -> double { return std::atof(argv[++i]); };
        auto next_str = [&]() -> std::string { return argv[++i]; };

        if      (a == "--founders")   n_founders_req = next_int();
        else if (a == "--snps")       n_snps         = next_int();
        else if (a == "--dup")        n_dup          = next_int();
        else if (a == "--po")         n_po           = next_int();
        else if (a == "--fs")         n_fs           = next_int();
        else if (a == "--hs")         n_hs           = next_int();
        else if (a == "--af-min")     af_min         = next_dbl();
        else if (a == "--af-max")     af_max         = next_dbl();
        else if (a == "--seed")       seed           = next_int();
        else if (a == "--prefix")     prefix         = next_str();
        else if (a == "--output-dir") out_dir        = next_str();
        else if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else { std::cerr << "Unknown option: " << a << "\n"; usage(argv[0]); return 1; }
    }

    if (af_min <= 0 || af_max >= 1 || af_min >= af_max) {
        std::cerr << "Error: invalid AF range [" << af_min << ", " << af_max << "]\n";
        return 1;
    }

    // ── Calculate founders needed ────────────────────────────────────
    //   extra PO  : 2 founders each
    //   FS family : 2 founders each  (also generates PO & avuncular)
    //   HS pair   : 3 founders each
    int founders_for_po = 2 * n_po;
    int founders_for_fs = 2 * n_fs;
    int founders_for_hs = 3 * n_hs;
    int min_founders    = founders_for_po + founders_for_fs + founders_for_hs;
    int n_founders      = std::max(n_founders_req, min_founders);

    mkdir(out_dir.c_str(), 0755);
    std::string full_prefix = out_dir + "/" + prefix;

    std::mt19937 rng(seed);

    // ── 1. Sample names ─────────────────────────────────────────────
    auto pad = [](int n) -> std::string {
        std::string s = std::to_string(n);
        while (s.size() < 3) s = "0" + s;   // 001, 002, …
        return s;
    };
    auto f_name = [&](int idx) { return "F" + pad(idx + 1); };

    std::vector<std::string> names;
    for (int i = 0; i < n_founders; ++i) names.push_back(f_name(i));

    // ── 2. SNP metadata ─────────────────────────────────────────────
    std::vector<int>         positions(n_snps);
    std::vector<double>      true_af(n_snps);
    std::vector<std::string> ref_alleles(n_snps), alt_alleles(n_snps);
    std::vector<std::string> chroms(n_snps);

    std::uniform_real_distribution<double> af_dist(af_min, af_max);
    const char bases[] = {'A', 'C', 'G', 'T'};
    int n_chroms = std::max(1, n_snps / 250);        // ~250 SNPs per chromosome

    for (int s = 0; s < n_snps; ++s) {
        chroms[s]     = "chr" + std::to_string(1 + s / ((n_snps + n_chroms - 1) / n_chroms));
        positions[s]  = 10000 + (s % 250) * 4000;
        true_af[s]    = af_dist(rng);
        int ri = rng() % 4, ai;
        do { ai = rng() % 4; } while (ai == ri);
        ref_alleles[s] = std::string(1, bases[ri]);
        alt_alleles[s] = std::string(1, bases[ai]);
    }

    // ── 3. Founder genotypes (HWE) ──────────────────────────────────
    std::vector<std::vector<int>> all_geno;
    all_geno.reserve(n_founders + n_dup + 2*n_fs + 2*n_hs + 2*n_po + 20);
    for (int i = 0; i < n_founders; ++i) {
        std::vector<int> g;
        // Need to generate per-SNP; store transposed
        g.resize(n_snps);
        all_geno.push_back(std::move(g));
    }
    {
        std::uniform_real_distribution<double> U(0.0, 1.0);
        for (int s = 0; s < n_snps; ++s) {
            double p = true_af[s];
            double pAA = (1-p)*(1-p), pAB = 2*p*(1-p);
            for (int i = 0; i < n_founders; ++i) {
                double u = U(rng);
                all_geno[i][s] = (u < pAA) ? 0 : (u < pAA + pAB) ? 1 : 2;
            }
        }
    }

    std::vector<GT> gt;   // ground-truth relationships

    // Helper: allocate a new child from two parents
    auto add_child = [&](int p1_idx, int p2_idx, const std::string& child_name) -> int {
        std::vector<int> child(n_snps);
        for (int s = 0; s < n_snps; ++s)
            child[s] = make_child(all_geno[p1_idx][s], all_geno[p2_idx][s], rng);
        int idx = static_cast<int>(all_geno.size());
        names.push_back(child_name);
        all_geno.push_back(std::move(child));
        return idx;
    };

    int fi = 0;   // founder allocation cursor

    // ── 4a. Duplicates ──────────────────────────────────────────────
    for (int d = 0; d < n_dup; ++d) {
        int src = d % n_founders;   // cycle through founders
        std::string dup_name = "Dup_" + names[src];
        int idx = static_cast<int>(all_geno.size());
        names.push_back(dup_name);
        all_geno.push_back(all_geno[src]);   // exact copy
        gt.push_back({names[src], dup_name, "Duplicate/MZ", 0, 0, 1, 1});
    }

    // ── 4b. Extra Parent-Offspring pairs ────────────────────────────
    for (int p = 0; p < n_po; ++p) {
        int father = fi++;
        int mother = fi++;
        std::string child_name = "Child_" + f_name(father) + "_" + f_name(mother);
        int child = add_child(father, mother, child_name);
        gt.push_back({names[father], child_name, "Parent-Offspring", 0, 1, 0, 0.5});
        gt.push_back({names[mother], child_name, "Parent-Offspring", 0, 1, 0, 0.5});
    }

    // ── 4c. Full-Sibling families (each: 2 parents → 2 children) ────
    for (int f = 0; f < n_fs; ++f) {
        int father = fi++;
        int mother = fi++;
        std::string tag = std::to_string(f + 1);
        std::string c1_name = "FS" + tag + "_1";
        std::string c2_name = "FS" + tag + "_2";
        int c1 = add_child(father, mother, c1_name);
        int c2 = add_child(father, mother, c2_name);

        // Full-sibling pair
        gt.push_back({c1_name, c2_name, "Full-Sibling", 0.25, 0.5, 0.25, 0.5});
        // Parent-offspring (4 pairs)
        gt.push_back({names[father], c1_name, "Parent-Offspring", 0, 1, 0, 0.5});
        gt.push_back({names[father], c2_name, "Parent-Offspring", 0, 1, 0, 0.5});
        gt.push_back({names[mother], c1_name, "Parent-Offspring", 0, 1, 0, 0.5});
        gt.push_back({names[mother], c2_name, "Parent-Offspring", 0, 1, 0, 0.5});

        // Avuncular: if there's a spare founder, mate it with one child
        // to create an uncle/aunt–nephew/niece relationship with the other parent.
        if (fi < n_founders) {
            int spouse = fi++;
            std::string nephew_name = "Avunc_" + tag;
            int nephew = add_child(c1, spouse, nephew_name);
            // Father (c1) ↔ nephew : PO
            gt.push_back({c1_name, nephew_name, "Parent-Offspring", 0, 1, 0, 0.5});
            gt.push_back({names[spouse], nephew_name, "Parent-Offspring", 0, 1, 0, 0.5});
            // Father's parents ↔ nephew : Grandparent-Grandchild
            // (same IBD coefficients as avuncular: k0=0.5, k1=0.5, k2=0)
            gt.push_back({names[father], nephew_name, "Grandparent", 0.5, 0.5, 0, 0.25});
            gt.push_back({names[mother], nephew_name, "Grandparent", 0.5, 0.5, 0, 0.25});
            // Father's full-sibling ↔ nephew : Avuncular (uncle/aunt–nephew/niece)
            gt.push_back({c2_name, nephew_name, "Avuncular", 0.5, 0.5, 0, 0.25});
        }
    }

    // ── 4d. Half-Sibling pairs (each: shared father + 2 mothers) ─────
    for (int h = 0; h < n_hs; ++h) {
        int father  = fi++;
        int mother1 = fi++;
        int mother2 = fi++;
        std::string tag = std::to_string(h + 1);
        std::string c1_name = "HS" + tag + "_1";
        std::string c2_name = "HS" + tag + "_2";
        int c1 = add_child(father, mother1, c1_name);
        int c2 = add_child(father, mother2, c2_name);

        // Half-sibling pair
        gt.push_back({c1_name, c2_name, "Half-Sibling", 0.5, 0.5, 0, 0.25});
        // Parent-offspring (3 pairs)
        gt.push_back({names[father],  c1_name, "Parent-Offspring", 0, 1, 0, 0.5});
        gt.push_back({names[father],  c2_name, "Parent-Offspring", 0, 1, 0, 0.5});
        gt.push_back({names[mother1], c1_name, "Parent-Offspring", 0, 1, 0, 0.5});
        gt.push_back({names[mother2], c2_name, "Parent-Offspring", 0, 1, 0, 0.5});
    }

    // ── 4e. Unrelated pairs (sample from remaining founders) ────────
    {
        int n_unrelated = std::min(n_founders / 4, 10);
        n_unrelated = std::max(n_unrelated, std::min(4, n_founders));
        // Pick pairs spaced apart to minimize accidental LD
        for (int u = 0; u < n_unrelated && (2 * u + 1) < n_founders; ++u) {
            int i = (u * 2) % n_founders;
            int j = (u * 2 + 1 + n_founders / 3) % n_founders;
            if (i == j) j = (j + 1) % n_founders;
            gt.push_back({names[i], names[j], "Unrelated", 1, 0, 0, 0});
        }
    }

    int n_total = static_cast<int>(all_geno.size());

    // ── 5. Write .fam ────────────────────────────────────────────────
    {
        std::ofstream ofs(full_prefix + ".fam");
        for (int i = 0; i < n_total; ++i)
            ofs << "FAM" << (i + 1) << "\t" << names[i]
                << "\t0\t0\t0\t-9\n";
    }

    // ── 6. Write .bim ────────────────────────────────────────────────
    {
        std::ofstream ofs(full_prefix + ".bim");
        for (int s = 0; s < n_snps; ++s) {
            std::string snp_id = chroms[s] + "_" + std::to_string(positions[s]);
            ofs << chroms[s] << "\t" << snp_id << "\t0\t"
                << positions[s] << "\t"
                << alt_alleles[s] << "\t" << ref_alleles[s] << "\n";
        }
    }

    // ── 7. Write .bed ────────────────────────────────────────────────
    {
        std::ofstream ofs(full_prefix + ".bed", std::ios::binary);
        unsigned char magic[] = {0x6c, 0x1b, 0x01};
        ofs.write(reinterpret_cast<const char*>(magic), 3);

        int bps = (n_total + 3) / 4;
        std::vector<uint8_t> buf(bps, 0);
        for (int s = 0; s < n_snps; ++s) {
            std::fill(buf.begin(), buf.end(), 0);
            for (int j = 0; j < n_total; ++j) {
                uint8_t code;
                switch (all_geno[j][s]) {
                case 0: code = 0; break;   // hom_ref → PLINK code 0
                case 1: code = 2; break;   // het → PLINK code 2
                case 2: code = 3; break;   // hom_alt → PLINK code 3
                default: code = 1;
                }
                buf[j / 4] |= (code << ((j % 4) * 2));
            }
            ofs.write(reinterpret_cast<const char*>(buf.data()), bps);
        }
    }

    // ── 8. Write .frq ────────────────────────────────────────────────
    {
        std::ofstream ofs(full_prefix + ".frq");
        ofs << "CHR\tSNP\tA1\tA2\tMAF\tNCHROBS\n";
        for (int s = 0; s < n_snps; ++s) {
            std::string snp_id = chroms[s] + "_" + std::to_string(positions[s]);
            double af = true_af[s];
            double maf = std::min(af, 1.0 - af);
            std::string a1 = (af <= 0.5) ? alt_alleles[s] : ref_alleles[s];
            std::string a2 = (af <= 0.5) ? ref_alleles[s] : alt_alleles[s];
            ofs << chroms[s] << "\t" << snp_id << "\t" << a1 << "\t" << a2
                << "\t" << maf << "\t" << (n_total * 2) << "\n";
        }
    }

    // ── 9. Write ground_truth.tsv ────────────────────────────────────
    {
        std::ofstream ofs(out_dir + "/ground_truth.tsv");
        ofs << "Ind1\tInd2\tRelationship\tExpected_k0\tExpected_k1\t"
               "Expected_k2\tExpected_PI_HAT\n";
        for (const auto& g : gt)
            ofs << g.ind1 << "\t" << g.ind2 << "\t" << g.rel << "\t"
                << g.k0 << "\t" << g.k1 << "\t" << g.k2 << "\t"
                << g.pihat << "\n";
    }

    // ── 10. Summary ──────────────────────────────────────────────────
    // Count by type
    int cnt[7] = {};
    const char* type_names[] = {
        "Duplicate/MZ", "Parent-Offspring", "Full-Sibling",
        "Half-Sibling", "Avuncular", "Grandparent", "Unrelated"
    };
    for (const auto& g : gt) {
        for (int t = 0; t < 7; ++t)
            if (g.rel == type_names[t]) { ++cnt[t]; break; }
    }

    std::cerr << "[generate_known_relationships] Done.\n"
              << "  " << n_total << " samples (" << n_founders
              << " founders + " << (n_total - n_founders) << " derived)\n"
              << "  " << n_snps << " SNPs, AF ∈ [" << af_min << ", "
              << af_max << "], seed=" << seed << "\n"
              << "  " << gt.size() << " ground-truth pairs:\n";
    for (int t = 0; t < 7; ++t)
        if (cnt[t] > 0)
            std::cerr << "    " << type_names[t] << ": " << cnt[t] << "\n";
    std::cerr << "  Output: " << full_prefix << ".{bed,bim,fam,frq}\n"
              << "          " << out_dir << "/ground_truth.tsv\n";

    return 0;
}
