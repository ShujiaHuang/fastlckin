/**
 * @file generate_test_data.cpp
 * @brief Generate synthetic PLINK .bed/.bim/.fam and VCF test data
 *
 * This program creates minimal test datasets for unit testing.
 * Generates: tests/data/test.bed, test.bim, test.fam, test.vcf, test.frq
 *
 * Run: generate_test_data
 */

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <sys/stat.h>

// Generate synthetic test data
int main() {
    const std::string dir = "tests/data";
    mkdir(dir.c_str(), 0755);

    const int n_samples = 20;
    const int n_snps = 50;
    const std::string chrom = "chr1";

    std::vector<std::string> sample_names;
    for (int i = 0; i < n_samples; ++i) {
        sample_names.push_back("Sample" + std::to_string(i + 1));
    }

    // SNP positions (evenly spaced)
    std::vector<int> positions(n_snps);
    std::vector<std::string> ref_alleles(n_snps);
    std::vector<std::string> alt_alleles(n_snps);

    std::mt19937 rng(12345);
    std::uniform_int_distribution<int> pos_dist(10000, 1000000);
    const char bases[] = {'A', 'C', 'G', 'T'};

    for (int s = 0; s < n_snps; ++s) {
        positions[s] = 10000 + s * 1000;
        int ri = rng() % 4;
        int ai;
        do { ai = rng() % 4; } while (ai == ri);
        ref_alleles[s] = std::string(1, bases[ri]);
        alt_alleles[s] = std::string(1, bases[ai]);
    }

    // Generate random genotypes (0,1,2) with allele frequency ~0.3
    std::vector<std::vector<int>> genotypes(n_samples, std::vector<int>(n_snps));
    std::uniform_real_distribution<double> uni(0.0, 1.0);

    double target_af = 0.3;
    std::vector<double> true_af(n_snps);
    for (int s = 0; s < n_snps; ++s) {
        // Vary AF per SNP
        true_af[s] = 0.1 + (s % 8) * 0.1;
        double p = true_af[s];
        for (int i = 0; i < n_samples; ++i) {
            double u = uni(rng);
            if (u < (1.0 - p) * (1.0 - p)) {
                genotypes[i][s] = 0;  // hom ref
            } else if (u < (1.0 - p) * (1.0 - p) + 2.0 * p * (1.0 - p)) {
                genotypes[i][s] = 1;  // het
            } else {
                genotypes[i][s] = 2;  // hom alt
            }
        }
    }

    // ── Write .fam ──
    {
        std::ofstream ofs(dir + "/test.fam");
        for (int i = 0; i < n_samples; ++i) {
            ofs << "FAM" << (i + 1) << "\t" << sample_names[i]
                << "\t0\t0\t0\t-9\n";
        }
    }

    // ── Write .bim ──
    {
        std::ofstream ofs(dir + "/test.bim");
        for (int s = 0; s < n_snps; ++s) {
            std::string snp_id = chrom + "_" + std::to_string(positions[s]);
            // PLINK .bim: chr snp_id cm pos A1(minor) A2(major)
            ofs << chrom << "\t" << snp_id << "\t0\t" << positions[s]
                << "\t" << alt_alleles[s] << "\t" << ref_alleles[s] << "\n";
        }
    }

    // ── Write .bed (SNP-major mode) ──
    {
        std::ofstream ofs(dir + "/test.bed", std::ios::binary);
        // Magic bytes for SNP-major
        unsigned char magic[] = {0x6c, 0x1b, 0x01};
        ofs.write(reinterpret_cast<const char*>(magic), 3);

        int bytes_per_snp = (n_samples + 3) / 4;
        std::vector<uint8_t> buffer(bytes_per_snp, 0);

        for (int s = 0; s < n_snps; ++s) {
            std::fill(buffer.begin(), buffer.end(), 0);
            for (int j = 0; j < n_samples; ++j) {
                int byte_idx = j / 4;
                int bit_idx = (j % 4) * 2;

                uint8_t code;
                switch (genotypes[j][s]) {
                case 0: code = 0; break;   // hom_ref → PLINK code 0
                case 1: code = 2; break;   // het → PLINK code 2
                case 2: code = 3; break;   // hom_alt → PLINK code 3
                default: code = 1; break;  // missing
                }

                buffer[byte_idx] |= (code << bit_idx);
            }
            ofs.write(reinterpret_cast<const char*>(buffer.data()), bytes_per_snp);
        }
    }

    // ── Write .frq ──
    {
        std::ofstream ofs(dir + "/test.frq");
        ofs << "CHR\tSNP\tA1\tA2\tMAF\tNCHROBS\n";
        for (int s = 0; s < n_snps; ++s) {
            std::string snp_id = chrom + "_" + std::to_string(positions[s]);
            double af = true_af[s];
            double maf = std::min(af, 1.0 - af);
            std::string a1 = (af <= 0.5) ? alt_alleles[s] : ref_alleles[s];
            std::string a2 = (af <= 0.5) ? ref_alleles[s] : alt_alleles[s];
            ofs << chrom << "\t" << snp_id << "\t" << a1 << "\t" << a2
                << "\t" << maf << "\t" << (n_samples * 2) << "\n";
        }
    }

    // ── Write VCF ──
    {
        std::ofstream ofs(dir + "/test.vcf");
        // Header
        ofs << "##fileformat=VCFv4.2\n";
        ofs << "##contig=<ID=" << chrom << ",length=1000000>\n";
        ofs << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
        ofs << "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">\n";
        ofs << "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
        ofs << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
        for (const auto& s : sample_names) ofs << "\t" << s;
        ofs << "\n";

        std::uniform_int_distribution<int> pl_dist(0, 50);

        for (int s = 0; s < n_snps; ++s) {
            std::string snp_id = chrom + "_" + std::to_string(positions[s]);
            ofs << chrom << "\t" << positions[s] << "\t" << snp_id
                << "\t" << ref_alleles[s] << "\t" << alt_alleles[s]
                << "\t100\tPASS\t.\tGT:PL:GQ";

            for (int i = 0; i < n_samples; ++i) {
                int g = genotypes[i][s];
                // GT
                std::string gt;
                switch (g) {
                case 0: gt = "0/0"; break;
                case 1: gt = "0/1"; break;
                case 2: gt = "1/1"; break;
                default: gt = "./."; break;
                }

                // PL: simulate realistic Phred-scaled likelihoods
                int pl0, pl1, pl2;
                if (g == 0) {
                    pl0 = 0; pl1 = 20 + pl_dist(rng); pl2 = 40 + pl_dist(rng);
                } else if (g == 1) {
                    pl0 = 20 + pl_dist(rng); pl1 = 0; pl2 = 20 + pl_dist(rng);
                } else {
                    pl0 = 40 + pl_dist(rng); pl1 = 20 + pl_dist(rng); pl2 = 0;
                }

                // GQ: min of non-best PL values
                int gq = std::min({pl0, pl1, pl2});
                int sorted[3] = {pl0, pl1, pl2};
                std::sort(sorted, sorted + 3);
                gq = sorted[1] - sorted[0];

                ofs << "\t" << gt << ":" << pl0 << "," << pl1 << "," << pl2
                    << ":" << gq;
            }
            ofs << "\n";
        }
    }

    std::cerr << "[generate_test_data] Generated test data in " << dir << "/\n";
    std::cerr << "  " << n_samples << " samples, " << n_snps << " SNPs\n";
    std::cerr << "  Files: test.bed, test.bim, test.fam, test.vcf, test.frq\n";

    return 0;
}
