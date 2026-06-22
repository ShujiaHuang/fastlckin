/**
 * @file frequency.cpp
 * @brief Allele frequency computation, .bed/.bim/.frq I/O
 * @author Shujia Huang
 * @date 2025-06-23
 */

#include "frequency.h"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <map>
#include <cstring>

namespace fastlckin {

std::vector<SNPInfo> load_bim_snps(const std::string& bim_path) {
    std::ifstream ifs(bim_path);
    if (!ifs) throw std::runtime_error("Cannot open .bim file: " + bim_path);

    std::vector<SNPInfo> snps;
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        SNPInfo s;
        int cm;
        std::string a1, a2;
        iss >> s.chrom >> s.id >> cm >> s.pos >> a1 >> a2;

        // In PLINK .bim: A1 is minor, A2 is major
        // We store ref=A2, alt=A1 (alt = the one whose freq we track)
        s.ref = a2;
        s.alt = a1;
        s.af = 0.0;
        s.af_masked = false;
        snps.push_back(std::move(s));
    }
    return snps;
}

// Read .fam file to get sample count
static int read_fam_count(const std::string& fam_path) {
    std::ifstream ifs(fam_path);
    if (!ifs) throw std::runtime_error("Cannot open .fam file: " + fam_path);
    int count = 0;
    std::string line;
    while (std::getline(ifs, line)) {
        if (!line.empty()) ++count;
    }
    return count;
}

std::vector<double> read_plink_freq(
    const std::string& frq_path,
    const std::vector<SNPInfo>& snp_infos)
{
    std::ifstream ifs(frq_path);
    if (!ifs) throw std::runtime_error("Cannot open .frq file: " + frq_path);

    // Build lookup: CHR_POS -> index in snp_infos
    std::map<std::string, size_t> snp_index_map;
    for (size_t i = 0; i < snp_infos.size(); ++i) {
        std::string key = snp_infos[i].chrom + "_" + std::to_string(snp_infos[i].pos);
        snp_index_map[key] = i;
    }

    std::vector<double> freqs(snp_infos.size(), -1.0);

    std::string line;
    // Skip header line
    std::getline(ifs, line);

    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        std::istringstream iss(line);
        std::string chr, snp_id, a1, a2;
        double maf;
        int nchroms;
        iss >> chr >> snp_id >> a1 >> a2 >> maf >> nchroms;

        // Try to match by position using snp_id format CHR_POS
        // Or try to match by the snp_id directly
        size_t idx = snp_infos.size();

        // Try CHR_POS key
        // The snp_id in .frq might be different from our internal ID
        // We try matching by checking all snp_infos
        for (size_t i = 0; i < snp_infos.size(); ++i) {
            if (snp_infos[i].id == snp_id ||
                (snp_infos[i].chrom == chr &&
                 snp_infos[i].pos == std::stol(snp_id.substr(snp_id.find('_') + 1)))) {
                idx = i;
                break;
            }
        }

        if (idx >= snp_infos.size()) {
            // Also try position-based match
            for (size_t i = 0; i < snp_infos.size(); ++i) {
                if (snp_infos[i].chrom == chr) {
                    // Try to extract position from snp_id
                    auto upos = snp_id.find('_');
                    if (upos != std::string::npos) {
                        try {
                            long snp_pos = std::stol(snp_id.substr(upos + 1));
                            if (static_cast<long>(snp_infos[i].pos) == snp_pos) {
                                idx = i;
                                break;
                            }
                        } catch (...) {}
                    }
                }
            }
        }

        if (idx >= snp_infos.size()) continue;

        const auto& si = snp_infos[idx];

        // Match REF/ALT direction
        if (si.alt == a1 && si.ref == a2) {
            // A1 in .frq matches our ALT → MAF is alt allele freq
            freqs[idx] = maf;
        } else if (si.alt == a2 && si.ref == a1) {
            // A1 in .frq matches our REF → need 1-MAF
            freqs[idx] = 1.0 - maf;
        } else {
            std::cerr << "[fastlckin] Warning: allele mismatch at " << si.id
                      << " (VCF ref=" << si.ref << " alt=" << si.alt
                      << ", FRQ A1=" << a1 << " A2=" << a2 << "). Masking.\n";
            freqs[idx] = -1.0;
        }
    }

    return freqs;
}

std::vector<std::vector<int8_t>> read_bed_genotypes(
    const std::string& bed_path,
    const std::string& bim_path,
    const std::vector<int>& snp_indices)
{
    // Get sample count from .fam
    std::string fam_path = bed_path.substr(0, bed_path.size() - 4) + ".fam";
    int n_samples = read_fam_count(fam_path);

    // Total SNPs from .bim
    auto all_snps = load_bim_snps(bim_path);
    int n_total_snps = static_cast<int>(all_snps.size());

    int bytes_per_snp = (n_samples + 3) / 4;

    std::ifstream ifs(bed_path, std::ios::binary);
    if (!ifs) throw std::runtime_error("Cannot open .bed file: " + bed_path);

    // Check magic number
    unsigned char magic[3];
    ifs.read(reinterpret_cast<char*>(magic), 3);
    if (magic[0] != 0x6c || magic[1] != 0x1b || magic[2] != 0x01) {
        throw std::runtime_error("Invalid .bed magic number. Expected SNP-major mode.");
    }

    int n_snps = static_cast<int>(snp_indices.size());
    std::vector<std::vector<int8_t>> genotypes(n_samples, std::vector<int8_t>(n_snps, -1));

    std::vector<uint8_t> buffer(bytes_per_snp);

    for (int si = 0; si < n_snps; ++si) {
        int snp_idx = snp_indices[si];
        if (snp_idx < 0 || snp_idx >= n_total_snps) continue;

        // Seek to SNP data
        ifs.seekg(3 + static_cast<std::streamoff>(snp_idx) * bytes_per_snp);
        ifs.read(reinterpret_cast<char*>(buffer.data()), bytes_per_snp);

        for (int j = 0; j < n_samples; ++j) {
            int byte_idx = j / 4;
            int bit_idx = (j % 4) * 2;
            uint8_t code = (buffer[byte_idx] >> bit_idx) & 0x03;

            // PLINK encoding: 00=hom_first(0), 01=missing(-1), 10=het(1), 11=hom_second(2)
            switch (code) {
            case 0: genotypes[j][si] = 0; break;  // hom A1A1
            case 1: genotypes[j][si] = -1; break;  // missing
            case 2: genotypes[j][si] = 1; break;  // het
            case 3: genotypes[j][si] = 2; break;  // hom A2A2
            }
        }
    }

    return genotypes;
}

std::vector<double> compute_allele_frequencies(
    const std::string& bed_path,
    const std::string& bim_path,
    const std::vector<int>& snp_indices)
{
    auto genotypes = read_bed_genotypes(bed_path, bim_path, snp_indices);
    int n_samples = static_cast<int>(genotypes.size());
    int n_snps = static_cast<int>(snp_indices.size());

    std::vector<double> freqs(n_snps, 0.0);

    for (int s = 0; s < n_snps; ++s) {
        int allele_count = 0;
        int valid_count = 0;
        for (int i = 0; i < n_samples; ++i) {
            int8_t g = genotypes[i][s];
            if (g < 0) continue;
            // In PLINK .bed: 0=A1A1, 1=A1A2, 2=A2A2
            // We want alt(A1) frequency. A1 is .bim A1 column.
            // genotype value 0 means hom_A1A1 → 2 A1 alleles
            // genotype value 1 means het → 1 A1 allele
            // genotype value 2 means hom_A2A2 → 0 A1 alleles
            allele_count += (2 - g);
            valid_count += 2;
        }
        freqs[s] = (valid_count > 0) ? static_cast<double>(allele_count) / valid_count : 0.0;
    }

    return freqs;
}

}  // namespace fastlckin
