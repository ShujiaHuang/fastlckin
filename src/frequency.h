#ifndef _FASTLCKIN_FREQUENCY_H_
#define _FASTLCKIN_FREQUENCY_H_

/**
 * @file frequency.h
 * @brief Allele frequency computation and file I/O
 * @author Shujia Huang
 * @date 2025-06-23
 */

#include <vector>
#include <string>
#include "ibd_model.h"  // SNPInfo

namespace fastlckin {

/// Load SNP metadata from PLINK .bim file
/// @param bim_path  Path to .bim file
/// @return Vector of SNPInfo (chrom, pos, id, ref, alt; af not yet filled)
std::vector<SNPInfo> load_bim_snps(const std::string& bim_path);

/// Read allele frequencies from PLINK .frq file
/// Matches REF/ALT direction and warns on mismatches.
/// @param frq_path   Path to .frq file
/// @param snp_infos  SNP metadata (for matching REF/ALT)
/// @return Frequency vector aligned to snp_infos order; -1.0 for unmatched
std::vector<double> read_plink_freq(
    const std::string& frq_path,
    const std::vector<SNPInfo>& snp_infos
);

/// Compute allele frequencies from PLINK .bed file
/// @param bed_path     Path to .bed file
/// @param bim_path     Path to .bim file
/// @param snp_indices  Indices of SNPs to compute (into .bim order)
/// @return Frequency vector (alt allele freq) for each requested SNP
std::vector<double> compute_allele_frequencies(
    const std::string& bed_path,
    const std::string& bim_path,
    const std::vector<int>& snp_indices
);

/// Read genotypes from PLINK .bed file for specified SNPs
/// @param bed_path     Path to .bed file
/// @param bim_path     Path to .bim file
/// @param snp_indices  Indices of SNPs to extract
/// @return Genotype matrix [n_samples][n_snps], values 0/1/2, -1=missing
std::vector<std::vector<int8_t>> read_bed_genotypes(
    const std::string& bed_path,
    const std::string& bim_path,
    const std::vector<int>& snp_indices
);

}  // namespace fastlckin

#endif
