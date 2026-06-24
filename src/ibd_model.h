#ifndef _FASTLCKIN_IBD_MODEL_H_
#define _FASTLCKIN_IBD_MODEL_H_

/**
 * @file ibd_model.h
 * @brief Anderson & Weir (2007) IBS|IBD conditional probability model
 * @author Shujia Huang
 * @date 2026-06-23
 */

#include <vector>
#include <array>
#include <string>

#include <htslib/hts.h>

namespace fastlckin {

/// SNP metadata
struct SNPInfo {
    std::string chrom;      ///< Chromosome name
    hts_pos_t pos = 0;      ///< 0-based position
    std::string id;         ///< CHR_POS format ID
    std::string ref;        ///< Reference allele
    std::string alt;        ///< Alternate allele (biallelic only)
    double af = 0.0;        ///< ALT allele frequency (freq of gl[2] allele)
    bool af_masked = false; ///< AF filter flag
};

/// 9 genotype combination encoding
enum GenotypeCombo : int {
    PPQQ = 0, QQPP = 1,                        // IBS=0
    PPPQ = 2, PQPP = 3, PQQQ = 4, QQPQ = 5,   // IBS=1
    PPPP = 6, PQPQ = 7, QQQQ = 8                // IBS=2
};

/// IBS|IBD probability matrix: [combo][snp][ibd_state]
/// ibd_state: 0=IBD0, 1=IBD1, 2=IBD2
using IBS_IBD_Matrix = std::vector<std::vector<std::array<double, 3>>>;

/// Mij function: FST-corrected allele frequency
/// M(p, FST, i) = (1 - FST) * p + i * FST
inline double Mij(double p, double fst, int i) {
    return (1.0 - fst) * p + i * fst;
}

/// Compute 9 genotype combos × 3 IBD states for a single SNP
/// @param af  ALT allele frequency (P = gl[2] allele)
/// @param fst FST value
/// @return array[combo][ibd_state]
std::array<std::array<double, 3>, 9> compute_ibs_ibd_probs(double af, double fst);

/// Batch pre-compute IBS|IBD probability matrix for all SNPs
/// @param snp_infos  SNP metadata list (with AF)
/// @param fst  FST value
/// @return IBS_IBD_Matrix [9][N][3]
IBS_IBD_Matrix precompute_ibs_ibd_matrix(const std::vector<SNPInfo>& snp_infos, double fst);

/// Batch pre-compute with per-SNP FST vector (v0.4.0)
/// @param snp_infos  SNP metadata list (with AF)
/// @param fst_vector  Per-SNP FST values (size must match snp_infos)
/// @return IBS_IBD_Matrix [9][N][3]
IBS_IBD_Matrix precompute_ibs_ibd_matrix_fst_vector(
    const std::vector<SNPInfo>& snp_infos,
    const std::vector<double>& fst_vector
);

}  // namespace fastlckin

#endif
