#ifndef _FASTLCKIN_GENOTYPE_LIKELIHOOD_H_
#define _FASTLCKIN_GENOTYPE_LIKELIHOOD_H_

/**
 * @file genotype_likelihood.h
 * @brief GL/PL extraction and quality control from VCF records
 * @author Shujia Huang
 * @date 2025-06-23
 */

#include <vector>
#include <string>
#include <cmath>

#include "io/vcf.h"
#include "io/vcf_record.h"
#include "io/vcf_header.h"

namespace fastlckin {

/// Per-individual per-SNP genotype likelihoods
/// Index: 0 = pp (REF/REF), 1 = pq (REF/ALT), 2 = qq (ALT/ALT)
struct GenotypeLikelihoods {
    double gl[3] = {1.0, 1.0, 1.0};  ///< Linear-scale likelihoods P(Data|G)
    int gq = 0;        ///< Genotype quality (Phred scale)
    bool masked = false; ///< QC failed flag
};

/// All samples × all SNPs likelihood matrix
/// Storage: [sample_idx][snp_idx] -> GenotypeLikelihoods
using GLMatrix = std::vector<std::vector<GenotypeLikelihoods>>;

/// Extract genotype likelihoods for all samples from a single VCF record.
/// PL is preferred over GL. If neither exists, all samples are masked.
/// @param rec  VCF record (must have unpacked BCF_UN_FMT)
/// @param hdr  VCF header
/// @param n_samples  Number of samples
/// @param gq_min  Minimum GQ threshold (samples below this are masked)
/// @return Vector of GenotypeLikelihoods, length = n_samples
std::vector<GenotypeLikelihoods> extract_genotype_likelihoods(
    const ngslib::VCFRecord& rec,
    const ngslib::VCFHeader& hdr,
    int n_samples,
    int gq_min = 1
);

}  // namespace fastlckin

#endif
