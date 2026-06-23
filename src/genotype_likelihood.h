#ifndef _FASTLCKIN_GENOTYPE_LIKELIHOOD_H_
#define _FASTLCKIN_GENOTYPE_LIKELIHOOD_H_

/**
 * @file genotype_likelihood.h
 * @brief Genotype likelihood extraction and quality control from VCF records
 *
 * "Genotype likelihood" here refers to the general concept P(Data|Genotype),
 * regardless of the VCF FORMAT field used. By default, the PL (Phred-scaled)
 * field is preferred because it is more commonly available (e.g., GATK,
 * DeepVariant). The GL (log10-scaled) field is used as a fallback.
 * Use the `pl_field` parameter to specify a custom field name.
 */

#include <vector>
#include <string>
#include <cmath>

#include "io/vcf.h"
#include "io/vcf_record.h"
#include "io/vcf_header.h"

namespace fastlckin {

/// Per-individual per-SNP genotype likelihoods.
/// Index: 0 = pp (REF/REF), 1 = pq (REF/ALT), 2 = qq (ALT/ALT)
struct GenotypeLikelihood {
    double gl[3] = {1.0, 1.0, 1.0};  ///< Linear-scale likelihoods P(Data|G)
    int gq = 0;        ///< Genotype quality (Phred scale)
    bool masked = false; ///< QC failed flag
};

/// All samples × all SNPs likelihood matrix
/// Storage: [sample_idx][snp_idx] -> GenotypeLikelihood
using LikelihoodMatrix = std::vector<std::vector<GenotypeLikelihood>>;

/// Extract genotype likelihoods for all samples from a single VCF record.
///
/// By default, PL (Phred-scaled) is preferred over GL (log10-scaled).
/// Use `pl_field` to specify a custom FORMAT field name for Phred-scaled
/// genotype likelihoods (default: "PL"). If the specified field is not
/// found, GL is tried as a fallback.
///
/// @param rec         VCF record (must have unpacked BCF_UN_FMT)
/// @param hdr         VCF header
/// @param n_samples   Number of samples
/// @param gq_min      Minimum GQ threshold (samples below this are masked)
/// @param pl_field    VCF FORMAT field name for Phred-scaled GL (default: "PL")
/// @return Vector of GenotypeLikelihood, length = n_samples
std::vector<GenotypeLikelihood> extract_genotype_likelihoods(
    const ngslib::VCFRecord& rec,
    const ngslib::VCFHeader& hdr,
    int n_samples,
    int gq_min = 1,
    const std::string& pl_field = "PL"
);

}  // namespace fastlckin

#endif
