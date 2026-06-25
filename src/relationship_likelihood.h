#ifndef _FASTLCKIN_RELATIONSHIP_LIKELIHOOD_H_
#define _FASTLCKIN_RELATIONSHIP_LIKELIHOOD_H_

/**
 * @file relationship_likelihood.h
 * @brief Exact per-relationship genotype joint probability P(G1, G2 | R)
 *
 * For each standard relationship type, computes the exact probability of
 * observing a pair of genotypes at a biallelic SNP, given the relationship.
 *
 * Each relationship's likelihood is derived from Mendelian transmission
 * probabilities summed over shared ancestors, replacing the Mij model's
 * population-averaged IBD=1 approximation with exact pedigree-based formulas.
 *
 * Supported relationship types (9):
 *   UNRELATED, DUPLICATE/MZ, PARENT-OFFSPRING, FULL-SIBLING,
 *   HALF-SIBLING, AVUNCULAR, GRANDPARENT, FIRST-COUSIN,
 *   GREAT-GRANDPARENT
 *
 * @author Shujia Huang
 * @date 2026-06-24
 */

#include <vector>
#include <array>

#include "ibd_model.h"

namespace fastlckin {

/// Number of supported relationship types
constexpr int NUM_RELATIONSHIP_TYPES = 9;

/// Number of IBD equivalence classes (relationships with identical likelihoods)
constexpr int NUM_IBD_CLASSES = 6;

/// Relationship type enumeration
enum class RelationshipType : int {
    UNRELATED = 0,
    DUPLICATE_MZ,
    PARENT_OFFSPRING,
    FULL_SIBLING,
    HALF_SIBLING,
    AVUNCULAR,
    GRANDPARENT,
    FIRST_COUSIN,
    GREAT_GRANDPARENT
};

/// IBD equivalence class enumeration
/// Due to the "random spouse" assumption in the likelihood model, relationships
/// with the same IBD coefficients (k0, k1, k2) produce identical likelihoods.
/// We group them into equivalence classes for model selection.
enum class IBDEquivalenceClass : int {
    UNRELATED = 0,       // k = (1.0, 0.0, 0.0)
    DUPLICATE_MZ,        // k = (0.0, 0.0, 1.0)
    PARENT_OFFSPRING,    // k = (0.0, 1.0, 0.0)
    FULL_SIBLING,        // k = (0.25, 0.5, 0.25)
    SECOND_DEGREE,       // k = (0.5, 0.5, 0.0): HS, GP, AV
    THIRD_DEGREE         // k = (0.75, 0.25, 0.0): FC, GGP
};

/// Map relationship type to IBD equivalence class
inline IBDEquivalenceClass relationship_to_ibd_class(RelationshipType t) {
    switch (t) {
        case RelationshipType::UNRELATED:         return IBDEquivalenceClass::UNRELATED;
        case RelationshipType::DUPLICATE_MZ:      return IBDEquivalenceClass::DUPLICATE_MZ;
        case RelationshipType::PARENT_OFFSPRING:  return IBDEquivalenceClass::PARENT_OFFSPRING;
        case RelationshipType::FULL_SIBLING:      return IBDEquivalenceClass::FULL_SIBLING;
        case RelationshipType::HALF_SIBLING:
        case RelationshipType::AVUNCULAR:
        case RelationshipType::GRANDPARENT:       return IBDEquivalenceClass::SECOND_DEGREE;
        case RelationshipType::FIRST_COUSIN:
        case RelationshipType::GREAT_GRANDPARENT: return IBDEquivalenceClass::THIRD_DEGREE;
        default:                                  return IBDEquivalenceClass::UNRELATED;
    }
}

/// Return human-readable name for an IBD equivalence class
inline const char* ibd_class_name(IBDEquivalenceClass c) {
    switch (c) {
        case IBDEquivalenceClass::UNRELATED:       return "Unrelated";
        case IBDEquivalenceClass::DUPLICATE_MZ:    return "Duplicate/MZ";
        case IBDEquivalenceClass::PARENT_OFFSPRING:return "Parent-Offspring";
        case IBDEquivalenceClass::FULL_SIBLING:    return "Full-Sibling";
        case IBDEquivalenceClass::SECOND_DEGREE:   return "2nd-Degree";
        case IBDEquivalenceClass::THIRD_DEGREE:    return "3rd-Degree";
        default:                                   return "Unknown";
    }
}

/// Exact relationship genotype probability matrix: [rel_type][combo][snp]
/// Stores P(G1=g1, G2=g2 | R) for each of the 9 genotype combinations,
/// each SNP, and each of the 9 relationship types.
using RelationshipLikelihoodMatrix =
    std::vector<std::vector<std::vector<double>>>;

/// Compute exact P(G1, G2 | R) for all 9 relationship types at a single SNP.
///
/// Uses Mendelian transmission probabilities summed over shared ancestors.
/// For FST-corrected allele frequencies, uses:
///   p_0 = (1-F)*p      (het path, no coancestry correction)
///   p_1 = (1-F)*p + F  (hom path, with coancestry correction)
///
/// @param af  Alternate allele frequency
/// @param fst FST value
/// @return array[relationship_type][9 genotype combos]
std::array<std::array<double, NUM_RELATIONSHIP_TYPES>, 9>
compute_relationship_probs_single_snp(double af, double fst);

/// Batch pre-compute exact relationship likelihoods for all SNPs.
///
/// @param snp_infos  SNP metadata list (with AF)
/// @param fst        FST value
/// @return RelationshipLikelihoodMatrix [NUM_RELATIONSHIP_TYPES][9][N]
RelationshipLikelihoodMatrix precompute_relationship_likelihoods(
    const std::vector<SNPInfo>& snp_infos, double fst);

/// Compute log-likelihood for a single pair under a specific relationship type.
///
/// log L(R) = sum_s log( sum_c PIBS[c] * P(G1,G2 | R, SNP_s) )
///
/// @param rel_matrix   Pre-computed relationship likelihood matrix
/// @param pibs_per_snp PIBS values per SNP for this pair
/// @param snp_indices  LD-pruned (or all unmasked) SNP indices
/// @param type         Relationship type to evaluate
/// @return Log-likelihood value
double compute_relationship_log_likelihood(
    const RelationshipLikelihoodMatrix& rel_matrix,
    const std::vector<std::array<double, 9>>& pibs_per_snp,
    const std::vector<int>& snp_indices,
    RelationshipType type);

}  // namespace fastlckin

#endif
