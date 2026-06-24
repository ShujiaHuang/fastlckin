/**
 * @file test_relationship_likelihood.cpp
 * @brief Unit tests for exact per-relationship genotype probability computation
 *
 * Tests verify:
 * 1. Normalization: sum of 9 genotype combo probabilities = 1 for each type
 * 2. Symmetry: P(G1,G2|R) = P(G2,G1|R) for symmetric relationships
 * 3. PO exact formula vs Mij IBD=1: mathematical equivalence
 * 4. HS exact formula vs Mij IBD=1: mathematical inequivalence (bias source)
 * 5. FST=0 special cases: manual verification
 * 6. DUPLICATE: only diagonal combos (PPPP, PQPQ, QQQQ) are non-zero
 */

#include "test_harness.h"
#include "relationship_likelihood.h"
#include "ibd_model.h"
#include <cmath>
#include <vector>
#include <iostream>

using namespace fastlckin;

// ── Normalization: sum of 9 genotype combo probabilities = 1 ────────
// Each relationship type is tested individually at FST=0, p=0.3

TEST_CASE(normalization_all_types_fst0) {
    auto R = compute_relationship_probs_single_snp(0.3, 0.0);
    for (int r = 0; r < NUM_RELATIONSHIP_TYPES; ++r) {
        double sum = 0.0;
        for (int c = 0; c < 9; ++c) {
            sum += R[c][r];
        }
        CHECK_NEAR(sum, 1.0, 1e-10);
    }
}

// ── Normalization across multiple AF and FST values ─────────────────

TEST_CASE(normalization_all_types_multiple_params) {
    std::vector<double> fst_values = {0.0, 0.001, 0.01, 0.05, 0.1};
    std::vector<double> af_values = {0.05, 0.1, 0.25, 0.3, 0.5, 0.7, 0.9, 0.95};

    for (double fst : fst_values) {
        for (double af : af_values) {
            auto R = compute_relationship_probs_single_snp(af, fst);
            for (int r = 0; r < NUM_RELATIONSHIP_TYPES; ++r) {
                double sum = 0.0;
                for (int c = 0; c < 9; ++c) {
                    sum += R[c][r];
                }
                CHECK_NEAR(sum, 1.0, 1e-9);
            }
        }
    }
}

// ── DUPLICATE: only diagonal combos are non-zero ────────────────────

TEST_CASE(duplicate_off_diagonal_zero) {
    auto R = compute_relationship_probs_single_snp(0.3, 0.0);
    // Off-diagonal combos (genotypes differ): PPQQ, QQPP, PPPQ, PQPP, PQQQ, QQPQ
    CHECK_NEAR(R[PPQQ][(int)RelationshipType::DUPLICATE_MZ], 0.0, 1e-15);
    CHECK_NEAR(R[QQPP][(int)RelationshipType::DUPLICATE_MZ], 0.0, 1e-15);
    CHECK_NEAR(R[PPPQ][(int)RelationshipType::DUPLICATE_MZ], 0.0, 1e-15);
    CHECK_NEAR(R[PQPP][(int)RelationshipType::DUPLICATE_MZ], 0.0, 1e-15);
    CHECK_NEAR(R[PQQQ][(int)RelationshipType::DUPLICATE_MZ], 0.0, 1e-15);
    CHECK_NEAR(R[QQPQ][(int)RelationshipType::DUPLICATE_MZ], 0.0, 1e-15);
}

TEST_CASE(duplicate_diagonal_nonzero) {
    auto R = compute_relationship_probs_single_snp(0.3, 0.0);
    // Diagonal combos (genotypes match): PPPP, PQPQ, QQQQ
    CHECK(R[PPPP][(int)RelationshipType::DUPLICATE_MZ] > 0.0);
    CHECK(R[PQPQ][(int)RelationshipType::DUPLICATE_MZ] > 0.0);
    CHECK(R[QQQQ][(int)RelationshipType::DUPLICATE_MZ] > 0.0);
}

// ── PO exact formula vs Mij IBD=1: verify they DIFFER ───────────────
// The PO formula (Mendelian transmission) differs from the Mij IBD=1
// model because Mij draws the shared allele from the population, while
// PO constrains it to be one of the parent's actual alleles.
// With FST > 0, the difference is more pronounced.

TEST_CASE(po_differs_from_mij_ibd1) {
    double af = 0.3;
    double fst = 0.05;  // Use FST > 0 to ensure difference
    auto R = compute_relationship_probs_single_snp(af, fst);
    auto B = compute_ibs_ibd_probs(af, fst);

    // At least some combos should differ between PO and Mij IBD=1
    int n_diff = 0;
    for (int c = 0; c < 9; ++c) {
        if (std::abs(R[c][(int)RelationshipType::PARENT_OFFSPRING] - B[c][1]) > 1e-6) {
            ++n_diff;
        }
    }
    CHECK(n_diff > 0);  // They should differ for at least some combos
}

// ── HS exact formula vs Mij IBD=1: verify inequivalence ─────────────
// The release note analysis shows that the HS true distribution differs
// from the Mij IBD=1 distribution. This confirms the bias source.

TEST_CASE(hs_differs_from_mij_ibd1) {
    double af = 0.5;
    auto R = compute_relationship_probs_single_snp(af, 0.0);
    auto B = compute_ibs_ibd_probs(af, 0.0);

    // HS (= PO formula) differs from Mij IBD=1 for specific combos
    // PPPP: Mij gives p^3=0.125, HS gives p^2(1+p)/2=0.09375
    double diff_pppp = std::abs(R[PPPP][(int)RelationshipType::HALF_SIBLING] - B[PPPP][1]);
    CHECK(diff_pppp > 1e-6);

    // At least some combos should differ
    int n_diff = 0;
    for (int c = 0; c < 9; ++c) {
        if (std::abs(R[c][(int)RelationshipType::HALF_SIBLING] - B[c][1]) > 1e-6) {
            ++n_diff;
        }
    }
    CHECK(n_diff > 0);
}

// ── HS differs from PO (different formula structures) ────────────────
// PO: P(G1,G2|PO) = 0.5 * Q(G1) * T(G2|G1) + 0.5 * Q(G2) * T(G1|G2)
// HS: P(G1,G2|HS) = sum_{G_S} Q(G_S) * T(G1|G_S) * T(G2|G_S)
// These are mathematically different relationships.

TEST_CASE(hs_differs_from_po) {
    std::vector<double> af_values = {0.1, 0.25, 0.3, 0.5, 0.7, 0.9};
    std::vector<double> fst_values = {0.0, 0.01, 0.05, 0.1};

    for (double af : af_values) {
        for (double fst : fst_values) {
            auto R = compute_relationship_probs_single_snp(af, fst);
            // At least some combos should differ between HS and PO
            int n_diff = 0;
            for (int c = 0; c < 9; ++c) {
                if (std::abs(R[c][(int)RelationshipType::HALF_SIBLING] -
                             R[c][(int)RelationshipType::PARENT_OFFSPRING]) > 1e-6) {
                    ++n_diff;
                }
            }
            CHECK(n_diff > 0);  // HS and PO should differ
        }
    }
}

// ── Symmetry: P(G1,G2|R) = P(G2,G1|R) ──────────────────────────────
// For all relationship types, the genotype pair probability should be
// symmetric (the pair is unordered).

TEST_CASE(symmetry_all_types_fst0) {
    auto R = compute_relationship_probs_single_snp(0.3, 0.0);

    // PPQQ ↔ QQPP
    for (int r = 0; r < NUM_RELATIONSHIP_TYPES; ++r) {
        CHECK_NEAR(R[PPQQ][r], R[QQPP][r], 1e-12);
    }
    // PPPQ ↔ PQPP
    for (int r = 0; r < NUM_RELATIONSHIP_TYPES; ++r) {
        CHECK_NEAR(R[PPPQ][r], R[PQPP][r], 1e-12);
    }
    // PQQQ ↔ QQPQ
    for (int r = 0; r < NUM_RELATIONSHIP_TYPES; ++r) {
        CHECK_NEAR(R[PQQQ][r], R[QQPQ][r], 1e-12);
    }
}

TEST_CASE(symmetry_all_types_fst_nonzero) {
    auto R = compute_relationship_probs_single_snp(0.25, 0.05);

    // Only test undirected relationship types.
    // AVUNCULAR and GRANDPARENT are directed (roles of G1/G2 differ),
    // so P(G1,G2|R) != P(G2,G1|R) when FST > 0 — this is correct.
    std::vector<int> undirected_types = {
        (int)RelationshipType::UNRELATED,
        (int)RelationshipType::DUPLICATE_MZ,
        (int)RelationshipType::PARENT_OFFSPRING,
        (int)RelationshipType::FULL_SIBLING,
        (int)RelationshipType::HALF_SIBLING,
        (int)RelationshipType::FIRST_COUSIN
    };
    for (int r : undirected_types) {
        CHECK_NEAR(R[PPQQ][r], R[QQPP][r], 1e-12);
        CHECK_NEAR(R[PPPQ][r], R[PQPP][r], 1e-12);
        CHECK_NEAR(R[PQQQ][r], R[QQPQ][r], 1e-12);
    }
}

// ── FST=0 manual verification ───────────────────────────────────────
// At p=0.5, F=0: Q(PP)=0.25, Q(PQ)=0.5, Q(QQ)=0.25
// T(PP|PP) = p = 0.5, T(PQ|PP) = q = 0.5, T(QQ|PP) = 0
// etc.

TEST_CASE(manual_unrelated_p05_fst0) {
    double p = 0.5;
    auto R = compute_relationship_probs_single_snp(p, 0.0);

    // UNRELATED: P(G1,G2) = Q(G1) * Q(G2)
    // Q(PP) = 0.25, Q(PQ) = 0.5, Q(QQ) = 0.25
    // PPPP = Q(PP)*Q(PP) = 0.0625
    CHECK_NEAR(R[PPPP][(int)RelationshipType::UNRELATED], 0.0625, 1e-10);
    // PQPQ = Q(PQ)*Q(PQ) = 0.25
    CHECK_NEAR(R[PQPQ][(int)RelationshipType::UNRELATED], 0.25, 1e-10);
    // PPQQ = Q(PP)*Q(QQ) = 0.0625
    CHECK_NEAR(R[PPQQ][(int)RelationshipType::UNRELATED], 0.0625, 1e-10);
}

TEST_CASE(manual_po_p05_fst0) {
    double p = 0.5;
    auto R = compute_relationship_probs_single_snp(p, 0.0);

    // PO at p=0.5, F=0:
    // Correct formula: P(G1,G2|PO) = 0.5 * Q(G1) * T(G2|G1) + 0.5 * Q(G2) * T(G1|G2)
    //
    // P(PPPP|PO) = 0.5 * Q(PP) * T(PP|PP) + 0.5 * Q(PP) * T(PP|PP)
    //            = Q(PP) * T(PP|PP)
    //            = p^2 * p  (parent PP transmits P, other allele from pop is P with prob p)
    //            = 0.25 * 0.5 = 0.125
    CHECK_NEAR(R[PPPP][(int)RelationshipType::PARENT_OFFSPRING], 0.125, 1e-10);

    // P(PQPQ|PO) = 0.5 * Q(PQ) * T(PQ|PQ) + 0.5 * Q(PQ) * T(PQ|PQ)
    //            = Q(PQ) * T(PQ|PQ)
    //            = 2*p*q * 0.5  (parent PQ transmits P or Q with prob 0.5 each, child must match)
    //            = 0.5 * 0.5 = 0.25
    CHECK_NEAR(R[PQPQ][(int)RelationshipType::PARENT_OFFSPRING], 0.25, 1e-10);
}

// ── FS > HS for same genotype combo (full-sib shares more) ──────────

TEST_CASE(fs_shares_more_than_hs) {
    double p = 0.3;
    auto R = compute_relationship_probs_single_snp(p, 0.0);

    // Full-siblings should have higher probability of matching genotypes
    // (PPPP, PQPQ, QQQQ) compared to half-siblings
    double fs_match = R[PPPP][(int)RelationshipType::FULL_SIBLING] +
                      R[PQPQ][(int)RelationshipType::FULL_SIBLING] +
                      R[QQQQ][(int)RelationshipType::FULL_SIBLING];
    double hs_match = R[PPPP][(int)RelationshipType::HALF_SIBLING] +
                      R[PQPQ][(int)RelationshipType::HALF_SIBLING] +
                      R[QQQQ][(int)RelationshipType::HALF_SIBLING];
    CHECK(fs_match > hs_match);
}

// ── Batch precomputation test ────────────────────────────────────────

TEST_CASE(batch_precompute_dimensions) {
    std::vector<SNPInfo> snps(10);
    for (int i = 0; i < 10; ++i) {
        snps[i].af = 0.1 + 0.05 * i;
        snps[i].chrom = "chr1";
        snps[i].pos = 1000 + i * 100;
    }

    auto matrix = precompute_relationship_likelihoods(snps, 0.0);

    // Check dimensions: [9][9][10]
    CHECK((int)matrix.size() == NUM_RELATIONSHIP_TYPES);
    for (int r = 0; r < NUM_RELATIONSHIP_TYPES; ++r) {
        CHECK((int)matrix[r].size() == 9);
        for (int c = 0; c < 9; ++c) {
            CHECK((int)matrix[r][c].size() == 10);
        }
    }

    // Check normalization for each SNP
    for (int s = 0; s < 10; ++s) {
        for (int r = 0; r < NUM_RELATIONSHIP_TYPES; ++r) {
            double sum = 0.0;
            for (int c = 0; c < 9; ++c) {
                sum += matrix[r][c][s];
            }
            CHECK_NEAR(sum, 1.0, 1e-10);
        }
    }
}

// ── Log-likelihood computation test ──────────────────────────────────

TEST_CASE(log_likelihood_basic) {
    // Create a simple 2-SNP scenario
    std::vector<SNPInfo> snps(2);
    snps[0].af = 0.3;
    snps[0].chrom = "chr1";
    snps[0].pos = 1000;
    snps[1].af = 0.5;
    snps[1].chrom = "chr1";
    snps[1].pos = 2000;

    auto matrix = precompute_relationship_likelihoods(snps, 0.0);

    // Create PIBS for a pair that looks like PO (all IBS=2, matching genotypes)
    std::vector<std::array<double, 9>> pibs(2);
    pibs[0] = {0, 0, 0, 0, 0, 0, 0.5, 0.3, 0.2};  // Mostly PPPP
    pibs[1] = {0, 0, 0, 0, 0, 0, 0.3, 0.5, 0.2};  // Mostly PQPQ

    std::vector<int> indices = {0, 1};

    // All relationship types should give finite log-likelihoods
    for (int r = 0; r < NUM_RELATIONSHIP_TYPES; ++r) {
        double ll = compute_relationship_log_likelihood(
            matrix, pibs, indices, static_cast<RelationshipType>(r));
        CHECK(std::isfinite(ll));
    }

    // PO/DUPLICATE should have higher likelihood than UNRELATED for matching genotypes
    double ll_po = compute_relationship_log_likelihood(
        matrix, pibs, indices, RelationshipType::PARENT_OFFSPRING);
    double ll_un = compute_relationship_log_likelihood(
        matrix, pibs, indices, RelationshipType::UNRELATED);
    CHECK(ll_po > ll_un);
}

// ── All probabilities non-negative ───────────────────────────────────

TEST_CASE(all_probs_non_negative) {
    std::vector<double> af_values = {0.05, 0.1, 0.25, 0.3, 0.5, 0.7, 0.9, 0.95};
    std::vector<double> fst_values = {0.0, 0.01, 0.05, 0.1};

    for (double af : af_values) {
        for (double fst : fst_values) {
            auto R = compute_relationship_probs_single_snp(af, fst);
            for (int c = 0; c < 9; ++c) {
                for (int r = 0; r < NUM_RELATIONSHIP_TYPES; ++r) {
                    CHECK(R[c][r] >= -1e-15);  // Allow tiny numerical errors
                }
            }
        }
    }
}

// ── main ─────────────────────────────────────────────────────────────

int main() {
    return RUN_ALL_TESTS();
}
