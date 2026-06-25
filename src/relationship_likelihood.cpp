/**
 * @file relationship_likelihood.cpp
 * @brief Exact per-relationship genotype joint probability P(G1, G2 | R)
 *
 * Implements Mendelian transmission-based likelihood computation for 8 standard
 * relationship types. Each type's probability is derived by summing over shared
 * ancestor genotypes, replacing the Mij model's population-averaged IBD=1
 * approximation with exact pedigree-based formulas.
 *
 * Key insight: The Mij model's IBD=1 distribution is exact for parent-offspring
 * but systematically biased for half-siblings, avuncular, and grandparent
 * relationships. This module computes the exact formulas for each type.
 *
 * @author Shujia Huang
 * @date 2026-06-24
 */

#include "relationship_likelihood.h"
#include <cmath>

namespace fastlckin {

// ────────────────────────────────────────────────────────────────────
// Helper: Population genotype probability Q(G | p, F)
// ────────────────────────────────────────────────────────────────────
// Q(PP) = p * p_1,  Q(PQ) = 2*p*q_0,  Q(QQ) = q * q_1
// where p_0 = (1-F)*p, q_0 = (1-F)*q, p_1 = (1-F)*p + F, q_1 = (1-F)*q + F

inline double Q_func(int g, double p, double q, double p0, double q0,
                     double p1, double q1) {
    switch (g) {
        case 0: return p * p1;       // PP
        case 1: return 2.0 * p * q0; // PQ
        case 2: return q * q1;       // QQ
        default: return 0.0;
    }
}

// ────────────────────────────────────────────────────────────────────
// Helper: Mendelian transmission probability T(G_child | G_parent)
// ────────────────────────────────────────────────────────────────────
// One parent has known genotype G_parent; the other parent is drawn from
// the population (allele frequencies p, q with FST correction F).
//
// T(G_child | G_parent) = sum_{a transmitted from parent} P(a | G_parent)
//                          * P(other allele from pop = b)
// where G_child = {a, b}

inline double T_func(int g_child, int g_parent, double p, double q,
                     double p0, double q0, double p1, double q1) {
    switch (g_parent) {
        case 0: // PP: always transmits P
            // Child = {P, X}: X from population
            switch (g_child) {
                case 0: return p1;        // PP: X=P
                case 1: return q0;        // PQ: X=Q
                case 2: return 0.0;       // QQ: impossible
            }
            break;
        case 1: // PQ: transmits P with prob 1/2, Q with prob 1/2
            switch (g_child) {
                case 0: return 0.5 * p1;  // PP: transmit P, pop gives P
                case 1: return 0.5;       // PQ: (transmit P, pop gives Q) + (transmit Q, pop gives P) = 0.5*q0 + 0.5*p0 = 0.5
                case 2: return 0.5 * q0;  // QQ: transmit Q, pop gives Q
            }
            break;
        case 2: // QQ: always transmits Q
            switch (g_child) {
                case 0: return 0.0;       // PP: impossible
                case 1: return p0;        // PQ: pop gives P
                case 2: return q1;        // QQ: pop gives Q
            }
            break;
    }
    return 0.0;
}

// ────────────────────────────────────────────────────────────────────
// Compute all 8 relationship types for a single SNP
// ────────────────────────────────────────────────────────────────────

std::array<std::array<double, NUM_RELATIONSHIP_TYPES>, 9>
compute_relationship_probs_single_snp(double af, double fst) {
    std::array<std::array<double, NUM_RELATIONSHIP_TYPES>, 9> R{};

    double p = af;
    double q = 1.0 - p;
    double F = fst;

    // FST-corrected allele frequencies
    double p0 = (1.0 - F) * p;         // het path, no coancestry
    double q0 = (1.0 - F) * q;
    double p1 = (1.0 - F) * p + F;     // hom path, with coancestry
    double q1 = (1.0 - F) * q + F;

    // Pre-compute Q(g) for g = 0,1,2
    double Qg[3];
    for (int g = 0; g < 3; ++g) Qg[g] = Q_func(g, p, q, p0, q0, p1, q1);

    // Pre-compute T(g_child | g_parent) for all g_child, g_parent in {0,1,2}
    double Tg[3][3];
    for (int gc = 0; gc < 3; ++gc)
        for (int gp = 0; gp < 3; ++gp)
            Tg[gc][gp] = T_func(gc, gp, p, q, p0, q0, p1, q1);

    // ── (a) UNRELATED: P(G1,G2) = Q(G1) * Q(G2) ───────────────────
    // Combo (g1, g2): PPQQ=(0,2), QQPP=(2,0), PPPQ=(0,1), PQPP=(1,0),
    //                 PQQQ=(1,2), QQPQ=(2,1), PPPP=(0,0), PQPQ=(1,1), QQQQ=(2,2)
    for (int g1 = 0; g1 < 3; ++g1) {
        for (int g2 = 0; g2 < 3; ++g2) {
            int combo = g1 * 3 + g2;
            // Map from sequential (0-8) to GenotypeCombo enum
            // Sequential: 0=PPPP, 1=PPPQ, 2=PPQQ, 3=PQPP, 4=PQPQ, 5=PQQQ, 6=QQPP, 7=QQPQ, 8=QQQQ
            // But we use the GenotypeCombo enum mapping:
            // PPQQ=0, QQPP=1, PPPQ=2, PQPP=3, PQQQ=4, QQPQ=5, PPPP=6, PQPQ=7, QQQQ=8
            // We'll use a direct mapping function below
        }
    }

    // Helper: genotype pair (g1, g2) → GenotypeCombo index
    // g1, g2 in {0=PP, 1=PQ, 2=QQ}
    auto combo_idx = [](int g1, int g2) -> int {
        // PPQQ=0, QQPP=1, PPPQ=2, PQPP=3, PQQQ=4, QQPQ=5, PPPP=6, PQPQ=7, QQQQ=8
        static const int table[3][3] = {
            // g2: PP(0)  PQ(1)  QQ(2)
            {6, 2, 0},   // g1=PP(0): PPPP=6, PPPQ=2, PPQQ=0
            {3, 7, 4},   // g1=PQ(1): PQPP=3, PQPQ=7, PQQQ=4
            {1, 5, 8}    // g1=QQ(2): QQPP=1, QQPQ=5, QQQQ=8
        };
        return table[g1][g2];
    };

    // (a) UNRELATED
    for (int g1 = 0; g1 < 3; ++g1)
        for (int g2 = 0; g2 < 3; ++g2)
            R[combo_idx(g1, g2)][(int)RelationshipType::UNRELATED] = Qg[g1] * Qg[g2];

    // (b) DUPLICATE/MZ: P(G1,G2) = Q(G1) * delta(G1==G2)
    for (int g1 = 0; g1 < 3; ++g1)
        for (int g2 = 0; g2 < 3; ++g2)
            R[combo_idx(g1, g2)][(int)RelationshipType::DUPLICATE_MZ] =
                (g1 == g2) ? Qg[g1] : 0.0;

    // (c) PARENT-OFFSPRING: G1 is direct parent of G2 (or vice versa)
    // P(G1,G2|PO) = 0.5 * Q(G1) * T(G2|G1) + 0.5 * Q(G2) * T(G1|G2)
    // Since we don't know who is parent and who is child, we average both directions.
    for (int g1 = 0; g1 < 3; ++g1) {
        for (int g2 = 0; g2 < 3; ++g2) {
            double prob = 0.5 * Qg[g1] * Tg[g2][g1] + 0.5 * Qg[g2] * Tg[g1][g2];
            R[combo_idx(g1, g2)][(int)RelationshipType::PARENT_OFFSPRING] = prob;
        }
    }

    // (d) HALF-SIBLING: G1 and G2 share one common ancestor S
    // P(G1,G2|HS) = sum_{G_S} Q(G_S) * T(G1|G_S) * T(G2|G_S)
    for (int g1 = 0; g1 < 3; ++g1) {
        for (int g2 = 0; g2 < 3; ++g2) {
            double prob = 0.0;
            for (int gs = 0; gs < 3; ++gs) {
                prob += Qg[gs] * Tg[g1][gs] * Tg[g2][gs];
            }
            R[combo_idx(g1, g2)][(int)RelationshipType::HALF_SIBLING] = prob;
        }
    }

    // (e) AVUNCULAR: G1 = uncle/aunt, G2 = nephew/niece
    // Shared grandparent (GP), intermediate parent (P = sibling of G1)
    // P(G1,G2|AV) = sum_{G_GP} Q(G_GP) * [sum_{G_P} T(G_P|G_GP)*Q(G_P)] * T(G1|G_P) * T(G2|G_GP)
    //
    // Here: G1 is the uncle/aunt (child of GP + GP_spouse)
    //       G2 is the nephew/niece (child of P + pop_spouse)
    //       P is the parent of G2 and sibling of G1 (child of GP + pop_spouse)
    //
    // Actually, let me re-derive:
    // G1 and P are full siblings (share both parents GP1, GP2)
    // G2 is child of P
    // But for avuncular, G1 and P are half-siblings (share one parent GP)
    // Wait, avuncular means G1 is the uncle/aunt of G2.
    // G2's parent (call them P) is G1's nephew/niece's parent.
    // G1 and P are siblings. If they are full siblings, then G1 and G2 are
    // full-avuncular. If half-siblings, half-avuncular.
    //
    // For standard "avuncular" (full-avuncular):
    // G1 and P are full siblings (share both parents GP1, GP2)
    // G2 is child of P + pop_spouse
    //
    // But the release note analysis treats avuncular as k=(0.5, 0.5, 0),
    // which is the same as half-sibling. This suggests the "standard"
    // avuncular in genetics is actually half-avuncular (sharing one ancestor).
    //
    // Let me use the half-avuncular model (one shared ancestor GP):
    // G1 is child of GP + pop
    // P is child of GP + pop (half-sibling of G1)
    // G2 is child of P + pop
    //
    // P(G1,G2|AV) = sum_{G_GP} Q(G_GP) *
    //   [sum_{G_P} T(G_P|G_GP) * Q(G_P)] * T(G1|G_GP) * T(G2|G_P)
    //
    // Wait, I need to be more careful. Let me re-derive:
    // GP is the shared ancestor.
    // G1's lineage: GP → G1 (one meiosis, G1 also gets allele from pop)
    // G2's lineage: GP → P → G2 (two meioses)
    //   P is child of GP + pop_spouse
    //   G2 is child of P + pop_spouse
    //
    // P(G1, G2 | AV) = sum_{G_GP} Q(G_GP) *
    //   T(G1 | G_GP) *   [G1 is child of GP]
    //   sum_{G_P} T(G_P | G_GP) * Q(G_P) *   [P is child of GP + pop]
    //   T(G2 | G_P)   [G2 is child of P + pop]
    //
    // Hmm, but this doesn't account for G1's other parent. Let me think again.
    //
    // G1 is the uncle/aunt. G1's parent is GP. G1's other parent is from pop.
    // So G1 is a child of GP: P(G1 | G_GP) = T(G1 | G_GP)
    //
    // P (the parent of G2) is also a child of GP. P's other parent is from pop.
    // So P(G_P | G_GP) = T(G_P | G_GP)
    //
    // G2 is a child of P. G2's other parent is from pop.
    // So P(G2 | G_P) = T(G2 | G_P)
    //
    // G1 and P are half-siblings (share GP).
    // Given GP, G1 and P are conditionally independent.
    //
    // P(G1, G2 | AV) = sum_{G_GP} Q(G_GP) * T(G1 | G_GP) *
    //   [sum_{G_P} T(G_P | G_GP) * T(G2 | G_P)]
    //
    // This is the half-avuncular model.

    for (int g1 = 0; g1 < 3; ++g1) {
        for (int g2 = 0; g2 < 3; ++g2) {
            double prob = 0.0;
            for (int ggp = 0; ggp < 3; ++ggp) {
                // Inner sum: sum_{G_P} T(G_P | G_GP) * T(G2 | G_P)
                double inner = 0.0;
                for (int gp = 0; gp < 3; ++gp) {
                    inner += Tg[gp][ggp] * Tg[g2][gp];
                }
                prob += Qg[ggp] * Tg[g1][ggp] * inner;
            }
            R[combo_idx(g1, g2)][(int)RelationshipType::AVUNCULAR] = prob;
        }
    }

    // (f) GRANDPARENT: G1 = grandparent, G2 = grandchild
    // G1 is the grandparent.
    // P (intermediate) is child of G1 + pop: P(G_P | G1) = T(G_P | G1)
    // G2 is child of P + pop: P(G2 | G_P) = T(G2 | G_P)
    //
    // P(G1, G2 | GP) = sum_{G1} Q(G1) * T(G1 | G1) *
    //   [sum_{G_P} T(G_P | G1) * T(G2 | G_P)]
    //
    // Wait, that's circular. Let me re-derive:
    //
    // P(G1, G2 | GP_rel) = sum_{G_GP} Q(G_GP) *
    //   T(G1 | G_GP) *   [G1 is child of GP + pop]
    //   sum_{G_P} T(G_P | G_GP) *   [P is child of GP + pop]
    //   T(G2 | G_P)   [G2 is child of P + pop]
    //
    // This is the same formula structure as avuncular! But the roles are different:
    // - In avuncular: G1 is the uncle/aunt (1 meiosis from GP), G2 is the
    //   nephew/niece (2 meioses from GP via P)
    // - In grandparent: G1 is the grandparent (0 meioses from GP = GP themselves),
    //   G2 is the grandchild (2 meioses from GP via P)
    //
    // For grandparent, G1 IS the grandparent, so there's no T(G1|G_GP):
    //
    // P(G1, G2 | GP_rel) = Q(G1) * sum_{G_P} T(G_P | G1) * T(G2 | G_P)
    //
    // where Q(G1) is the population frequency of G1,
    // T(G_P | G1) is the probability of intermediate parent given G1,
    // T(G2 | G_P) is the probability of grandchild given intermediate parent.

    for (int g1 = 0; g1 < 3; ++g1) {
        for (int g2 = 0; g2 < 3; ++g2) {
            double prob = 0.0;
            for (int gp = 0; gp < 3; ++gp) {
                prob += Tg[gp][g1] * Tg[g2][gp];
            }
            R[combo_idx(g1, g2)][(int)RelationshipType::GRANDPARENT] = Qg[g1] * prob;
        }
    }

    // ── Helper: probability of transmitting allele A from parent with genotype G
    // A: 0=P, 1=Q
    auto transmit = [](int allele, int g_parent) -> double {
        switch (g_parent) {
            case 0: return (allele == 0) ? 1.0 : 0.0;  // PP
            case 1: return 0.5;                          // PQ
            case 2: return (allele == 0) ? 0.0 : 1.0;  // QQ
            default: return 0.0;
        }
    };

    // ── Helper: P(G_child | G_father, G_mother) - both parents known
    // For unordered genotypes, we must sum over both allele assignments:
    //   P(PQ | gf, gm) = P(P←fa,Q←ma) + P(Q←fa,P←ma)
    // This symmetrizes over father/mother allele roles.
    auto T2_func = [&](int g_child, int g_father, int g_mother) -> double {
        if (g_child == 0) {
            // PP: need P from father AND P from mother
            return transmit(0, g_father) * transmit(0, g_mother);
        } else if (g_child == 2) {
            // QQ: need Q from father AND Q from mother
            return transmit(1, g_father) * transmit(1, g_mother);
        } else {
            // PQ (het): P from father + Q from mother, OR Q from father + P from mother
            return transmit(0, g_father) * transmit(1, g_mother)
                 + transmit(1, g_father) * transmit(0, g_mother);
        }
    };

    // (g) FULL-SIBLING: G1 and G2 share both parents (Father F, Mother M)
    //
    // P(G1, G2 | FS) = sum_{G_F, G_M} Q(G_F) * Q(G_M) *
    //                    T2(G1 | G_F, G_M) * T2(G2 | G_F, G_M)
    //
    // where T2(g | g_f, g_m) = P(child genotype g | father=g_f, mother=g_m)
    // = sum over (allele from father, allele from mother) that produce g.
    //
    // CRITICAL: Cannot use T(g|gf)*T(g|gm) because that gives BOTH alleles
    // from the same parent. T2 correctly assigns one allele from each parent.

    for (int g1 = 0; g1 < 3; ++g1) {
        for (int g2 = 0; g2 < 3; ++g2) {
            double prob = 0.0;
            for (int gf = 0; gf < 3; ++gf) {
                for (int gm = 0; gm < 3; ++gm) {
                    prob += Qg[gf] * Qg[gm] * T2_func(g1, gf, gm) * T2_func(g2, gf, gm);
                }
            }
            R[combo_idx(g1, g2)][(int)RelationshipType::FULL_SIBLING] = prob;
        }
    }

    // (h) FIRST-COUSIN: G1 and G2 are first cousins
    //
    // Pedigree:
    //   Shared grandparents: GP1 + GP2
    //   P1 = child of GP1+GP2 (parent of G1)
    //   P2 = child of GP1+GP2 (parent of G2)
    //   P1 and P2 are full siblings
    //   G1 = child of P1 + spouse1
    //   G2 = child of P2 + spouse2
    //
    // P(G1, G2 | FC) = sum_{G_GP1, G_GP2} Q(G_GP1) * Q(G_GP2) *
    //   [sum_{G_P1, G_P2} T(G_P1|G_GP1)*T(G_P1|G_GP2) *
    //                       T(G_P2|G_GP1)*T(G_P2|G_GP2) *
    //                       T(G1|G_P1)*Q(G1_spouse) *
    //                       T(G2|G_P2)*Q(G2_spouse)]
    //
    // where Q(G1_spouse) and Q(G2_spouse) are population frequencies of the
    // spouses' genotypes (who marry into the family). But we need to sum over
    // the spouses' genotypes too... actually no. The spouses contribute to
    // G1 and G2 through T(G1|G_P1) which already accounts for the other parent
    // being from the population.
    //
    // Wait, T(G1|G_P1) = P(G1 | one parent = G_P1, other from pop).
    // This already includes the spouse's contribution. So:
    //
    // P(G1, G2 | FC) = sum_{G_GP1, G_GP2} Q(G_GP1) * Q(G_GP2) *
    //   [sum_{G_P1} T(G_P1|G_GP1)*T(G_P1|G_GP2) * T(G1|G_P1)] *
    //   [sum_{G_P2} T(G_P2|G_GP1)*T(G_P2|G_GP2) * T(G2|G_P2)]
    //
    // where T(G_P1|G_GP1)*T(G_P1|G_GP2) is NOT right. P1 gets one allele from
    // GP1 and one from GP2. So:
    //
    // P(G_P1 | G_GP1, G_GP2) = T(G_P1 | G_GP1, G_GP2)
    //   = sum over which allele comes from which parent
    //
    // Actually, T(G_child | G_parent) already handles this: one parent is known,
    // the other is from population. But here BOTH parents are known (GP1 and GP2).
    //
    // I need a different transmission function for when both parents are known:
    // T2(G_child | G_father, G_mother) = P(G_child | father=G_f, mother=G_m)
    //
    // For unordered genotype G_child = {a, b}:
    // P(G_child | G_f, G_m) = sum over (allele from father, allele from mother)
    //   that produces G_child
    //
    // T2(PP | G_f, G_m) = P(P from father) * P(P from mother)
    // T2(PQ | G_f, G_m) = P(P from father)*P(Q from mother) + P(Q from father)*P(P from mother)
    // T2(QQ | G_f, G_m) = P(Q from father) * P(Q from mother)
    //
    // P(A from parent with genotype G):
    // G=PP: P(A=P)=1, P(A=Q)=0
    // G=PQ: P(A=P)=1/2, P(A=Q)=1/2
    // G=QQ: P(A=P)=0, P(A=Q)=1

    // First-cousin computation:
    // P(G1, G2 | FC) = sum_{G_GP1, G_GP2} Q(G_GP1) * Q(G_GP2) *
    //   [sum_{G_P1} T2(G_P1 | G_GP1, G_GP2) * T(G1 | G_P1)] *
    //   [sum_{G_P2} T2(G_P2 | G_GP1, G_GP2) * T(G2 | G_P2)]

    for (int g1 = 0; g1 < 3; ++g1) {
        for (int g2 = 0; g2 < 3; ++g2) {
            double prob = 0.0;
            for (int ggp1 = 0; ggp1 < 3; ++ggp1) {
                for (int ggp2 = 0; ggp2 < 3; ++ggp2) {
                    double gp_prob = Qg[ggp1] * Qg[ggp2];

                    // Sum over P1 genotypes
                    double sum_p1 = 0.0;
                    for (int gp1 = 0; gp1 < 3; ++gp1) {
                        sum_p1 += T2_func(gp1, ggp1, ggp2) * Tg[g1][gp1];
                    }

                    // Sum over P2 genotypes
                    double sum_p2 = 0.0;
                    for (int gp2 = 0; gp2 < 3; ++gp2) {
                        sum_p2 += T2_func(gp2, ggp1, ggp2) * Tg[g2][gp2];
                    }

                    prob += gp_prob * sum_p1 * sum_p2;
                }
            }
            R[combo_idx(g1, g2)][(int)RelationshipType::FIRST_COUSIN] = prob;
        }
    }

    // (i) GREAT-GRANDPARENT: G1 = great-grandparent, G2 = great-grandchild
    //
    // Linear chain of 3 transmissions:
    //   G1 → G_A (intermediate parent 1) → G_B (intermediate parent 2) → G2
    //
    // P(G1, G2 | GGP) = Q(G1) * sum_{G_A} T(G_A | G1) *
    //                    sum_{G_B} T(G_B | G_A) * T(G2 | G_B)
    //
    // Note: This is a directed relationship (G1 = ancestor, G2 = descendant).

    for (int g1 = 0; g1 < 3; ++g1) {
        for (int g2 = 0; g2 < 3; ++g2) {
            double prob = 0.0;
            for (int ga = 0; ga < 3; ++ga) {
                double t_ga = Tg[ga][g1];  // T(G_A | G1)
                if (t_ga == 0.0) continue;
                double inner = 0.0;
                for (int gb = 0; gb < 3; ++gb) {
                    inner += Tg[gb][ga] * Tg[g2][gb];  // T(G_B | G_A) * T(G2 | G_B)
                }
                prob += t_ga * inner;
            }
            R[combo_idx(g1, g2)][(int)RelationshipType::GREAT_GRANDPARENT] = Qg[g1] * prob;
        }
    }

    return R;
}

// ────────────────────────────────────────────────────────────────────
// Batch pre-compute for all SNPs
// ────────────────────────────────────────────────────────────────────

RelationshipLikelihoodMatrix precompute_relationship_likelihoods(
    const std::vector<SNPInfo>& snp_infos, double fst
) {
    size_t n_snps = snp_infos.size();

    // Allocate: [NUM_RELATIONSHIP_TYPES][9][n_snps]
    RelationshipLikelihoodMatrix matrix(
        NUM_RELATIONSHIP_TYPES,
        std::vector<std::vector<double>>(
            9, std::vector<double>(n_snps, 0.0)
        )
    );

    for (size_t s = 0; s < n_snps; ++s) {
        auto probs = compute_relationship_probs_single_snp(snp_infos[s].af, fst);
        for (int c = 0; c < 9; ++c) {
            for (int r = 0; r < NUM_RELATIONSHIP_TYPES; ++r) {
                matrix[r][c][s] = probs[c][r];
            }
        }
    }

    return matrix;
}

// ────────────────────────────────────────────────────────────────────
// Compute log-likelihood for a pair under a relationship type
// ────────────────────────────────────────────────────────────────────

double compute_relationship_log_likelihood(
    const RelationshipLikelihoodMatrix& rel_matrix,
    const std::vector<std::array<double, 9>>& pibs_per_snp,
    const std::vector<int>& snp_indices,
    RelationshipType type
) {
    int r = static_cast<int>(type);
    double log_lik = 0.0;

    for (size_t ri = 0; ri < snp_indices.size(); ++ri) {
        int s = snp_indices[ri];
        double site_prob = 0.0;
        for (int c = 0; c < 9; ++c) {
            site_prob += pibs_per_snp[ri][c] * rel_matrix[r][c][s];
        }
        if (site_prob <= 0.0 || std::isinf(site_prob)) {
            return -1e10;  // Effectively zero likelihood
        }
        log_lik += std::log(site_prob);
    }

    return log_lik;
}

}  // namespace fastlckin
