/**
 * @file ibd_model.cpp
 * @brief IBS|IBD probability model using Mij sequential-draw framework
 *
 * Computes P(IBS=c | IBD=z) for 9 genotype combos x 3 IBD states.
 *
 * Uses the Mij sequential-draw model from Anderson & Weir (2007):
 *   M(p, F, i) = (1-F)*p + i*F
 *
 * Each IBD state's numerators are computed as products of Mij terms,
 * then normalized by dividing by their sum, ensuring exact normalization:
 *   sum_c B[c][z] = 1  for all z and all allele frequencies.
 *
 * IBD=0: Four independent sequential allele draws (two unrelated individuals)
 * IBD=1: One shared allele + two non-shared (three sequential draws)
 * IBD=2: Both alleles shared (genotypes must match; two sequential draws)
 *
 * @author Shujia Huang
 * @date 2026-06-23
 */

#include "ibd_model.h"

namespace fastlckin {

std::array<std::array<double, 3>, 9> compute_ibs_ibd_probs(double af, double fst) {
    std::array<std::array<double, 3>, 9> B{};

    double p = af;         // ALT allele frequency (P = gl[2])
    double q = 1.0 - p;    // REF allele frequency (Q = gl[0])
    double F = fst;

    // Mij shorthand: M(freq, FST, i) = (1-FST)*freq + i*FST
    auto M = [&](double freq, int i) -> double { return Mij(freq, F, i); };

    // ── IBS=0 combos (IBD=1 and IBD=2 are impossible) ────────
    // Four sequential allele draws: two from ind1, two from ind2

    // PPQQ: G1=PP, G2=QQ
    B[PPQQ][0] = M(p, 1) * M(p, 0) * M(q, 1) * M(q, 0);
    B[PPQQ][1] = 0.0;
    B[PPQQ][2] = 0.0;

    // QQPP: G1=QQ, G2=PP (symmetric)
    B[QQPP][0] = M(q, 1) * M(q, 0) * M(p, 1) * M(p, 0);
    B[QQPP][1] = 0.0;
    B[QQPP][2] = 0.0;

    // ── IBS=1 combos ──────────────────────────────────────────
    // IBD=0: Four sequential draws (standard Mij model)
    // IBD=1: One shared allele S + two independent unshared alleles.
    //         Genotypes are UNORDERED pairs {allele1, allele2}.
    //         P(G1,G2|IBD=1) = sum over valid (S,U1,U2) assignments
    //         Key: het genotypes use M(freq,0)=(1-F)*freq (no coancestry
    //         correction for the first draw from population), while hom
    //         genotypes use M(freq,1)=(1-F)*freq+F.
    // IBD=2: impossible (genotypes must match for IBS=1)

    double p1 = M(p, 1);   // (1-F)*p + F  (hom: one match already counted)
    double q1 = M(q, 1);   // (1-F)*q + F
    double p0 = M(p, 0);   // (1-F)*p      (het: first draw, no match)
    double q0 = M(q, 0);   // (1-F)*q

    // PPPQ: G1=PP, G2=PQ → S=P, U1=P, U2=Q
    B[PPPQ][0] = 2.0 * M(p, 2) * M(p, 1) * M(p, 0) * M(q, 0);
    B[PPPQ][1] = p * p1 * q0;
    B[PPPQ][2] = 0.0;

    // PQPP: G1=PQ, G2=PP → S=P, U1=Q, U2=P
    B[PQPP][0] = 2.0 * M(p, 2) * M(p, 1) * M(p, 0) * M(q, 0);
    B[PQPP][1] = p * q0 * p1;
    B[PQPP][2] = 0.0;

    // PQQQ: G1=PQ, G2=QQ → S=Q, U1=P, U2=Q
    B[PQQQ][0] = 2.0 * M(q, 2) * M(q, 1) * M(q, 0) * M(p, 0);
    B[PQQQ][1] = q * p0 * q1;
    B[PQQQ][2] = 0.0;

    // QQPQ: G1=QQ, G2=PQ → S=Q, U1=Q, U2=P
    B[QQPQ][0] = 2.0 * M(q, 2) * M(q, 1) * M(q, 0) * M(p, 0);
    B[QQPQ][1] = q * q1 * p0;
    B[QQPQ][2] = 0.0;

    // ── IBS=2 combos ──────────────────────────────────────────
    // IBD=0: Four sequential draws
    // IBD=1: One shared + two independent (PPPP/QQQQ possible when F>0)
    // IBD=2: Both alleles shared (two sequential draws)

    // PPPP: G1=PP, G2=PP → S=P, U1=P, U2=P
    B[PPPP][0] = M(p, 3) * M(p, 2) * M(p, 1) * M(p, 0);
    B[PPPP][1] = p * p1 * p1;
    B[PPPP][2] = M(p, 1) * M(p, 0);

    // PQPQ: G1=PQ, G2=PQ → S=P,U1=Q,U2=Q  OR  S=Q,U1=P,U2=P
    B[PQPQ][0] = 4.0 * M(p, 1) * M(p, 0) * M(q, 1) * M(q, 0);
    B[PQPQ][1] = p * q0 * q0 + q * p0 * p0;
    B[PQPQ][2] = 2.0 * M(p, 0) * M(q, 1);

    // QQQQ: G1=QQ, G2=QQ → S=Q, U1=Q, U2=Q
    B[QQQQ][0] = M(q, 3) * M(q, 2) * M(q, 1) * M(q, 0);
    B[QQQQ][1] = q * q1 * q1;
    B[QQQQ][2] = M(q, 1) * M(q, 0);

    // ── Dynamic normalization ─────────────────────────────────
    // Sum numerators for each IBD state, then divide.
    // This guarantees sum_c B[c][z] = 1 for each z, regardless of
    // allele frequency or FST value.
    for (int z = 0; z < 3; ++z) {
        double sum = 0.0;
        for (int c = 0; c < 9; ++c) {
            sum += B[c][z];
        }
        if (sum > 1e-15) {
            for (int c = 0; c < 9; ++c) {
                B[c][z] /= sum;
            }
        }
    }

    return B;
}

IBS_IBD_Matrix precompute_ibs_ibd_matrix(const std::vector<SNPInfo>& snp_infos, double fst) {
    size_t n_snps = snp_infos.size();

    // Allocate: [9][n_snps][3]
    IBS_IBD_Matrix matrix(9, std::vector<std::array<double, 3>>(n_snps, {0.0, 0.0, 0.0}));

    for (size_t s = 0; s < n_snps; ++s) {
        auto B = compute_ibs_ibd_probs(snp_infos[s].af, fst);
        for (int c = 0; c < 9; ++c) {
            matrix[c][s] = B[c];
        }
    }

    return matrix;
}

IBS_IBD_Matrix precompute_ibs_ibd_matrix_fst_vector(
    const std::vector<SNPInfo>& snp_infos,
    const std::vector<double>& fst_vector
) {
    size_t n_snps = snp_infos.size();
    IBS_IBD_Matrix matrix(9, std::vector<std::array<double, 3>>(n_snps, {0.0, 0.0, 0.0}));

    for (size_t s = 0; s < n_snps; ++s) {
        // Use per-SNP FST if available, otherwise default to 0
        double fst = (s < fst_vector.size()) ? fst_vector[s] : 0.0;
        auto B = compute_ibs_ibd_probs(snp_infos[s].af, fst);
        for (int c = 0; c < 9; ++c) {
            matrix[c][s] = B[c];
        }
    }

    return matrix;
}

}  // namespace fastlckin

