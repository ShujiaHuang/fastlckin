/**
 * @file ibd_model.cpp
 * @brief Anderson & Weir (2007) IBS|IBD probability model implementation
 * @author Shujia Huang
 * @date 2025-06-23
 */

#include "ibd_model.h"

namespace fastlckin {

std::array<std::array<double, 3>, 9> compute_ibs_ibd_probs(double af, double fst) {
    std::array<std::array<double, 3>, 9> B{};

    double p = af;         // alt allele frequency
    double q = 1.0 - p;    // ref allele frequency
    double F = fst;

    // Denominators
    double D0 = (1.0 + 2.0 * F) * (1.0 + F) * (1.0 - F);
    double D1 = (1.0 + F) * (1.0 - F);
    double D2 = (1.0 - F);

    // Guard against zero denominators when FST is close to 1
    if (D0 < 1e-15) D0 = 1e-15;
    if (D1 < 1e-15) D1 = 1e-15;
    if (D2 < 1e-15) D2 = 1e-15;

    // Mij shorthand: M(allele_freq, FST, count)
    auto M = [&](double freq, int i) -> double { return Mij(freq, F, i); };

    // ── IBS=0 combos ──────────────────────────────────────────

    // PPQQ: ind1=pp(REF/REF), ind2=qq(ALT/ALT)
    B[PPQQ][0] = M(p, 1) * M(p, 0) * M(q, 1) * M(q, 0) / D0;
    B[PPQQ][1] = 0.0;
    B[PPQQ][2] = 0.0;

    // QQPP: ind1=qq, ind2=pp (symmetric)
    B[QQPP][0] = M(q, 1) * M(q, 0) * M(p, 1) * M(p, 0) / D0;
    B[QQPP][1] = 0.0;
    B[QQPP][2] = 0.0;

    // ── IBS=1 combos ──────────────────────────────────────────

    // PPPQ: ind1=pp, ind2=pq
    B[PPPQ][0] = 2.0 * M(p, 2) * M(p, 1) * M(p, 0) * M(q, 0) / D0;
    B[PPPQ][1] = M(p, 1) * M(p, 0) * M(q, 0) / D1;
    B[PPPQ][2] = 0.0;

    // PQPP: ind1=pq, ind2=pp
    B[PQPP][0] = 2.0 * M(p, 2) * M(p, 1) * M(p, 0) * M(q, 0) / D0;
    B[PQPP][1] = M(p, 1) * M(p, 0) * M(q, 0) / D1;
    B[PQPP][2] = 0.0;

    // PQQQ: ind1=pq, ind2=qq
    B[PQQQ][0] = 2.0 * M(q, 2) * M(q, 1) * M(q, 0) * M(p, 0) / D0;
    B[PQQQ][1] = M(q, 1) * M(q, 0) * M(p, 0) / D1;
    B[PQQQ][2] = 0.0;

    // QQPQ: ind1=qq, ind2=pq
    B[QQPQ][0] = 2.0 * M(q, 2) * M(q, 1) * M(q, 0) * M(p, 0) / D0;
    B[QQPQ][1] = M(q, 1) * M(q, 0) * M(p, 0) / D1;
    B[QQPQ][2] = 0.0;

    // ── IBS=2 combos ──────────────────────────────────────────

    // PPPP: ind1=pp, ind2=pp
    B[PPPP][0] = M(p, 3) * M(p, 2) * M(p, 1) * M(p, 0) / D0;
    B[PPPP][1] = M(p, 2) * M(p, 1) * M(p, 0) / D1;
    B[PPPP][2] = M(p, 1) * M(p, 0) / D2;

    // PQPQ: ind1=pq, ind2=pq
    B[PQPQ][0] = 4.0 * M(p, 1) * M(p, 0) * M(q, 1) * M(q, 0) / D0;
    B[PQPQ][1] = M(p, 0) * M(q, 0) * (M(p, 1) + M(q, 1)) / D1;
    B[PQPQ][2] = 2.0 * M(p, 0) * M(q, 0) / D2;

    // QQQQ: ind1=qq, ind2=qq
    B[QQQQ][0] = M(q, 3) * M(q, 2) * M(q, 1) * M(q, 0) / D0;
    B[QQQQ][1] = M(q, 2) * M(q, 1) * M(q, 0) / D1;
    B[QQQQ][2] = M(q, 1) * M(q, 0) / D2;

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

}  // namespace fastlckin
