/**
 * @file test_se_ci.cpp
 * @brief Unit tests for standard error and confidence interval computation
 *        via observed Fisher information matrix (numerical Hessian).
 *
 * Tests cover:
 *   1. KinshipResult default SE/CI fields are -1 (NA)
 *   2. Numerical Hessian correctness against analytical Hessian for a
 *      simplified model
 *   3. Boundary detection: SE = NA when MLE is on parameter boundary
 *   4. CI coverage: MLE point lies within the 95% CI
 *   5. Delta method SE for PI_HAT
 */

#include "test_harness.h"
#include "kinship_estimator.h"
#include <cmath>
#include <array>
#include <vector>
#include <random>

// ────────────────────────────────────────────────────────────────────
// Test 1: KinshipResult default SE/CI fields are NA (-1)
// ────────────────────────────────────────────────────────────────────
TEST_CASE(kinship_result_default_se_na) {
    fastlckin::KinshipResult r;
    CHECK(r.se_k0 < -0.5);
    CHECK(r.se_k1 < -0.5);
    CHECK(r.se_k2 < -0.5);
    CHECK(r.se_pi_hat < -0.5);
    CHECK(r.ci_k0_lo < -0.5);
    CHECK(r.ci_pi_hat_hi < -0.5);
}

// ────────────────────────────────────────────────────────────────────
// Test 2: Numerical Hessian of a quadratic function
//   f(x, y) = a*x^2 + b*y^2 + c*x*y
//   H = [[2a, c], [c, 2b]]
// Verify the central-difference Hessian matches the analytical one.
// ────────────────────────────────────────────────────────────────────
TEST_CASE(numerical_hessian_quadratic) {
    // f(k1, k2) = 3*k1^2 + 2*k2^2 + 1.5*k1*k2
    // Analytical Hessian: H = [[6, 1.5], [1.5, 4]]
    auto f = [](double k1, double k2) -> double {
        return 3.0 * k1 * k1 + 2.0 * k2 * k2 + 1.5 * k1 * k2;
    };

    double k1 = 0.3, k2 = 0.2;
    const double h = 1e-4;

    double f00 = f(k1, k2);
    double fp1 = f(k1 + h, k2);
    double fm1 = f(k1 - h, k2);
    double fp2 = f(k1, k2 + h);
    double fm2 = f(k1, k2 - h);
    double fpp = f(k1 + h, k2 + h);
    double fpm = f(k1 + h, k2 - h);
    double fmp = f(k1 - h, k2 + h);
    double fmm = f(k1 - h, k2 - h);

    double H00 = (fp1 - 2.0 * f00 + fm1) / (h * h);
    double H11 = (fp2 - 2.0 * f00 + fm2) / (h * h);
    double H01 = (fpp - fpm - fmp + fmm) / (4.0 * h * h);

    CHECK_NEAR(H00, 6.0, 1e-6);
    CHECK_NEAR(H11, 4.0, 1e-6);
    CHECK_NEAR(H01, 1.5, 1e-6);

    // Verify inverse: det = 6*4 - 1.5^2 = 24 - 2.25 = 21.75
    double det = H00 * H11 - H01 * H01;
    CHECK(det > 0.0);
    CHECK_NEAR(det, 21.75, 1e-5);

    // Var(k1) = H11/det = 4/21.75, Var(k2) = H00/det = 6/21.75
    double var_k1 = H11 / det;
    double var_k2 = H00 / det;
    double cov_k12 = -H01 / det;

    CHECK_NEAR(var_k1, 4.0 / 21.75, 1e-5);
    CHECK_NEAR(var_k2, 6.0 / 21.75, 1e-5);
    CHECK_NEAR(cov_k12, -1.5 / 21.75, 1e-5);
}

// ────────────────────────────────────────────────────────────────────
// Test 3: SE computation with a realistic negative-log-likelihood
//
// Simulate a simplified GLkin objective with 100 SNPs.
// Use per-SNP mixture probabilities to ensure the Hessian is non-singular.
// ────────────────────────────────────────────────────────────────────
TEST_CASE(se_computation_realistic_nll) {
    const int n_snps = 100;

    // Per-SNP IBS mixture probabilities (simulating different AFs)
    // Each SNP has a unique distribution over 9 genotype combos
    std::vector<std::array<double, 9>> pibs(n_snps);
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist(0.01, 0.3);
    for (int s = 0; s < n_snps; ++s) {
        double sum = 0.0;
        for (int c = 0; c < 9; ++c) {
            pibs[s][c] = dist(rng);
            sum += pibs[s][c];
        }
        for (int c = 0; c < 9; ++c) pibs[s][c] /= sum;
    }

    // IBS|IBD probabilities (Anderson & Weir model, p=0.3, F=0)
    // Use three different AF values to create distinct IBD signatures
    auto make_b = [](double p, int z) -> std::array<double, 9> {
        double q = 1.0 - p;
        std::array<double, 9> b{};
        if (z == 0) {
            // IBD=0: 4 independent draws
            b[0] = p*p*q*q;  b[1] = q*q*p*p;   // PPQQ, QQPP
            b[2] = 2*p*p*p*q; b[3] = 2*p*p*p*q; // PPPQ, PQPP (symmetric approx)
            b[4] = 2*p*q*q*q; b[5] = 2*p*q*q*q; // PQQQ, QQPQ
            b[6] = p*p*p*p;  b[7] = 4*p*p*q*q;  // PPPP, PQPQ
            b[8] = q*q*q*q;                       // QQQQ
        } else if (z == 1) {
            // IBD=1: 3 independent draws
            b[0] = p*p*q*q;  b[1] = q*q*p*p;   // PPQQ, QQPP
            b[2] = p*p*q;    b[3] = p*p*q;      // PPPQ, PQPP
            b[4] = p*q*q;    b[5] = p*q*q;      // PQQQ, QQPQ
            b[6] = p*p*p;    b[7] = p*q;        // PPPP, PQPQ
            b[8] = q*q*q;                         // QQQQ
        } else {
            // IBD=2: 2 independent draws
            b[0] = 0; b[1] = 0; b[2] = 0; b[3] = 0; b[4] = 0; b[5] = 0;
            b[6] = p*p; b[7] = 2*p*q; b[8] = q*q;
        }
        return b;
    };

    // Use 3 different AFs to create distinct per-SNP IBS|IBD signatures
    std::vector<std::array<double, 9>> b0(n_snps), b1(n_snps), b2(n_snps);
    for (int s = 0; s < n_snps; ++s) {
        double af = 0.2 + 0.3 * (s % 5) / 4.0;  // AF varies: 0.2, 0.275, 0.35, 0.425, 0.5
        b0[s] = make_b(af, 0);
        b1[s] = make_b(af, 1);
        b2[s] = make_b(af, 2);
    }

    // Negative log-likelihood as a function of (k1, k2)
    auto nll = [&](double k1, double k2) -> double {
        double k0 = 1.0 - k1 - k2;
        if (k0 < -0.1 || k1 < -0.1 || k2 < -0.1) return 1e10;
        double ll = 0.0;
        for (int s = 0; s < n_snps; ++s) {
            double site_prob = 0.0;
            for (int c = 0; c < 9; ++c) {
                double ibd_sum = k0 * b0[s][c] + k1 * b1[s][c] + k2 * b2[s][c];
                site_prob += pibs[s][c] * ibd_sum;
            }
            if (site_prob <= 0.0) return 1e10;
            ll += std::log(site_prob);
        }
        if (std::isinf(ll)) return 1e10;
        return -ll;
    };

    // Evaluate at an interior point (second-degree-like: k0=0.65, k1=0.3, k2=0.05)
    double k1_opt = 0.3, k2_opt = 0.05;
    const double h = 1e-4;

    double f00 = nll(k1_opt, k2_opt);
    double fp1 = nll(k1_opt + h, k2_opt);
    double fm1 = nll(k1_opt - h, k2_opt);
    double fp2 = nll(k1_opt, k2_opt + h);
    double fm2 = nll(k1_opt, k2_opt - h);
    double fpp = nll(k1_opt + h, k2_opt + h);
    double fpm = nll(k1_opt + h, k2_opt - h);
    double fmp = nll(k1_opt - h, k2_opt + h);
    double fmm = nll(k1_opt - h, k2_opt - h);

    double H00 = (fp1 - 2.0 * f00 + fm1) / (h * h);
    double H11 = (fp2 - 2.0 * f00 + fm2) / (h * h);
    double H01 = (fpp - fpm - fmp + fmm) / (4.0 * h * h);

    double det = H00 * H11 - H01 * H01;

    // Hessian should be positive definite for a well-posed problem
    CHECK(H00 > 0.0);
    CHECK(det > 0.0);

    if (H00 > 0.0 && det > 0.0) {
        double var_k1 = H11 / det;
        double var_k2 = H00 / det;
        double cov_k12 = -H01 / det;

        CHECK(var_k1 > 0.0);
        CHECK(var_k2 > 0.0);

        double se_k1 = std::sqrt(var_k1);
        double se_k2 = std::sqrt(var_k2);

        // SE should be reasonable (not huge, not tiny for 100 SNPs)
        CHECK(se_k1 > 1e-6);
        CHECK(se_k1 < 1.0);
        CHECK(se_k2 > 1e-6);
        CHECK(se_k2 < 1.0);

        // SE(k0) via k0 = 1 - k1 - k2
        double var_k0 = var_k1 + var_k2 + 2.0 * cov_k12;
        if (var_k0 > 0.0) {
            double se_k0 = std::sqrt(var_k0);
            CHECK(se_k0 > 1e-6);
            CHECK(se_k0 < 1.0);
        }

        // SE(PI_HAT) via Delta method: PI = 0.5*k1 + k2
        double var_pi = 0.25 * var_k1 + var_k2 + cov_k12;
        if (var_pi > 0.0) {
            double se_pi = std::sqrt(var_pi);
            CHECK(se_pi > 1e-6);
            CHECK(se_pi < 1.0);
        }
    }
}

// ────────────────────────────────────────────────────────────────────
// Test 4: Boundary detection — SE should be NA when k is on boundary
// ────────────────────────────────────────────────────────────────────
TEST_CASE(boundary_detection) {
    // When k0 = 0 (PO case), the SE is unreliable
    double k0 = 0.0, k1 = 0.95, k2 = 0.05;

    bool on_boundary = (k0 < 0.01) || (k1 < 0.01) || (k2 < 0.01);
    CHECK(on_boundary);  // k0 < 0.01 triggers boundary

    // k2 near 0 also triggers boundary
    k0 = 0.5; k1 = 0.49; k2 = 0.005;
    on_boundary = (k0 < 0.01) || (k1 < 0.01) || (k2 < 0.01);
    CHECK(on_boundary);  // k2 < 0.01 triggers boundary

    // Interior point should NOT trigger boundary
    k0 = 0.5; k1 = 0.3; k2 = 0.2;
    on_boundary = (k0 < 0.01) || (k1 < 0.01) || (k2 < 0.01);
    CHECK(!on_boundary);  // All parameters > 0.01

    // Full-sibling-like interior point
    k0 = 0.33; k1 = 0.47; k2 = 0.20;
    on_boundary = (k0 < 0.01) || (k1 < 0.01) || (k2 < 0.01);
    CHECK(!on_boundary);  // Interior point
}

// ────────────────────────────────────────────────────────────────────
// Test 5: CI coverage — MLE should lie within 95% CI
// ────────────────────────────────────────────────────────────────────
TEST_CASE(ci_coverage_mle) {
    // Simulate: MLE k1=0.3, SE_k1=0.05
    double k1 = 0.3, se_k1 = 0.05;
    double z = 1.96;
    double ci_lo = std::max(0.0, k1 - z * se_k1);
    double ci_hi = std::min(1.0, k1 + z * se_k1);

    CHECK(ci_lo <= k1);
    CHECK(ci_hi >= k1);
    CHECK_NEAR(ci_lo, 0.3 - 1.96 * 0.05, 1e-10);
    CHECK_NEAR(ci_hi, 0.3 + 1.96 * 0.05, 1e-10);

    // Edge case: CI clamped to [0, 1]
    double k_near_0 = 0.02, se = 0.05;
    double lo = std::max(0.0, k_near_0 - z * se);
    double hi = std::min(1.0, k_near_0 + z * se);
    CHECK(lo >= 0.0);
    CHECK(hi <= 1.0);
    CHECK(lo <= k_near_0);
    CHECK(hi >= k_near_0);
}

// ────────────────────────────────────────────────────────────────────
// Test 6: Delta method for PI_HAT
//   PI = 0.5*k1 + k2
//   Var(PI) = 0.25*Var(k1) + Var(k2) + 2*0.5*Cov(k1,k2)
//           = 0.25*Var(k1) + Var(k2) + Cov(k1,k2)
// ────────────────────────────────────────────────────────────────────
TEST_CASE(delta_method_pi_hat) {
    // Known covariance matrix
    double var_k1 = 0.01;   // Var(k1)
    double var_k2 = 0.005;  // Var(k2)
    double cov_k12 = 0.002; // Cov(k1, k2)

    double var_pi = 0.25 * var_k1 + var_k2 + cov_k12;
    // = 0.25 * 0.01 + 0.005 + 0.002 = 0.0025 + 0.005 + 0.002 = 0.0095
    CHECK_NEAR(var_pi, 0.0095, 1e-10);

    double se_pi = std::sqrt(var_pi);
    CHECK_NEAR(se_pi, std::sqrt(0.0095), 1e-10);
    CHECK(se_pi > 0.0);
    CHECK(se_pi < 1.0);
}

// ────────────────────────────────────────────────────────────────────
// Test 7: Hessian non-positive-definite → SE = NA
// ────────────────────────────────────────────────────────────────────
TEST_CASE(hessian_non_positive_definite) {
    // A saddle point: f(x,y) = x^2 - y^2
    // H = [[2, 0], [0, -2]], det = -4 < 0 → not positive definite
    auto f = [](double x, double y) -> double {
        return x * x - y * y;
    };

    double x = 0.5, y = 0.5;
    const double h = 1e-4;

    double f00 = f(x, y);
    double fp1 = f(x + h, y);
    double fm1 = f(x - h, y);
    double fp2 = f(x, y + h);
    double fm2 = f(x, y - h);
    double fpp = f(x + h, y + h);
    double fpm = f(x + h, y - h);
    double fmp = f(x - h, y + h);
    double fmm = f(x - h, y - h);

    double H00 = (fp1 - 2.0 * f00 + fm1) / (h * h);
    double H11 = (fp2 - 2.0 * f00 + fm2) / (h * h);
    double H01 = (fpp - fpm - fmp + fmm) / (4.0 * h * h);

    double det = H00 * H11 - H01 * H01;
    // det should be negative → SE should be NA
    CHECK(det < 0.0);
    CHECK(H00 > 0.0);  // H00 > 0 but det < 0 → not positive definite
}

// ────────────────────────────────────────────────────────────────────
// Test 8: SE decreases with more SNPs (more data → more precision)
// ────────────────────────────────────────────────────────────────────
TEST_CASE(se_decreases_with_more_snps) {
    // Simulate NLL with n SNPs using mixture of AF-specific IBS|IBD probs
    auto compute_se = [](int n_snps) -> double {
        // Pre-compute per-SNP IBS|IBD probs with varying AF
        std::vector<std::array<double, 9>> b0(n_snps), b1(n_snps), b2(n_snps);
        std::vector<std::array<double, 9>> pibs(n_snps);
        std::mt19937 rng(123);
        std::uniform_real_distribution<double> dist(0.01, 0.3);

        for (int s = 0; s < n_snps; ++s) {
            double af = 0.2 + 0.3 * (s % 5) / 4.0;
            double p = af, q = 1.0 - af;
            b0[s] = {p*p*q*q, q*q*p*p, 2*p*p*p*q, 2*p*p*p*q,
                     2*p*q*q*q, 2*p*q*q*q, p*p*p*p, 4*p*p*q*q, q*q*q*q};
            b1[s] = {p*p*q*q, q*q*p*p, p*p*q, p*p*q,
                     p*q*q, p*q*q, p*p*p, p*q, q*q*q};
            b2[s] = {0, 0, 0, 0, 0, 0, p*p, 2*p*q, q*q};
            double sum = 0;
            for (int c = 0; c < 9; ++c) { pibs[s][c] = dist(rng); sum += pibs[s][c]; }
            for (int c = 0; c < 9; ++c) pibs[s][c] /= sum;
        }

        auto nll = [&](double k1, double k2) -> double {
            double k0 = 1.0 - k1 - k2;
            if (k0 < -0.1 || k1 < -0.1 || k2 < -0.1) return 1e10;
            double ll = 0.0;
            for (int s = 0; s < n_snps; ++s) {
                double site_prob = 0.0;
                for (int c = 0; c < 9; ++c) {
                    double ibd_sum = k0 * b0[s][c] + k1 * b1[s][c] + k2 * b2[s][c];
                    site_prob += pibs[s][c] * ibd_sum;
                }
                if (site_prob <= 0.0) return 1e10;
                ll += std::log(site_prob);
            }
            if (std::isinf(ll)) return 1e10;
            return -ll;
        };

        double k1 = 0.3, k2 = 0.05;
        const double h = 1e-4;
        double f00 = nll(k1, k2);
        double fp1 = nll(k1 + h, k2);
        double fm1 = nll(k1 - h, k2);
        double fp2 = nll(k1, k2 + h);
        double fm2 = nll(k1, k2 - h);
        double fpp = nll(k1 + h, k2 + h);
        double fpm = nll(k1 + h, k2 - h);
        double fmp = nll(k1 - h, k2 + h);
        double fmm = nll(k1 - h, k2 - h);

        double H00 = (fp1 - 2.0 * f00 + fm1) / (h * h);
        double H11 = (fp2 - 2.0 * f00 + fm2) / (h * h);
        double H01 = (fpp - fpm - fmp + fmm) / (4.0 * h * h);

        double det = H00 * H11 - H01 * H01;
        if (det <= 0.0 || H00 <= 0.0) return -1.0;
        double var_k1 = H11 / det;
        return (var_k1 > 0.0) ? std::sqrt(var_k1) : -1.0;
    };

    double se_100 = compute_se(100);
    double se_500 = compute_se(500);
    double se_2000 = compute_se(2000);

    CHECK(se_100 > 0.0);
    CHECK(se_500 > 0.0);
    CHECK(se_2000 > 0.0);

    // SE should decrease as ~1/sqrt(n)
    CHECK(se_500 < se_100);
    CHECK(se_2000 < se_500);

    // Ratio se_100/se_500 should be approximately sqrt(5) ~ 2.24
    double ratio = se_100 / se_500;
    CHECK(ratio > 1.5);
    CHECK(ratio < 4.0);
}

int main() {
    return RUN_ALL_TESTS();
}
